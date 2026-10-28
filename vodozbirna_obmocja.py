# -*- coding: utf-8 -*-
"""
Hidrološki profili in vodozbirna območja za izbrano območje HGO III.

Potek:
  1. vnos parametrov (input s privzetimi vrednostmi),
  2. branje HGO III in vodne mreže HIDRO5 (EPSG:3794),
  3. prenos listov DMR (LAZ, CLSS) za HGO + buffer (vir in pravila: CLAUDE.md),
  4. začasni DMR: gridiranje (celica < 2 m: TIN; >= 2 m: povprečje v celici + zapolnitev
     vrzeli), vžiganje HIDRO5, breach in zapolnitev depresij, D8 smer in akumulacija,
  5. smer toka linij (smer digitalizacije, preverjena z DMR) in topologija med HGO,
  6. profili na vsakih d m vzdolž vodotokov (stacionaža od izliva vodotoka navzgor),
     le do minimalne površine zaledja,
  7. vodozbirna območja (DMR del + poligoni gorvodnih HGO III), topološka poenostavitev
     Douglas-Peucker (toleranca = velikost celice),
  8. zapis v <izhod>/<HGOIII_ID>.gpkg (sloji z oznako HGO, npr. 111_profili; tudi uporabljena
     vodna mreža od izliva do zadnjega profila),
     brisanje začasnih podatkov (LAZ, DMR, rastri).

Obdelava DMR (gridiranje z zapolnitvijo vrzeli, breach, vžiganje) je izrecno zahtevana
in zabeležena v tabeli 'parametri'.

Zagon (64-bitni Python, knjižnice iz requirements.txt):
    .venv64\\Scripts\\python vodozbirna_obmocja.py
"""
import datetime as dt
import hashlib
import math
import os
import re
import shutil
import sys
import tempfile
import time
import warnings
from collections import defaultdict, deque

# PROJ_LIB/GDAL_DATA drugih namestitev (npr. PostGIS) imajo nezdružljiv proj.db;
# brez njih rasterio, pyogrio in pyproj uporabijo svoje podatke (velja le za ta proces).
for _spr in ("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"):
    os.environ.pop(_spr, None)

import numpy as np
import pandas as pd
import geopandas as gpd
import pyogrio
import rasterio
import rasterio.crs
import rasterio.features
import rasterio.fill
import rasterio.transform
from rasterio.windows import Window
import requests
import shapely
from shapely.geometry import box
import shapely.ops
from shapely.ops import unary_union
from shapely.prepared import prep
import laspy
from pyproj import CRS

# --- konstante -------------------------------------------------------------------------

EPSG = 3794
BUFFER_M = 500.0            # buffer okoli HGO za DMR
BURN_M = 5.0                # globina vžiga linij HIDRO5 v DMR
BREACH_M = 100.0            # največja dolžina preboja depresij
FLAT_M = 0.001              # naklon ravnin po zapolnitvi depresij [m na celico]
PRAG_HITRO_M = 2.0          # celica >= prag: povprečje točk v celici namesto TIN
ZAPOLNI_M = 1000.0          # največja razdalja zapolnjevanja vrzeli pri povprečenju [m]
NODATA = -32768.0
SIMPLIFY_CELIC = 1.0        # toleranca poenostavitve meja zaledij [x velikost celice DMR]
SNAP_CELIC = 1              # radij pripenjanja profila na celico z max. akumulacijo
PREMOSTITEV_M = 150.0       # največja vrzel v HIDRO5, premoščena v sloju vodne mreže
ZAMIK_IZLIV_CELIC = 10     # profil ob izlivu: vsaj toliko celic (in >= 50 m) gorvodno
PRAG_SMERI_M = 0.5          # min. višinska razlika za odločitev o smeri linije
ZAOKROZITEV_VOZLISC = 2     # decimalke pri združevanju vozlišč (1 cm)

# vir DMR (CLAUDE.md)
WFS_URL = "https://lift.clss.si/geoserver/wfs"
WFS_SLOJ = "tileindex_clss_zdruzeno"
WFS_MAX_FEATURES = 10000
PRENOS_URL = "https://assets.flycom.si/"
PREMOR_S = 3.0
MAX_PONOVITEV = 5          # ponovitve prenosa (HTTP 429, prekinjena povezava)

POLJE_HGO = "HGOIII_ID"
POLJE_ID = "ID_REG_HLI"                 # ID_REG_HLINPV, v SHP skrajšano
# VODE_ID je hierarhična šifra vodotoka: "1" = Sava, "1-31552" = pritok Save, ...
# (VTOK_ID je vrsta toka, npr. 1 = glavni vodni tok, in ni ID vodotoka)
POLJE_VODE = "VODE_ID"
POLJA_SKUPINE = [POLJE_VODE]            # prvo neprazno določa vodotok
POLJA_PRENOS = [
    "ID_REG_HLI", "HLINPV_ID", "HLINPV_IM", "VODE_ID", "HMZ_ID", "IME_ID", "IME",
    "VRSTA_ID", "VRSTA_IME", "VTOK_ID", "VTOK_IM", "TIPTV_ID", "TIPTV_IM",
    "POTEK_ID", "POTEK_IM", "STALN_ID", "VERZIJA",
]

# WhiteboxTools D8 (esri_pntr=False): vrednost -> (d_vrstica, d_stolpec)
D8 = {1: (-1, 1), 2: (0, 1), 4: (1, 1), 8: (1, 0),
      16: (1, -1), 32: (0, -1), 64: (-1, -1), 128: (-1, 0)}

WBT = None


# --- vnos ------------------------------------------------------------------------------

def vprasaj(besedilo, privzeto=None, tip=str):
    while True:
        poziv = "{} [{}]: ".format(besedilo, privzeto) if privzeto is not None else besedilo + ": "
        odg = input(poziv).strip().strip('"')
        if not odg:
            if privzeto is None:
                print("  Vnos je obvezen.")
                continue
            odg = str(privzeto)
        if tip in (int, float):
            try:
                v = tip(odg.replace(",", "."))
            except ValueError:
                print("  Neveljavno število.")
                continue
            if v <= 0:
                print("  Vrednost mora biti večja od 0.")
                continue
            return v
        return odg


def vprasaj_datoteko(besedilo_mapa, privzeta_mapa, besedilo_ime, privzeto_ime):
    while True:
        mapa = vprasaj(besedilo_mapa, privzeta_mapa)
        ime = vprasaj(besedilo_ime, privzeto_ime)
        if not os.path.splitext(ime)[1]:
            ime += ".shp"
        pot = os.path.join(mapa, ime)
        if os.path.isfile(pot):
            return mapa, pot
        print("  Datoteka ne obstaja: {}".format(pot))


def vnos_parametrov():
    mapa_skripte = os.path.dirname(os.path.abspath(__file__))
    p = {}
    p["hidro_mapa"], p["hidro_pot"] = vprasaj_datoteko(
        "Lokacija sloja vodne mreže HIDRO5", os.path.join(mapa_skripte, "SHP"),
        "Ime sloja HIDRO5", "HIDRO5_LIN_PV_TIPTV1_2.shp")
    _, p["hgo_pot"] = vprasaj_datoteko(
        "Lokacija sloja HGO", p["hidro_mapa"],
        "Ime datoteke HGO", "DRSV_HGO_III.shp")
    p["hgo_id"] = vprasaj("ID obravnavanega območja ({})".format(POLJE_HGO))
    p["razdalja"] = vprasaj("Razdalja med hidrološkimi profili [m]", 2500.0, float)
    p["celica"] = vprasaj("Velikost celice DMR [m]", 5.0, float)
    p["a_min"] = vprasaj("Minimalna površina vodozbirnega območja [km2]", 10.0, float)
    p["max_listov"] = vprasaj("Največje število listov LAZ na zagon", 500, int)
    p["izhod_mapa"] = vprasaj("Lokacija izhodnega sloja (GPKG)", os.path.join(mapa_skripte, "IZHOD"))
    ime = re.sub(r'[\\/:*?"<>|]', "_", p["hgo_id"])
    p["izhod_pot"] = os.path.join(p["izhod_mapa"], ime + ".gpkg")
    return p


# --- pomožne ---------------------------------------------------------------------------

ST_KORAKOV = 7
_korak = [0]
_vrstica = [""]


def korak(besedilo):
    konec_napredka()
    _korak[0] += 1
    print("\n[{}] {}/{} {}".format(dt.datetime.now().strftime("%H:%M:%S"),
                                   _korak[0], ST_KORAKOV, besedilo), flush=True)


def napredek(besedilo):
    """Sprotni izpis v isti vrstici terminala."""
    if besedilo == _vrstica[0]:
        return
    sys.stdout.write("\r" + besedilo.ljust(max(len(_vrstica[0]), len(besedilo))))
    sys.stdout.flush()
    _vrstica[0] = besedilo


def konec_napredka():
    if _vrstica[0]:
        sys.stdout.write("\n")
        sys.stdout.flush()
        _vrstica[0] = ""


def cas(sekunde):
    m, s = divmod(int(sekunde), 60)
    return "{}h {:02d}min".format(m // 60, m % 60) if m >= 60 else "{}min {:02d}s".format(m, s)


def stevec(besedilo, i, n, korak_izpisa=25):
    if i == n or i % korak_izpisa == 0:
        napredek("  {}: {}/{} ({:.0%})".format(besedilo, i, n, i / n if n else 1))
        if i == n:
            konec_napredka()


def preveri_crs(gdf, ime):
    if gdf.crs is None:
        sys.exit("Sloj {} nima koordinatnega sistema.".format(ime))
    crs = CRS(gdf.crs)
    if crs.to_epsg(min_confidence=20) != EPSG and \
            not crs.equals(CRS.from_epsg(EPSG), ignore_axis_order=True):
        sys.exit("Sloj {} ni v EPSG:{} ({}).".format(ime, EPSG, crs.name))
    return gdf.set_crs(EPSG, allow_override=True)


def wbt(ime, ponovitev=True, **kw):
    """Zažene orodje WBT; napako prepozna tudi iz izpisa (panic vrne kodo 0).

    WBT se občasno sesuje pri zapisu izhoda (npr. ko datoteko zaklene protivirusni
    program); takrat se orodje enkrat ponovi.
    """
    try:
        _wbt(ime, **kw)
    except RuntimeError as e:
        if not ponovitev:
            raise
        konec_napredka()
        print("  {} – ponavljam.".format(e))
        if kw.get("output") and os.path.exists(kw["output"]):
            os.remove(kw["output"])
        time.sleep(5)
        _wbt(ime, **kw)


def _wbt(ime, **kw):
    izpis = []
    zacetek = time.time()

    def odziv(vrstica):
        izpis.append(vrstica)
        v = vrstica.strip()
        if v.endswith("%") or v.endswith("..."):
            napredek("  WBT {}: {} [{}]".format(ime, v, cas(time.time() - zacetek)))

    napredek("  WBT {} ...".format(ime))
    koda = getattr(WBT, ime)(callback=odziv, **kw)
    napredek("  WBT {}: končano [{}]".format(ime, cas(time.time() - zacetek)))
    konec_napredka()
    napake = [v for v in izpis if "panicked" in v or v.strip().lower().startswith("error")]
    izhod = kw.get("output")
    if koda != 0 or napake or (izhod and not os.path.exists(izhod)):
        raise RuntimeError("WhiteboxTools '{}' ni uspel (koda {}): {}".format(
            ime, koda, " | ".join(napake or izpis[-5:])))


def niz(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v).strip()


class Rastr:
    """Odprt raster z branjem posameznih celic in oken."""

    def __init__(self, pot):
        self.src = rasterio.open(pot)
        self.nodata = self.src.nodata

    def close(self):
        self.src.close()

    def okno(self, r, c, rad):
        r0, c0 = max(0, r - rad), max(0, c - rad)
        r1, c1 = min(self.src.height, r + rad + 1), min(self.src.width, c + rad + 1)
        if r0 >= r1 or c0 >= c1:
            return None, r0, c0
        a = self.src.read(1, window=Window(c0, r0, c1 - c0, r1 - r0)).astype("float64")
        if self.nodata is not None:
            a[a == self.nodata] = np.nan
        return a, r0, c0

    def pripni(self, x, y, rad):
        """Celica z največjo vrednostjo v radiju rad okoli (x, y): (r, c, vrednost)."""
        r, c = self.src.index(x, y)
        a, r0, c0 = self.okno(r, c, rad)
        if a is None or np.all(np.isnan(a)):
            return None
        k = int(np.nanargmax(a))
        rr, cc = np.unravel_index(k, a.shape)
        return r0 + int(rr), c0 + int(cc), float(a[rr, cc])

    def xy(self, r, c):
        return self.src.xy(r, c)


# --- 1. vektorski podatki --------------------------------------------------------------

def kodiranje(pot):
    """.cpg lahko navaja UTF-8, besedila pa so v cp1250 – preveri na vzorcu."""
    vzorec = pyogrio.read_dataframe(pot, read_geometry=False, max_features=5000)
    besedila = vzorec.select_dtypes(include=["object", "string"]).to_numpy().ravel()
    if any(isinstance(v, str) and "�" in v for v in besedila):
        print("  {}: besedila niso UTF-8, berem kot cp1250.".format(os.path.basename(pot)))
        return "cp1250"
    return None


def preberi_hgo(pot, hgo_id):
    hgo = preveri_crs(pyogrio.read_dataframe(pot, encoding=kodiranje(pot)), pot)
    print("  Polja HGO: {}".format(", ".join(c for c in hgo.columns if c != "geometry")))
    if POLJE_HGO not in hgo.columns:
        sys.exit("Sloj HGO nima polja {}.".format(POLJE_HGO))
    hgo[POLJE_HGO] = hgo[POLJE_HGO].map(niz)
    hgo = hgo[hgo.geometry.notna() & ~hgo.geometry.is_empty].reset_index(drop=True)
    izbor = hgo[hgo[POLJE_HGO] == hgo_id.strip()]
    if izbor.empty:
        print("  Razpoložljive vrednosti {}:".format(POLJE_HGO))
        print("  " + ", ".join(sorted(hgo[POLJE_HGO].unique())))
        sys.exit("HGO {} ne obstaja.".format(hgo_id))
    if len(izbor) > 1:
        print("  OPOZORILO: {} poligonov z ID {} – združeni.".format(len(izbor), hgo_id))
    hgo_poligoni = hgo.dissolve(by=POLJE_HGO).geometry  # en (multi)poligon na HGO
    return hgo_poligoni, hgo_poligoni.loc[hgo_id.strip()]


def preberi_hidro(pot, obmocje):
    polja = pyogrio.read_info(pot)["fields"].tolist()
    if POLJE_ID not in polja:
        sys.exit("HIDRO5 nima polja {}. Polja: {}".format(POLJE_ID, ", ".join(polja)))
    kod = kodiranje(pot)
    # (a) vsa polja v območju DMR
    atributi = preveri_crs(pyogrio.read_dataframe(pot, bbox=tuple(obmocje.bounds), encoding=kod), pot)
    atributi = atributi.explode(index_parts=False).reset_index(drop=True)
    atributi["geometry"] = shapely.force_2d(atributi.geometry.to_numpy())
    print("  HIDRO5 v območju: {} linij".format(len(atributi)))
    # (b) celotna mreža, le ključna polja – za smer, topologijo in stacionažo
    stolpci = [POLJE_ID] + [p for p in POLJA_SKUPINE if p in polja]
    mreza = preveri_crs(pyogrio.read_dataframe(pot, columns=stolpci, encoding=kod), pot)
    mreza = mreza.explode(index_parts=False).reset_index(drop=True)
    mreza = mreza[mreza.geometry.notna() & ~mreza.geometry.is_empty]
    mreza = mreza[mreza.geom_type == "LineString"].reset_index(drop=True)
    mreza["geometry"] = shapely.force_2d(mreza.geometry.to_numpy())
    print("  HIDRO5 celotna mreža: {} linij".format(len(mreza)))
    return atributi, mreza


# --- 3. prenos DMR ---------------------------------------------------------------------

def zadnji_datum(besedilo):
    datumi = re.findall(r"(\d{2})\.(\d{2})\.(\d{4})", besedilo or "")
    if not datumi:
        return None
    return max(dt.date(int(y), int(m), int(d)) for d, m, y in datumi)


def seznam_listov(obmocje, max_listov):
    minx, miny, maxx, maxy = obmocje.bounds
    pricakovani = {}
    for xk in range(int(math.floor(minx / 1000)), int(math.ceil(maxx / 1000))):
        for yk in range(int(math.floor(miny / 1000)), int(math.ceil(maxy / 1000))):
            b = box(xk * 1000, yk * 1000, (xk + 1) * 1000, (yk + 1) * 1000)
            if b.intersects(obmocje):
                pricakovani[(xk, yk)] = b

    params = {
        "service": "WFS", "version": "1.1.0", "request": "GetFeature",
        "typeName": WFS_SLOJ, "outputFormat": "application/json",
        "maxFeatures": WFS_MAX_FEATURES,
        "bbox": "{:.0f},{:.0f},{:.0f},{:.0f},EPSG:{}".format(minx, miny, maxx, maxy, EPSG),
    }
    odg = requests.get(WFS_URL, params=params, timeout=120)
    odg.raise_for_status()
    zapisi = odg.json().get("features", [])
    print("  WFS vrnil {} zapisov.".format(len(zapisi)))
    if len(zapisi) >= WFS_MAX_FEATURES:
        print("  OPOZORILO: doseženo maxFeatures – seznam je lahko nepopoln.")

    listi = {}
    izven = 0
    podvojeni = defaultdict(int)
    for z in zapisi:
        lastnosti = z.get("properties", {})
        if "path_dmr" not in lastnosti or not lastnosti["path_dmr"]:
            print("  Atributi zapisa: {}".format(", ".join(sorted(lastnosti))))
            sys.exit("Zapis {} nima atributa path_dmr.".format(lastnosti.get("ti_name")))
        kljuc = (int(float(lastnosti["x_ime"])), int(float(lastnosti["y_ime"])))
        if kljuc not in pricakovani:
            izven += 1
            continue
        if kljuc in listi:
            obstojeci = listi[kljuc]
            d_nov = zadnji_datum(lastnosti.get("datum_za_1"))
            d_obs = zadnji_datum(obstojeci.get("datum_za_1"))
            podvojeni[tuple(sorted((niz(obstojeci.get("ime_ob")), niz(lastnosti.get("ime_ob")))))] += 1
            if d_nov and (d_obs is None or d_nov > d_obs):
                listi[kljuc] = lastnosti
            continue
        listi[kljuc] = lastnosti
    if izven:
        print("  {} zapisov izven območja izločenih (preverjanje bbox).".format(izven))
    for (a, b), n in podvojeni.items():
        print("  OPOZORILO: {} listov v dveh območjih zajema ({} / {}) – uporabljen novejši zajem.".format(
            n, a, b))

    manjkajoci = [k for k in pricakovani if k not in listi]
    if manjkajoci:
        print("  OPOZORILO: za {} listov ni podatkov DMR: {}".format(
            len(manjkajoci), ", ".join("{}_{}".format(*k) for k in manjkajoci[:30])))
    if len(listi) > max_listov:
        sys.exit("Potrebnih {} listov LAZ, omejitev je {}. Ustavljeno.".format(len(listi), max_listov))
    print("  Za prenos: {} listov.".format(len(listi)))
    return listi, [pricakovani[k] for k in manjkajoci]


def preveri_laz(pot, xk, yk):
    with laspy.open(pot) as f:
        glava = f.header
        try:
            crs = glava.parse_crs()
        except Exception as e:  # noqa: BLE001 – poškodovan VLR ne sme ustaviti prenosa
            print("  {}: CRS iz glave ni berljiv ({}).".format(os.path.basename(pot), e))
            crs = None
        mins, maxs, st = glava.mins, glava.maxs, glava.point_count
    vis = "ni naveden"
    if crs is None:
        crs_opis = "ni naveden"
        konec_napredka()
        print("  {}: CRS v glavi ni naveden – ni pripisan.".format(os.path.basename(pot)))
    else:
        horiz = crs.sub_crs_list[0] if crs.is_compound else crs
        if crs.is_compound and len(crs.sub_crs_list) > 1:
            vis = crs.sub_crs_list[1].name
        if horiz.to_epsg(min_confidence=20) != EPSG and \
                not horiz.equals(CRS.from_epsg(EPSG), ignore_axis_order=True):
            sys.exit("{}: CRS {} ni EPSG:{}.".format(pot, horiz.name, EPSG))
        crs_opis = crs.name
    tol = 1.0
    if mins[0] < xk * 1000 - tol or maxs[0] > (xk + 1) * 1000 + tol or \
            mins[1] < yk * 1000 - tol or maxs[1] > (yk + 1) * 1000 + tol:
        sys.exit("{}: koordinate točk ({:.0f}-{:.0f}, {:.0f}-{:.0f}) niso v listu {}_{}.".format(
            pot, mins[0], maxs[0], mins[1], maxs[1], xk, yk))
    return crs_opis, vis, int(st)


def prenesi_datoteko(seja, url, pot):
    """Prenos z upoštevanjem HTTP 429 (Retry-After) in ponovitvami ob prekinjeni povezavi.

    Zapisuje v .part in preimenuje šele po uspešnem prenosu. Vrne (velikost, sha256).
    """
    zacasna = pot + ".part"
    for poskus in range(MAX_PONOVITEV + 1):
        try:
            odg = seja.get(url, stream=True, timeout=300)
            if odg.status_code == 429 and poskus < MAX_PONOVITEV:
                cakaj = float(odg.headers.get("Retry-After", 60))
                konec_napredka()
                print("  HTTP 429 – čakam {:.0f} s.".format(cakaj))
                odg.close()
                time.sleep(cakaj)
                continue
            odg.raise_for_status()
            sha, velikost = hashlib.sha256(), 0
            with open(zacasna, "wb") as f:
                for kos in odg.iter_content(1 << 20):
                    f.write(kos)
                    sha.update(kos)
                    velikost += len(kos)
            os.replace(zacasna, pot)
            return velikost, sha.hexdigest()
        except (requests.ConnectionError, requests.Timeout,
                requests.exceptions.ChunkedEncodingError) as e:
            if poskus == MAX_PONOVITEV:
                raise
            cakaj = 10 * 2 ** poskus
            konec_napredka()
            print("  Prekinjena povezava ({}), ponovni poskus čez {} s.".format(
                type(e).__name__, cakaj))
            time.sleep(cakaj)
    raise RuntimeError("Prenos {} ni uspel.".format(url))


def prenesi_liste(listi, mapa):
    os.makedirs(mapa, exist_ok=True)
    manifest = []
    seja = requests.Session()
    zacetek, skupaj_b = time.time(), 0
    for i, ((xk, yk), l) in enumerate(sorted(listi.items())):
        if i:
            time.sleep(PREMOR_S)
        url = PRENOS_URL + l["path_dmr"].lstrip("/")
        pot = os.path.join(mapa, "DMR_{}_{}.laz".format(xk, yk))
        velikost, sha = prenesi_datoteko(seja, url, pot)
        cas_utc = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        crs_opis, vis, st = preveri_laz(pot, xk, yk)
        porocilo = re.search(r'href="?([^"\s>]+)', l.get("tehnicno_porocilo") or "")
        manifest.append({
            "LIST": l.get("ti_name", "{}_{}".format(xk, yk)), "URL": url,
            "VELIKOST_B": velikost, "SHA256": sha, "CAS_UTC": cas_utc,
            "OBMOCJE": l.get("ime_ob"), "DATUM_ZAJEMA": l.get("datum_za_1"),
            "CRS_GLAVA": crs_opis, "VIS_SISTEM": vis, "ST_TOCK": st,
            "TEH_POROCILO": porocilo.group(1) if porocilo else None,
        })
        skupaj_b += velikost
        porabljeno = time.time() - zacetek
        preostalo = porabljeno / (i + 1) * (len(listi) - i - 1)
        napredek("  Preneseno: {}/{} listov ({:.0%}), {:.0f} MB, zadnji {}, preostalo ~{}".format(
            i + 1, len(listi), (i + 1) / len(listi), skupaj_b / 1e6, manifest[-1]["LIST"], cas(preostalo)))
    konec_napredka()
    print("  Prenos končan: {} listov, {:.0f} MB, {}.".format(
        len(listi), skupaj_b / 1e6, cas(time.time() - zacetek)))
    return pd.DataFrame(manifest)


# --- 4. DMR ----------------------------------------------------------------------------

BLOK_VRSTIC = 2048


def profil_rastra(src, **spremembe):
    prof = src.profile.copy()
    prof.update(driver="GTiff", BIGTIFF="IF_SAFER", **spremembe)
    return prof


def vzgi_linije(dmr_pot, linije, izhod):
    """DMR, znižan za BURN_M v celicah, ki jih sekajo linije HIDRO5."""
    with rasterio.open(dmr_pot) as src:
        # all_touched: neprekinjena struga, sicer tok na diagonalah uide v vzporedne celice
        maska = rasterio.features.rasterize(
            ((g, 1) for g in linije if g is not None and not g.is_empty),
            out_shape=(src.height, src.width), transform=src.transform, fill=0, dtype="uint8",
            all_touched=True)
        nd = src.nodata
        with rasterio.open(izhod, "w", **profil_rastra(src)) as dst:
            for r0 in range(0, src.height, BLOK_VRSTIC):
                stevec("Vžiganje HIDRO5 (vrstice)", r0, src.height, BLOK_VRSTIC)
                okno = Window(0, r0, src.width, min(BLOK_VRSTIC, src.height - r0))
                a = src.read(1, window=okno)
                vzig = maska[r0:r0 + okno.height] == 1
                if nd is not None:
                    vzig &= a != nd
                a[vzig] -= BURN_M
                dst.write(a, 1, window=okno)
            stevec("Vžiganje HIDRO5 (vrstice)", src.height, src.height)


def rastr_iztokov(vzorec_pot, celice, izhod):
    """Raster iztočnih točk za WBT watershed: vrednost = zaporedna številka profila (1..n)."""
    po_vrsticah = defaultdict(list)
    for n, (r, c) in enumerate(celice, start=1):
        po_vrsticah[r].append((c, n))
    with rasterio.open(vzorec_pot) as src:
        prof = profil_rastra(src, dtype="int32", nodata=-32768)
        with rasterio.open(izhod, "w", **prof) as dst:
            for r0 in range(0, src.height, BLOK_VRSTIC):
                okno = Window(0, r0, src.width, min(BLOK_VRSTIC, src.height - r0))
                a = np.zeros((okno.height, src.width), dtype="int32")
                for r in range(r0, r0 + okno.height):
                    for c, n in po_vrsticah.get(r, ()):
                        a[r - r0, c] = n
                dst.write(a, 1, window=okno)


def izdelaj_dmr(mapa_laz, mapa, obmocje_shp, celica, linije):
    dmr = gridiraj_dmr(mapa_laz, mapa, obmocje_shp, celica)
    return hidroloski_dmr(dmr, mapa, celica, linije)


def gridiraj_dmr(mapa_laz, mapa, obmocje_shp, celica):
    """LAZ -> raster DMR, obrezan na HGO + buffer. Vrne pot do DMR.

    Celica >= PRAG_HITRO_M: povprečje točk v celici + zapolnitev vrzeli (hitro);
    sicer TIN po listih (WBT), ki pri gostoti ~0,5 m ohrani podrobnosti.
    """
    if celica >= PRAG_HITRO_M:
        return gridiraj_povprecje(mapa_laz, mapa, obmocje_shp, celica)
    return gridiraj_tin(mapa_laz, mapa, obmocje_shp, celica)


def gridiraj_povprecje(mapa_laz, mapa, obmocje_shp, celica):
    obmocje = pyogrio.read_dataframe(obmocje_shp).geometry.iloc[0]
    minx, miny, maxx, maxy = obmocje.bounds
    x0, y1 = math.floor(minx / celica) * celica, math.ceil(maxy / celica) * celica
    sirina = int(math.ceil((maxx - x0) / celica))
    visina = int(math.ceil((y1 - miny) / celica))
    transform = rasterio.transform.from_origin(x0, y1, celica, celica)
    vsota = np.zeros((visina, sirina), dtype="float64")
    stevilo = np.zeros((visina, sirina), dtype="uint32")

    listi = sorted(f for f in os.listdir(mapa_laz) if f.lower().endswith(".laz"))
    pokritost = []
    for st, f in enumerate(listi, start=1):
        stevec("Gridiranje LAZ (povprečje v celici)", st, len(listi), 5)
        xk, yk = (int(v) for v in re.findall(r"(\d+)_(\d+)", f)[-1])
        pokritost.append(box(xk * 1000, yk * 1000, (xk + 1) * 1000, (yk + 1) * 1000))
        tocke = laspy.read(os.path.join(mapa_laz, f), laz_backend=laspy.LazBackend.LazrsParallel)
        c = np.floor((np.asarray(tocke.x) - x0) / celica).astype("int64")
        r = np.floor((y1 - np.asarray(tocke.y)) / celica).astype("int64")
        v = (c >= 0) & (c < sirina) & (r >= 0) & (r < visina)
        if not v.any():
            continue
        c, r, zt = c[v], r[v], np.asarray(tocke.z)[v]
        del tocke, v
        # kopičenje v oknu lista (bincount čez celo mrežo bi bil prepočasen)
        c0, r0 = c.min(), r.min()
        w, h = c.max() - c0 + 1, r.max() - r0 + 1
        idx = (r - r0) * w + (c - c0)
        vsota[r0:r0 + h, c0:c0 + w] += np.bincount(idx, weights=zt, minlength=w * h).reshape(h, w)
        stevilo[r0:r0 + h, c0:c0 + w] += np.bincount(idx, minlength=w * h).reshape(h, w).astype("uint32")
        del c, r, zt, idx

    ima = stevilo > 0
    z = np.full((visina, sirina), NODATA, dtype="float32")
    z[ima] = (vsota[ima] / stevilo[ima]).astype("float32")
    del vsota, stevilo
    print("  Celic s točkami: {:.1%}".format(ima.mean()))

    napredek("  Zapolnjevanje vrzeli (GDAL fillnodata) ...")
    z = rasterio.fill.fillnodata(z, mask=ima.astype("uint8"),
                                 max_search_distance=max(10.0, ZAPOLNI_M / celica))
    del ima
    # brez podatkov ostane: izven HGO + buffer in izven prenesenih listov (tujina, manjkajoči)
    veljavno = rasterio.features.geometry_mask(
        [obmocje.intersection(unary_union(pokritost))], (visina, sirina), transform, invert=True)
    z[~veljavno | ~np.isfinite(z)] = NODATA
    konec_napredka()

    dmr = os.path.join(mapa, "dmr.tif")
    with rasterio.open(dmr, "w", driver="GTiff", width=sirina, height=visina, count=1,
                       dtype="float32", nodata=NODATA, transform=transform,
                       crs=rasterio.crs.CRS.from_epsg(EPSG), BIGTIFF="IF_SAFER") as dst:
        dst.write(z, 1)
    return dmr


def gridiraj_tin(mapa_laz, mapa, obmocje_shp, celica):
    WBT.set_working_dir(mapa_laz)
    wbt("lidar_tin_gridding", parameter="elevation", returns="all", resolution=celica)
    tif = sorted(f for f in os.listdir(mapa_laz) if f.lower().endswith(".tif"))
    if not tif:
        # starejši WBT ne bere LAZ – pretvorba v LAS
        print("  WBT ne bere LAZ, pretvarjam v LAS ...")
        mapa_las = os.path.join(mapa, "las")
        os.makedirs(mapa_las, exist_ok=True)
        for f in os.listdir(mapa_laz):
            if f.lower().endswith(".laz"):
                laspy.read(os.path.join(mapa_laz, f)).write(
                    os.path.join(mapa_las, f[:-4] + ".las"))
        shutil.rmtree(mapa_laz, ignore_errors=True)
        mapa_laz = mapa_las
        WBT.set_working_dir(mapa_laz)
        wbt("lidar_tin_gridding", parameter="elevation", returns="all", resolution=celica)
        tif = sorted(f for f in os.listdir(mapa_laz) if f.lower().endswith(".tif"))
        if not tif:
            raise RuntimeError("TIN gridiranje ni ustvarilo rastrov.")

    mozaik, dmr = os.path.join(mapa, "mozaik.tif"), os.path.join(mapa, "dmr.tif")
    wbt("mosaic", inputs=";".join(tif), output=mozaik, method="nn")
    for f in tif:
        os.remove(os.path.join(mapa_laz, f))
    wbt("clip_raster_to_polygon", i=mozaik, polygons=obmocje_shp, output=dmr,
        maintain_dimensions=False)
    os.remove(mozaik)
    return dmr


def hidroloski_dmr(dmr, mapa, celica, linije):
    """Vžig HIDRO5, breach depresij, D8 smer in akumulacija."""
    p = {k: os.path.join(mapa, k + ".tif") for k in
         ("dmr_vzig", "dmr_breach", "dmr_hidro", "pntr", "akum")}
    p["dmr"] = dmr
    vzgi_linije(p["dmr"], linije, p["dmr_vzig"])
    dist = max(10, int(round(BREACH_M / celica)))
    wbt("breach_depressions_least_cost", dem=p["dmr_vzig"], output=p["dmr_breach"],
        dist=dist, fill=False)
    os.remove(p["dmr_vzig"])
    # preostale depresije zapolni, ravnine dobijo naklon FLAT_M na celico (sicer celice brez
    # smeri); privzeti prirastek WBT je manjši od natančnosti float32 na višinah > 1000 m
    wbt("fill_depressions", dem=p["dmr_breach"], output=p["dmr_hidro"], fix_flats=True,
        flat_increment=FLAT_M)
    os.remove(p["dmr_breach"])
    wbt("d8_pointer", dem=p["dmr_hidro"], output=p["pntr"], esri_pntr=False)
    wbt("d8_flow_accumulation", i=p["pntr"], output=p["akum"], out_type="cells",
        pntr=True, esri_pntr=False)
    return p


def preberi_kazalec(pot):
    with rasterio.open(pot) as src:
        a = src.read(1, out_dtype="int32")
        nd = src.nodata
    if nd is not None:
        a[a == int(nd)] = 0
    a[(a < 0) | (a > 128)] = 0
    return a.astype("uint8")


def sledi(kazalec, r, c, cilji, max_korakov=50_000_000):
    """Sledi D8 dolvodno od (r, c) do prve celice iz cilji; vrne njeno vrednost ali None."""
    nr, nc = kazalec.shape
    for _ in range(max_korakov):
        smer = D8.get(int(kazalec[r, c]))
        if smer is None:
            return None
        r += smer[0]
        c += smer[1]
        if not (0 <= r < nr and 0 <= c < nc):
            return None
        cilj = cilji.get((r, c))
        if cilj is not None:
            return cilj
    return None


# --- 5. smer toka in topologija --------------------------------------------------------

def doloci_smer(mreza, obmocje, dmr_pot):
    """Smer digitalizacije, preverjena z DMR.

    Vrne (obrni, neskladje): obrni = skupna smer za vso mrežo (obrne se le, če je
    večina linij digitalizirana proti toku); neskladje = linije, pri katerih DMR
    kaže obratno smer (mostovi, prehodi, ravnine) – le oznaka, linija se ne obrne,
    da se veriga vodotoka ne prekine.
    """
    zacetki = shapely.get_point(mreza.geometry.to_numpy(), 0)
    konci = shapely.get_point(mreza.geometry.to_numpy(), -1)
    v_obm = mreza.sindex.query(obmocje, predicate="intersects")
    xy_z = shapely.get_coordinates(zacetki[v_obm])
    xy_k = shapely.get_coordinates(konci[v_obm])
    with rasterio.open(dmr_pot) as src:
        nd = src.nodata
        z_z = np.array([v[0] for v in src.sample(xy_z)], dtype="float64")
        z_k = np.array([v[0] for v in src.sample(xy_k)], dtype="float64")
    if nd is not None:
        z_z[z_z == nd] = np.nan
        z_k[z_k == nd] = np.nan
    obrni_dmr = z_z < z_k - PRAG_SMERI_M
    ohrani_dmr = z_z > z_k + PRAG_SMERI_M
    odloceno = obrni_dmr | ohrani_dmr
    delez = obrni_dmr[odloceno].mean() if odloceno.any() else 0.0
    globalno = bool(delez > 0.5)
    print("  Preverjenih {} linij, odločenih {}, proti DMR {:.1%}.".format(
        len(v_obm), int(odloceno.sum()), delez))
    if globalno:
        print("  Večina linij je digitalizirana proti toku – privzeta obratna smer za vso mrežo.")
    obrni = np.full(len(mreza), globalno)
    neskladje = np.zeros(len(mreza), dtype=bool)
    neskladje[v_obm[odloceno]] = obrni_dmr[odloceno] != globalno
    return obrni, neskladje


def usmeri(mreza, obrni, neskladje):
    geom = mreza.geometry.to_numpy().copy()
    geom[obrni] = shapely.reverse(geom[obrni])
    mreza = mreza.copy()
    mreza["geometry"] = gpd.GeoSeries(geom, index=mreza.index, crs=EPSG)
    mreza["OBRNJENA"] = obrni.astype(int)
    mreza["NESKLADJE"] = neskladje.astype(int)
    z = shapely.get_coordinates(shapely.get_point(geom, 0))
    k = shapely.get_coordinates(shapely.get_point(geom, -1))
    mreza["XZ"], mreza["YZ"], mreza["XK"], mreza["YK"] = z[:, 0], z[:, 1], k[:, 0], k[:, 1]
    mreza["DOLZINA"] = shapely.length(geom)
    return mreza


def topologija_hgo(mreza, hgo_poligoni):
    """Povezave g -> h (linija začne v g, konča v h) iz usmerjene mreže."""
    ids = list(hgo_poligoni.index)
    drevo = shapely.STRtree(np.array(list(hgo_poligoni), dtype=object))

    def kandidati(x, y):
        pi, hi = drevo.query(shapely.points(x, y), predicate="intersects")
        k = defaultdict(list)
        for p, h in zip(pi, hi):
            k[int(p)].append(ids[int(h)])
        return k

    kz = kandidati(mreza["XZ"].values, mreza["YZ"].values)
    kk = kandidati(mreza["XK"].values, mreza["YK"].values)
    povezave = []   # (indeks linije, hgo začetka ali None, hgo konca)
    for i in set(kz) | set(kk):
        cz, ck = kz.get(i, []), kk.get(i, [])
        if not ck:
            continue
        konec = ([h for h in ck if h not in cz] or ck)[0]
        zac = [h for h in cz if h != konec]
        if zac:
            povezave.append((i, zac[0], konec))
        elif not cz:
            povezave.append((i, None, konec))   # začetek izven vseh HGO (npr. tujina)

    # Glavni odtok HGO g = vodotok z največjo dolžino znotraj g med tistimi, ki iz g izstopajo.
    # Ostala prečkanja meje (jarki, drobni potoki čez razvodnico) ne prenašajo zaledja g.
    vode = kljuci_skupin(mreza).to_numpy()
    dolzina = mreza["DOLZINA"].to_numpy()
    dolzina_v_hgo = defaultdict(float)
    for i, hs in kz.items():
        if len(hs) == 1 and vode[i] is not None:
            dolzina_v_hgo[(hs[0], vode[i])] += dolzina[i]
    glavni = {}
    for i, g, _ in povezave:
        if g is not None and vode[i] is not None:
            if dolzina_v_hgo[(g, vode[i])] > dolzina_v_hgo.get((g, glavni.get(g)), -1.0):
                glavni[g] = vode[i]
    sprejete = [(i, g, h) for i, g, h in povezave if g is None or vode[i] == glavni.get(g)]
    zavrnjene = [(i, g, h) for i, g, h in povezave if not (g is None or vode[i] == glavni.get(g))]
    return sprejete, zavrnjene


def gorvodni_hgo(povezave, cilj):
    predniki = defaultdict(set)
    for _, g, h in povezave:
        if g is not None:
            predniki[h].add(g)

    def zaprtje(g):
        mn, vrsta = {g}, deque([g])
        while vrsta:
            for p in predniki[vrsta.popleft()]:
                if p != cilj and p not in mn:
                    mn.add(p)
                    vrsta.append(p)
        return mn

    vtoki = []
    for i, g, h in povezave:
        if h != cilj:
            continue
        if g is None:
            vtoki.append((i, None, frozenset(), True))
        else:
            vtoki.append((i, g, frozenset(zaprtje(g)), False))
    return vtoki


# --- 6. profili ------------------------------------------------------------------------

def kljuci_skupin(mreza):
    """Ključ vodotoka za vsako linijo: prvo neprazno polje iz POLJA_SKUPINE."""
    kljuci = pd.Series([None] * len(mreza), index=mreza.index, dtype=object)
    for p in POLJA_SKUPINE:
        if p in mreza.columns:
            v = mreza[p].map(niz)
            m = kljuci.isna() & (v != "")
            kljuci[m] = p + ":" + v[m]
    return kljuci


def debla(indeksi, mreza):
    """Glavni tokovi skupine linij: seznam (komponenta, [indeksi linij od izliva navzgor])."""
    vozel = lambda x, y: (round(x, ZAOKROZITEV_VOZLISC), round(y, ZAOKROZITEV_VOZLISC))  # noqa: E731
    zac, kon, dolz = {}, {}, {}
    vhodne = defaultdict(list)
    izhodne = defaultdict(int)
    for i in indeksi:
        v = mreza.loc[i]
        zac[i], kon[i], dolz[i] = vozel(v.XZ, v.YZ), vozel(v.XK, v.YK), v.DOLZINA
        vhodne[kon[i]].append(i)
        izhodne[zac[i]] += 1
    ponori = [n for n in vhodne if izhodne[n] == 0]

    def najdaljsa(ponor):
        stanje, najb = {}, {}
        sklad = [(ponor, False)]
        while sklad:
            n, konc = sklad.pop()
            if konc:
                bl, be = 0.0, None
                for e in vhodne[n]:
                    u = zac[e]
                    if stanje.get(u) == 2 and najb[u][0] + dolz[e] > bl:
                        bl, be = najb[u][0] + dolz[e], e
                najb[n] = (bl, be)
                stanje[n] = 2
                continue
            if stanje.get(n):
                continue
            stanje[n] = 1
            sklad.append((n, True))
            for e in vhodne[n]:
                if not stanje.get(zac[e]):
                    sklad.append((zac[e], False))
        pot, n = [], ponor
        while najb[n][1] is not None:
            e = najb[n][1]
            pot.append(e)
            n = zac[e]
        return najb[ponor][0], pot

    kandidati = sorted((najdaljsa(p) for p in ponori), key=lambda t: -t[0])
    uporabljene, rezultat = set(), []
    for _, pot in kandidati:
        # dodatni ponori (razcepi, prekinitve): deblo do stika z že obdelanim
        prvi_upor = next((j for j, e in enumerate(pot) if e in uporabljene), len(pot))
        pot = pot[:prvi_upor]
        if not pot:
            continue
        uporabljene.update(pot)
        rezultat.append((len(rezultat), pot))
    stranske = len(indeksi) - len(uporabljene)
    return rezultat, stranske


def tvori_profile(mreza, cilj, razdalja, celica):
    cilj_p = prep(cilj)
    zamik0 = max(ZAMIK_IZLIV_CELIC * celica, 50.0)
    v_hgo = mreza.sindex.query(cilj, predicate="intersects")
    kljuci = kljuci_skupin(mreza)
    brez = int(kljuci.iloc[v_hgo].isna().sum())
    if brez:
        print("  {} linij v HGO nima {} – izpuščene.".format(brez, "/".join(POLJA_SKUPINE)))
    skupine = set(kljuci.iloc[v_hgo].dropna())
    po_skupinah = defaultdict(list)
    for i, k in kljuci.items():
        if k in skupine:
            po_skupinah[k].append(i)

    profili, s_stranskimi, vsa_debla = [], 0, {}
    for st, k in enumerate(sorted(skupine), start=1):
        stevec("Vodotoki (debla in stacionaže)", st, len(skupine))
        deb, stranske = debla(po_skupinah[k], mreza)
        s_stranskimi += bool(stranske)
        for komp, pot in deb:
            dolzine = [mreza.at[e, "DOLZINA"] for e in pot]
            kum = np.concatenate([[0.0], np.cumsum(dolzine)])
            vsa_debla[(k, komp)] = (pot, kum)
            skupaj = kum[-1]
            for n in range(int(skupaj // razdalja) + 1):
                if n == 0 and skupaj <= zamik0:
                    continue    # prekratek vodotok: točka bi se pripela na sprejemnik
                s = n * razdalja if n else zamik0
                j = min(int(np.searchsorted(kum, s, side="right")) - 1, len(pot) - 1)
                e = pot[j]
                lok = s - kum[j]                          # od dolvodnega konca linije
                tocka = mreza.at[e, "geometry"].interpolate(max(dolzine[j] - lok, 0.0))
                if not cilj_p.contains(tocka):
                    continue
                profili.append({"SKUPINA": k, "KOMPONENTA": komp, "STAC_M": int(round(s)),
                                "LINIJA": e, "geometry": tocka})
    if s_stranskimi:
        print("  {} vodotokov ima stranske krake (razcepi), ti so izpuščeni.".format(s_stranskimi))
    print("  Kandidatov za profile v HGO: {}".format(len(profili)))
    return profili, vsa_debla


def uporabljena_mreza(profili, vsa_debla, mreza, cilj):
    """Vodna mreža, ki povezuje profile: od vsakega profila po usmerjeni mreži HIDRO5 navzdol
    (po vseh linijah, ne le po deblu vodotoka) do prvega naslednjega profila oz. do izstopa
    iz HGO. Skupaj torej od najbolj dolvodne točke v HGO do najbolj gorvodnih profilov.

    Vrne seznam (indeks linije, STAC_OD, STAC_DO, geometrija, SKUPINA, KOMPONENTA);
    stacionaža in komponenta sta znani le za linije na deblu vodotoka.
    """
    vozel = lambda x, y: (round(x, ZAOKROZITEV_VOZLISC), round(y, ZAOKROZITEV_VOZLISC))  # noqa: E731
    izhodne, vhodne = defaultdict(list), defaultdict(list)
    for i, x, y in zip(mreza.index, mreza["XZ"].to_numpy(), mreza["YZ"].to_numpy()):
        izhodne[vozel(x, y)].append(i)
    for i, x, y in zip(mreza.index, mreza["XK"].to_numpy(), mreza["YK"].to_numpy()):
        vhodne[vozel(x, y)].append(i)
    kljuci = kljuci_skupin(mreza)
    v_deblu = {}
    for (k, komp), (pot, kum) in vsa_debla.items():
        for j, e in enumerate(pot):
            v_deblu.setdefault(e, (komp, kum[j], kum[j + 1] - kum[j]))
    lege_profilov = defaultdict(list)
    for kp in profili:
        lege_profilov[kp["LINIJA"]].append(mreza.at[kp["LINIJA"], "geometry"].project(kp["geometry"]))
    cilj_p = prep(cilj)

    def gorvodne(v):
        """Vse linije, ki (posredno) tečejo v vozlišče v."""
        najdene, sklad = set(), [v]
        while sklad:
            for c in vhodne.get(sklad.pop(), ()):
                if c not in najdene:
                    najdene.add(c)
                    sklad.append(vozel(mreza.at[c, "XZ"], mreza.at[c, "YZ"]))
        return najdene

    def premosti(e, obiskane):
        """Slep konec linije e: najbližja druga linija v PREMOSTITEV_M, ki ne teče v isti slep
        konec (ponor) -> (linija, točka)."""
        konec = shapely.Point(mreza.at[e, "XK"], mreza.at[e, "YK"])
        kandidati = [c for c in mreza.sindex.query(konec.buffer(PREMOSTITEV_M))
                     if c != e and c not in obiskane]
        if kandidati:
            gorv = gorvodne(vozel(konec.x, konec.y))
            kandidati = [c for c in kandidati if c not in gorv]
        if not kandidati:
            return None, konec
        c = min(kandidati, key=lambda c: mreza.at[c, "geometry"].distance(konec))
        if mreza.at[c, "geometry"].distance(konec) > PREMOSTITEV_M:
            return None, konec
        return c, konec

    intervali = defaultdict(list)       # linija -> [(od, do)] v dolžini vzdolž usmerjene linije
    premostitve, slepi = [], []
    for kp in profili:
        e = kp["LINIJA"]
        zac = mreza.at[e, "geometry"].project(kp["geometry"])
        obiskane = set()
        while e is not None and e not in obiskane:
            obiskane.add(e)
            dolzina = mreza.at[e, "DOLZINA"]
            naprej = [s for s in lege_profilov.get(e, ()) if s > zac + 1e-6]
            if naprej:                                   # naslednji profil na isti liniji
                intervali[e].append((zac, min(naprej)))
                break
            intervali[e].append((zac, dolzina))
            if not cilj_p.contains(shapely.Point(mreza.at[e, "XK"], mreza.at[e, "YK"])):
                break                                    # izstop iz HGO
            nasl = izhodne.get(vozel(mreza.at[e, "XK"], mreza.at[e, "YK"]), [])
            if nasl:
                # na razcepu ostane na istem vodotoku, sicer prva linija
                e = next((c for c in nasl if kljuci[c] == kljuci[e]), nasl[0])
                zac = 0.0
                continue
            # slep konec (vrzel v HIDRO5): premostitev na najbližjo linijo
            c, konec = premosti(e, obiskane)
            if c is None:
                slepi.append((kp["PROFIL_ID"], mreza.at[e, POLJE_ID]))
                break
            zac = mreza.at[c, "geometry"].project(konec)
            most = shapely.LineString([konec, mreza.at[c, "geometry"].interpolate(zac)])
            if most.length > 0:
                premostitve.append(most)
            e = c

    odseki = []
    for e, ivs in intervali.items():
        linija = mreza.at[e, "geometry"]
        zdruzeni = []
        for a, b in sorted(ivs):
            if zdruzeni and a <= zdruzeni[-1][1] + 1e-6:
                zdruzeni[-1][1] = max(zdruzeni[-1][1], b)
            else:
                zdruzeni.append([a, b])
        for a, b in zdruzeni:
            g = shapely.ops.substring(linija, a, b).intersection(cilj)
            if g.is_empty or g.length == 0:
                continue
            lege = [linija.project(shapely.Point(xy)) for xy in shapely.get_coordinates(g)]
            komp, s_od, s_do = None, None, None
            if e in v_deblu:                             # stacionaža od izliva vodotoka
                komp, kum_j, dolzina = v_deblu[e]
                s_od, s_do = kum_j + dolzina - max(lege), kum_j + dolzina - min(lege)
            odseki.append((e, s_od, s_do, shapely.line_merge(g), kljuci[e], komp))
    for most in premostitve:
        g = most.intersection(cilj)
        if not g.is_empty and g.length > 0:
            odseki.append((None, None, None, g, "PREMOSTITEV", None))
    if premostitve:
        print("  Vrzeli v HIDRO5 premoščene (<= {:g} m): {}".format(PREMOSTITEV_M, len(premostitve)))
    for pid, id_linije in slepi:
        print("  OPOZORILO: pot od profila {} se konča na liniji {} (vrzel > {:g} m).".format(
            pid, id_linije, PREMOSTITEV_M))
    return odseki


# --- 7. poenostavitev ------------------------------------------------------------------

def poenostavi_pokritost(poligoni, toleranca):
    """Topološka poenostavitev Douglas-Peucker za pokritost brez prekrivanj (med-območja).

    Meje se razbijejo na loke med stičišči; vsak lok (tudi skupna meja dveh območij) se
    poenostavi enkrat, zato sosednja območja po poenostavitvi nimajo rež ali prekrivanj.
    Iz poenostavljenih lokov se ploskve ponovno sestavijo (polygonize) in dodelijo območju
    z največjim prekrivanjem. Vrne {ključ: poenostavljena geometrija}.
    """
    kljuci = list(poligoni)
    izvirni = [poligoni[k] for k in kljuci]
    meje = shapely.union_all(shapely.boundary(np.array(izvirni, dtype=object)))
    loki = shapely.line_merge(meje)
    poenostavljeni = []
    for lok in getattr(loki, "geoms", [loki]):
        # zaprt lok (otok brez stičišč) ne sme propasti v točko
        s = shapely.simplify(lok, toleranca, preserve_topology=bool(lok.is_ring))
        if not s.is_empty:
            poenostavljeni.append(s)
    vozlana = shapely.union_all(np.array(poenostavljeni, dtype=object))   # ponovno vozlanje
    ploskve = shapely.polygonize(list(getattr(vozlana, "geoms", [vozlana])))
    drevo = shapely.STRtree(np.array(izvirni, dtype=object))
    deli = defaultdict(list)
    for f in getattr(ploskve, "geoms", [ploskve]):
        kandidati = drevo.query(f, predicate="intersects")
        if len(kandidati) == 0:
            continue
        prekrivanja = [izvirni[j].intersection(f).area for j in kandidati]
        j = int(np.argmax(prekrivanja))
        if prekrivanja[j] > 0.5 * f.area:      # ploskve v režah med območji se izpustijo
            deli[kljuci[kandidati[j]]].append(f)
    return {k: shapely.union_all(np.array(v, dtype=object)) for k, v in deli.items()}


# --- glavni potek ----------------------------------------------------------------------

def zapisi_tabelo(df, pot, sloj):
    try:
        pyogrio.write_dataframe(df, pot, layer=sloj, driver="GPKG")
    except Exception as e:  # noqa: BLE001
        csv = os.path.splitext(pot)[0] + "_{}.csv".format(sloj)
        df.to_csv(csv, index=False, encoding="utf-8")
        print("  Tabela {} ni zapisana v GPKG ({}), zapisana v {}.".format(sloj, e, csv))


def main():
    global WBT
    print("Hidrološki profili in vodozbirna območja (HIDRO5, HGO III, DMR CLSS)\n")
    p = vnos_parametrov()
    os.makedirs(p["izhod_mapa"], exist_ok=True)
    if os.path.exists(p["izhod_pot"]):
        if input("Datoteka {} obstaja. Prepišem? [d/N]: ".format(p["izhod_pot"])).strip().lower() != "d":
            sys.exit("Prekinjeno.")
        # stara datoteka ostane do uspešnega konca; tu le preverimo, da ni zaklenjena
        try:
            open(p["izhod_pot"], "r+b").close()
        except PermissionError:
            sys.exit("Datoteke {} ni mogoče prepisati – je odprta v drugem programu (npr. QGIS)? "
                     "Zaprite jo in ponovite.".format(p["izhod_pot"]))
    celica, hgo_id = p["celica"], p["hgo_id"].strip()
    zacetek = dt.datetime.now()

    korak("Branje HGO in HIDRO5")
    hgo_poligoni, cilj = preberi_hgo(p["hgo_pot"], hgo_id)
    obmocje = cilj.buffer(BUFFER_M)
    print("  HGO {}: {:.2f} km2".format(hgo_id, cilj.area / 1e6))
    atributi, mreza = preberi_hidro(p["hidro_pot"], obmocje)
    atributi_id = atributi.drop_duplicates(POLJE_ID).set_index(POLJE_ID)

    import whitebox
    WBT = whitebox.WhiteboxTools()
    WBT.set_verbose_mode(True)   # izpis gre v callback funkcije wbt(), ne na zaslon

    zacasna = tempfile.mkdtemp(prefix="_tmp_vzo_", dir=p["izhod_mapa"])
    akum = None
    try:
        korak("Seznam in prenos listov DMR")
        listi, manjkajoci = seznam_listov(obmocje, p["max_listov"])
        if not listi:
            sys.exit("Za območje ni listov DMR.")
        manifest = prenesi_liste(listi, os.path.join(zacasna, "laz"))

        korak("Izdelava začasnega DMR ({} m)".format(celica))
        obmocje_shp = os.path.join(zacasna, "obmocje.shp")
        pyogrio.write_dataframe(gpd.GeoDataFrame({"ID": [1]}, geometry=[obmocje], crs=EPSG),
                                obmocje_shp)
        r = izdelaj_dmr(os.path.join(zacasna, "laz"), zacasna, obmocje_shp, celica,
                        atributi.geometry.to_numpy())

        korak("Smer toka in topologija HGO")
        mreza = usmeri(mreza, *doloci_smer(mreza, obmocje, r["dmr"]))
        povezave, zavrnjene = topologija_hgo(mreza, hgo_poligoni)
        vtoki = gorvodni_hgo(povezave, hgo_id)
        # prečkanja meje, ki niso glavni odtok sosednjega HGO: le zabeležena, brez zaledja
        stranski = [(i, g, frozenset(), False) for i, g, h in zavrnjene if h == hgo_id]
        print("  Vtokov v HGO: {} glavnih odtokov gorvodnih HGO, {} iz tujine, "
              "{} drugih prečkanj meje (brez zaledja HGO)".format(
                  sum(1 for v in vtoki if not v[3]), sum(1 for v in vtoki if v[3]), len(stranski)))
        ploscine_hgo = hgo_poligoni.area

        korak("Profili")
        kandidati, vsa_debla = tvori_profile(mreza, cilj, p["razdalja"], celica)
        akum = Rastr(r["akum"])
        izbrani, po_celici, podvojeni = [], {}, 0
        for st, kp in enumerate(kandidati, start=1):
            stevec("Pripenjanje profilov na DMR", st, len(kandidati))
            pr = akum.pripni(kp["geometry"].x, kp["geometry"].y, SNAP_CELIC)
            if pr is None:
                continue
            celica_rc = (pr[0], pr[1])
            if celica_rc in po_celici:
                podvojeni += 1
                continue
            kp["CELICA"], kp["AKUM"] = celica_rc, pr[2]
            po_celici[celica_rc] = len(izbrani)
            izbrani.append(kp)
        kandidati = izbrani
        if podvojeni:
            print("  {} profilov pade v isto celico kot že obstoječi – izpuščeni.".format(podvojeni))

        napredek("  Branje rastra smeri toka ...")
        kazalec = preberi_kazalec(r["pntr"])
        dolv = []
        for st, kp in enumerate(kandidati, start=1):
            stevec("Sledenje toka do dolvodnega profila", st, len(kandidati))
            dolv.append(sledi(kazalec, kp["CELICA"][0], kp["CELICA"][1], po_celici))

        gorv = [set() for _ in kandidati]
        izven = [False] * len(kandidati)
        zapisi_vtokov = []
        vsi_vtoki = [v + (True,) for v in vtoki] + [v + (False,) for v in stranski]
        for st, (i, iz_hgo, mn, tujina, glavni) in enumerate(vsi_vtoki, start=1):
            stevec("Vtoki iz gorvodnih HGO", st, len(vsi_vtoki), 1)
            linija = mreza.at[i, "geometry"]
            presek = shapely.get_coordinates(linija.intersection(cilj.boundary))
            # najbolj dolvodno prečkanje meje HGO
            proj = max((linija.project(shapely.Point(xy)) for xy in presek), default=0.0)
            vstop = linija.interpolate(min(proj + 2 * celica, linija.length))
            pr = akum.pripni(vstop.x, vstop.y, SNAP_CELIC)
            prvi = None
            if pr is not None:
                prvi = po_celici.get((pr[0], pr[1]))
                if prvi is None:
                    prvi = sledi(kazalec, pr[0], pr[1], po_celici)
            n = prvi
            while n is not None:
                gorv[n] |= mn
                izven[n] = izven[n] or tujina
                n = dolv[n]
            zapisi_vtokov.append({
                POLJE_ID: mreza.at[i, POLJE_ID], "VODE_ID": niz(mreza.at[i, POLJE_VODE]),
                "IZ_HGO": iz_hgo, "GLAVNI_ODTOK": int(glavni), "GORV_HGO": ";".join(sorted(mn)),
                "IZVEN_HGO": int(tujina), "A_GORV_KM2": sum(ploscine_hgo[h] for h in mn) / 1e6,
                "PRVI_PROFIL": None if prvi is None else prvi, "geometry": vstop})

        a_dmr = [kp["AKUM"] * celica ** 2 / 1e6 for kp in kandidati]
        a_gorv = [sum(ploscine_hgo[h] for h in g) / 1e6 for g in gorv]
        ohrani = [n for n in range(len(kandidati)) if a_dmr[n] + a_gorv[n] >= p["a_min"]]
        print("  Profilov z zaledjem >= {} km2: {} od {}".format(p["a_min"], len(ohrani), len(kandidati)))
        if not ohrani:
            sys.exit("Noben profil ne doseže minimalne površine.")

        # PROFIL_ID in ponovno drevo samo med ohranjenimi
        profili = []
        for n in ohrani:
            kp = kandidati[n]
            _, vred = kp["SKUPINA"].split(":", 1)
            pid = "{}_{}_{:06d}".format(hgo_id, vred, kp["STAC_M"])
            if kp["KOMPONENTA"]:
                pid += "_K{}".format(kp["KOMPONENTA"])
            kp.update({"PROFIL_ID": pid, "A_DMR": a_dmr[n], "A_GORV": a_gorv[n],
                       "GORV": gorv[n], "IZVEN": izven[n], "KAND": n})
            profili.append(kp)
        cilji = {kp["CELICA"]: m for m, kp in enumerate(profili)}
        for st, kp in enumerate(profili, start=1):
            stevec("Drevo profilov", st, len(profili))
            nd = sledi(kazalec, kp["CELICA"][0], kp["CELICA"][1], cilji)
            kp["DOLV"] = nd
        for z in zapisi_vtokov:
            # prvi ohranjeni profil dolvodno od vtoka
            n = z["PRVI_PROFIL"]
            while n is not None and "PROFIL_ID" not in kandidati[n]:
                n = dolv[n]
            z["PRVI_PROFIL"] = kandidati[n]["PROFIL_ID"] if n is not None else None
        del kazalec

        korak("Vodozbirna območja")
        xy = [akum.xy(*kp["CELICA"]) for kp in profili]
        akum.close()
        akum = None
        iztoki_tif = os.path.join(zacasna, "iztoki.tif")
        rastr_iztokov(r["pntr"], [kp["CELICA"] for kp in profili], iztoki_tif)
        ws_tif = os.path.join(zacasna, "porecja.tif")
        wbt("watershed", d8_pntr=r["pntr"], pour_pts=iztoki_tif, output=ws_tif, esri_pntr=False)
        ws_shp = os.path.join(zacasna, "porecja.shp")
        wbt("raster_to_vector_polygons", i=ws_tif, output=ws_shp)
        # koordinate rastrov izhajajo iz LAZ, preverjenih v EPSG:3794 (preveri_laz)
        with warnings.catch_warnings():
            # WBT zapiše obroče v obratnem vrstnem redu; GDAL jih ob branju popravi
            warnings.filterwarnings("ignore", message=".*winding order.*")
            pol = pyogrio.read_dataframe(ws_shp).set_crs(EPSG, allow_override=True)
        polje_v = next((c for c in pol.columns if c.upper() == "VALUE"), None)
        if polje_v is None:
            raise RuntimeError("Poligoni porečij nimajo polja VALUE: {}".format(list(pol.columns)))
        # vrednost porečja = zaporedna številka profila v rastru iztokov (1..n)
        pol = pol[(pol[polje_v] >= 1) & (pol[polje_v] <= len(profili))]
        pol = pol.dissolve(by=polje_v)
        inkrement = {int(v) - 1: g for v, g in pol.geometry.items()}
        napredek("  Topološka poenostavitev (Douglas-Peucker, toleranca {:g} m) ...".format(
            celica * SIMPLIFY_CELIC))
        pred = sum(g.area for g in inkrement.values())
        tock_pred = sum(shapely.get_num_coordinates(g) for g in inkrement.values())
        inkrement = poenostavi_pokritost(inkrement, celica * SIMPLIFY_CELIC)
        konec_napredka()
        print("  Poenostavitev: točk {} -> {}, sprememba ploščine {:+.3%}".format(
            tock_pred, sum(shapely.get_num_coordinates(g) for g in inkrement.values()),
            sum(g.area for g in inkrement.values()) / pred - 1))

        otroci = defaultdict(list)
        for m, kp in enumerate(profili):
            if kp["DOLV"] is not None:
                otroci[kp["DOLV"]].append(m)
        rob = unary_union([obmocje.boundary] + [b.boundary for b in manjkajoci]).buffer(2 * celica)
        rob_p = prep(rob)
        skupno_dmr, hgo_unija = {}, {}
        for st, m in enumerate(sorted(range(len(profili)), key=lambda m: profili[m]["AKUM"]), start=1):
            stevec("Združevanje zaledij (od izvira navzdol)", st, len(profili))
            deli = [inkrement[m]] if m in inkrement else []
            deli += [skupno_dmr[o] for o in otroci[m] if o in skupno_dmr]
            if deli:
                skupno_dmr[m] = unary_union(deli)
        profili_vrstice, obmocja_vrstice, medobm_vrstice = [], [], []
        danes = dt.date.today().isoformat()
        for m, kp in enumerate(profili):
            stevec("Gorvodni HGO in atributi", m + 1, len(profili))
            g = frozenset(kp["GORV"])
            if g and g not in hgo_unija:
                hgo_unija[g] = unary_union([hgo_poligoni[h] for h in g])
            deli = [x for x in (skupno_dmr.get(m), hgo_unija.get(g)) if x is not None]
            skupaj = unary_union(deli) if deli else None
            dotik = m in skupno_dmr and rob_p.intersects(skupno_dmr[m])
            linija = mreza.loc[kp["LINIJA"]]
            vir = atributi_id.loc[linija[POLJE_ID]] if linija[POLJE_ID] in atributi_id.index else None
            x_s, y_s = xy[m]
            zapis = {"PROFIL_ID": kp["PROFIL_ID"], POLJE_HGO: hgo_id,
                     "SKUPINA": kp["SKUPINA"], "KOMPONENTA": kp["KOMPONENTA"], "STAC_M": kp["STAC_M"]}
            for f in POLJA_PRENOS:
                if f == POLJE_ID:
                    zapis[f] = linija[POLJE_ID]
                elif vir is not None and f in vir.index:
                    zapis[f] = vir[f]
            vode = niz(zapis.get(POLJE_VODE))
            zapis["VODE_RED"] = vode.count("-") + 1 if vode else None
            zapis["VODE_ID_NADR"] = vode.rsplit("-", 1)[0] if "-" in vode else None
            zapis.update({
                "A_DMR_KM2": round(kp["A_DMR"], 4), "A_GORV_KM2": round(kp["A_GORV"], 4),
                "A_KM2": round(skupaj.area / 1e6, 4) if skupaj is not None else None,
                "GORV_HGO": ";".join(sorted(g)), "GORV_IZVEN": int(kp["IZVEN"]),
                "DOLV_PROFIL": profili[kp["DOLV"]]["PROFIL_ID"] if kp["DOLV"] is not None else None,
                "ROB_DMR": (2 if g or kp["IZVEN"] else 1) if dotik else 0,
                "SMER_OBRNJENA": int(linija["OBRNJENA"]),
                "SMER_DMR_NESKLADJE": int(linija["NESKLADJE"]), "DMR_CELICA_M": celica,
                "DAT_IZRACUNA": danes,
            })
            profili_vrstice.append(dict(zapis, X=round(kp["geometry"].x, 2), Y=round(kp["geometry"].y, 2),
                                        X_SNAP=round(x_s, 2), Y_SNAP=round(y_s, 2),
                                        geometry=kp["geometry"]))
            if skupaj is not None:
                obmocja_vrstice.append(dict(zapis, OBSEG_KM=round(skupaj.length / 1000, 3),
                                            geometry=skupaj))
            if m in inkrement:
                medobm_vrstice.append({"PROFIL_ID": kp["PROFIL_ID"], "DOLV_PROFIL": zapis["DOLV_PROFIL"],
                                       "A_KM2": round(inkrement[m].area / 1e6, 4),
                                       "geometry": inkrement[m]})

        korak("Zapis {}".format(p["izhod_pot"]))
        pot = os.path.join(zacasna, "izhod.gpkg")    # na cilj se premakne šele na koncu
        sloj = (hgo_id + "_{}").format    # oznaka HGO III v imenu sloja, npr. 111_profili
        pyogrio.write_dataframe(gpd.GeoDataFrame(profili_vrstice, crs=EPSG), pot,
                                layer=sloj("profili"), driver="GPKG")
        pyogrio.write_dataframe(gpd.GeoDataFrame(obmocja_vrstice, crs=EPSG), pot,
                                layer=sloj("vodozbirna_obmocja"), driver="GPKG", promote_to_multi=True)
        pyogrio.write_dataframe(gpd.GeoDataFrame(medobm_vrstice, crs=EPSG), pot,
                                layer=sloj("medobmocja"), driver="GPKG", promote_to_multi=True)
        mreza_vrstice = []
        for e, s_od, s_do, g, k, komp in uporabljena_mreza(profili, vsa_debla, mreza, cilj):
            id_linije = None if e is None else mreza.at[e, POLJE_ID]     # None = premostitev
            vir = atributi_id.loc[id_linije] if id_linije in atributi_id.index else None
            zapis = {POLJE_HGO: hgo_id, "SKUPINA": k, "KOMPONENTA": komp,
                     "STAC_OD": None if s_od is None else round(s_od, 1),
                     "STAC_DO": None if s_do is None else round(s_do, 1)}
            for f in POLJA_PRENOS:
                zapis[f] = id_linije if f == POLJE_ID else (
                    vir[f] if vir is not None and f in vir.index else None)
            zapis.update({"SMER_OBRNJENA": None if e is None else int(mreza.at[e, "OBRNJENA"]),
                          "SMER_DMR_NESKLADJE": None if e is None else int(mreza.at[e, "NESKLADJE"]),
                          "VIR": "HIDRO5" if e is not None else
                                 "premostitev vrzeli v HIDRO5 (<= {:g} m)".format(PREMOSTITEV_M),
                          "DOLZINA_M": round(g.length, 1), "geometry": g})
            mreza_vrstice.append(zapis)
        if mreza_vrstice:
            pyogrio.write_dataframe(gpd.GeoDataFrame(mreza_vrstice, crs=EPSG), pot,
                                    layer=sloj("vodna_mreza"), driver="GPKG", promote_to_multi=True)
        if zapisi_vtokov:
            pyogrio.write_dataframe(gpd.GeoDataFrame(zapisi_vtokov, crs=EPSG), pot,
                                    layer=sloj("vtoki"), driver="GPKG")
        zapisi_tabelo(manifest, pot, sloj("manifest_dmr"))
        parametri = {
            "HIDRO5": p["hidro_pot"], "HGO": p["hgo_pot"], POLJE_HGO: hgo_id,
            "RAZDALJA_M": p["razdalja"], "DMR_CELICA_M": celica, "A_MIN_KM2": p["a_min"],
            "MAX_LISTOV": p["max_listov"], "BUFFER_M": BUFFER_M, "BURN_M": BURN_M,
            "BREACH_M": BREACH_M, "FLAT_M": FLAT_M, "PRAG_HITRO_M": PRAG_HITRO_M,
            "ZAPOLNI_M": ZAPOLNI_M, "SNAP_CELIC": SNAP_CELIC, "ZAMIK_IZLIV_CELIC": ZAMIK_IZLIV_CELIC,
            "PRAG_SMERI_M": PRAG_SMERI_M,
            "POENOSTAVITEV": "Douglas-Peucker, topološko (skupne meje med-območij poenostavljene "
                             "enkrat), toleranca {:g} m = {} x celica DMR; poligoni gorvodnih HGO "
                             "nespremenjeni".format(celica * SIMPLIFY_CELIC, SIMPLIFY_CELIC),
            "VIR_DMR": "{} ({}), {}".format(WFS_URL, WFS_SLOJ, PRENOS_URL),
            "OBDELAVA_DMR": ("povprečje točk LAZ v celici + GDAL fillnodata (do {} m)".format(ZAPOLNI_M)
                             if celica >= PRAG_HITRO_M else
                             "WBT lidar_tin_gridding (zapolnitev vrzeli)") + ", clip na HGO+buffer, "
                            "vžig HIDRO5 ({} m, all_touched), breach_depressions_least_cost, "
                            "fill_depressions (fix_flats, {} m), d8_pointer, d8_flow_accumulation, "
                            "watershed".format(BURN_M, FLAT_M),
            "VISINSKI_SISTEM": "ni privzet; za smer toka so uporabljene le relativne višine",
            "STACIONAZA": "od dolvodnega konca vodotoka (HIDRO5, skupina {}) navzgor".format(
                "/".join(POLJA_SKUPINE)),
            "SMER_TOKA": "smer digitalizacije HIDRO5; DMR določi le skupno smer mreže, neskladja "
                         "posameznih linij (prag {} m) so označena v SMER_DMR_NESKLADJE".format(PRAG_SMERI_M),
            "A_KM2": "ploščina geometrije: DMR zaledje v HGO+buffer ∪ poligoni gorvodnih HGO",
            "ROB_DMR": "0 = ne seže do roba DMR, 1 = seže do roba (zaledje morda nepopolno), "
                       "2 = seže do roba, gorvodni HGO dodani",
            "LISTI_MANJKAJO": ";".join("{:.0f}_{:.0f}".format(b.bounds[0] / 1000, b.bounds[1] / 1000)
                                       for b in manjkajoci),
            "ZACETEK": zacetek.isoformat(timespec="seconds"),
            "KONEC": dt.datetime.now().isoformat(timespec="seconds"),
            "VERZIJE": "python {}, geopandas {}, shapely {}, pyogrio {}, rasterio {}, laspy {}, {}".format(
                sys.version.split()[0], gpd.__version__, shapely.__version__, pyogrio.__version__,
                rasterio.__version__, laspy.__version__, WBT.version().splitlines()[0]),
        }
        zapisi_tabelo(pd.DataFrame({"PARAMETER": list(parametri),
                                    "VREDNOST": [str(v) for v in parametri.values()]}),
                      pot, sloj("parametri"))

        cilj_pot = p["izhod_pot"]
        try:
            os.replace(pot, cilj_pot)
        except PermissionError:
            cilj_pot = os.path.splitext(cilj_pot)[0] + "_{}.gpkg".format(
                dt.datetime.now().strftime("%Y%m%d_%H%M%S"))
            shutil.move(pot, cilj_pot)
            print("  OPOZORILO: {} je zaklenjena, rezultat zapisan v {}.".format(p["izhod_pot"], cilj_pot))
        print("\nKončano: {}".format(cilj_pot))
        print("  profilov: {}, največje zaledje: {:.2f} km2".format(
            len(profili_vrstice), max((z["A_KM2"] or 0) for z in profili_vrstice)))
        n_rob = sum(1 for z in profili_vrstice if z["ROB_DMR"] == 1)
        if n_rob:
            print("  OPOZORILO: {} zaledij seže do roba DMR brez gorvodnih HGO (ROB_DMR=1).".format(n_rob))
    finally:
        if akum is not None:
            akum.close()
        WBT.set_working_dir(p["izhod_mapa"])
        shutil.rmtree(zacasna, ignore_errors=True)
        if os.path.exists(zacasna):
            print("OPOZORILO: začasne mape ni bilo mogoče v celoti izbrisati: {}".format(zacasna))
        else:
            print("Začasni podatki (LAZ, DMR, rastri) izbrisani.")


if __name__ == "__main__":
    main()
