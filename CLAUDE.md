# Vir podatkov: DMR CLSS 2023–2025 (GURS)

Status: delovna verzija, 27. 10. 2026. Storitev ni uradno dokumentirana kot API GURS.

## Seznam listov – WFS

- Naslov: `https://lift.clss.si/geoserver/wfs`
- Sloj (`typeName`): `tileindex_clss_zdruzeno`
- Preizkušena zahteva (deluje, 27. 10. 2026):
  `service=WFS&version=1.1.0&request=GetFeature&typeName=tileindex_clss_zdruzeno&maxFeatures=1&outputFormat=application/json`
- Izbor po območju: parameter `bbox=xmin,ymin,xmax,ymax,EPSG:3794` (iz javne kode github.com/larchman01/strmina, NI preizkušeno).
- Koordinatni sistem geometrije v odgovoru: `urn:ogc:def:crs:EPSG::3794` (D96/TM) – potrjeno iz odgovora storitve.
- En zapis = en list 1 × 1 km; geometrija MultiPolygon; `bbox` lista na celih kilometrih.

## Atributi zapisa (potrjeno iz odgovora, list 430_135)

| Atribut | Pomen | Primer |
|---|---|---|
| `ti_name` | oznaka lista (x_y v km) | `430_135` |
| `x_ime`, `y_ime` | levi spodnji vogal v km | `430`, `135` |
| `ime_ob`, `orig_ime` | območje zajema | `Jesenice`, `04-jesenice` |
| `path_dmr` | relativna pot do DMR | `clss/raw/04-jesenice/zls/dmr/DMR_430_135.laz` |
| `size_dmr` | velikost DMR v MB | `4.232` |
| `datum_za_1` | obdobje zajema (besedilo) | `od 07.04.2025 do 09.04.2025` |
| `min_visina`, `max_visina` | višine lista – ENOTA IN VIŠINSKI SISTEM NISTA NAVEDENA | `423.828`, `662.087` |
| `tehnicno_porocilo` | HTML povezava na tehnično poročilo območja | `.../clss/tehnicna_porocila/Jesenice_porocilo.pdf` |
| `max_datum_zajema` | lahko `null` | – |

Druge poti po istem vzorcu: `path_gkot`, `path_dmp`, `path_ndmp`, `path_pas`, `path_pof`, `path_pofi`.
Atributi `link_*` vsebujejo HTML (`<a href=...>`), za prenos uporabljaj `path_*`.

## Prenos datoteke

- URL datoteke = `https://assets.flycom.si/` + `path_dmr`
- Oblika DMR: `*.laz` (točke v mreži, ne raster). Ločljivost 0,5 m po e-prostor.gov.si – na datoteki še NI preverjeno.
- Prenos datoteke iz kode še NI preizkušen.

## Pravila za kodo

- Med prenosi premor (≥ 3 s); ob HTTP 429 upoštevaj `Retry-After`. (Vir: javna koda tretje osebe, ni uradno.)
- Omejitev števila listov na zagon; večje količine naroči pri GURS: gic.gu@gov.si (vir: e-prostor.gov.si, Daljinsko zaznavanje).
- Koordinatnega in višinskega sistema LAZ ne privzemaj – preberi ga iz glave datoteke; če ga ni, to izpiši in ne pripisuj.
- Brez interpolacije in dopolnjevanja podatkov, razen če je izrecno zahtevano.
- Vsak prenos zapiši v manifest: list, URL, velikost, SHA-256, čas (UTC).
- Če atributa `path_dmr` ni, se ustavi in izpiši dejanske atribute.
- V kodo ne vpisuj gesel, ključev ali povezovalnih nizov do internih baz IzVRS.

## Odprto

- Enota višin in višinski sistem: preveriti v tehničnem poročilu območja (`tehnicno_porocilo`).
- Ali se listi različnih območij zajema prekrivajo (isti `ti_name` večkrat).
- Pogoji uporabe samodejnega prenosa: preveriti pri GURS.
