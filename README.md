# Vodozbirna območja

Python skripta `vodozbirna_obmocja.py` za izbrano hidrografsko območje HGO III vzdolž vodne mreže HIDRO5 postavi hidrološke profile na izbrani razdalji in za vsak profil izračuna vodozbirno območje. Profili se ustvarijo le, dokler je zaledje večje od minimalne površine.

DMR se sproti ustvari iz listov LAZ CLSS (GURS), vir in pravila so v `CLAUDE.md`. Po obdelavi se izbriše.

## Namestitev

Potreben je 64-bitni Python 3.10 ali novejši:

```
python -m venv .venv64
.venv64\Scripts\pip install -r requirements.txt
```

## Zagon

```
.venv64\Scripts\python vodozbirna_obmocja.py
```

Skripta ob zagonu vpraša za:
- mapo in ime sloja HIDRO5;
- mapo in ime sloja HGO;
- `HGOIII_ID`;
- razdaljo med profili (privzeto 2500 m);
- velikost celice DMR (privzeto 5 m);
- minimalno zaledje (privzeto 10 km²);
- največje število listov LAZ (privzeto 500);
- mapo izhoda.

## Izhod

Rezultat je `<izhod>/<HGOIII_ID>.gpkg` v EPSG:3794. Vsebuje sloje s predpono HGO:
- `profili`;
- `vodozbirna_obmocja`;
- `medobmocja`;
- `vodna_mreza`;
- `vtoki`;
- `manifest_dmr`;
- `parametri`.

## Podatki

- `SHP/DRSV_HGO_III.*`: hidrografska območja III. ravni.
- `SHP/HIDRO5_LIN_PV_TIPTV1_2.*`: vodna mreža HIDRO5. Zaradi velikosti (~2,1 GB) ni v repozitoriju, priložena so le metapodatki (`.shp.xml`, `.prj`, `.cpg`).
- `izhod/111.gpkg`, `112.gpkg`, `113.gpkg`: testni rezultati (DMR 5 m, razmik 2500 m).
- `komunikacija_claude.txt`: zapis razvoja skripte.
