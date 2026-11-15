# Xenos

## Problem

A field team planning a lunar geology exercise needs to choose a specific patch of ground on Earth. Two sites can resemble the same crater in shape while exposing very different rock types. Choosing between them requires both elevation and geological evidence at the scale of the exercise.

## Solution

Xenos searches selected Earth regions for circular footprints that resemble a planetary target such as McGetchin crater on the Moon. It compares the elevation difference between the center and outer ring, then calculates how much of each footprint falls within each mapped rock type. The dashboard maps the candidates and shows separate terrain and geology scores so a team can see why a site ranked highly. Xenos withholds the combined score when coverage or resolution falls below its thresholds.

## Data

- **Lunar reference:** The LROC team's McGetchin’s location, 222 m diameter, and appearance before and after impact. [McGetchin crater overview](https://lroc.im-ldi.com/images/1501)
- **Lunar imagery:** Orbital images and their observation records provide the visual reference for the crater. [LROC Planetary Data System archive](https://pds.lroc.im-ldi.com/data/)
- **Lunar elevation:** When installed, the 3 m McGetchin digital terrain model lets Xenos measure center-to-ring relief using the same footprint geometry as on Earth. [Published McGetchin terrain dataset](https://doi.org/10.5281/zenodo.21481275)
- **Earth elevation:** Copernicus GLO-30 surface elevation tiles provide the heights used to search for candidate footprints. [Copernicus DEM on AWS](https://registry.opendata.aws/copernicus-dem/)
- **Earth geology:** Xenos intersects each footprint with Global Lithological Map (GLiM) polygons to calculate the area covered by each rock class. [GLiM database paper](https://doi.org/10.1029/2012GC004370) · [CCGM printed map overview](https://www.ccgm.org/en/product/lithological-map-of-the-world/)


## Run locally

1. Download `lunar.zip` and `earth.zip` from [Xenos GitHub Releases](https://github.com/eyedecay/xenos/releases).
2. Unzip both archives, then move the extracted `lunar/` and `earth/` folders into the project's `data/` directory. The final paths should be `data/lunar/` and `data/earth/`.

The calculation suite was run with Python 3.9.25. From the project directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
streamlit run app.py
```

## Code layout

`app.py` is the small Streamlit entrypoint. The application lives under `src/xenos/` and separates the scientific pipeline from the dashboard:

- `analysis/terrain.py`, `analysis/geology.py`, `analysis/scoring.py` and `analysis/discovery.py` contain terrain metrics, area-weighted GLiM overlays, coverage-gated PSI and regional candidate search.
- `references.py`, `regions.py`, `validation.py` and `outputs.py` handle cited target/region catalogs, planetary DTM measurements, independent-site checks, and result serialization. `terramechanics.py` processes soil and rover scenarios.
- `ui/page.py` coordinates the dashboard state and scan; `ui/reference_panel.py`, `ui/candidate_panel.py`, `ui/briefing.py`, `ui/validation.py`, `ui/exports.py`, `ui/mobility.py`, `ui/provenance.py` and `ui/layout.py` render separate page features.
- The root `xenos.py` and `terramechanics.py` files remain as compatibility imports for existing scripts and tests. New code should import from `src.xenos`.
