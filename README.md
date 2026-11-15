# Xenos

Xenos is a local-first Streamlit app for discovering Earth terrain candidates that resemble cited lunar or Martian crater targets. It processes the local Copernicus DEM and GLiM GeoPackage, LROC observations, and—when installed—the published McGetchin LROC NAC DTM. Other lunar/Martian target profiles and separately mapped Earth regions can be added locally. Sentinel-2 and terramechanics remain optional inputs; no soil or spectral measurements are fabricated. Original datasets remain unchanged.

## Run locally

Python 3.10 or newer is required. From the project directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
streamlit run app.py
```

The **Run regional analog scan** uses the bundled Columbia DEM and GLiM data, user-supplied georeferenced DEM GeoTIFFs and GLiM-compatible polygons, or a newly prepared area from public Copernicus GLO-30 AWS tiles plus the local global GLiM database. The on-demand region builder is limited to four 1° tiles per request, records unavailable tiles, and writes a new ignored region pack without editing source files. The scan constructs a region-centered equal-area projection, calculates Riley TRI, circular annular relief, flank slope distributions and area-weighted lithology. User uploads are staged in a temporary directory and never replace project data. Annular relief is an outer-ring minus inner-disk elevation contrast, not crater depth or d/D; Earth flank slope remains descriptive only.

For each retained center, the application creates a circular footprint in a region-centered Lambert azimuthal equal-area projection and intersects it with the geology polygons. Coverage is the union area of mapped polygons divided by total footprint area; overlapping polygon area is measured separately. Area-weighted class shares are calculated from the intersected `Litho` records. GLiM class `vb` is basic volcanic, usually basalt-type. Since the GLiM class scheme has no anorthosite class, other classes are counted only as generic mapped contrast, not as a highland analog.

The terrain similarity score is `max(0, 100 × (1 − |Earth_annular_relief_ratio − target_comparison_ratio| / target_comparison_ratio))`. When the McGetchin NAC DTM is installed, Xenos measures the same circular inner-disk/outer-annulus relief ratio on both bodies. Otherwise it falls back to the published crater depth/diameter ratio as an unlike-measure screening proxy, labeled as such. The McGetchin paper reports a median interior wall slope measured on a 12 m baseline; the 30 m Earth DEM flank-band mean is not methodologically comparable and is excluded from PSI.


Candidate sampling is stratified across the full terrain-score range per tile, and GLiM composition is computed for that pool before PSI ranking. This reduces terrain-only preselection bias. The default McGetchin geology score is `0.5 × (100 × vb_share) + 0.5 × (200 × min(vb_share, non_vb_share))`, with a maximum of 75. Uploaded cited Moon/Mars crater profiles may provide target-relevant GLiM class codes; if they do not, geology composition remains visible but geology score and PSI are withheld. The 0–75 geology score is normalized to 0–100 for combination. Without Sentinel-2, PSI is `0.65 × terrain_proxy + 0.35 × normalized_geology`. With an eligible Sentinel-2 chip, selected-candidate PSI is `0.55 × terrain + 0.25 × normalized_geology + 0.20 × spectral_heterogeneity`. PSI is withheld if mapped GLiM coverage is below 90%, valid DEM coverage below 95%, footprint width below 6 cells, Sentinel coverage (when used) below 90%, or class areas are ambiguous from >1% polygon overlap. A score that passes these gates still carries warnings: the DEM spans about 7 cells across the 222 m target footprint, and GLiM source-map scale is absent from the clipped GeoPackage. These are transparent screening rules, not a calibrated uncertainty model.

## Adding target profiles and Earth regions

Target profiles stored as cited JSON files under `data/targets/` appear in the target selector. The same menu supports a one-off JSON upload. For the default McGetchin target, Xenos automatically uses `data/lunar/reference_dtms/NAC_DTM_McGetchin_3m_crop.tif` if present; fetch this small, published NAC stereo DTM with `python scripts/fetch_mgetchin_dtm.py`. The fetch script verifies the source archive checksum and refuses to overwrite an existing file. The DTM data are from Robinson et al., “A New 222-m Diameter Lunar Crater,” supplemental data, [Zenodo DOI 10.5281/zenodo.21481275](https://doi.org/10.5281/zenodo.21481275). The source archive includes a 3 m/px NAC DTM cropped around McGetchin; verify reuse/license terms from the source record before redistribution. Xenos computes the same circular relief metric and only uses it when coverage and resolution gates pass.

Earth regions can be scanned together from the installed-region catalog. Add separately sourced data packs with `data/earth/regions/<id>/region.json`; see [the region-pack instructions](docs/region-packs/README.md). Only complete packs with real DEM and geological polygon data are shown. At present, Columbia River Basalts is the only installed region in this workspace; other environments are not simulated or represented by invented measurements. Local cited crater profiles can be added as described in [target profiles](docs/target-profiles.md).

A custom crater target is a source-cited JSON object with `name`, `body` (`Moon` or `Mars`), `diameter_m`, `depth_m`, and `source_citation`. Add `relevant_glim_codes` as a list of evidence-based GLiM two-letter classes to enable target geology scoring. Without relevant codes, composition is still shown, while geology score and PSI are withheld. This custom profile mode supports crater/depression analog screening only; it does not claim to analyze caves, ice, aridity, temperature, or other environments without corresponding planetary and Earth datasets and scoring methods.

For a new Earth region, upload one or more georeferenced DEM GeoTIFFs and a GLiM-compatible GeoPackage or GeoJSON containing `Litho` and `xx` fields. The target-centered local equal-area projection is calculated from the DEM region. Check that the geology layer actually covers the DEM extent and that its class codes are meaningful for the selected target. Source files are copied to temporary storage for processing and are not modified.

## PRD feature audit

| PRD feature | Implementation | Data needed / interpretation |
| --- | --- | --- |
| Target profiles | Bundled McGetchin, local cited crater profile, or one-off cited upload for Moon/Mars | Custom JSON provides name, body, diameter_m, depth_m, source_citation; optional relevant_glim_codes control target geology score |
| Regional and custom Earth scan | Bundled Columbia, uploaded DEM + geology, or on-demand Copernicus GLO-30 plus a spatial clip of local global GLiM | On-demand area is limited to four 1° tiles; GLO-30 public coverage can be incomplete |
| TRI and slope distribution | Riley TRI and 5° flank-band histogram per candidate | Lunar 24–41° band is context only; Earth mean slope is not the same measurement |
| Scale-matched candidates | Circular DEM footprint, center/flank/outer annulus, spacing suppression | Candidates are not confirmed craters; GLO-30 resolution is coarse at 222 m |
| Geology composition | Polygon-union coverage and area-weighted `Litho`; default basalt (`vb`) relevance or profile-specified relevant GLiM classes plus mapped contrast | GLiM has no explicit anorthosite class and local GeoPackage lacks source-map scale |
| Sentinel-2 L2A | B02/B03/B04/B11/B12 chips align to B02; B04/B03/B02 true color, iron proxy, SWIR contrast, contact-gradient score | Upload georeferenced reflectance GeoTIFF chips for the selected candidate; no Sentinel data is bundled. Optional SCL masks clouds/shadows |
| Planetary elevation comparison | Samples a mapped lunar DTM profile and measures the same circular annular-relief ratio as on Earth | McGetchin's published 3 m/px LROC NAC DTM is automatically used when locally installed; custom Moon/Mars DTM uploads require cited target coordinates and a matching planetary CRS |
| Terramechanics | PSD CSV quantiles, Bekker sinkage curve, simplified shear traction/efficiency and risk heuristic | Requires uploaded PSD plus all 11 user-entered scenario parameters; no soil defaults or measurements are fabricated |
| Viewer and analytics | Planetary imagery/DTM reference, Earth DEM/Sentinel, PSI, geology and optional mobility scenario | DTM is auto-used when installed or can be uploaded; Sentinel and soil panels require user inputs |
| Geographic candidate exploration | Interactive regional coordinate map, candidate footprint outlines, class/coverage hover details, PSI-ranked candidate selector and advanced multi-candidate comparison | Candidate outlines are geographic approximations; dashed box is the tile bounding extent, not an authoritative region boundary; no basemap is downloaded |
| Candidate-pool and ranking diagnostics | Terrain-score-stratified sampling; cross-region results with per-region representation; terrain/geology ablations and variable weights | Known-site validation accepts a user-supplied held-out CSV; no independent validation set is bundled and no accuracy claim is made |
| Reproducible exports | Candidate CSV, GeoJSON point features, and JSON run manifest with target citation, source signatures, thresholds, and scan summary | Source rasters and pixel-level arrays are not embedded; local input files are required to reproduce |
| Multi-region catalog | Selects and scans multiple independently mapped region packs, then merges candidates | Additional packs can be prepared from public Copernicus DEM and local GLiM or supplied with other cited products |
| Target profile catalog | McGetchin plus cited Moon/Mars crater profiles under `data/targets/`, or a one-off upload | Scores crater profiles only; other environments need relevant target and Earth measurements |
| Offline regional scan | Local DEM and geology datasets are used directly | Installed-pack scans require no remote data request; preparing a Copernicus area needs internet |

Sentinel chips must be projected and co-registered, no larger than 2048 pixels per side, and named with their band (`B02`, `B03`, `B04`, `B11`, `B12`; optional `SCL`). Upload reflectance values in 0–1 with scale 1, or apply the product's documented scale/offset before upload. The displayed iron and SWIR indices are exploratory proxies; they do not identify mineral species. The uploaded chip modifies PSI for the selected candidate only; the regional candidate list remains ordered by terrain + geology until spectral coverage is available for every candidate.

The planetary DTM processor validates the raster's lunar or Mars body radius and reference coordinates before deriving terrain. For McGetchin, the published 3 m/px NAC DTM is preferred over a coarser global grid and compared with Earth using the same annular geometry. Without that local DTM file, Xenos clearly falls back to the cited depth/diameter proxy. See [the data-source record](https://doi.org/10.5281/zenodo.21481275).

Terramechanics inputs use a Bekker pressure-sinkage curve and a simplified uniform-contact shear approximation. The risk index is a stated heuristic (`70% sinkage/radius + 30% (1 − drawbar efficiency)`) and is not a rover certification model. PSD alone does not determine Bekker coefficients.

Run the calculation tests with:

```bash
python -m unittest discover -s tests -v
```

## Source data and limits

The scan processes local files without downloading source datasets during the session. Running the optional McGetchin DTM setup script downloads only the published 4.6 MB supplemental archive, verifies its MD5 checksum, and extracts the DTM to a new ignored local file. Custom region uploads are processed from temporary copies and are removed when the scan returns. The interface reports processed DEM and GLiM counts and tile extent. Provenance fields absent from supplied manifests or clipped products are explicitly marked unavailable; product identifiers are read from LROC PDS labels where present. The supplied Copernicus and GLiM bundles do not include a complete release/license manifest. The Zenodo landing page supplies the DTM archive and citation; confirm terms there before redistributing its contents.

- `data/lunar/crater.json` and `data/lunar/lunar_reference.json` supply the crater reference values, coordinates, verification note and PDS product names.
- `data/lunar/pds/` contains the LROC EDR images. The app reads their PDS labels and sampled brightness statistics. These are digital numbers, not elevation. `Pre-impact.png` and `Post-impact.png` are used for the lunar viewer.
- `data/earth/elevation/columbia/` contains Copernicus GLO-30 GeoTIFFs. A 222 m footprint spans about 7 pixels across its coarser grid direction. Annular relief and flank slope are DEM-derived proxies; local windows are candidate terrain patches, not confirmed impact craters.
- `data/earth/geology/glim/columbia_glim_clipped.gpkg` contains 3,314 polygons in WGS 84 with `Litho` and `xx` fields. The GeoPackage lacks the GLiM source-map scale attribute; full geometric coverage does not imply fine-scale map accuracy.
- McGetchin's published d/D uses median rim elevation minus maximum crater depth. Its 24° interior slope is a median wall slope measured on a 12 m baseline. The Earth proxies do not reproduce those measurement methods.
- The bundled crater reference JSON retains its “To be verified” status. Its measurements and highland-mare setting were cross-checked against the published study; the source JSON was not modified. The downloaded NAC DTM is stored separately under `data/lunar/reference_dtms/` and is ignored by Git. GLiM cannot identify anorthosite directly.
- The earlier Bekker plot used assumed soil constants. It is now labeled illustrative and not shown as a computed measurement. The supplied data do not contain calibrated soil mechanics measurements.

Earth gravity, atmosphere, weathering, vegetation and DEM resolution differ from lunar conditions. PSI is an explainable screening index, not a probability, validated analog determination or rover engineering prediction.

## Code layout

`app.py` is the small Streamlit entrypoint. The application lives under `src/xenos/` and separates the scientific pipeline from the dashboard:

- `analysis/terrain.py`, `analysis/geology.py`, `analysis/scoring.py` and `analysis/discovery.py` contain terrain metrics, area-weighted GLiM overlays, coverage-gated PSI and regional candidate search.
- `references.py`, `regions.py`, `validation.py` and `outputs.py` handle cited target/region catalogs, planetary DTM measurements, independent-site checks, and result serialization. `sentinel.py` and `terramechanics.py` process optional uploaded chips and explicitly entered scenarios.
- `ui/page.py` coordinates the dashboard state and scan; `ui/reference_panel.py`, `ui/candidate_panel.py`, `ui/briefing.py`, `ui/validation.py`, `ui/exports.py`, `ui/mobility.py`, `ui/provenance.py` and `ui/layout.py` render separate page features.
- The root `xenos.py`, `sentinel.py` and `terramechanics.py` files remain as compatibility imports for existing scripts and tests. New code should import from `src.xenos`.
