# Add a cited planetary crater target

Xenos supports a target library at `data/targets/*.json`. Every target must be a source-cited crater profile. Use these fields:

- `name`, `body` (`Moon` or `Mars`), `diameter_m`, `depth_m`, and `source_citation` are required.
- `relevant_glim_codes` may contain evidence-based GLiM two-letter classes; without them, Xenos withholds geology score and PSI.
- `latitude_deg` and `longitude_deg` are needed to compare a planetary DTM against the Earth DEM at the crater location.

Xenos currently scores crater morphology. Do not label non-crater environments (caves, ice deposits, arid regions) as supported until the app has matching planetary and Earth measurements and a separately documented method. Local profiles remain on disk and are not bundled into a Git commit.
