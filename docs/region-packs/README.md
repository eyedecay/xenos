# Add a local Earth region pack

Xenos discovers regional datasets under `data/earth/regions/<region-id>/region.json`. Each region directory must contain the DEM GeoTIFFs and one compatible GLiM GeoPackage/GeoJSON; a region is offered in the app only when every referenced file exists.

Manifest example:

```json
{
  "id": "iceland-basalt",
  "name": "Iceland basalt fields",
  "dem_globs": ["dem/*.tif"],
  "geology_file": "geology/units.gpkg",
  "citation": "Dataset citation and version URL"
}
```

Paths are relative to the region directory and cannot escape it. DEMs must have a CRS. The geology layer must provide `Litho` and `xx` fields and cover the DEM footprint. Keep the original downloaded products untouched; put derived clips in this region directory and cite the source, resolution, license, and preprocessing in the manifest. The application does not treat a region as scientifically characterized just because files exist: each selected scan reads and intersects the actual rasters and polygons.
