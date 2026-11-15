import unittest

import numpy as np

from xenos import (
    annular_relief_proxy_score,
    circular_window_metrics,
    circular_zone_masks,
    coverage_aware_psi,
    geology_footprint_metrics,
    geology_relevance_score,
    lunar_dtm_profile,
    parse_crater_target_profile,
    planetary_dtm_terrain_metrics,
    riley_terrain_ruggedness,
    scan_uploaded_region,
    slope_distribution,
    stratified_candidate_indices,
)
from sentinel import contact_gradient_score, spectral_indices
from sentinel import process_sentinel2_uploads, sentinel_candidate_metrics
from terramechanics import bekker_wong_scenario, particle_size_summary


class CircularTerrainTests(unittest.TestCase):
    def test_masks_follow_metric_circle_for_anisotropic_pixels(self):
        zones = circular_zone_masks(200, pixel_x_m=20, pixel_y_m=40)
        footprint = zones["footprint"]
        self.assertEqual(footprint.shape, (7, 11))
        self.assertTrue(footprint[3, 10])  # 100 m east is on the 100 m radius.
        self.assertFalse(footprint[0, 5])  # 120 m north is outside the circle.
        self.assertLess(int(footprint.sum()), footprint.size)
        self.assertTrue(zones["center"].any())
        self.assertTrue(zones["flank"].any())
        self.assertTrue(zones["outer_annulus"].any())

    def test_relief_is_reported_as_annular_proxy_not_crater_depth(self):
        zones = circular_zone_masks(400, pixel_x_m=20, pixel_y_m=20)
        footprint = zones["footprint"]
        elevation = np.full(footprint.shape, 20.0)
        elevation[zones["center"]] = 0.0
        metrics = circular_window_metrics(
            elevation, np.ones(footprint.shape, dtype=bool), 20, 20, 400,
            slope_degrees=np.full(footprint.shape, 12.0),
        )
        self.assertAlmostEqual(metrics["annular_relief_m"], 20.0)
        self.assertAlmostEqual(metrics["annular_relief_ratio"], 0.05)
        self.assertNotIn("crater_depth_m", metrics)
        self.assertAlmostEqual(metrics["flank_mean_slope_deg"], 12.0)
        self.assertAlmostEqual(metrics["dem_coverage_fraction"], 1.0)

    def test_proxy_score_has_transparent_linear_endpoints(self):
        self.assertEqual(annular_relief_proxy_score(0.19, 0.19), 100.0)
        self.assertEqual(annular_relief_proxy_score(0.38, 0.19), 0.0)
        self.assertEqual(annular_relief_proxy_score(-0.1, 0.19), 0.0)

    def test_planetary_dtm_is_average_resampled_to_earth_spacing(self):
        import rasterio
        from rasterio.io import MemoryFile
        from rasterio.transform import from_origin

        n, resolution = 501, 3.0
        yy, xx = np.mgrid[:n, :n]
        radius = np.hypot((xx - n / 2) * resolution, (yy - n / 2) * resolution)
        elevation = np.where(radius < 30, 0.0, 40.0).astype("float32")
        lunar_crs = "+proj=aeqd +lat_0=1.3536 +lon_0=67.1765 +a=1737400 +b=1737400 +units=m +no_defs +type=crs"
        with MemoryFile() as memory:
            with memory.open(driver="GTiff", height=n, width=n, count=1,
                             dtype="float32", crs=lunar_crs,
                             transform=from_origin(-n * resolution / 2,
                                                   n * resolution / 2,
                                                   resolution, resolution)) as dataset:
                dataset.write(elevation, 1)
            measured = planetary_dtm_terrain_metrics(
                memory.read(), 240, comparison_pixel_size_m=30)
        self.assertAlmostEqual(measured["native_pixel_size_m"], 3.0, places=3)
        self.assertAlmostEqual(measured["analysis_pixel_size_m"], 30.0, places=3)
        self.assertAlmostEqual(measured["cells_across"], 8.0, places=3)
        self.assertTrue(np.isfinite(measured["annular_relief_ratio"]))


class GeologicalScoreTests(unittest.TestCase):
    def test_unrelated_balanced_classes_do_not_score_as_target_mix(self):
        result = geology_relevance_score({"ss": 0.5, "pa": 0.5})
        self.assertEqual(result["score"], 0.0)

    def test_basalt_relevance_and_contrast_are_both_required_for_maximum(self):
        mixed = geology_relevance_score({"vb": 0.5, "ss": 0.5})
        pure_basalt = geology_relevance_score({"vb": 1.0})
        self.assertEqual(mixed["score"], 75.0)
        self.assertEqual(pure_basalt["score"], 50.0)
        self.assertFalse(mixed["highland_anorthosite_supported"])

    def test_custom_target_classes_drive_geology_relevance_and_can_be_withheld(self):
        custom = geology_relevance_score({"ss": 0.5, "pa": 0.5}, relevant_codes=("ss",))
        self.assertEqual(custom["score"], 75.0)
        self.assertEqual(custom["target_share"], 0.5)
        self.assertIsNone(geology_relevance_score({"ss": 1.0}, relevant_codes=())["score"])

    def test_area_weighted_overlay_and_coverage_use_polygon_union(self):
        import geopandas as gpd
        from shapely.geometry import Point, box

        circle = Point(0, 0).buffer(100, quad_segs=32)
        # Split a circular footprint into equal left/right mapped units.
        frame = gpd.GeoDataFrame(
            {"Litho": ["vb____", "ss____"], "xx": ["vb", "ss"]},
            geometry=[box(-100, -100, 0, 100), box(0, -100, 100, 100)],
            crs="EPSG:5070",
        )
        result = geology_footprint_metrics(circle, frame)
        self.assertAlmostEqual(result["coverage_fraction"], 1.0, places=3)
        self.assertAlmostEqual(result["class_shares"]["vb"], 0.5, places=3)
        self.assertAlmostEqual(result["score"], 75.0, places=2)
        self.assertAlmostEqual(result["overlap_fraction"], 0.0, places=5)

    def test_partial_and_overlapping_polygons_are_reported(self):
        import geopandas as gpd
        from shapely.geometry import Point, box

        circle = Point(0, 0).buffer(100, quad_segs=32)
        half = gpd.GeoDataFrame(
            {"Litho": ["vb____"], "xx": ["vb"]},
            geometry=[box(0, -100, 100, 100)], crs="EPSG:5070",
        )
        partial = geology_footprint_metrics(circle, half)
        self.assertAlmostEqual(partial["coverage_fraction"], 0.5, places=2)
        self.assertIsNone(coverage_aware_psi(80, partial["score"], partial["coverage_fraction"], 1.0, 10)["psi"])

        overlapping = gpd.GeoDataFrame(
            {"Litho": ["vb____", "ss____"], "xx": ["vb", "ss"]},
            geometry=[circle, circle], crs="EPSG:5070",
        )
        duplicate = geology_footprint_metrics(circle, overlapping)
        self.assertAlmostEqual(duplicate["coverage_fraction"], 1.0, places=4)
        self.assertAlmostEqual(duplicate["overlap_fraction"], 1.0, places=4)
        self.assertIsNone(coverage_aware_psi(80, duplicate["score"], 1.0, 1.0, 10,
                                             duplicate["overlap_fraction"])["psi"])

    def test_psi_is_withheld_below_coverage_and_resolution_gates(self):
        low_geo = coverage_aware_psi(80, 60, 0.89, 1.0, 10)
        low_dem = coverage_aware_psi(80, 60, 1.0, 0.94, 10)
        too_coarse = coverage_aware_psi(80, 60, 1.0, 1.0, 5.9)
        accepted = coverage_aware_psi(80, 60, 0.95, 1.0, 10)
        with_spectral = coverage_aware_psi(80, 60, 0.95, 1.0, 10,
                                           sentinel_score=70, sentinel_coverage=0.95,
                                           sentinel_available=True)
        low_spectral = coverage_aware_psi(80, 60, 0.95, 1.0, 10,
                                          sentinel_score=70, sentinel_coverage=0.5,
                                          sentinel_available=True)
        self.assertIsNone(low_geo["psi"])
        self.assertIsNone(low_dem["psi"])
        self.assertIsNone(too_coarse["psi"])
        self.assertAlmostEqual(accepted["psi"], 80.0)
        self.assertTrue(accepted["reasons"])  # coarse DEM / missing GLiM scale cautions
        self.assertAlmostEqual(with_spectral["psi"], 78.0)
        self.assertIsNone(low_spectral["psi"])

    def test_psi_normalizes_geology_to_same_scale_as_terrain(self):
        result = coverage_aware_psi(100, 75, 1.0, 1.0, 12)
        self.assertAlmostEqual(result["psi"], 100.0)


class DiscoveryFrameworkTests(unittest.TestCase):
    def test_stratified_sampling_spans_lower_and_higher_terrain_scores(self):
        scores = np.arange(1000, dtype=float)
        selected = stratified_candidate_indices(scores, 10)
        self.assertEqual(len(selected), 10)
        self.assertEqual(scores[selected[0]], 0)
        self.assertEqual(scores[selected[-1]], 999)

    def test_custom_crater_profile_requires_citable_valid_metrics(self):
        profile = parse_crater_target_profile({
            "name": "Example crater", "body": "Mars", "diameter_m": 1000,
            "depth_m": 100, "source_citation": "DOI: example",
        })
        self.assertEqual(profile["body"], "Mars")
        self.assertAlmostEqual(profile["depth_diameter_ratio"], 0.1)
        with self.assertRaises(ValueError):
            parse_crater_target_profile({"name": "uncited", "body": "Mars",
                                         "diameter_m": 100, "depth_m": 200,
                                         "source_citation": ""})

    def test_uploaded_region_processes_real_geotiff_and_polygon_bytes(self):
        import geopandas as gpd
        import rasterio
        from shapely.geometry import box
        from rasterio.transform import from_origin
        from tempfile import TemporaryDirectory
        from pathlib import Path

        with TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            dem_path, geology_path = root / "small_dem.tif", root / "geology.gpkg"
            n, res, lon0, lat0 = 100, 0.0001, -120.01, 45.005
            yy, xx = np.mgrid[:n, :n]
            radius_m = np.hypot((xx - n/2) * 7.86, (yy - n/2) * 11.13)
            elevation = 100 - 20 * np.clip(1 - radius_m / 140, 0, 1)
            with rasterio.open(dem_path, "w", driver="GTiff", height=n, width=n,
                               count=1, dtype="float32", crs="EPSG:4326",
                               transform=from_origin(lon0, lat0, res, res)) as dataset:
                dataset.write(elevation.astype("float32"), 1)
            polygons = gpd.GeoDataFrame(
                {"Litho": ["vb____"], "xx": ["vb"]},
                geometry=[box(lon0, lat0 - n*res, lon0 + n*res, lat0)],
                crs="EPSG:4326",
            )
            polygons.to_file(geology_path, layer="units", driver="GPKG")
            candidates, summary = scan_uploaded_region(
                200, 0.19, [(dem_path.name, dem_path.read_bytes())],
                geology_path.name, geology_path.read_bytes(), region_name="Test region",
                top_per_tile=10, max_candidates=5,
            )
        self.assertGreater(summary["candidate_pool_size"], 0)
        self.assertEqual(summary["region_name"], "Test region")
        self.assertGreater(summary["tiles_read"], 0)
        self.assertTrue(candidates)
        self.assertIn("geology", candidates[0])


class RegionBuilderTests(unittest.TestCase):
    def test_copernicus_tile_names_and_small_bounds(self):
        from src.xenos.region_builder import copernicus_tile_name, copernicus_tiles_for_bounds

        self.assertEqual(copernicus_tile_name(45, -121),
                         "Copernicus_DSM_COG_10_N45_00_W121_00_DEM")
        self.assertEqual(copernicus_tile_name(-1, 0),
                         "Copernicus_DSM_COG_10_S01_00_E000_00_DEM")
        tiles = copernicus_tiles_for_bounds((-121, 45, -120, 46))
        self.assertEqual(len(tiles), 1)
        self.assertEqual(tiles[0], "Copernicus_DSM_COG_10_N45_00_W121_00_DEM")
        self.assertEqual(len(copernicus_tiles_for_bounds((-121, 45, -119, 47))), 4)
        self.assertEqual(len(copernicus_tiles_for_bounds((-121, 45, -118, 47))), 6)

    def test_build_pack_clips_local_glim_and_preserves_source(self):
        import geopandas as gpd
        import rasterio
        from rasterio.transform import from_origin
        from shapely.geometry import box
        from tempfile import TemporaryDirectory
        from pathlib import Path
        from unittest.mock import patch
        from src.xenos import region_builder

        with TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "global.gpkg"
            global_frame = gpd.GeoDataFrame(
                {"Litho": ["vb____", "ss____"], "xx": ["vb", "ss"]},
                geometry=[box(-2, -2, 2, 2), box(2, 2, 3, 3)], crs="EPSG:4326",
            )
            global_frame.to_file(source, layer="GLiM_export", driver="GPKG")

            def mock_download(tile, destination):
                destination = Path(destination)
                destination.parent.mkdir(parents=True, exist_ok=True)
                with rasterio.open(destination, "w", driver="GTiff", height=4, width=4,
                                   count=1, dtype="float32", crs="EPSG:4326",
                                   transform=from_origin(-1, 1, 0.25, 0.25)) as ds:
                    ds.write(np.ones((4, 4), dtype="float32"), 1)
                return True

            with patch.object(region_builder, "_download_tile", side_effect=mock_download):
                result = region_builder.build_copernicus_region_pack(
                    "Test Area", (-1, 0, 0, 1), output_root=root / "regions",
                    glim_source=source,
                )
            clipped = gpd.read_file(result["path"] / "geology/glim_region.gpkg")
            original = gpd.read_file(source, layer="GLiM_export")
            manifest_exists = (result["path"] / "region.json").is_file()
        self.assertEqual(result["tiles_downloaded"], 1)
        self.assertEqual(result["glim_features"], 1)
        self.assertEqual(len(clipped), 1)
        self.assertEqual(set(clipped["xx"]), {"vb"})
        self.assertEqual(len(original), 2)
        self.assertTrue(manifest_exists)


class PRDModuleTests(unittest.TestCase):
    def test_riley_tri_is_sum_of_eight_neighbor_differences(self):
        elevation = np.zeros((5, 5), dtype=float)
        elevation[2, 2] = 10
        tri = riley_terrain_ruggedness(elevation)
        self.assertEqual(tri[2, 2], 80.0)
        self.assertTrue(np.isnan(tri[0, 0]))

    def test_slope_distribution_retains_histogram_counts(self):
        values = np.array([[2.0, 7.0], [24.0, 41.0]])
        result = slope_distribution(values, np.ones_like(values, dtype=bool))
        self.assertEqual(result["sample_count"], 4)
        self.assertEqual(sum(result["counts"]), 4)
        self.assertEqual(result["bin_width_deg"], 5)

    def test_spectral_contact_index_uses_explicit_band_inputs(self):
        bands = {"B02": np.full((5, 5), 0.1), "B04": np.full((5, 5), 0.2),
                 "B11": np.tile(np.linspace(0.1, 0.3, 5), (5, 1)),
                 "B12": np.full((5, 5), 0.2)}
        iron, swir = spectral_indices(bands)
        self.assertTrue(np.allclose(iron, 1 / 3))
        result = contact_gradient_score(iron, swir, np.ones((5, 5), dtype=bool), 0.01)
        self.assertGreater(result["score"], 0)
        self.assertEqual(result["sample_count"], 25)

    def test_sentinel_geotiff_bands_align_and_measure_candidate_coverage(self):
        import rasterio
        from pyproj import Transformer
        from rasterio.io import MemoryFile
        from rasterio.transform import from_origin

        class Upload:
            def __init__(self, name, content):
                self.name = name
                self.content = content

            def getvalue(self):
                return self.content

        crs = "EPSG:32610"
        transform = from_origin(500000, 5100000, 10, 10)
        uploads = []
        for band, value in [("B02", 1000), ("B03", 1500), ("B04", 2000),
                            ("B11", 3000), ("B12", 2500)]:
            array = np.full((20, 20), value, dtype="uint16")
            if band == "B11":
                array[:, 10:] += 500
            with MemoryFile() as mem:
                with mem.open(driver="GTiff", height=20, width=20, count=1,
                              dtype="uint16", crs=crs, transform=transform,
                              nodata=0) as ds:
                    ds.write(array, 1)
                uploads.append(Upload(f"scene_{band}_10m.tif", mem.read()))
        product = process_sentinel2_uploads(uploads, reflectance_scale=0.0001)
        center_lon, center_lat = Transformer.from_crs(crs, "EPSG:4326", always_xy=True).transform(500100, 5099900)
        result = sentinel_candidate_metrics(product, center_lon, center_lat, 100, 0.005)
        self.assertEqual(product["shape"], (20, 20))
        self.assertAlmostEqual(float(product["bands"]["B02"][0, 0]), 0.1, places=5)
        self.assertGreaterEqual(result["coverage_fraction"], 0.9)
        self.assertIsNotNone(result["score"])

    def test_psd_quantiles_and_bekker_curve_use_supplied_values(self):
        import pandas as pd

        psd = particle_size_summary(pd.DataFrame({
            "particle_size_mm": [2.0, 1.0], "mass_fraction": [0.5, 0.5],
        }))
        self.assertEqual(psd["quantiles"]["D10_mm"], 1.0)
        self.assertEqual(psd["quantiles"]["D90_mm"], 2.0)
        simulation = bekker_wong_scenario(
            kc=1, kphi=1_000_000, exponent=1, wheel_width_m=0.2,
            contact_length_m=0.3, wheel_load_n=60, wheel_radius_m=0.1,
            cohesion_pa=1000, internal_friction_deg=30,
            shear_modulus_m=0.02, max_shear_displacement_m=0.1,
        )
        self.assertAlmostEqual(simulation["normal_pressure_pa"], 1000)
        self.assertGreater(simulation["sinkage_m"], 0)
        self.assertEqual(len(simulation["pressure_pa"]), 100)
        self.assertGreaterEqual(simulation["rover_mobility_risk_index"], 0)
        self.assertLessEqual(simulation["rover_mobility_risk_index"], 100)

    def test_lunar_profile_requires_and_samples_lunar_georeferencing(self):
        import rasterio
        from pyproj import CRS
        from rasterio.io import MemoryFile
        from rasterio.transform import from_origin

        moon = CRS.from_proj4("+proj=longlat +a=1737400 +b=1737400 +no_defs +type=crs")
        transform = from_origin(67.1565, 1.3736, 0.001, 0.001)
        with MemoryFile() as memory:
            with memory.open(driver="GTiff", width=40, height=40, count=1,
                             dtype="float32", crs=moon, transform=transform,
                             nodata=-9999) as dataset:
                dataset.write(np.tile(np.arange(40, dtype="float32"), (40, 1)), 1)
            result = lunar_dtm_profile(memory.read(), half_width_m=50, samples=11)
        self.assertEqual(len(result["distance_m"]), 11)
        self.assertEqual(result["valid_fraction"], 1.0)
        self.assertIn("1737400", result["crs"])


if __name__ == "__main__":
    unittest.main()
