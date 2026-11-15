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
from terramechanics import bekker_wong_scenario, particle_size_summary
from src.xenos.validation import validate_known_sites


class CandidateValidationTests(unittest.TestCase):
    def test_top_k_retrieval_compares_ranked_hits_with_pool_baseline(self):
        pool = [
            {"lat": 0.0, "lon": 1.0, "psi": 95.0, "terrain_score": 95.0},
            {"lat": 0.0, "lon": 0.01, "psi": 90.0, "terrain_score": 90.0},
            {"lat": 0.0, "lon": 2.0, "psi": 80.0, "terrain_score": 80.0},
        ]
        data = b"site_name,latitude,longitude,tolerance_m\nheld-out,0,0,2000\n"
        rows, summary = validate_known_sites(data, pool, top_k=(1, 2))
        self.assertEqual(rows[0]["first_matching_rank"], 2)
        self.assertEqual(summary["hit_at_1"], 0.0)
        self.assertEqual(summary["hit_at_2"], 1.0)
        self.assertAlmostEqual(summary["random_pool_hit_fraction"], 1 / 3)
        self.assertGreater(summary["enrichment_at_2"], 1.0)

    def test_withheld_score_cannot_count_as_retrieved(self):
        pool = [{"lat": 0, "lon": 0, "psi": None}, {"lat": 1, "lon": 1, "psi": 0}]
        rows, summary = validate_known_sites(b"site_name,latitude,longitude,tolerance_m\na,0,0,222", pool)
        self.assertIsNone(rows[0]["first_matching_rank"])
        self.assertEqual(summary["eligible_pool_size"], 1)
        self.assertEqual(summary["random_baseline_at_10"], 0)

    def test_sensitivity_compares_common_pool_and_flags_cutoff_tie(self):
        from src.xenos.ranking import compare_rankings
        pool = [{"terrain_score": t, "geology_score": 0, "psi": psi}
                for t, psi in [(10, 6.5), (50, 32.5), (90, None)]]
        rows, _ = compare_rankings(pool, 0.65, k=1)
        self.assertEqual(rows[-1]["Shared with baseline"], 1)
        self.assertEqual(rows[-1]["Samples compared"], 2)
        self.assertTrue(rows[2]["Tie at cutoff"])

    def test_validation_rejects_invalid_coordinates(self):
        data = b"site_name,latitude,longitude\nbad,91,0\n"
        with self.assertRaisesRegex(ValueError, "out of range"):
            validate_known_sites(data, [{"lat": 0, "lon": 0, "psi": 1}])


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
        self.assertIsNone(low_geo["psi"])
        self.assertIsNone(low_dem["psi"])
        self.assertIsNone(too_coarse["psi"])
        self.assertAlmostEqual(accepted["psi"], 80.0)
        self.assertTrue(accepted["reasons"])  # coarse DEM / missing GLiM scale cautions

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
                    "Test Area", (-0.75, 0.25, -0.25, 0.75), output_root=root / "regions",
                    glim_source=source,
                )
            clipped = gpd.read_file(result["path"] / "geology/glim_region.gpkg")
            original = gpd.read_file(source, layer="GLiM_export")
            manifest_exists = (result["path"] / "region.json").is_file()
            with rasterio.open(next((result["path"] / "dem").glob("*.tif"))) as analysis_dem:
                analysis_bounds = analysis_dem.bounds
                analysis_shape = analysis_dem.shape
            with rasterio.open(next((result["path"] / "source_dem").glob("*.tif"))) as raw_dem:
                source_shape = raw_dem.shape
        self.assertEqual(result["tiles_downloaded"], 1)
        self.assertEqual(result["glim_features"], 1)
        self.assertEqual(len(clipped), 1)
        self.assertEqual(set(clipped["xx"]), {"vb"})
        self.assertEqual(len(original), 2)
        self.assertTrue(manifest_exists)
        self.assertEqual(analysis_shape, (2, 2))
        self.assertGreater(source_shape[0], analysis_shape[0])
        self.assertGreater(source_shape[1], analysis_shape[1])
        self.assertAlmostEqual(analysis_bounds.left, -0.75)
        self.assertAlmostEqual(analysis_bounds.right, -0.25)

    def test_place_search_results_and_automatic_square_bounds(self):
        from unittest.mock import patch
        from types import SimpleNamespace
        from src.xenos.place_search import search_named_places, square_bounds_around_place

        payload = {"features": [{
            "geometry": {"coordinates": [-21.5, 63.5]},
            "properties": {"name": "Iceland", "country": "Iceland", "osm_type": "R", "osm_id": 9},
        }]}
        response = SimpleNamespace(raise_for_status=lambda: None, json=lambda: payload)
        with patch("requests.get", return_value=response) as get:
            results = search_named_places("Iceland")
        get.assert_called_once()
        self.assertEqual(results[0]["label"], "Iceland")
        bounds = square_bounds_around_place(-21.5, 63.5, 50)
        self.assertAlmostEqual((bounds[0] + bounds[2]) / 2, -21.5)
        self.assertAlmostEqual((bounds[1] + bounds[3]) / 2, 63.5)

    def test_region_catalog_merges_same_area_aliases_for_scoring(self):
        import json
        from pathlib import Path
        from tempfile import TemporaryDirectory
        from unittest.mock import patch
        from src.xenos import regions

        bounds = [-22, 63, -21, 64]
        tile = "Copernicus_DSM_COG_10_N63_00_W022_00_DEM"
        with TemporaryDirectory() as temp:
            root = Path(temp)
            for index, name in enumerate(("canada", "Iceland volcanic field", "p")):
                pack = root / f"pack{index}"
                (pack / "dem").mkdir(parents=True)
                (pack / "geology").mkdir()
                (pack / "dem" / f"{tile}.tif").write_bytes(b"same source tile")
                (pack / "geology" / "glim.gpkg").write_bytes(b"same geology")
                (pack / "region.json").write_text(json.dumps({
                    "id": f"id{index}", "name": name, "dem_globs": ["dem/*.tif"],
                    "geology_file": "geology/glim.gpkg", "bounds_wgs84": bounds,
                    "copernicus_tiles": [tile], "relevant_glim_codes": [],
                }))
            with patch.object(regions, "REGION_PACKS_DIR", root), \
                 patch.object(regions, "DEM_DIR", root / "no-columbia"), \
                 patch.object(regions, "GEOLOGY", root / "no-columbia.gpkg"):
                profiles = regions.discover_region_profiles()
        self.assertEqual(len(profiles), 1)
        self.assertEqual(profiles[0]["name"], "Iceland volcanic field")
        self.assertEqual(set(profiles[0]["duplicate_aliases"]), {"canada", "p"})


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
