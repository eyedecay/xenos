"""Place-name lookup used to help users choose an Earth search area."""
from __future__ import annotations


def search_named_places(query, limit=8):
    """Find OSM place-name matches through Photon; this is not Copernicus metadata."""
    import requests

    query = str(query).strip()
    if len(query) < 2:
        raise ValueError("Enter at least two characters to search for a place")
    response = requests.get(
        "https://photon.komoot.io/api/",
        params={"q": query, "limit": min(max(int(limit), 1), 10), "lang": "en"},
        headers={"User-Agent": "Xenos-Space-Apps/1.0 (place search)"},
        timeout=(5, 15),
    )
    response.raise_for_status()
    features = response.json().get("features", [])
    results = []
    seen = set()
    for feature in features:
        coordinates = feature.get("geometry", {}).get("coordinates", [])
        if len(coordinates) < 2:
            continue
        try:
            longitude, latitude = float(coordinates[0]), float(coordinates[1])
        except (TypeError, ValueError):
            continue
        if not (-180 <= longitude <= 180 and -90 <= latitude <= 90):
            continue
        properties = feature.get("properties", {})
        name = str(properties.get("name") or query).strip()
        qualifiers = []
        for key in ("city", "county", "state", "country"):
            value = properties.get(key)
            if value and str(value).strip().casefold() != name.casefold():
                label = str(value).strip()
                if label.casefold() not in {item.casefold() for item in qualifiers}:
                    qualifiers.append(label)
        label = " · ".join([name] + qualifiers)
        osm_id = properties.get("osm_id")
        osm_type = properties.get("osm_type", "")
        key = (osm_type, str(osm_id), round(latitude, 6), round(longitude, 6))
        if key in seen:
            continue
        seen.add(key)
        results.append({
            "id": ":".join(map(str, key)), "name": name, "label": label,
            "latitude": latitude, "longitude": longitude,
            "osm_type": osm_type, "osm_id": osm_id,
        })
    return results


def square_bounds_around_place(longitude, latitude, side_km=50):
    """Create a WGS84 square around a geocoder point, for a small DEM search area."""
    import math

    longitude, latitude, side_km = float(longitude), float(latitude), float(side_km)
    if not (-180 <= longitude <= 180 and -90 <= latitude <= 90 and side_km > 0):
        raise ValueError("Place coordinates or area width are invalid")
    half_latitude = side_km / (2 * 111.32)
    cos_latitude = max(abs(math.cos(math.radians(latitude))), 1e-4)
    half_longitude = side_km / (2 * 111.32 * cos_latitude)
    west = max(-180.0, longitude - half_longitude)
    east = min(180.0, longitude + half_longitude)
    south = max(-90.0, latitude - half_latitude)
    north = min(90.0, latitude + half_latitude)
    if west >= east or south >= north:
        raise ValueError("Could not construct an area around that place")
    return tuple(round(value, 6) for value in (west, south, east, north))
