#!/usr/bin/env python3
"""Fetch only the published McGetchin LROC NAC DTM into the local ignored data folder."""
from __future__ import annotations
import hashlib
import io
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile

URL = "https://zenodo.org/records/21481275/files/DTM_data.zip?download=1"
ARCHIVE_MD5 = "d8c6172a4940ec054f0730eb91426619"
MEMBER = "NAC_DTM_NEW200CRAT3_E014N0671_crop.tif"
DESTINATION = Path(__file__).resolve().parents[1] / "data/lunar/reference_dtms/NAC_DTM_McGetchin_3m_crop.tif"


def main():
    if DESTINATION.exists():
        raise SystemExit(f"Refusing to overwrite existing data: {DESTINATION}")
    with urlopen(URL, timeout=60) as response:
        archive_bytes = response.read()
    digest = hashlib.md5(archive_bytes).hexdigest()
    if digest != ARCHIVE_MD5:
        raise SystemExit(f"Archive checksum mismatch: expected {ARCHIVE_MD5}, got {digest}")
    with ZipFile(io.BytesIO(archive_bytes)) as archive:
        if MEMBER not in archive.namelist():
            raise SystemExit(f"Expected DTM file {MEMBER!r} is missing from source archive")
        dtm_bytes = archive.read(MEMBER)
    DESTINATION.parent.mkdir(parents=True, exist_ok=True)
    DESTINATION.write_bytes(dtm_bytes)
    print(f"Saved {len(dtm_bytes):,} bytes to {DESTINATION}")
    print("Source: Robinson et al. (2026), Zenodo DOI 10.5281/zenodo.21481275")


if __name__ == "__main__":
    main()
