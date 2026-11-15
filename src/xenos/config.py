"""Project paths shared by Xenos modules."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
DEM_DIR = DATA_DIR / "earth/elevation/columbia"
GEOLOGY = DATA_DIR / "earth/geology/glim/columbia_glim_clipped.gpkg"
LUNAR = DATA_DIR / "lunar"

# Explicit evidence gates used by PSI; kept in one shared configuration module.
MIN_GEOLOGY_COVERAGE = 0.90
MIN_DEM_COVERAGE = 0.95
MIN_SPECTRAL_COVERAGE = 0.90
MIN_CELLS_ACROSS = 6.0
MAX_GEOLOGY_OVERLAP = 0.01

MGETCHIN_NAC_DTM = LUNAR / "reference_dtms/NAC_DTM_McGetchin_3m_crop.tif"
