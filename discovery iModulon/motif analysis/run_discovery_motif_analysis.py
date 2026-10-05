"""Run discovery MEME/TOMTOM and AME analyses."""
from pathlib import Path
import sys
LOCATION = Path(__file__).resolve().parent
sys.path.insert(0, str(LOCATION.parents[1]))
from scripts.motif_analysis_core import run_cli
if __name__ == "__main__":
    run_cli("discovery", LOCATION)
