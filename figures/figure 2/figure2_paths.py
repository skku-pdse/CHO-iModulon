"""Shared input locations; transient calculation files stay outside the repository."""
from pathlib import Path
import os
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "figures/data"
RESULTS = Path(os.environ["IMODULON_FIGURE_WORK"])
RESULTS.mkdir(parents=True, exist_ok=True)
