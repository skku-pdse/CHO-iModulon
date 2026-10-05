"""Generate one figure group, keeping PNG panels only by default.

Calculation tables and temporary compatibility files are not left in the repository.
Use --format pdf to save vector versions instead of PNGs.
"""
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def main():
    groups = json.loads((HERE / "figure_groups.json").read_text())
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("group", choices=list(groups) + ["figure 5"])
    parser.add_argument("--output-dir", type=Path, help="Output directory for generated panels")
    parser.add_argument("--format", choices=["png", "pdf"], default="png")
    parser.add_argument("--panel", help="Run a single script from the selected group")
    args = parser.parse_args()
    if args.panel and args.panel not in groups.get(args.group, []):
        parser.error("--panel must name a script in the selected group")
    output = HERE / args.group
    with tempfile.TemporaryDirectory(prefix="imodulon_figure_") as directory:
        work = Path(directory)
        env = dict(os.environ, IMODULON_FIGURE_WORK=str(work), IMODULON_COMPARISON_WORK=str(work), PYTHONIOENCODING="utf-8")
        if args.group == "figure 5":
            comparison = ROOT / "expanded iModulon/iModulon comparison"
            scripts = [comparison / "compare_components.py", comparison / "analyze_activities.py", output / "make_figure5.py"]
        else:
            scripts = [output / (name + ".py") for name in ([args.panel] if args.panel else groups[args.group])]
        for script in scripts:
            subprocess.run([sys.executable, str(script)], env=env, check=True)
        images = list(work.rglob("*." + args.format))
        if args.group == "figure 5":
            images = [p for p in images if p.stem == "figure5_compendium_comparison"]
        if not images:
            raise RuntimeError("No figures generated")
        for image in images:
            destination = args.output_dir or output
            destination.mkdir(parents=True, exist_ok=True)
            shutil.copy2(image, destination / image.name)
    print(f"Saved {len(images)} {args.format.upper()} files in {args.output_dir or output}")


if __name__ == "__main__":
    main()
