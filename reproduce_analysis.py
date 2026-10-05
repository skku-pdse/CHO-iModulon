"""Reproduce CHO iModulon analyses without replacing packaged publication data."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
STEPS = ["validate", "membership", "genomic", "association", "projection", "comparison", "figures", "motifs", "string", "ica"]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", nargs="+", choices=STEPS, default=STEPS[:8])
    parser.add_argument("--dataset", choices=["discovery", "expanded", "both"], default="both")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output == ROOT or ROOT in output.parents and output.parts[len(ROOT.parts)] not in ("outputs",):
        parser.error("Use outputs/ or a directory outside the repository to preserve packaged files.")
    output.mkdir(parents=True, exist_ok=True)
    cohorts = ["discovery", "expanded"] if args.dataset == "both" else [args.dataset]
    tasks = []
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONIOENCODING="utf-8",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

    def run(script, *arguments):
        command = [sys.executable, str(ROOT / script), *map(str, arguments)]
        log = output / f"{len(tasks)+1:02d}_{Path(script).stem}.log"
        print(f"Running {script}", flush=True)
        started = time.perf_counter()
        with log.open("w", encoding="utf-8") as stream:
            result = subprocess.run(command, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
        tasks.append({"script": script, "arguments": list(map(str, arguments)),
                      "seconds": round(time.perf_counter()-started, 3), "exit_code": result.returncode,
                      "log": log.name})
        (output / "analysis_run.json").write_text(json.dumps(tasks, indent=2) + "\n", encoding="utf-8")
        if result.returncode:
            raise RuntimeError(f"Step failed; see {log}")

    if "validate" in args.steps:
        manifests = [ROOT / "ICA/TPM matrix/tpm_input_manifest.json",
                     ROOT / "projection/data/validation_tpm_manifest.json"]
        for manifest in manifests:
            records = json.loads(manifest.read_text())
            for record in records if isinstance(records, list) else [records]:
                path = manifest.parent / record["file"]
                if digest(path) != record["sha256"]:
                    raise ValueError(f"Input checksum mismatch: {path}")
        run("validation/validate_package.py")
        for cohort in cohorts:
            run(f"{cohort} iModulon/functional enrichment/analyze_{cohort}_string_functional_enrichment.py", "--check-inputs")
    for cohort in cohorts:
        data = ROOT / f"{cohort} iModulon/data"
        destination = output / cohort
        destination.mkdir(exist_ok=True)
        if "motifs" in args.steps:
            run(f"{cohort} iModulon/motif analysis/summarize_{cohort}_motif_analysis.py",
                "--output-dir", destination / "motifs")
        if "ica" in args.steps:
            run("ICA/code/run_robust_ica.py", "--dataset", cohort, "--tpm",
                ROOT / f"ICA/TPM matrix/{cohort}_tpm.tsv.gz", "--output-dir", destination / "ica")
        if "membership" in args.steps:
            run("ICA/code/determine_thresholds.py", data / f"{cohort}_gene_weight_matrix.csv",
                "--output-dir", destination, "--prefix", f"{cohort}_membership")
        for step, folder, analysis in [("genomic", "genomic", "chromosome_enrichment"),
                                       ("association", "conditional association", "metadata_association")]:
            if step in args.steps:
                run(f"{cohort} iModulon/{folder}/analyze_{cohort}_{analysis}.py",
                    "--output", destination / f"{cohort}_{analysis}_results.csv")
        if "string" in args.steps:
            run(f"{cohort} iModulon/functional enrichment/analyze_{cohort}_string_functional_enrichment.py",
                "--output", destination / f"{cohort}_string_functional_enrichment_results.csv")
    if "projection" in args.steps:
        run("projection/project.py", "--tpm", ROOT / "projection/data/validation_tpm.tsv.gz",
            "--samples", ROOT / "projection/data/projected_sample_metadata.csv",
            "--output", output / "projection/reproduced_projected_activities.csv")
    if "comparison" in args.steps:
        run("expanded iModulon/iModulon comparison/compare_discovery_expanded_imodulons.py",
            "--output-dir", output / "comparison")
    if "figures" in args.steps:
        for group in [f"figure {n}" for n in range(1, 6)] + ["supplementary"]:
            run("figures/generate_manuscript_figures.py", group, "--output-dir", output / "figures" / group)
    if set(args.steps) & {"membership", "genomic", "association", "projection", "comparison"}:
        run("validation/compare_reference_outputs.py", "--output-dir", output)
        print((output / "reference_comparison.csv").read_text(), flush=True)
    print(f"Completed {len(tasks)} tasks. Outputs and logs: {output}")


if __name__ == "__main__":
    main()
