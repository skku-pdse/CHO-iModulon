"""Small numerical checks that do not run the full research dataset."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "ICA/code/run_robust_ica.py"
spec = importlib.util.spec_from_file_location("robust_ica", SCRIPT)
workflow = importlib.util.module_from_spec(spec)
spec.loader.exec_module(workflow)


class ICAPipelineTests(unittest.TestCase):
    def test_sign_alignment_and_signed_rms(self):
        weights = [pd.DataFrame({0: [-3., 1., 0.]}), pd.DataFrame({0: [-6., 2., 0.]})]
        mixing = [pd.DataFrame({0: [1., -2.]}), pd.DataFrame({0: [3., -4.]})]
        m, a, stats = workflow.aggregate_clusters(weights, mixing, np.array([0, 0]))
        np.testing.assert_allclose(m[0], [4.5, -1.5, 0.])
        np.testing.assert_allclose(a.loc[0], [-np.sqrt(5), np.sqrt(10)])
        self.assertEqual(stats.loc[0, "count"], 2)

    def test_preprocessing_padding_and_missing_values(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test_tpm.tsv"
            columns = pd.MultiIndex.from_tuples([("p", "a"), ("p", "b"), ("", "")])
            frame = pd.DataFrame([[2, 3, np.nan], [4, 8, np.nan], [np.nan, np.nan, 0]],
                                 index=["gene1", "gene2", None], columns=columns)
            frame.to_csv(path, sep="\t")
            tpm, transformed = workflow.preprocess_tpm(path)
            self.assertEqual(tpm.shape, (2, 2))
            np.testing.assert_allclose(transformed.mean(axis=0), 0, atol=1e-12)
            np.testing.assert_allclose(transformed.std(axis=0), 1)
            frame.iloc[0, 0] = np.nan
            frame.to_csv(path, sep="\t")
            with self.assertRaises(ValueError):
                workflow.preprocess_tpm(path)

    def test_full_small_run_and_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            rng = np.random.default_rng(123)
            values = np.exp(np.clip(rng.laplace(size=(600, 3)) @ rng.normal(size=(3, 12))/3 + 2, -5, 8))
            columns = pd.MultiIndex.from_tuples([("project", f"sample{i}") for i in range(12)])
            pd.DataFrame(values, index=[f"gene{i}" for i in range(600)], columns=columns).to_csv(folder / "input.tsv", sep="\t")
            (folder / "seeds.txt").write_text("".join(f"Run {i}: Seed {i*10}\n" for i in range(1, 6)))
            command = [sys.executable, str(SCRIPT), "--tpm", str(folder / "input.tsv"),
                       "--seeds", str(folder / "seeds.txt"), "--dataset", "discovery",
                       "--components", "3", "--output-dir", str(folder / "output")]
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            matrix = folder / "output/discovery_ica_gene_weights.csv"
            before = matrix.read_bytes()
            self.assertEqual(pd.read_csv(matrix, index_col=0).shape, (600, 3))
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            self.assertEqual(before, matrix.read_bytes())


if __name__ == "__main__":
    unittest.main()
