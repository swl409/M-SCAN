"""Regression tests for reference checks, RNA inversion, and CLI precedence."""

import argparse
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from mscan.cli import add_common_options, build_overrides
from mscan.workflow import load_workflow


class ReferenceChecks(unittest.TestCase):
    def test_reference_checks(self):
        for kind in ("vmax", "pfba"):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as directory:
                directory = Path(directory)
                if kind == "vmax":
                    ref = pd.DataFrame({"id": ["R1", "R2"], "S1": [1., 2.], "S2": [3., 4.]})
                    key, cols = "id", ["S1", "S2"]
                else:
                    ref = pd.DataFrame({"sample": ["S1", "S2"], "R1": [1., 2.], "R2": [3., 4.],
                                        "growth_rate": [.1, .2], "total_flux_sum": [4., 6.]})
                    key, cols = "sample", ["R1", "R2"]
                reference = directory / "reference.csv"
                observed = directory / "observed.csv"

                def check(obs, expected, ref_table=ref, extra=()):
                    obs.to_csv(observed, index=False)
                    ref_table.to_csv(reference, index=False)
                    command = [sys.executable, str(ROOT / f"scripts/mscan_check_{kind}.py"),
                               "--observed", str(observed), "--reference", str(reference), *extra]
                    proc = subprocess.run(command, capture_output=True, text=True)
                    self.assertEqual(proc.returncode == 0, expected, proc.stdout + proc.stderr)
                    if expected:
                        self.assertIn("reference check passed", proc.stdout)

                cases = {"clean": (ref, True), "reordered": (ref.iloc[::-1], True),
                         "missing_row": (ref.iloc[:1], False),
                         "missing_column": (ref.drop(columns=cols[0]), False),
                         "one_sample_one_reaction": (ref.iloc[:1].drop(columns=cols[0]), False),
                         "empty": (ref.iloc[:0], False),
                         "duplicate_id": (pd.concat([ref, ref.iloc[:1]], ignore_index=True), False)}
                for name, value in [("nan", np.nan), ("inf", np.inf), ("negative_inf", -np.inf), ("text", "invalid")]:
                    table = ref.copy()
                    table[cols[0]] = table[cols[0]].astype(object)
                    table.at[0, cols[0]] = value
                    cases[name] = (table, False)
                all_nan = ref.copy()
                all_nan[[c for c in ref if c != key]] = np.nan
                cases["all_nan"] = (all_nan, False)
                wrong = ref.copy()
                wrong.at[0, cols[0]] += 1.
                cases["wrong_finite"] = (wrong, False)
                missing_id = ref.copy()
                missing_id.at[0, key] = np.nan
                cases["missing_id"] = (missing_id, False)
                extra_row = ref.iloc[:1].copy()
                extra_row[key] = "extra"
                extra_obs = pd.concat([ref, extra_row], ignore_index=True)
                extra_obs["extra_column"] = 9.
                cases["extra_observed_data"] = (extra_obs, True)
                for name, (table, expected) in cases.items():
                    with self.subTest(kind=kind, case=name):
                        check(table, expected)
                for tol in ("nan", "inf", "-1"):
                    with self.subTest(kind=kind, tolerance=tol):
                        check(ref, False, extra=["--tolerance", tol])
                with self.subTest(kind=kind, case="invalid_reference"):
                    check(ref, False, ref_table=all_nan)
                if kind == "pfba":
                    reactions = directory / "reactions.csv"
                    pd.DataFrame({"reaction_id": ["R1", "R2"]}).to_csv(reactions, index=False)
                    args = ["--reaction-list", str(reactions)]
                    check(ref, True, extra=args)
                    check(ref.drop(columns="R2"), False, extra=args)
                    pd.DataFrame({"reaction_id": ["R1", "R404"]}).to_csv(reactions, index=False)
                    check(ref, False, extra=args)
                    for column in ("growth_rate", "total_flux_sum"):
                        bad = ref.copy()
                        bad.at[0, column] = np.nan
                        check(bad, False)


class RNAAndOverrides(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("gene_matrix", ROOT / "scripts/mscan_prepare_gene_matrix.py")
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def test_rna_scales(self):
        genes = {"ENSG000001", "ENSG000002"}
        cases = [("log2", [1., 2.], [250000., 750000.], False),
                 ("log2", [0., 0.], [0., 0.], False),
                 ("log2", [0., 1.], [0., 1e6], False),
                 ("log2", [0., 1.], [1e6, 1e6], True),
                 ("linear", [1., 3.], [250000., 750000.], False),
                 ("already_normalized", [1., 3.], [1., 3.], False)]
        for scale, values, expected, replace_zero in cases:
            with self.subTest(scale=scale, values=values, replace_zero=replace_zero):
                rna = pd.DataFrame({"gene_id": sorted(genes), "sample": values})
                actual, _ = self.module.prepare_rna(rna, genes, None, scale, 1e6, replace_zero)
                np.testing.assert_allclose(actual["sample"], expected)

    def test_cli_mutation_precedence(self):
        parser = argparse.ArgumentParser()
        add_common_options(parser)
        for value in (None, "0.10", "0"):
            with self.subTest(value=value):
                args = parser.parse_args([] if value is None else ["--mutation-factor", value])
                overrides = build_overrides(args)
                workflow = load_workflow(ROOT, Path("configs/demo_ccle.yaml"), overrides)
                expected = .05 if value is None else float(value)
                self.assertEqual(workflow.mutation_factor, expected)
                self.assertEqual(workflow.config["severe_factor"], expected)
                self.assertEqual(workflow.config["lof_factor"], expected)
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / "config.yaml"
            config.write_text("mutation_factor: 0.05\nsevere_factor: 0.02\nlof_factor: 0.08\n")
            workflow = load_workflow(ROOT, config, build_overrides(parser.parse_args([])))
            self.assertEqual(workflow.config["severe_factor"], .02)
            self.assertEqual(workflow.config["lof_factor"], .08)


if __name__ == "__main__":
    unittest.main()
