#!/usr/bin/env python3
"""Compare M-SCAN pFBA output with a stored reference subset."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from mscan_reference_utils import finite_numeric, require_coverage, validate_ids, validate_tolerance


META_COLS = {"sample", "total_flux_sum", "growth_rate", "status", "error_msg"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--reaction-list", type=Path, default=None)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()
    validate_tolerance(args.tolerance)

    observed = pd.read_csv(args.observed, low_memory=False)
    reference = pd.read_csv(args.reference, low_memory=False)

    validate_ids(observed, "sample", "Observed pFBA table")
    validate_ids(reference, "sample", "Reference pFBA table")
    samples = sorted(reference["sample"])
    require_coverage(samples, observed["sample"], "reference samples in observed pFBA table")
    reaction_cols = sorted(set(reference.columns) - META_COLS)
    if args.reaction_list is not None:
        requested = pd.read_csv(args.reaction_list)
        validate_ids(requested, "reaction_id", "Reaction list")
        keep = requested["reaction_id"].astype(str).tolist()
        require_coverage(keep, reaction_cols, "requested reactions in reference pFBA table")
        reaction_cols = keep
    if not reaction_cols:
        raise ValueError("No reactions selected for pFBA reference comparison")
    numeric_cols = [*reaction_cols, "growth_rate", "total_flux_sum"]
    require_coverage(numeric_cols, observed.columns, "required columns in observed pFBA table")
    require_coverage(numeric_cols, reference.columns, "required columns in reference pFBA table")
    observed[numeric_cols] = finite_numeric(observed[numeric_cols], "Observed pFBA values")
    reference[numeric_cols] = finite_numeric(reference[numeric_cols], "Reference pFBA values")

    failed = False
    for sample in samples:
        obs_row = observed.loc[observed["sample"].eq(sample)].iloc[0]
        ref_row = reference.loc[reference["sample"].eq(sample)].iloc[0]
        flux_diff = (
            pd.to_numeric(obs_row[reaction_cols], errors="coerce")
            - pd.to_numeric(ref_row[reaction_cols], errors="coerce")
        ).abs()
        growth_diff = abs(float(obs_row["growth_rate"]) - float(ref_row["growth_rate"]))
        total_diff = abs(float(obs_row["total_flux_sum"]) - float(ref_row["total_flux_sum"]))
        max_flux_diff = float(flux_diff.max())
        n_flux_diff = int((flux_diff > args.tolerance).sum())
        print(
            f"{sample}: compared_reactions={len(reaction_cols)} max_flux_abs_diff={max_flux_diff:.6g} "
            f"growth_abs_diff={growth_diff:.6g} total_flux_abs_diff={total_diff:.6g} "
            f"n_flux_diff_gt_tol={n_flux_diff}",
            flush=True,
        )
        failed = failed or max_flux_diff > args.tolerance or growth_diff > args.tolerance or total_diff > args.tolerance

    if failed:
        raise SystemExit("pFBA output differs from the stored reference subset")
    print("M-SCAN pFBA reference check passed", flush=True)


if __name__ == "__main__":
    main()
