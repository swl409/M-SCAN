#!/usr/bin/env python3
"""Compare generated Vmax constraints with a stored reference subset."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from mscan_reference_utils import finite_numeric, require_coverage, validate_ids, validate_tolerance


META_COLS = {"id", "name", "direction", "original_lb", "original_ub"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()
    validate_tolerance(args.tolerance)

    observed = pd.read_csv(args.observed, low_memory=False)
    reference = pd.read_csv(args.reference, low_memory=False)
    validate_ids(observed, "id", "Observed Vmax table")
    validate_ids(reference, "id", "Reference Vmax table")
    samples = sorted(set(reference.columns) - META_COLS)
    if not samples:
        raise ValueError("Reference Vmax table contains no sample columns")
    require_coverage(samples, observed.columns, "reference sample columns in observed Vmax table")
    require_coverage(reference["id"], observed["id"], "reference reaction ids in observed Vmax table")
    observed[samples] = finite_numeric(observed[samples], "Observed Vmax values")
    reference[samples] = finite_numeric(reference[samples], "Reference Vmax values")

    merged = observed[["id", *samples]].merge(
        reference[["id", *samples]],
        on="id",
        suffixes=("_observed", "_reference"),
        how="inner",
        validate="one_to_one",
    )
    if merged.empty:
        raise ValueError("No shared reaction ids between observed and reference Vmax tables")

    failed = False
    for sample in samples:
        diff = (
            pd.to_numeric(merged[f"{sample}_observed"], errors="coerce")
            - pd.to_numeric(merged[f"{sample}_reference"], errors="coerce")
        ).abs()
        max_diff = float(diff.max())
        n_diff = int((diff > args.tolerance).sum())
        print(
            f"{sample}: compared_reactions={len(merged)} "
            f"max_vmax_abs_diff={max_diff:.6g} n_diff_gt_tol={n_diff}",
            flush=True,
        )
        failed = failed or max_diff > args.tolerance

    if failed:
        raise SystemExit("Vmax output differs from the stored reference subset")
    print("M-SCAN Vmax reference check passed", flush=True)


if __name__ == "__main__":
    main()
