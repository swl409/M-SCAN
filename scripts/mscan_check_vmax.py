#!/usr/bin/env python3
"""Compare generated Vmax constraints with a stored reference subset."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


META_COLS = {"id", "name", "direction", "original_lb", "original_ub"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()

    observed = pd.read_csv(args.observed, low_memory=False)
    reference = pd.read_csv(args.reference, low_memory=False)
    if "id" not in observed.columns or "id" not in reference.columns:
        raise ValueError("Both Vmax tables must contain an 'id' column")

    samples = sorted((set(observed.columns) - META_COLS).intersection(set(reference.columns) - META_COLS))
    if not samples:
        raise ValueError("No shared sample columns between observed and reference Vmax tables")

    merged = observed[["id", *samples]].merge(
        reference[["id", *samples]],
        on="id",
        suffixes=("_observed", "_reference"),
        how="inner",
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
