#!/usr/bin/env python3
"""Compare M-SCAN pFBA output with a stored reference subset."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


META_COLS = {"sample", "total_flux_sum", "growth_rate", "status", "error_msg"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--reaction-list", type=Path, default=None)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()

    observed = pd.read_csv(args.observed, low_memory=False)
    reference = pd.read_csv(args.reference, low_memory=False)

    samples = sorted(set(observed["sample"]).intersection(reference["sample"]))
    reaction_cols = sorted((set(observed.columns) - META_COLS).intersection(set(reference.columns) - META_COLS))
    if args.reaction_list is not None:
        keep = pd.read_csv(args.reaction_list)["reaction_id"].dropna().astype(str).tolist()
        reaction_cols = [reaction_id for reaction_id in keep if reaction_id in reaction_cols]
    if not samples:
        raise ValueError("No shared samples between observed and reference pFBA tables")
    if not reaction_cols:
        raise ValueError("No shared reaction columns between observed and reference pFBA tables")

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
