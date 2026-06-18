#!/usr/bin/env python3
"""Apply optional mutation attenuation to a prepared gene-by-sample RNA table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REQUIRED_MUTATION_COLS = {"sample", "gene_id"}


def apply_mutations(rna: pd.DataFrame, mutations: pd.DataFrame, default_factor: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    if "gene_id" not in rna.columns:
        raise ValueError("RNA table must contain a 'gene_id' column")
    missing = REQUIRED_MUTATION_COLS - set(mutations.columns)
    if missing:
        raise ValueError(f"mutation table is missing required columns: {sorted(missing)}")

    out = rna.copy()
    sample_cols = [col for col in out.columns if col != "gene_id"]
    out = out.set_index("gene_id")

    applied_rows = []
    for _, mut in mutations.iterrows():
        sample = str(mut["sample"])
        gene_id = str(mut["gene_id"])
        if sample not in sample_cols or gene_id not in out.index:
            continue

        factor = float(mut["mutation_factor"]) if "mutation_factor" in mutations.columns and pd.notna(mut["mutation_factor"]) else default_factor
        original = float(out.at[gene_id, sample])
        adjusted = original * factor
        out.at[gene_id, sample] = adjusted
        applied_rows.append(
            {
                "sample": sample,
                "gene_id": gene_id,
                "mutation_factor": factor,
                "rna_original": original,
                "rna_adjusted": adjusted,
            }
        )

    return out.reset_index(), pd.DataFrame(applied_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna", type=Path, required=True)
    parser.add_argument("--mutation-table", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--applied-out", type=Path, required=True)
    parser.add_argument("--default-factor", type=float, default=0.05)
    args = parser.parse_args()

    rna = pd.read_csv(args.rna)
    mutations = pd.read_csv(args.mutation_table)
    adjusted, applied = apply_mutations(rna, mutations, args.default_factor)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.applied_out.parent.mkdir(parents=True, exist_ok=True)
    adjusted.to_csv(args.out, index=False)
    applied.to_csv(args.applied_out, index=False)

    print(f"applied mutations: {len(applied)}")
    print(f"saved adjusted RNA: {args.out}")
    print(f"saved applied mutation log: {args.applied_out}")


if __name__ == "__main__":
    main()
