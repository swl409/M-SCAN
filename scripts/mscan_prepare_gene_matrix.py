#!/usr/bin/env python3
"""Prepare RNA expression for Human-GEM Vmax construction."""

from __future__ import annotations

import argparse
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


ENSG_RE = re.compile(r"ENSG[0-9]+")


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def load_gene_list(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    genes = {
        line.strip().split(".")[0]
        for line in path.read_text().splitlines()
        if line.strip()
    }
    return genes or None


def load_humangem_genes(human_gem_xlsx: Path) -> set[str]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        rxns = pd.read_excel(human_gem_xlsx, sheet_name="RXNS")
    genes: set[str] = set()
    for gpr in rxns["GENE ASSOCIATION"].dropna():
        genes.update(ENSG_RE.findall(str(gpr)))
    return genes


def prepare_rna(
    rna: pd.DataFrame,
    human_gem_genes: set[str],
    protein_coding_genes: set[str] | None,
    input_scale: str,
    target_sum: float,
    replace_zero_min: bool,
) -> tuple[pd.DataFrame, dict[str, int]]:
    if rna.empty:
        raise ValueError("RNA table is empty")

    gene_col = rna.columns[0]
    rna = rna.rename(columns={gene_col: "gene_id"}).copy()
    rna["gene_id"] = rna["gene_id"].astype(str).str.split(".").str[0]

    sample_cols = [col for col in rna.columns if col != "gene_id"]
    rna[sample_cols] = rna[sample_cols].apply(pd.to_numeric, errors="coerce")

    n_input = len(rna)
    if protein_coding_genes is not None:
        rna = rna[rna["gene_id"].isin(protein_coding_genes)].copy()
    n_protein_coding = len(rna)

    rna = rna.drop_duplicates(subset=["gene_id"], keep="first").set_index("gene_id")

    if input_scale == "log2":
        values = np.exp2(rna[sample_cols]) - 0.0001
        values[values < 0] = 0
    elif input_scale == "linear":
        values = rna[sample_cols].clip(lower=0)
    elif input_scale == "already_normalized":
        values = rna[sample_cols].clip(lower=0)
    else:
        raise ValueError(f"Unsupported input scale: {input_scale}")

    if input_scale != "already_normalized":
        sample_sums = values.sum(axis=0)
        sample_sums[sample_sums == 0] = 1.0
        values = values.div(sample_sums, axis=1) * target_sum

    values = values.reset_index()
    values = values[values["gene_id"].isin(human_gem_genes)].copy()

    if replace_zero_min:
        for col in sample_cols:
            vals = pd.to_numeric(values[col], errors="coerce")
            positive = vals[vals > 0]
            if len(positive) == 0:
                continue
            values.loc[vals == 0, col] = positive.min()

    values = values.reset_index(drop=True)
    stats = {
        "n_input_rows": int(n_input),
        "n_protein_coding_rows": int(n_protein_coding),
        "n_humangem_rows": int(len(values)),
        "n_samples": int(len(sample_cols)),
    }
    return values, stats


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna", type=Path, required=True)
    parser.add_argument("--human-gem-xlsx", type=Path, required=True)
    parser.add_argument("--protein-coding-genes", type=Path, default=None)
    parser.add_argument(
        "--input-scale",
        choices=["log2", "linear", "already_normalized"],
        default="log2",
        help="Scale of the RNA input values. CCLE/DepMap RNA values are usually log2(TPM+1).",
    )
    parser.add_argument("--target-sum", type=float, default=1_000_000.0)
    parser.add_argument("--no-replace-zero-min", action="store_true")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    rna = read_table(args.rna)
    protein_coding = load_gene_list(args.protein_coding_genes)
    humangem_genes = load_humangem_genes(args.human_gem_xlsx)
    prepared, stats = prepare_rna(
        rna,
        human_gem_genes=humangem_genes,
        protein_coding_genes=protein_coding,
        input_scale=args.input_scale,
        target_sum=args.target_sum,
        replace_zero_min=not args.no_replace_zero_min,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    prepared.to_csv(args.out, index=False)
    print(
        "prepared RNA: "
        f"input_rows={stats['n_input_rows']} "
        f"protein_coding_rows={stats['n_protein_coding_rows']} "
        f"humangem_rows={stats['n_humangem_rows']} "
        f"samples={stats['n_samples']}",
        flush=True,
    )
    print(f"saved prepared RNA: {args.out}", flush=True)


if __name__ == "__main__":
    main()
