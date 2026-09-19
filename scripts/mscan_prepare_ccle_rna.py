#!/usr/bin/env python3
"""Prepare CCLE RNA for the M-SCAN Human-GEM workflow."""

from __future__ import annotations

import argparse
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


ENSG_RE = re.compile(r"ENSG[0-9]+")
METADATA_COLS = [
    "depmap_id",
    "lineage_1",
    "lineage_2",
    "lineage_3",
    "lineage_6",
    "lineage_4",
    "Unnamed: 0",
]
DEFAULT_SAMPLE_COLUMN = "cell_line_display_name"
ROWNAME_SAMPLE_COLUMNS = {"", "Unnamed: 0", "index"}
ROWNAMES_LOWER = {col.lower() for col in ROWNAME_SAMPLE_COLUMNS}
GENE_ID_COLUMNS = {"gene_id", "gene", "genes", "ensembl_gene_id", "Unnamed: 0", "index"}
GENE_ID_COLUMNS_LOWER = {col.lower() for col in GENE_ID_COLUMNS}
SEVERE_TERMS = [
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_lost",
]


def load_gene_list(path: Path) -> set[str]:
    return {
        line.strip().split(".")[0]
        for line in path.read_text().splitlines()
        if line.strip()
    }


def load_humangem_genes(human_gem_xlsx: Path) -> set[str]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        rxns = pd.read_excel(human_gem_xlsx, sheet_name="RXNS")

    genes: set[str] = set()
    for gpr in rxns["GENE ASSOCIATION"].dropna():
        genes.update(ENSG_RE.findall(str(gpr)))
    return genes


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def standardize_sample_column(rna: pd.DataFrame, sample_column: str) -> pd.DataFrame:
    """Accept explicit sample columns or rowname-style first columns."""
    if sample_column in rna.columns:
        if sample_column == DEFAULT_SAMPLE_COLUMN:
            return rna
        return rna.rename(columns={sample_column: DEFAULT_SAMPLE_COLUMN})

    for col in rna.columns:
        col_name = str(col)
        if col_name in ROWNAME_SAMPLE_COLUMNS or col_name.lower() in ROWNAMES_LOWER:
            return rna.rename(columns={col: DEFAULT_SAMPLE_COLUMN})

    first_col = str(rna.columns[0])
    if first_col in ROWNAME_SAMPLE_COLUMNS or first_col.lower() in ROWNAMES_LOWER:
        return rna.rename(columns={rna.columns[0]: DEFAULT_SAMPLE_COLUMN})

    raise ValueError(
        "CCLE RNA input must contain a sample column named "
        f"'{sample_column}' or a rowname-style first column such as 'Unnamed: 0'."
    )


def standardize_gene_id_column(rna: pd.DataFrame) -> pd.DataFrame:
    """Accept explicit gene_id columns or rowname-style gene columns."""
    if "gene_id" in rna.columns:
        out = rna.copy()
    else:
        first_col = rna.columns[0]
        first_name = str(first_col)
        if first_name in GENE_ID_COLUMNS or first_name.lower() in GENE_ID_COLUMNS_LOWER:
            out = rna.rename(columns={first_col: "gene_id"}).copy()
        else:
            raise ValueError(
                "Gene-row RNA input must contain a gene ID first column, such as "
                "'gene_id' or rowname-style 'Unnamed: 0'."
            )
    out["gene_id"] = out["gene_id"].astype(str).str.split(".").str[0]
    return out


def fraction_ensembl_like(values: pd.Series) -> float:
    if values.empty:
        return 0.0
    text = values.dropna().astype(str).str.split(".").str[0]
    if text.empty:
        return 0.0
    return float(text.str.match(r"^ENSG[0-9]+$").mean())


def infer_rna_layout(rna: pd.DataFrame, sample_column: str, requested_layout: str) -> str:
    if requested_layout != "auto":
        return requested_layout

    if sample_column in rna.columns or DEFAULT_SAMPLE_COLUMN in rna.columns:
        return "sample_rows"

    for col in rna.columns:
        col_name = str(col)
        if col_name in GENE_ID_COLUMNS or col_name.lower() in GENE_ID_COLUMNS_LOWER:
            if fraction_ensembl_like(rna[col]) >= 0.5:
                return "gene_rows"

    first_col = rna.columns[0]
    if fraction_ensembl_like(rna[first_col]) >= 0.5:
        return "gene_rows"

    return "sample_rows"


def preprocess_rna(
    rna: pd.DataFrame,
    protein_coding_genes: set[str],
    humangem_genes: set[str],
    target_sum: float,
    sample_column: str = DEFAULT_SAMPLE_COLUMN,
    rna_layout: str = "auto",
) -> tuple[pd.DataFrame, dict[str, int]]:
    layout = infer_rna_layout(rna, sample_column, rna_layout)
    if layout == "gene_rows":
        return preprocess_gene_row_rna(
            rna,
            protein_coding_genes=protein_coding_genes,
            humangem_genes=humangem_genes,
            target_sum=target_sum,
            input_layout=layout,
        )
    if layout != "sample_rows":
        raise ValueError("rna_layout must be 'auto', 'sample_rows', or 'gene_rows'")

    rna = standardize_sample_column(rna, sample_column)

    existing_meta = [col for col in METADATA_COLS if col in rna.columns]
    meta_df = rna[existing_meta].copy()
    meta_df[DEFAULT_SAMPLE_COLUMN] = rna[DEFAULT_SAMPLE_COLUMN].astype(str).values

    rna = rna.set_index(DEFAULT_SAMPLE_COLUMN)
    data_numeric = rna.drop(columns=existing_meta, errors="ignore").apply(pd.to_numeric, errors="coerce")

    # The M-SCAN CCLE workflow uses DepMap log2(TPM+1) values converted back to linear TPM-like values.
    data_linear = np.exp2(data_numeric) - 1.0
    data_linear[data_linear < 0] = 0.0

    clean_cols = [col.split(".")[0] for col in data_linear.columns]
    valid_indices = [i for i, gene_id in enumerate(clean_cols) if gene_id in protein_coding_genes]
    data_coding = data_linear.iloc[:, valid_indices].copy()
    data_coding.columns = [data_coding.columns[i] for i in range(data_coding.shape[1])]

    current_sums = data_coding.sum(axis=1)
    current_sums[current_sums == 0] = 1.0
    data_normalized = data_coding.div(current_sums, axis=0) * target_sum

    rna_norm = data_normalized.reset_index()
    for col in meta_df.columns:
        if col != "cell_line_display_name":
            rna_norm[col] = meta_df[col].values

    metadata_cols2 = [DEFAULT_SAMPLE_COLUMN, *METADATA_COLS]
    existing_meta2 = [col for col in metadata_cols2 if col in rna_norm.columns]
    gene_cols = [col for col in rna_norm.columns if col not in existing_meta2]
    gene_col_map = {col: col.split(".")[0] for col in gene_cols}
    valid_gene_cols = [col for col in gene_cols if gene_col_map[col] in humangem_genes]

    meta_df2 = rna_norm[existing_meta2].copy()
    gene_df = rna_norm[valid_gene_cols].copy()
    gene_df.columns = [col.split(".")[0] for col in gene_df.columns]
    gene_df = gene_df.T.groupby(level=0).mean().T
    rna_filtered = pd.concat([meta_df2, gene_df], axis=1)

    gene_only_cols = [col for col in rna_filtered.columns if col not in existing_meta2]
    for idx in rna_filtered.index:
        vals = pd.to_numeric(rna_filtered.loc[idx, gene_only_cols], errors="coerce")
        positive_vals = vals[vals > 0]
        if len(positive_vals) == 0:
            continue
        min_positive = positive_vals.min()
        rna_filtered.loc[idx, gene_only_cols] = np.where(vals <= 0, min_positive, vals)

    gene = rna_filtered.set_index(DEFAULT_SAMPLE_COLUMN)
    gene_values = gene.drop(columns=[col for col in METADATA_COLS if col in gene.columns], errors="ignore")
    out = gene_values.T.reset_index().rename(columns={"index": "gene_id"})
    out["gene_id"] = out["gene_id"].astype(str).str.split(".").str[0]

    stats = {
        "input_layout": layout,
        "input_samples": int(len(rna)),
        "input_gene_columns": int(data_numeric.shape[1]),
        "protein_coding_columns": int(data_coding.shape[1]),
        "humangem_genes": int(len(out)),
        "samples": int(len([col for col in out.columns if col != "gene_id"])),
    }
    return out, stats


def preprocess_gene_row_rna(
    rna: pd.DataFrame,
    protein_coding_genes: set[str],
    humangem_genes: set[str],
    target_sum: float,
    input_layout: str,
) -> tuple[pd.DataFrame, dict[str, int]]:
    rna = standardize_gene_id_column(rna)
    sample_cols = [col for col in rna.columns if col != "gene_id"]
    if not sample_cols:
        raise ValueError("Gene-row RNA input must contain at least one sample column")

    rna[sample_cols] = rna[sample_cols].apply(pd.to_numeric, errors="coerce")
    n_input_genes = len(rna)
    rna = rna[rna["gene_id"].isin(protein_coding_genes)].copy()
    n_protein_coding = len(rna)
    rna = rna.set_index("gene_id")

    # CCLE-mode gene-row input is still interpreted as DepMap log2(TPM+1).
    values = np.exp2(rna[sample_cols]) - 1.0
    values[values < 0] = 0.0

    sample_sums = values.sum(axis=0)
    sample_sums[sample_sums == 0] = 1.0
    values = values.div(sample_sums, axis=1) * target_sum

    values = values.reset_index()
    values = values[values["gene_id"].isin(humangem_genes)].copy()
    values = values.groupby("gene_id", as_index=False)[sample_cols].mean()
    for col in sample_cols:
        vals = pd.to_numeric(values[col], errors="coerce")
        positive_vals = vals[vals > 0]
        if len(positive_vals) == 0:
            continue
        min_positive = positive_vals.min()
        values.loc[vals <= 0, col] = min_positive

    values = values.reset_index(drop=True)
    stats = {
        "input_layout": input_layout,
        "input_samples": int(len(sample_cols)),
        "input_gene_columns": int(n_input_genes),
        "protein_coding_columns": int(n_protein_coding),
        "humangem_genes": int(len(values)),
        "samples": int(len(sample_cols)),
    }
    return values, stats


def apply_ccle_mutations(
    gene: pd.DataFrame,
    mutation: pd.DataFrame,
    severe_factor: float,
    lof_factor: float,
    include_lof: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {"cell_line_display_name", "Gene", "ensembl_gene_id", "variant_info", "Consequence", "final_prediction"}
    missing = required - set(mutation.columns)
    if missing:
        raise ValueError(f"Mutation table is missing required CCLE annotation columns: {sorted(missing)}")

    gene = gene.copy()
    sample_cols = [col for col in gene.columns if col != "gene_id"]
    gene["gene_id"] = gene["gene_id"].astype(str).str.split(".").str[0]
    gene_idx_map = dict(zip(gene["gene_id"], gene.index))

    mutation = mutation.copy()
    mutation["mutation_gene_id"] = mutation["Gene"]
    missing_gene = mutation["mutation_gene_id"].isna()

    # Reference workflow compatibility: true missing Gene values were not
    # rescued from ensembl_gene_id, but literal strings such as "nan" were.
    literal_bad_gene = (
        mutation["mutation_gene_id"].notna()
        & mutation["mutation_gene_id"].astype(str).isin(["nan", "NaN", "None"])
    )
    mutation.loc[literal_bad_gene, "mutation_gene_id"] = mutation.loc[
        literal_bad_gene, "ensembl_gene_id"
    ].astype(str)
    mutation["mutation_gene_id"] = mutation["mutation_gene_id"].astype(str).str.split(".").str[0]
    mutation.loc[missing_gene, "mutation_gene_id"] = np.nan
    mutation["sample"] = mutation["cell_line_display_name"].astype(str)

    text = (
        mutation["variant_info"].fillna("").astype(str)
        + " "
        + mutation["Consequence"].fillna("").astype(str)
    ).str.lower()
    is_severe = pd.Series(False, index=mutation.index)
    for term in SEVERE_TERMS:
        is_severe = is_severe | text.str.contains(term, regex=False)
    is_lof = mutation["final_prediction"].fillna("").astype(str).str.strip().eq("LOF")

    target = mutation.loc[
        mutation["sample"].isin(sample_cols)
        & mutation["mutation_gene_id"].notna()
        & (is_severe | is_lof)
    ].copy()
    target["is_severe"] = is_severe.loc[target.index]
    target["is_lof"] = is_lof.loc[target.index]

    applied: dict[tuple[str, str], dict[str, object]] = {}
    for _, row in target.iterrows():
        sample = row["sample"]
        gene_id = row["mutation_gene_id"]
        if sample not in sample_cols or gene_id not in gene_idx_map:
            continue
        if bool(row["is_severe"]):
            factor = severe_factor
            reason = "severe"
        elif include_lof and bool(row["is_lof"]):
            factor = lof_factor
            reason = "LoGoFunc_LOF"
        else:
            continue

        key = (gene_id, sample)
        previous = applied.get(key)
        if previous is None or factor < float(previous["mutation_factor"]):
            applied[key] = {
                "gene_id": gene_id,
                "sample": sample,
                "mutation_factor": factor,
                "mutation_reason": reason,
                "n_matching_mutations": 1,
            }
        else:
            previous["n_matching_mutations"] = int(previous["n_matching_mutations"]) + 1

    for (gene_id, sample), info in applied.items():
        idx = gene_idx_map[gene_id]
        gene.at[idx, sample] = float(gene.at[idx, sample]) * float(info["mutation_factor"])

    applied_df = pd.DataFrame(applied.values())
    if not applied_df.empty:
        applied_df = applied_df.sort_values(["sample", "gene_id"]).reset_index(drop=True)
    return gene, applied_df


def validate_residual_factor(value: float, name: str) -> float:
    value = float(value)
    if value < 0.0 or value > 1.0:
        raise ValueError(f"{name} must be between 0 and 1, got {value}")
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna", type=Path, required=True)
    parser.add_argument("--mutation-table", type=Path, default=None)
    parser.add_argument("--human-gem-xlsx", type=Path, required=True)
    parser.add_argument("--protein-coding-genes", type=Path, required=True)
    parser.add_argument(
        "--sample-column",
        default=DEFAULT_SAMPLE_COLUMN,
        help=(
            "RNA sample/cell-line column name. If absent, a rowname-style first "
            "column such as 'Unnamed: 0' is also accepted."
        ),
    )
    parser.add_argument(
        "--rna-layout",
        choices=["auto", "sample_rows", "gene_rows"],
        default="auto",
        help=(
            "RNA input layout. sample_rows means rows are cell lines and columns "
            "are genes. gene_rows means rows are genes and columns are samples."
        ),
    )
    parser.add_argument("--target-sum", type=float, default=1_000_000.0)
    parser.add_argument(
        "--mutation-factor",
        type=float,
        default=None,
        help="Residual-capacity factor applied to all mutation-affected genes unless overridden.",
    )
    parser.add_argument("--severe-factor", type=float, default=None)
    parser.add_argument("--lof-factor", type=float, default=None)
    parser.add_argument("--include-lof", action="store_true")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--preprocessed-out", type=Path, default=None)
    parser.add_argument("--applied-out", type=Path, default=None)
    args = parser.parse_args()
    default_mutation_factor = 0.05 if args.mutation_factor is None else args.mutation_factor
    severe_factor = validate_residual_factor(
        default_mutation_factor if args.severe_factor is None else args.severe_factor,
        "severe_factor",
    )
    lof_factor = validate_residual_factor(
        default_mutation_factor if args.lof_factor is None else args.lof_factor,
        "lof_factor",
    )

    protein_coding = load_gene_list(args.protein_coding_genes)
    humangem_genes = load_humangem_genes(args.human_gem_xlsx)
    rna = read_table(args.rna)
    prepared, stats = preprocess_rna(
        rna,
        protein_coding_genes=protein_coding,
        humangem_genes=humangem_genes,
        target_sum=args.target_sum,
        sample_column=args.sample_column,
        rna_layout=args.rna_layout,
    )

    if args.preprocessed_out is not None:
        args.preprocessed_out.parent.mkdir(parents=True, exist_ok=True)
        prepared.to_csv(args.preprocessed_out, index=False)

    final = prepared
    applied = pd.DataFrame()
    if args.mutation_table is not None:
        mutation = read_table(args.mutation_table)
        final, applied = apply_ccle_mutations(
            prepared,
            mutation,
            severe_factor=severe_factor,
            lof_factor=lof_factor,
            include_lof=args.include_lof,
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(args.out, index=False)
    if args.applied_out is not None:
        args.applied_out.parent.mkdir(parents=True, exist_ok=True)
        applied.to_csv(args.applied_out, index=False)

    print(
        "prepared M-SCAN CCLE RNA: "
        f"input_layout={stats['input_layout']} "
        f"input_samples={stats['input_samples']} "
        f"input_gene_columns={stats['input_gene_columns']} "
        f"protein_coding_columns={stats['protein_coding_columns']} "
        f"humangem_genes={stats['humangem_genes']} "
        f"samples={stats['samples']} "
        f"severe_factor={severe_factor} "
        f"lof_factor={lof_factor} "
        f"applied_mutations={len(applied)}",
        flush=True,
    )
    print(f"saved RNA for Vmax: {args.out}", flush=True)


if __name__ == "__main__":
    main()
