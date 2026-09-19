#!/usr/bin/env python3
"""Build Human-GEM reaction Vmax constraints from RNA activity."""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from cobra.io import load_matlab_model

from mscan_gpr import collect_genes, eval_gpr, parse_gpr


META_COLS = ["id", "name", "direction", "original_lb", "original_ub"]
GENE_ID_ALIASES = {"", "Unnamed: 0", "index", "gene", "genes", "ensembl_gene_id"}
GENE_ID_ALIASES_LOWER = {col.lower() for col in GENE_ID_ALIASES}


def standardize_gene_id_column(rna_activity: pd.DataFrame) -> pd.DataFrame:
    """Accept gene_id columns or rowname-style first columns from CSV exports."""
    if "gene_id" in rna_activity.columns:
        out = rna_activity.copy()
    else:
        first_col = rna_activity.columns[0]
        first_name = str(first_col)
        if first_name in GENE_ID_ALIASES or first_name.lower() in GENE_ID_ALIASES_LOWER:
            out = rna_activity.rename(columns={first_col: "gene_id"}).copy()
        else:
            raise ValueError(
                "RNA table must contain a 'gene_id' column, or use a rowname-style "
                "first column such as 'Unnamed: 0'."
            )
    out["gene_id"] = out["gene_id"].astype(str).str.split(".").str[0]
    return out


def read_humangem_reactions(human_gem_xlsx: Path, human_gem_mat: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    model = load_matlab_model(str(human_gem_mat))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        gem = pd.read_excel(human_gem_xlsx, sheet_name="RXNS")
    gem_lookup = gem.set_index("ID").to_dict("index")

    expanded = []
    meta = []
    for reaction in model.reactions:
        if reaction.id in gem_lookup:
            row = gem_lookup[reaction.id].copy()
            row["ID"] = reaction.id
        else:
            row = {
                "ID": reaction.id,
                "NAME": reaction.name,
                "EQUATION": reaction.reaction,
                "GENE ASSOCIATION": reaction.gene_reaction_rule,
                "LOWER BOUND": reaction.lower_bound,
                "UPPER BOUND": reaction.upper_bound,
                "SUBSYSTEM": reaction.subsystem,
                "EC-NUMBER": np.nan,
                "OBJECTIVE": reaction.objective_coefficient,
            }
        row["direction"] = "forward"
        expanded.append(row)
        meta.append(
            {
                "id": reaction.id,
                "name": reaction.name,
                "direction": "forward",
                "original_lb": reaction.lower_bound,
                "original_ub": reaction.upper_bound,
            }
        )
        if reaction.lower_bound < 0:
            reverse_row = row.copy()
            reverse_row["ID"] = f"{reaction.id}_REV"
            reverse_row["NAME"] = f"{reaction.name}_REV"
            reverse_row["direction"] = "reverse"
            expanded.append(reverse_row)
            meta.append(
                {
                    "id": f"{reaction.id}_REV",
                    "name": f"{reaction.name}_REV",
                    "direction": "reverse",
                    "original_lb": reaction.lower_bound,
                    "original_ub": reaction.upper_bound,
                }
            )

    return pd.DataFrame(expanded), pd.DataFrame(meta)


def read_kcat_table(kcat_path: Path) -> tuple[pd.DataFrame, pd.Series, float]:
    kcat = pd.read_csv(kcat_path, low_memory=False)
    required = {"reaction_id", "gene_id", "predicted_kcat_s_inv"}
    missing = required - set(kcat.columns)
    if missing:
        raise ValueError(f"kcat table is missing required columns: {sorted(missing)}")

    kcat = kcat[["reaction_id", "gene_id", "predicted_kcat_s_inv"]].copy()
    kcat["reaction_id"] = kcat["reaction_id"].astype(str)
    kcat["gene_id"] = kcat["gene_id"].astype(str)
    kcat["predicted_kcat_s_inv"] = pd.to_numeric(kcat["predicted_kcat_s_inv"], errors="coerce")
    kcat = kcat.dropna(subset=["predicted_kcat_s_inv"])
    kcat_agg = (
        kcat.groupby(["reaction_id", "gene_id"], as_index=False)["predicted_kcat_s_inv"]
        .mean()
        .rename(columns={"predicted_kcat_s_inv": "kcat_s_inv"})
    )
    kcat_agg["kcat_s_inv"] = kcat_agg["kcat_s_inv"].clip(lower=0.0)
    reaction_median = kcat_agg.groupby("reaction_id")["kcat_s_inv"].median()
    global_median = float(kcat_agg["kcat_s_inv"].median()) if len(kcat_agg) else 1.0
    return kcat_agg, reaction_median, global_median


def build_vmax(
    rna_activity: pd.DataFrame,
    expanded_reactions: pd.DataFrame,
    meta: pd.DataFrame,
    kcat_bundle: tuple[pd.DataFrame, pd.Series, float] | None = None,
    fallback_kcat: float = 1.0,
) -> pd.DataFrame:
    rna_activity = standardize_gene_id_column(rna_activity)

    sample_cols = [col for col in rna_activity.columns if col != "gene_id"]
    if not sample_cols:
        raise ValueError("RNA table must contain at least one sample column")

    rna_activity = rna_activity.set_index("gene_id")
    gene = rna_activity.reset_index()
    gene[sample_cols] = gene[sample_cols].apply(pd.to_numeric, errors="coerce")
    first_col = "gene_id"

    if kcat_bundle is not None:
        kcat_agg, reaction_median, global_median = kcat_bundle
    else:
        kcat_agg = pd.DataFrame(columns=["reaction_id", "gene_id", "kcat_s_inv"])
        reaction_median = pd.Series(dtype=float)
        global_median = fallback_kcat

    valid = expanded_reactions[expanded_reactions["GENE ASSOCIATION"].notna()].copy()
    valid["has_or"] = valid["GENE ASSOCIATION"].str.contains(" or ", case=False, na=False)
    valid["has_and"] = valid["GENE ASSOCIATION"].str.contains(" and ", case=False, na=False)
    or_rxns = valid.loc[~valid["has_and"]].reset_index(drop=True)
    and_rxns = valid.loc[(~valid["has_or"]) & valid["has_and"]].reset_index(drop=True)
    mixed_rxns = valid.loc[valid["has_or"] & valid["has_and"]].reset_index(drop=True)

    def simple_vmax(rxns: pd.DataFrame, reducer: str) -> pd.DataFrame:
        rmap = rxns[["ID", "GENE ASSOCIATION"]].copy()
        rmap["gene_list"] = (
            rmap["GENE ASSOCIATION"]
            .fillna("")
            .str.replace(" or ", " ", regex=False)
            .str.replace(" and ", " ", regex=False)
            .str.split()
        )
        exploded = rmap.explode("gene_list")
        if kcat_bundle is not None:
            merged = pd.merge(
                exploded,
                kcat_agg,
                left_on=["ID", "gene_list"],
                right_on=["reaction_id", "gene_id"],
                how="left",
            )
            merged["reaction_median"] = pd.to_numeric(merged["ID"].map(reaction_median), errors="coerce")
            merged["kcat_filled"] = (
                pd.to_numeric(merged["kcat_s_inv"], errors="coerce")
                .combine_first(merged["reaction_median"])
                .fillna(global_median)
            )
        else:
            merged = exploded.copy()
            merged["kcat_filled"] = fallback_kcat
        final = pd.merge(merged, gene, left_on="gene_list", right_on=first_col, how="left")
        final[sample_cols] = final[sample_cols].mul(final["kcat_filled"], axis=0)
        grouped = final.groupby("ID")[sample_cols]
        out = grouped.sum(min_count=1) if reducer == "sum" else grouped.min()
        return out.reindex(rxns["ID"]).reset_index().rename(columns={"ID": "id"})

    or_vmax = simple_vmax(or_rxns, "sum")
    and_vmax = simple_vmax(and_rxns, "min")

    mixed_rows = []
    for _, reaction in mixed_rxns.iterrows():
        rid = reaction["ID"]
        ast = parse_gpr(reaction["GENE ASSOCIATION"])
        genes = sorted(collect_genes(ast))
        gdf = pd.merge(pd.DataFrame({first_col: genes}), gene, on=first_col, how="left")
        if kcat_bundle is not None:
            reaction_kcat = kcat_agg.loc[kcat_agg["reaction_id"].eq(rid), ["gene_id", "kcat_s_inv"]]
            tmp = pd.merge(gdf, reaction_kcat, left_on=first_col, right_on="gene_id", how="left")
            fallback = reaction_median.get(rid, np.nan)
            fallback = fallback if pd.notna(fallback) else global_median
            tmp["kcat_filled"] = pd.to_numeric(tmp["kcat_s_inv"], errors="coerce").fillna(fallback)
        else:
            tmp = gdf.copy()
            tmp["kcat_filled"] = fallback_kcat
        tmp[sample_cols] = tmp[sample_cols].mul(tmp["kcat_filled"], axis=0)
        tmp = tmp.set_index(first_col)
        row = {"id": rid}
        for sample in sample_cols:
            row[sample] = eval_gpr(ast, tmp[sample].to_dict())
        mixed_rows.append(row)
    mixed_vmax = pd.DataFrame(mixed_rows)

    vmax = pd.concat([or_vmax, and_vmax, mixed_vmax], ignore_index=True)
    meta = meta.copy()
    meta["id_upper"] = meta["id"].astype(str).str.upper()
    vmax["id_upper"] = vmax["id"].astype(str).str.upper()
    filled = pd.merge(meta, vmax, on="id_upper", how="left", suffixes=("", "_drop"))
    filled = filled.drop(columns=[col for col in filled.columns if col == "id_upper" or col.endswith("_drop")])

    for sample in sample_cols:
        filled[sample] = pd.to_numeric(filled[sample], errors="coerce")
        max_val = filled[sample].max()
        if pd.isna(max_val):
            max_val = 1000.0
        filled[sample] = filled[sample].fillna(max_val)

    scaled = filled.copy()
    scaled[sample_cols] = np.log2(scaled[sample_cols] + 1.0)
    return scaled


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--human-gem-xlsx", type=Path, required=True)
    parser.add_argument("--human-gem-mat", type=Path, required=True)
    parser.add_argument("--rna", type=Path, required=True)
    parser.add_argument("--kcat", type=Path, default=None)
    parser.add_argument("--fallback-kcat", type=float, default=1.0)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    rna_activity = pd.read_csv(args.rna)
    expanded_reactions, meta = read_humangem_reactions(args.human_gem_xlsx, args.human_gem_mat)
    kcat_bundle = read_kcat_table(args.kcat) if args.kcat is not None else None
    vmax = build_vmax(
        rna_activity,
        expanded_reactions,
        meta,
        kcat_bundle=kcat_bundle,
        fallback_kcat=args.fallback_kcat,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    vmax.to_csv(args.out, index=False)
    sample_cols = [col for col in vmax.columns if col not in META_COLS]
    kcat_msg = "with kcat" if kcat_bundle is not None else "without kcat"
    print(f"built Human-GEM Vmax: reactions={len(vmax)} samples={len(sample_cols)} {kcat_msg} scale=log2(vmax+1)")
    print(f"saved: {args.out}")


if __name__ == "__main__":
    main()
