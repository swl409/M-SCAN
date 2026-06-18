#!/usr/bin/env python3
"""Run Human-GEM pFBA for M-SCAN samples."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from cobra.flux_analysis import pfba
from cobra.io import load_matlab_model
from cobra.util import linear_reaction_coefficients


VMAX_META = {"id", "name", "direction", "original_lb", "original_ub", "reaction_id", "reaction_name", "gpr"}


def configure_solver(model, solver: str) -> str:
    if solver == "auto":
        for candidate in ("gurobi", "glpk"):
            try:
                model.solver = candidate
                solver = candidate
                break
            except Exception:
                continue
        else:
            raise RuntimeError("No supported COBRApy solver is available")
    else:
        try:
            model.solver = solver
        except Exception as exc:
            raise RuntimeError(
                f"Requested solver '{solver}' is not available. "
                "Install it with the provided environment.yml or rerun with SOLVER=glpk."
            ) from exc

    if solver == "gurobi":
        try:
            model.solver.problem.Params.OutputFlag = 0
        except Exception:
            pass
    return solver


def apply_vmax_bounds(model, vmax: pd.DataFrame, sample: str) -> int:
    id_column = "id" if "id" in vmax.columns else "reaction_id"
    vmax_values = dict(zip(vmax[id_column], pd.to_numeric(vmax[sample], errors="coerce").fillna(0.0)))
    exchange_ids = {reaction.id for reaction in model.exchanges}
    changed = 0
    for reaction in model.reactions:
        if reaction.id in exchange_ids:
            continue
        if reaction.id in vmax_values:
            reaction.upper_bound = max(float(vmax_values[reaction.id]), 0.0)
            changed += 1
        reverse_id = f"{reaction.id}_REV"
        if reverse_id in vmax_values:
            reaction.lower_bound = -max(float(vmax_values[reverse_id]), 0.0)
        elif reaction.lower_bound < 0:
            reaction.lower_bound = 0.0
    return changed


def apply_media_constraints(model, media_path: Path | None) -> int:
    if media_path is None:
        return 0

    media = pd.read_csv(media_path)
    if "HMR" not in media.columns:
        raise ValueError("media table must contain an 'HMR' column")
    media_ids = set(media["HMR"].dropna().astype(str))

    for reaction in model.exchanges:
        reaction.lower_bound = 0.0
        reaction.upper_bound = 1000.0

    opened = 0
    for reaction_id in media_ids:
        if reaction_id not in model.reactions:
            continue
        model.reactions.get_by_id(reaction_id).lower_bound = -1000.0
        opened += 1
    return opened


def run_pfba_for_samples(
    human_gem_mat: Path,
    vmax: pd.DataFrame,
    max_samples: int,
    media_path: Path | None = None,
    solver: str = "gurobi",
) -> pd.DataFrame:
    sample_cols = [col for col in vmax.columns if col not in VMAX_META]
    if max_samples > 0:
        sample_cols = sample_cols[:max_samples]

    base_model = load_matlab_model(str(human_gem_mat))
    solver = configure_solver(base_model, solver)
    print(f"COBRApy solver: {solver}", flush=True)
    objective_coefficients = linear_reaction_coefficients(base_model)
    objective_reaction = next(iter(objective_coefficients), None)
    objective_reaction_id = objective_reaction.id if objective_reaction is not None else None

    rows = []
    for sample in sample_cols:
        with base_model:
            n_constrained = apply_vmax_bounds(base_model, vmax, sample)
            n_media_opened = apply_media_constraints(base_model, media_path)
            solution = pfba(base_model, fraction_of_optimum=1.0)
            fluxes = solution.fluxes
            row = fluxes.to_dict()
            row["sample"] = sample
            row["total_flux_sum"] = float(solution.objective_value)
            row["growth_rate"] = float(fluxes[objective_reaction_id]) if objective_reaction_id else float("nan")
            row["status"] = solution.status
            rows.append(row)
            print(
                f"{sample}: status={row['status']} growth_rate={row['growth_rate']:.6g} "
                f"total_flux_sum={row['total_flux_sum']:.6g} constrained={n_constrained} media={n_media_opened}",
                flush=True,
            )

    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--human-gem-mat", type=Path, required=True)
    parser.add_argument("--vmax", type=Path, required=True)
    parser.add_argument("--media", type=Path, default=None)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-samples", type=int, default=3)
    parser.add_argument("--solver", choices=["gurobi", "glpk", "auto"], default="gurobi")
    args = parser.parse_args()

    vmax = pd.read_csv(args.vmax, low_memory=False)
    pfba_results = run_pfba_for_samples(
        args.human_gem_mat,
        vmax,
        args.max_samples,
        media_path=args.media,
        solver=args.solver,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    pfba_results.to_csv(args.out, index=False)
    print(f"saved: {args.out}")


if __name__ == "__main__":
    main()
