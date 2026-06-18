"""Command-line interface for M-SCAN."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from .workflow import load_workflow


def add_common_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--project-dir", type=Path, default=Path.cwd())
    parser.add_argument("--config", type=Path, default=Path("configs/demo_ccle.yaml"))
    parser.add_argument("--rna-mode", choices=["ccle", "gene_matrix"], default=None)
    parser.add_argument("--rna", type=str, default=None)
    parser.add_argument("--mutation-table", type=str, default=None)
    parser.add_argument("--mutation-factor", type=float, default=None)
    parser.add_argument("--solver", choices=["gurobi", "glpk", "auto"], default=None)
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--output-dir", type=str, default=None)
    parser.add_argument("--skip-mutation", action="store_true")
    parser.add_argument("--no-kcat", action="store_true")
    parser.add_argument("--no-media", action="store_true")
    parser.add_argument("--skip-reference-check", action="store_true")
    parser.add_argument("--dry-run", action="store_true")


def build_overrides(args: argparse.Namespace) -> dict[str, Any]:
    overrides: dict[str, Any] = {
        "rna_mode": args.rna_mode,
        "rna": args.rna,
        "mutation_table": args.mutation_table,
        "mutation_factor": args.mutation_factor,
        "solver": args.solver,
        "max_samples": args.max_samples,
        "output_dir": args.output_dir,
    }
    if args.skip_mutation:
        overrides["use_mutation"] = False
    if args.no_kcat:
        overrides["use_kcat"] = False
    if args.no_media:
        overrides["use_media"] = False
    if args.skip_reference_check:
        overrides["check_reference"] = False
    return overrides


def get_workflow(args: argparse.Namespace):
    return load_workflow(
        project_dir=args.project_dir,
        config_path=args.config,
        overrides=build_overrides(args),
        dry_run=args.dry_run,
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="mscan", description="M-SCAN RNA-to-pFBA command line interface.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    commands = {
        "prepare-rna": "Prepare RNA expression for Vmax construction.",
        "apply-mutations": "Apply generic sample-gene mutation attenuation after gene-matrix RNA preparation.",
        "build-vmax": "Build Human-GEM reaction-level Vmax constraints.",
        "run-pfba": "Run Human-GEM pFBA from a Vmax table.",
        "check-reference": "Compare demo outputs with bundled reference files.",
        "run-all": "Run prepare-rna, build-vmax, run-pfba, and optional reference checks.",
    }
    for name, help_text in commands.items():
        sub = subparsers.add_parser(name, help=help_text)
        add_common_options(sub)

    args = parser.parse_args(argv)
    workflow = get_workflow(args)

    if args.command == "prepare-rna":
        workflow.prepare_rna()
    elif args.command == "apply-mutations":
        workflow.apply_mutations()
    elif args.command == "build-vmax":
        workflow.build_vmax()
    elif args.command == "run-pfba":
        workflow.run_pfba()
    elif args.command == "check-reference":
        workflow.check_reference()
    elif args.command == "run-all":
        workflow.run_all()
    else:
        parser.error(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
