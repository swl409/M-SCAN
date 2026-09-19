"""Step-wise M-SCAN workflow command builders."""

from __future__ import annotations

import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .config import as_bool, load_simple_yaml, resolve_path, set_if_present


@dataclass
class Workflow:
    project_dir: Path
    config: dict[str, Any]
    dry_run: bool = False

    @property
    def scripts_dir(self) -> Path:
        return self.project_dir / "scripts"

    @property
    def output_dir(self) -> Path:
        path = resolve_path(self.project_dir, self.config.get("output_dir", "outputs"), required=True)
        assert path is not None
        return path

    @property
    def human_gem_xlsx(self) -> Path:
        path = resolve_path(self.project_dir, self.config.get("human_gem_xlsx"), required=True)
        assert path is not None
        return path

    @property
    def human_gem_mat(self) -> Path:
        path = resolve_path(self.project_dir, self.config.get("human_gem_mat"), required=True)
        assert path is not None
        return path

    @property
    def protein_coding(self) -> Path | None:
        return resolve_path(self.project_dir, self.config.get("protein_coding_genes"), required=False)

    @property
    def rna(self) -> Path:
        path = resolve_path(self.project_dir, self.config.get("rna"), required=True)
        assert path is not None
        return path

    @property
    def mutation_table(self) -> Path | None:
        return resolve_path(self.project_dir, self.config.get("mutation_table"), required=False)

    @property
    def kcat(self) -> Path | None:
        return resolve_path(self.project_dir, self.config.get("kcat"), required=False)

    @property
    def media(self) -> Path | None:
        return resolve_path(self.project_dir, self.config.get("media"), required=False)

    @property
    def rna_preprocessed(self) -> Path:
        path = resolve_path(self.project_dir, self.config.get("rna_preprocessed", self.output_dir / "rna_preprocessed.csv"))
        assert path is not None
        return path

    @property
    def rna_for_vmax(self) -> Path:
        path = resolve_path(self.project_dir, self.config.get("rna_for_vmax", self.output_dir / "rna_for_vmax.csv"))
        assert path is not None
        return path

    @property
    def applied_mutations(self) -> Path:
        path = resolve_path(self.project_dir, self.config.get("applied_mutations", self.output_dir / "applied_mutations.csv"))
        assert path is not None
        return path

    @property
    def vmax(self) -> Path:
        path = resolve_path(self.project_dir, self.config.get("vmax", self.output_dir / "vmax.csv"))
        assert path is not None
        return path

    @property
    def pfba(self) -> Path:
        path = resolve_path(self.project_dir, self.config.get("pfba", self.output_dir / "pfba_results.csv"))
        assert path is not None
        return path

    @property
    def use_mutation(self) -> bool:
        return as_bool(self.config, "use_mutation", True) and self.mutation_table is not None

    @property
    def use_kcat(self) -> bool:
        return as_bool(self.config, "use_kcat", True) and self.kcat is not None

    @property
    def use_media(self) -> bool:
        return as_bool(self.config, "use_media", True) and self.media is not None

    @property
    def mutation_factor(self) -> float:
        return float(self.config.get("mutation_factor", 0.05))

    def run_command(self, command: list[str]) -> None:
        printable = " ".join(shlex.quote(part) for part in command)
        print(f"[run] {printable}", flush=True)
        if self.dry_run:
            return
        subprocess.run(command, check=True)

    def prepare_rna(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        rna_mode = str(self.config.get("rna_mode", "ccle"))
        if rna_mode == "ccle":
            command = [
                sys.executable,
                str(self.scripts_dir / "mscan_prepare_ccle_rna.py"),
                "--rna",
                str(self.rna),
                "--human-gem-xlsx",
                str(self.human_gem_xlsx),
                "--protein-coding-genes",
                str(self.protein_coding),
                "--sample-column",
                str(self.config.get("sample_column", "cell_line_display_name")),
                "--rna-layout",
                str(self.config.get("rna_layout", "auto")),
                "--target-sum",
                str(self.config.get("target_sum", 1_000_000)),
                "--preprocessed-out",
                str(self.rna_preprocessed),
                "--out",
                str(self.rna_for_vmax),
            ]
            if self.use_mutation:
                command.extend(
                    [
                        "--mutation-table",
                        str(self.mutation_table),
                        "--applied-out",
                        str(self.applied_mutations),
                        "--mutation-factor",
                        str(self.mutation_factor),
                        "--severe-factor",
                        str(self.config.get("severe_factor", self.mutation_factor)),
                        "--lof-factor",
                        str(self.config.get("lof_factor", self.mutation_factor)),
                    ]
                )
                if as_bool(self.config, "include_lof", True):
                    command.append("--include-lof")
            self.run_command(command)
            return

        if rna_mode == "gene_matrix":
            prepared = self.rna_preprocessed if self.use_mutation else self.rna_for_vmax
            command = [
                sys.executable,
                str(self.scripts_dir / "mscan_prepare_gene_matrix.py"),
                "--rna",
                str(self.rna),
                "--human-gem-xlsx",
                str(self.human_gem_xlsx),
                "--input-scale",
                str(self.config.get("input_scale", "log2")),
                "--target-sum",
                str(self.config.get("target_sum", 1_000_000)),
                "--out",
                str(prepared),
            ]
            if self.protein_coding is not None:
                command.extend(["--protein-coding-genes", str(self.protein_coding)])
            if not as_bool(self.config, "replace_zero_min", True):
                command.append("--no-replace-zero-min")
            self.run_command(command)
            return

        raise ValueError("rna_mode must be 'ccle' or 'gene_matrix'")

    def apply_mutations(self) -> None:
        if not self.use_mutation:
            print("[skip] mutation table disabled or not provided", flush=True)
            return
        if str(self.config.get("rna_mode", "ccle")) == "ccle":
            print("[skip] CCLE mutation attenuation is handled by prepare-rna for reference compatibility", flush=True)
            return
        command = [
            sys.executable,
            str(self.scripts_dir / "mscan_apply_mutations.py"),
            "--rna",
            str(self.rna_preprocessed),
            "--mutation-table",
            str(self.mutation_table),
            "--default-factor",
            str(self.mutation_factor),
            "--out",
            str(self.rna_for_vmax),
            "--applied-out",
            str(self.applied_mutations),
        ]
        self.run_command(command)

    def build_vmax(self) -> None:
        command = [
            sys.executable,
            str(self.scripts_dir / "mscan_build_vmax.py"),
            "--human-gem-xlsx",
            str(self.human_gem_xlsx),
            "--human-gem-mat",
            str(self.human_gem_mat),
            "--rna",
            str(self.rna_for_vmax),
            "--fallback-kcat",
            str(self.config.get("fallback_kcat", 1.0)),
            "--out",
            str(self.vmax),
        ]
        if self.use_kcat:
            command.extend(["--kcat", str(self.kcat)])
        self.run_command(command)

    def run_pfba(self) -> None:
        command = [
            sys.executable,
            str(self.scripts_dir / "mscan_run_pfba.py"),
            "--human-gem-mat",
            str(self.human_gem_mat),
            "--vmax",
            str(self.vmax),
            "--out",
            str(self.pfba),
            "--max-samples",
            str(int(self.config.get("max_samples", 0))),
            "--solver",
            str(self.config.get("solver", "gurobi")),
        ]
        if self.use_media:
            command.extend(["--media", str(self.media)])
        self.run_command(command)

    def check_reference(self) -> None:
        if not as_bool(self.config, "check_reference", False):
            print("[skip] reference check disabled", flush=True)
            return
        reference_vmax = resolve_path(self.project_dir, self.config.get("reference_vmax"), required=False)
        reference_pfba = resolve_path(self.project_dir, self.config.get("reference_pfba"), required=False)
        reference_reactions = resolve_path(self.project_dir, self.config.get("reference_reactions"), required=False)
        if reference_vmax is not None:
            self.run_command(
                [
                    sys.executable,
                    str(self.scripts_dir / "mscan_check_vmax.py"),
                    "--observed",
                    str(self.vmax),
                    "--reference",
                    str(reference_vmax),
                ]
            )
        if reference_pfba is not None:
            command = [
                sys.executable,
                str(self.scripts_dir / "mscan_check_pfba.py"),
                "--observed",
                str(self.pfba),
                "--reference",
                str(reference_pfba),
            ]
            if reference_reactions is not None:
                command.extend(["--reaction-list", str(reference_reactions)])
            self.run_command(command)

    def run_all(self) -> None:
        self.prepare_rna()
        if str(self.config.get("rna_mode", "ccle")) != "ccle":
            self.apply_mutations()
        self.build_vmax()
        self.run_pfba()
        self.check_reference()
        print(f"M-SCAN pipeline complete. Outputs written under: {self.output_dir}", flush=True)


def load_workflow(
    project_dir: Path,
    config_path: Path,
    overrides: dict[str, Any] | None = None,
    dry_run: bool = False,
) -> Workflow:
    project_dir = project_dir.resolve()
    if not config_path.is_absolute():
        config_path = project_dir / config_path
    config = load_simple_yaml(config_path)
    overrides = overrides or {}
    for key, value in overrides.items():
        set_if_present(config, key, value)
    return Workflow(project_dir=project_dir, config=config, dry_run=dry_run)
