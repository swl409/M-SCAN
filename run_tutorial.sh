#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

CMD=(python -m mscan.cli run-all --config configs/demo_ccle.yaml)

if [[ "${USE_MUTATION:-1}" != "1" ]]; then
  CMD+=(--skip-mutation)
fi
if [[ "${USE_KCAT:-1}" != "1" ]]; then
  CMD+=(--no-kcat)
fi
if [[ "${USE_MEDIA:-1}" != "1" ]]; then
  CMD+=(--no-media)
fi
if [[ "${CHECK_FIGURE2_REFERENCE:-1}" != "1" ]]; then
  CMD+=(--skip-reference-check)
fi
if [[ -n "${RNA_INPUT:-}" ]]; then
  CMD+=(--rna "$RNA_INPUT")
fi
if [[ -n "${MUTATION_TABLE:-}" ]]; then
  CMD+=(--mutation-table "$MUTATION_TABLE")
fi
if [[ -n "${MUTATION_FACTOR:-}" ]]; then
  CMD+=(--mutation-factor "$MUTATION_FACTOR")
fi
if [[ -n "${SOLVER:-}" ]]; then
  CMD+=(--solver "$SOLVER")
fi

"${CMD[@]}"
