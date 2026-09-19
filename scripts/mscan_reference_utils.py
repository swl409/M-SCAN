"""Strict validation shared by the bundled reference checkers."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd


def validate_tolerance(tolerance: float) -> None:
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("Tolerance must be finite and non-negative")


def validate_ids(table: pd.DataFrame, key: str, label: str) -> None:
    if key not in table.columns:
        raise ValueError(f"{label} must contain a {key!r} column")
    if table.empty:
        raise ValueError(f"{label} is empty")
    if table[key].isna().any() or table[key].astype(str).str.strip().eq("").any():
        raise ValueError(f"{label} contains missing {key} values")
    if table[key].duplicated().any():
        raise ValueError(f"{label} contains duplicate {key} values")


def require_coverage(required, available, label: str) -> None:
    missing = sorted(set(required) - set(available))
    if missing:
        raise ValueError(f"Missing {label}: {missing[:10]} (total={len(missing)})")


def finite_numeric(values, label: str):
    if isinstance(values, pd.DataFrame):
        numeric = values.apply(pd.to_numeric, errors="raise")
    else:
        numeric = pd.to_numeric(values, errors="raise")
    if not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError(f"{label} contains non-finite values (NaN or infinity)")
    return numeric
