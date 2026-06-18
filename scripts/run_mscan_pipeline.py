#!/usr/bin/env python3
"""Backward-compatible wrapper for the M-SCAN CLI."""

from __future__ import annotations

import sys

from mscan.cli import main


if __name__ == "__main__":
    # Preserve the previous wrapper behavior: running this script without an
    # explicit subcommand executes the full pipeline.
    argv = sys.argv[1:]
    if not argv or argv[0].startswith("-"):
        argv = ["run-all", *argv]
    main(argv)
