"""Portable path helpers for optional full-regeneration scripts."""

from __future__ import annotations

import os
from pathlib import Path


def _env_path(name: str, default: Path) -> Path:
    value = os.environ.get(name)
    return Path(value).expanduser().resolve() if value else default.resolve()


PROJECT_ROOT = _env_path("GRN_PROJECT_ROOT", Path(__file__).resolve().parents[2])
MATCHING_ROOT = _env_path("GRN_MATCHING_ROOT", PROJECT_ROOT)
HALLMARK_GMT = _env_path(
    "GRN_HALLMARK_GMT",
    PROJECT_ROOT / "data" / "ground_truth" / "pathway_gene_sets" / "hallmark.gmt",
)


def resolve_project_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (PROJECT_ROOT / path).resolve()


def project_relative(path: str | Path) -> str:
    resolved = Path(path).expanduser().resolve()
    try:
        return resolved.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return str(resolved)
