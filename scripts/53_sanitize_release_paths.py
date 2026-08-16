#!/usr/bin/env python3
"""Remove machine-specific absolute paths from release inventory tables."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TARGETS = {
    ROOT
    / "results"
    / "extension_forward_visibility"
    / "data_inventory"
    / "zju_cellatlas_adult_celltype_contexts.csv": ["h5ad", "gene_summary"],
    ROOT
    / "results"
    / "extension_forward_visibility"
    / "data_inventory"
    / "zju_cellatlas_adult_tissue_contexts.csv": ["h5ad", "gene_summary"],
    ROOT / "results" / "extension_reverse_consistency" / "reverse_extension_unit_inventory.csv": ["metric_file"],
}


def release_path(value: str) -> str:
    normalized = value.replace("\\", "/")
    marker = "/grn_celltype_visibility/"
    if marker in normalized:
        return normalized.split(marker, 1)[1]
    root_text = ROOT.as_posix().rstrip("/") + "/"
    if normalized.startswith(root_text):
        return normalized[len(root_text) :]
    return value


def sanitize(path: Path, columns: list[str]) -> int:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames
    if not fieldnames:
        raise ValueError(f"Missing CSV header: {path}")
    missing = set(columns) - set(fieldnames)
    if missing:
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
    changed = 0
    for row in rows:
        for column in columns:
            revised = release_path(row[column])
            if revised != row[column]:
                row[column] = revised
                changed += 1
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return changed


def main() -> None:
    total = 0
    for path, columns in TARGETS.items():
        changed = sanitize(path, columns)
        total += changed
        print(f"[write] {path.relative_to(ROOT)} ({changed} paths updated)")
    print(f"[done] {total} machine-specific paths updated")


if __name__ == "__main__":
    main()
