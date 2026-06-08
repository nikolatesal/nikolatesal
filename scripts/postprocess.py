#!/usr/bin/env python3
"""Post-process droplet-particle jump diagnostics."""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="") as handle:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(handle)]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: postprocess.py diagnostics.csv", file=sys.stderr)
        return 2
    path = Path(sys.argv[1])
    rows = read_rows(path)
    if not rows:
        print("no diagnostics rows found", file=sys.stderr)
        return 1

    max_liquid_z = max(rows, key=lambda row: row["liquid_com_z"])
    max_particle_z = max(rows, key=lambda row: row["particle_z"])
    takeoff_rows = [row for row in rows if row["takeoff"] > 0.5]
    takeoff = takeoff_rows[0] if takeoff_rows else None
    max_lag = max(
        math.dist(
            (row["liquid_com_x"], row["liquid_com_y"], row["liquid_com_z"]),
            (row["particle_x"], row["particle_y"], row["particle_z"]),
        )
        for row in rows
    )

    print(f"rows: {len(rows)}")
    print(f"max liquid COM z: {max_liquid_z['liquid_com_z']:.6g} at step {int(max_liquid_z['step'])}")
    print(f"max particle z: {max_particle_z['particle_z']:.6g} at step {int(max_particle_z['step'])}")
    if takeoff:
        print(f"first take-off: step {int(takeoff['step'])}, time {takeoff['time']:.6g}")
    else:
        print("first take-off: not detected")
    print(f"max particle/liquid COM lag: {max_lag:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
