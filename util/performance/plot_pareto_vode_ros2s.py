#!/usr/bin/env python3
"""Plot VODE/ROS2S walltime versus rtol_spec on log-log axes."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot median walltime vs rtol_spec from pareto_vode_ros2s.py output."
    )
    parser.add_argument("csv", nargs="?", default="pareto_vode_ros2s.csv")
    parser.add_argument("-o", "--output", default="pareto_vode_ros2s.png")
    parser.add_argument(
        "--title",
        default="VODE vs ROS2S tolerance sweep",
    )
    parser.add_argument("--xlabel", default="rtol_spec")
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        required = {
            "integrator",
            "rtol_spec",
            "walltime_median_s",
            "error_l1_rel",
            "pareto_global",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing column(s): {', '.join(sorted(missing))}")
        return list(reader)


def main() -> int:
    args = parse_args()
    rows = read_rows(Path(args.csv))

    by_integrator: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if int(row["successful"]) != 1:
            continue
        rtol = float(row["rtol_spec"])
        walltime = float(row["walltime_median_s"])
        if rtol <= 0.0 or walltime <= 0.0:
            continue
        by_integrator[row["integrator"]].append(row)

    fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)

    markers = {"VODE": "o", "ROS2S": "s"}
    for integrator, group in sorted(by_integrator.items()):
        group = sorted(group, key=lambda row: float(row["rtol_spec"]))
        x = [float(row["rtol_spec"]) for row in group]
        y = [float(row["walltime_median_s"]) for row in group]
        ax.plot(
            x,
            y,
            marker=markers.get(integrator, "o"),
            linewidth=1.5,
            markersize=5,
            label=integrator,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel("median walltime [s]")
    ax.set_title(args.title)
    ax.grid(True, which="both", linestyle=":", linewidth=0.7)
    ax.legend()

    fig.savefig(args.output, dpi=200)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
