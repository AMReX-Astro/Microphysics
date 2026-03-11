#!/usr/bin/env python3

"""
Plot temperature and species fractions versus density from state_over_time.txt using Plotly.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots


def parse_header(path: Path) -> Sequence[str]:
    """Extract column names from the leading comment line."""
    with path.open("r", encoding="ascii") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                return stripped.lstrip("#").split()
            # If no comment header is present, fall back to whitespace split.
            return stripped.split()
    raise ValueError(f"Could not find a header row in {path}")


def load_state_data(path: Path) -> pd.DataFrame:
    columns = parse_header(path)
    # Skip the header we just parsed and load the remainder.
    return pd.read_csv(
        path,
        delim_whitespace=True,
        skiprows=1,
        names=columns,
        comment="#",
    )


def plot_state(path: Path, output: Path | None) -> None:
    data = load_state_data(path)

    required_columns = {"Density", "Temperature"}
    missing = required_columns - set(data.columns)
    if missing:
        missing_cols = ", ".join(sorted(missing))
        raise KeyError(f"Missing required column(s): {missing_cols}")

    density = data["Density"]
    temperature = data["Temperature"]
    species_cols = [
        col for col in data.columns if col not in {"Time", "Density", "Temperature"}
    ]

    if not species_cols:
        raise ValueError("No species fraction columns found in input data.")

    if (density <= 0).any():
        raise ValueError("Density contains non-positive values; cannot use log scale.")
    if (temperature <= 0).any():
        raise ValueError("Temperature contains non-positive values; cannot use log scale.")

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.08,
        subplot_titles=(
            "Temperature vs. Density",
            "Species Fractions vs. Density",
        ),
    )

    fig.add_trace(
        go.Scatter(
            x=density,
            y=temperature,
            mode="lines",
            name="Temperature",
            hovertemplate=(
                "Density: %{x:.3e}<br>Temperature: %{y:.3e}<extra>Temperature</extra>"
            ),
        ),
        row=1,
        col=1,
    )

    for column in species_cols:
        series = data[column]
        positive_mask = series > 0
        if positive_mask.any():
            fig.add_trace(
                go.Scatter(
                    x=density[positive_mask],
                    y=series[positive_mask],
                    mode="lines",
                    name=column,
                    hovertemplate=(
                        "Species: "
                        + column
                        + "<br>Density: %{x:.3e}"
                        + "<br>Fraction: %{y:.3e}<extra></extra>"
                    ),
                ),
                row=2,
                col=1,
            )

    fig.update_xaxes(type="log", row=1, col=1)
    fig.update_xaxes(type="log", title_text="Density", row=2, col=1)
    fig.update_yaxes(type="log", title_text="Temperature", row=1, col=1)
    fig.update_yaxes(type="log", title_text="Species Fraction", row=2, col=1)

    fig.update_layout(
        height=800,
        legend_title_text="Species",
        hovermode="x unified",
        margin=dict(l=70, r=200, t=80, b=60),
    )

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        pio.write_html(fig, file=str(output), include_plotlyjs="cdn", full_html=True)
    else:
        fig.show()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot temperature and species fractions versus density."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("state_over_time.txt"),
        help="Path to the state_over_time.txt file (default: state_over_time.txt).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output path for saving the Plotly figure as HTML.",
    )
    args = parser.parse_args()

    plot_state(args.input, args.output)


if __name__ == "__main__":
    main()
