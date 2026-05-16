#!/usr/bin/env python3
"""Run a VODE/ROS2S tolerance sweep using the internal collapse-loop timer."""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
import subprocess
import sys
from datetime import datetime
from pathlib import Path


FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?"
WALLTIME_RE = re.compile(rf"collapse loop walltime:\s*({FLOAT_RE})\s*s")
SUCCESS_RE = re.compile(r"successful\?\s+([01])")
FINAL_RE = {
    "final_T": re.compile(rf"T final\s*=\s*({FLOAT_RE})"),
    "final_Eint": re.compile(rf"Eint final\s*=\s*({FLOAT_RE})"),
    "final_rho": re.compile(rf"rho final\s*=\s*({FLOAT_RE})"),
}

TOLERANCE_KEYS = (
    "integrator.rtol_spec",
    "integrator.atol_spec",
    "integrator.rtol_enuc",
    "integrator.atol_enuc",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep rtol_spec for the VODE and ROS2S burn_cell executables. "
            "Other tolerances are scaled from their values in the inputs file."
        )
    )
    parser.add_argument("--input", default="inputs_primordial_chem")
    parser.add_argument("--vode-exe", default="./main1d.llvm.VODE.ex")
    parser.add_argument("--ros-exe", default="./main1d.llvm.ROS2S.ex")
    parser.add_argument("--rtol-min", type=float, default=1.0e-4)
    parser.add_argument("--rtol-max", type=float, default=1.0e-2)
    parser.add_argument("--num-rtols", type=int, default=9)
    parser.add_argument(
        "--rtols",
        help="Comma-separated rtol_spec values. Overrides --rtol-min/--rtol-max/--num-rtols.",
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", default="pareto_vode_ros2s.csv")
    parser.add_argument("--runs-output", default="pareto_vode_ros2s_runs.csv")
    parser.add_argument("--log-root", default="pareto_vode_ros2s_logs")
    parser.add_argument(
        "--reference-output",
        help=(
            "Existing burn_cell stdout to use as the accuracy reference. "
            "If omitted, the script tries reference_solution_vode_rtol_1e-7.out, "
            "then reference_solution.out, then the strictest successful VODE run."
        ),
    )
    parser.add_argument(
        "--no-auto-reference",
        action="store_true",
        help="Do not automatically use local reference_solution*.out files.",
    )
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the commands that would run without executing them.",
    )
    return parser.parse_args()


def parse_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def format_float(value: float) -> str:
    return f"{value:.16e}"


def slug_float(value: float) -> str:
    return f"{value:.3e}".replace("+", "").replace("-", "m").replace(".", "p")


def logspace(start: float, stop: float, count: int) -> list[float]:
    if count <= 0:
        raise ValueError("--num-rtols must be positive")
    if start <= 0.0 or stop <= 0.0:
        raise ValueError("--rtol-min and --rtol-max must be positive")
    if count == 1:
        return [start]
    log_start = math.log10(start)
    log_stop = math.log10(stop)
    return [
        10.0 ** (log_start + i * (log_stop - log_start) / (count - 1))
        for i in range(count)
    ]


def requested_rtols(args: argparse.Namespace) -> list[float]:
    if args.rtols:
        values = [parse_float(value.strip()) for value in args.rtols.split(",")]
        if not values:
            raise ValueError("--rtols did not contain any values")
        return values
    return logspace(args.rtol_min, args.rtol_max, args.num_rtols)


def read_base_tolerances(inputs: Path) -> dict[str, float]:
    text = inputs.read_text()
    values: dict[str, float] = {}
    for key in TOLERANCE_KEYS:
        pattern = re.compile(rf"^\s*{re.escape(key)}\s*=\s*({FLOAT_RE})", re.MULTILINE)
        match = pattern.search(text)
        if match is None:
            raise ValueError(f"{inputs} is missing {key}")
        values[key] = parse_float(match.group(1))
    return values


def parse_output(text: str) -> dict[str, object]:
    result: dict[str, object] = {}

    walltime_match = WALLTIME_RE.search(text)
    result["walltime_s"] = (
        parse_float(walltime_match.group(1)) if walltime_match is not None else math.nan
    )

    success_matches = list(SUCCESS_RE.finditer(text))
    result["successful"] = (
        int(success_matches[-1].group(1)) if success_matches else 0
    )

    for key, pattern in FINAL_RE.items():
        match = pattern.search(text)
        result[key] = parse_float(match.group(1)) if match is not None else math.nan

    result["number_densities"] = parse_number_densities(text)
    return result


def parse_number_densities(text: str) -> list[float]:
    values: list[float] = []
    in_block = False
    for line in text.splitlines():
        if "New number densities" in line:
            in_block = True
            continue
        if not in_block:
            continue
        stripped = line.strip()
        if not stripped:
            continue
        try:
            values.append(parse_float(stripped.split()[0]))
        except ValueError:
            if values:
                break
    return values


def final_vector(parsed: dict[str, object]) -> list[float]:
    vector = [
        float(parsed.get("final_T", math.nan)),
        float(parsed.get("final_Eint", math.nan)),
        float(parsed.get("final_rho", math.nan)),
    ]
    vector.extend(float(value) for value in parsed.get("number_densities", []))
    return vector


def symmetric_linf_relative_error(values: list[float], reference: list[float]) -> float:
    if not values or len(values) != len(reference):
        return math.nan
    error = 0.0
    for value, ref in zip(values, reference):
        if not math.isfinite(value) or not math.isfinite(ref):
            return math.nan
        denom = abs(value) + abs(ref) + 1.0e-300
        error = max(error, 2.0 * abs(value - ref) / denom)
    return error


def load_reference(args: argparse.Namespace, cwd: Path) -> tuple[list[float] | None, str]:
    candidates: list[Path] = []
    if args.reference_output:
        candidates.append(Path(args.reference_output))
    elif not args.no_auto_reference:
        candidates.extend(
            [
                cwd / "reference_solution_vode_rtol_1e-7.out",
                cwd / "reference_solution.out",
            ]
        )

    for path in candidates:
        if not path.exists():
            continue
        parsed = parse_output(path.read_text(errors="replace"))
        vector = final_vector(parsed)
        if vector and all(math.isfinite(value) for value in vector):
            return vector, str(path)
    return None, ""


def command_for_run(
    executable: Path,
    inputs: Path,
    tolerances: dict[str, float],
) -> list[str]:
    return [
        str(executable),
        str(inputs),
        *[f"{key}={format_float(value)}" for key, value in tolerances.items()],
    ]


def run_case(
    *,
    integrator: str,
    executable: Path,
    inputs: Path,
    tolerances: dict[str, float],
    rtol_spec: float,
    repeat: int,
    run_dir: Path,
    timeout: float,
) -> dict[str, object]:
    run_dir.mkdir(parents=True, exist_ok=True)
    cmd = command_for_run(executable, inputs, tolerances)
    (run_dir / "command.txt").write_text(" ".join(cmd) + "\n")

    completed = subprocess.run(
        cmd,
        cwd=run_dir,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    (run_dir / "stdout.txt").write_text(completed.stdout)
    (run_dir / "stderr.txt").write_text(completed.stderr)

    parsed = parse_output(completed.stdout + "\n" + completed.stderr)
    row: dict[str, object] = {
        "integrator": integrator,
        "repeat": repeat,
        "rtol_spec": rtol_spec,
        "atol_spec": tolerances["integrator.atol_spec"],
        "rtol_enuc": tolerances["integrator.rtol_enuc"],
        "atol_enuc": tolerances["integrator.atol_enuc"],
        "returncode": completed.returncode,
        "successful": int(parsed["successful"]) if completed.returncode == 0 else 0,
        "walltime_s": parsed["walltime_s"],
        "final_T": parsed["final_T"],
        "final_Eint": parsed["final_Eint"],
        "final_rho": parsed["final_rho"],
        "stdout": str(run_dir / "stdout.txt"),
        "stderr": str(run_dir / "stderr.txt"),
        "state_over_time": str(run_dir / "state_over_time.txt"),
    }
    for index, value in enumerate(parsed["number_densities"], start=1):
        row[f"number_density_{index}"] = value
    row["_final_vector"] = final_vector(parsed)
    return row


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def summarize(
    rows: list[dict[str, object]],
    repeats: int,
    reference_vector: list[float] | None,
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, float], list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault((str(row["integrator"]), float(row["rtol_spec"])), []).append(row)

    if reference_vector is None:
        reference_vector = best_sweep_reference(grouped)

    summary_rows: list[dict[str, object]] = []
    for (integrator, rtol_spec), group in sorted(grouped.items()):
        successful = [
            row
            for row in group
            if int(row["successful"]) == 1 and finite(row["walltime_s"])
        ]
        walltimes = [float(row["walltime_s"]) for row in successful]
        representative = successful[0] if successful else group[0]
        final_vector_values = representative.get("_final_vector", [])
        error = (
            symmetric_linf_relative_error(final_vector_values, reference_vector)
            if reference_vector is not None
            else math.nan
        )
        summary_rows.append(
            {
                "integrator": integrator,
                "rtol_spec": rtol_spec,
                "atol_spec": representative["atol_spec"],
                "rtol_enuc": representative["rtol_enuc"],
                "atol_enuc": representative["atol_enuc"],
                "repeats": repeats,
                "successful_runs": len(successful),
                "successful": int(len(successful) == repeats),
                "walltime_median_s": statistics.median(walltimes)
                if walltimes
                else math.nan,
                "walltime_mean_s": statistics.mean(walltimes)
                if walltimes
                else math.nan,
                "walltime_min_s": min(walltimes) if walltimes else math.nan,
                "walltime_max_s": max(walltimes) if walltimes else math.nan,
                "walltime_stdev_s": statistics.stdev(walltimes)
                if len(walltimes) > 1
                else 0.0,
                "error_linf_rel": error,
                "pareto_integrator": 0,
                "pareto_global": 0,
            }
        )

    mark_pareto(summary_rows, "pareto_global")
    for integrator in sorted({row["integrator"] for row in summary_rows}):
        mark_pareto(
            [row for row in summary_rows if row["integrator"] == integrator],
            "pareto_integrator",
        )
    return summary_rows


def best_sweep_reference(
    grouped: dict[tuple[str, float], list[dict[str, object]]]
) -> list[float] | None:
    candidates: list[tuple[int, float, list[float]]] = []
    for (integrator, rtol_spec), group in grouped.items():
        for row in group:
            if int(row["successful"]) != 1:
                continue
            priority = 0 if integrator == "VODE" else 1
            vector = row.get("_final_vector", [])
            if vector:
                candidates.append((priority, rtol_spec, vector))
                break
    if not candidates:
        return None
    candidates.sort(key=lambda item: (item[0], item[1]))
    return candidates[0][2]


def mark_pareto(rows: list[dict[str, object]], column: str) -> None:
    eligible = [
        row
        for row in rows
        if int(row["successful"]) == 1
        and finite(row["walltime_median_s"])
        and finite(row["error_linf_rel"])
    ]
    for row in eligible:
        time = float(row["walltime_median_s"])
        error = float(row["error_linf_rel"])
        dominated = False
        for other in eligible:
            if other is row:
                continue
            other_time = float(other["walltime_median_s"])
            other_error = float(other["error_linf_rel"])
            if (
                other_time <= time
                and other_error <= error
                and (other_time < time or other_error < error)
            ):
                dominated = True
                break
        row[column] = int(not dominated)


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: format_csv_value(row.get(key, ""))
                    for key in fieldnames
                }
            )


def format_csv_value(value: object) -> object:
    if isinstance(value, float):
        return format_float(value) if math.isfinite(value) else "nan"
    return value


def main() -> int:
    args = parse_args()
    cwd = Path.cwd()
    inputs = Path(args.input).resolve()
    executables = {
        "VODE": Path(args.vode_exe).resolve(),
        "ROS2S": Path(args.ros_exe).resolve(),
    }

    if not inputs.exists():
        raise FileNotFoundError(inputs)
    for executable in executables.values():
        if not executable.exists():
            raise FileNotFoundError(executable)
    if args.repeats <= 0:
        raise ValueError("--repeats must be positive")

    base_tolerances = read_base_tolerances(inputs)
    rtols = requested_rtols(args)
    base_rtol = base_tolerances["integrator.rtol_spec"]

    reference_vector, reference_source = load_reference(args, cwd)
    if reference_source:
        print(f"using accuracy reference: {reference_source}")
    else:
        print("using strictest successful sweep run as accuracy reference")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_root = Path(args.log_root).resolve() / timestamp

    dry_commands: list[list[str]] = []
    run_rows: list[dict[str, object]] = []

    for rtol_spec in rtols:
        scale = rtol_spec / base_rtol
        tolerances = {
            "integrator.rtol_spec": rtol_spec,
            "integrator.atol_spec": base_tolerances["integrator.atol_spec"] * scale,
            "integrator.rtol_enuc": base_tolerances["integrator.rtol_enuc"] * scale,
            "integrator.atol_enuc": base_tolerances["integrator.atol_enuc"] * scale,
        }
        for integrator, executable in executables.items():
            for repeat in range(1, args.repeats + 1):
                run_dir = (
                    log_root
                    / f"{integrator}_rtol_{slug_float(rtol_spec)}_run_{repeat:02d}"
                )
                if args.dry_run:
                    dry_commands.append(command_for_run(executable, inputs, tolerances))
                    continue
                print(
                    f"running {integrator} rtol_spec={rtol_spec:.3e} "
                    f"repeat {repeat}/{args.repeats}"
                )
                run_rows.append(
                    run_case(
                        integrator=integrator,
                        executable=executable,
                        inputs=inputs,
                        tolerances=tolerances,
                        rtol_spec=rtol_spec,
                        repeat=repeat,
                        run_dir=run_dir,
                        timeout=args.timeout,
                    )
                )

    if args.dry_run:
        for cmd in dry_commands:
            print(" ".join(cmd))
        return 0

    run_fieldnames = [
        "integrator",
        "repeat",
        "rtol_spec",
        "atol_spec",
        "rtol_enuc",
        "atol_enuc",
        "returncode",
        "successful",
        "walltime_s",
        "final_T",
        "final_Eint",
        "final_rho",
        *[f"number_density_{index}" for index in range(1, 15)],
        "stdout",
        "stderr",
        "state_over_time",
    ]
    summary_fieldnames = [
        "integrator",
        "rtol_spec",
        "atol_spec",
        "rtol_enuc",
        "atol_enuc",
        "repeats",
        "successful_runs",
        "successful",
        "walltime_median_s",
        "walltime_mean_s",
        "walltime_min_s",
        "walltime_max_s",
        "walltime_stdev_s",
        "error_linf_rel",
        "pareto_integrator",
        "pareto_global",
    ]

    summary_rows = summarize(run_rows, args.repeats, reference_vector)
    write_csv(Path(args.runs_output), run_rows, run_fieldnames)
    write_csv(Path(args.output), summary_rows, summary_fieldnames)

    print(f"wrote {args.output}")
    print(f"wrote {args.runs_output}")
    print(f"logs are under {log_root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
