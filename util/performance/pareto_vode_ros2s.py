#!/usr/bin/env python3
"""Run a VODE/Rosenbrock tolerance sweep using the internal collapse-loop timer."""

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
COUNTER_RE = {
    "rhs_evals": re.compile(r"collapse loop rhs evaluations:\s*(\d+)"),
    "jacobian_evals": re.compile(r"collapse loop jacobian evaluations:\s*(\d+)"),
    "internal_substeps": re.compile(r"collapse loop internal substeps:\s*(\d+)"),
}
TIMER_RE = {
    "rhs_walltime_s": re.compile(rf"collapse loop rhs walltime:\s*({FLOAT_RE})\s*s"),
    "jacobian_walltime_s": re.compile(
        rf"collapse loop jacobian walltime:\s*({FLOAT_RE})\s*s"
    ),
    "lu_walltime_s": re.compile(rf"collapse loop lu walltime:\s*({FLOAT_RE})\s*s"),
    "linear_solve_walltime_s": re.compile(
        rf"collapse loop linear solve walltime:\s*({FLOAT_RE})\s*s"
    ),
}
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

ROSENBROCK_TABLEAUX = {
    0: "ROS2S",
    1: "ROS2",
    2: "SanduA",
    3: "SanduB",
    4: "SanduC",
    5: "SanduD",
    6: "RosenbrockEuler",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep rtol_spec for VODE and all Rosenbrock tableaux. "
            "Other tolerances are scaled from their values in the inputs file."
        )
    )
    parser.add_argument("--input", default="inputs_primordial_chem")
    parser.add_argument("--vode-exe", default="./main1d.llvm.VODE.ex")
    parser.add_argument("--ros-exe", default="./main1d.llvm.ROS2S.ex")
    parser.add_argument(
        "--rosenbrock-tableaux",
        default="all",
        help=(
            "Comma-separated Rosenbrock tableau selectors or names, or 'all'. "
            "Supported: 0=ROS2S, 1=ROS2, 2=SanduA, 3=SanduB, 4=SanduC, 5=SanduD."
        ),
    )
    parser.add_argument(
        "--no-vode",
        action="store_true",
        help="Only run Rosenbrock tableaux, omitting the VODE baseline.",
    )
    parser.add_argument("--rtol-min", type=float, default=1.0e-4)
    parser.add_argument("--rtol-max", type=float, default=1.0e-2)
    parser.add_argument("--num-rtols", type=int, default=9)
    parser.add_argument(
        "--rtols",
        help="Comma-separated rtol_spec values. Overrides --rtol-min/--rtol-max/--num-rtols.",
    )
    parser.add_argument(
        "--sweep-key",
        default="integrator.rtol_spec",
        help=(
            "Runtime parameter to sweep. The output column remains rtol_spec for "
            "compatibility with the plotting scripts."
        ),
    )
    parser.add_argument(
        "--extra-override",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Additional runtime parameter override to apply to every run.",
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


def requested_rosenbrock_tableaux(args: argparse.Namespace) -> list[int]:
    text = args.rosenbrock_tableaux.strip()
    if text.lower() == "all":
        return sorted(ROSENBROCK_TABLEAUX)

    name_to_value = {
        name.lower(): value for value, name in ROSENBROCK_TABLEAUX.items()
    }
    values: list[int] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            value = int(item)
        except ValueError:
            value = name_to_value.get(item.lower(), -1)
        if value not in ROSENBROCK_TABLEAUX:
            supported = ", ".join(
                f"{value}={name}" for value, name in sorted(ROSENBROCK_TABLEAUX.items())
            )
            raise ValueError(f"unknown Rosenbrock tableau {item!r}; supported: {supported}")
        values.append(value)

    if not values:
        raise ValueError("--rosenbrock-tableaux did not contain any values")
    return values


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


def parse_extra_overrides(values: list[str]) -> dict[str, float | int]:
    overrides: dict[str, float | int] = {}
    for item in values:
        if "=" not in item:
            raise ValueError(f"--extra-override must be KEY=VALUE, got {item!r}")
        key, value_text = item.split("=", 1)
        key = key.strip()
        value_text = value_text.strip()
        if not key or not value_text:
            raise ValueError(f"--extra-override must be KEY=VALUE, got {item!r}")
        value = parse_float(value_text)
        overrides[key] = int(value) if value.is_integer() else value
    return overrides


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

    for key, pattern in COUNTER_RE.items():
        match = pattern.search(text)
        result[key] = int(match.group(1)) if match is not None else math.nan

    for key, pattern in TIMER_RE.items():
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


def relative_l1_error(values: list[float], reference: list[float]) -> float:
    if not values or len(values) != len(reference):
        return math.nan
    for value, ref in zip(values, reference):
        if not math.isfinite(value) or not math.isfinite(ref):
            return math.nan
    numerator = sum(abs(value - ref) for value, ref in zip(values, reference))
    denominator = sum(abs(ref) for ref in reference)
    if denominator == 0.0:
        return math.nan
    return numerator / denominator


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
    overrides: dict[str, float | int],
) -> list[str]:
    return [
        str(executable),
        str(inputs),
        *[f"{key}={format_override_value(value)}" for key, value in overrides.items()],
    ]


def format_override_value(value: float | int) -> str:
    if isinstance(value, int):
        return str(value)
    return format_float(value)


def run_case(
    *,
    integrator: str,
    executable: Path,
    inputs: Path,
    tolerances: dict[str, float],
    overrides: dict[str, float | int],
    rtol_spec: float,
    rosenbrock_tableau: int | None,
    tableau_name: str,
    repeat: int,
    run_dir: Path,
    timeout: float,
) -> dict[str, object]:
    run_dir.mkdir(parents=True, exist_ok=True)
    cmd = command_for_run(executable, inputs, overrides)
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
        "rosenbrock_tableau": rosenbrock_tableau,
        "tableau_name": tableau_name,
        "repeat": repeat,
        "rtol_spec": rtol_spec,
        "atol_spec": tolerances["integrator.atol_spec"],
        "rtol_enuc": tolerances["integrator.rtol_enuc"],
        "atol_enuc": tolerances["integrator.atol_enuc"],
        "returncode": completed.returncode,
        "successful": int(parsed["successful"]) if completed.returncode == 0 else 0,
        "walltime_s": parsed["walltime_s"],
        "rhs_evals": parsed["rhs_evals"],
        "jacobian_evals": parsed["jacobian_evals"],
        "internal_substeps": parsed["internal_substeps"],
        "rhs_walltime_s": parsed["rhs_walltime_s"],
        "jacobian_walltime_s": parsed["jacobian_walltime_s"],
        "lu_walltime_s": parsed["lu_walltime_s"],
        "linear_solve_walltime_s": parsed["linear_solve_walltime_s"],
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


def numeric_values(rows: list[dict[str, object]], key: str) -> list[float]:
    return [float(row[key]) for row in rows if finite(row.get(key))]


def add_stats(
    summary: dict[str, object],
    prefix: str,
    values: list[float],
) -> None:
    summary[f"{prefix}_median"] = statistics.median(values) if values else math.nan
    summary[f"{prefix}_mean"] = statistics.mean(values) if values else math.nan
    summary[f"{prefix}_min"] = min(values) if values else math.nan
    summary[f"{prefix}_max"] = max(values) if values else math.nan


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
            relative_l1_error(final_vector_values, reference_vector)
            if reference_vector is not None
            else math.nan
        )
        summary = {
            "integrator": integrator,
            "rosenbrock_tableau": representative["rosenbrock_tableau"],
            "tableau_name": representative["tableau_name"],
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
            "error_l1_rel": error,
            "pareto_integrator": 0,
            "pareto_global": 0,
        }
        add_stats(summary, "rhs_evals", numeric_values(successful, "rhs_evals"))
        add_stats(
            summary,
            "jacobian_evals",
            numeric_values(successful, "jacobian_evals"),
        )
        add_stats(
            summary,
            "internal_substeps",
            numeric_values(successful, "internal_substeps"),
        )
        add_stats(
            summary,
            "rhs_walltime_s",
            numeric_values(successful, "rhs_walltime_s"),
        )
        add_stats(
            summary,
            "jacobian_walltime_s",
            numeric_values(successful, "jacobian_walltime_s"),
        )
        add_stats(
            summary,
            "lu_walltime_s",
            numeric_values(successful, "lu_walltime_s"),
        )
        add_stats(
            summary,
            "linear_solve_walltime_s",
            numeric_values(successful, "linear_solve_walltime_s"),
        )
        summary_rows.append(summary)

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
        and finite(row["error_l1_rel"])
    ]
    for row in eligible:
        time = float(row["walltime_median_s"])
        error = float(row["error_l1_rel"])
        dominated = False
        for other in eligible:
            if other is row:
                continue
            other_time = float(other["walltime_median_s"])
            other_error = float(other["error_l1_rel"])
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
    vode_executable = Path(args.vode_exe).resolve()
    ros_executable = Path(args.ros_exe).resolve()

    if not inputs.exists():
        raise FileNotFoundError(inputs)
    if not args.no_vode and not vode_executable.exists():
        raise FileNotFoundError(vode_executable)
    if not ros_executable.exists():
        raise FileNotFoundError(ros_executable)
    if args.repeats <= 0:
        raise ValueError("--repeats must be positive")

    base_tolerances = read_base_tolerances(inputs)
    rtols = requested_rtols(args)
    rosenbrock_tableaux = requested_rosenbrock_tableaux(args)
    base_rtol = base_tolerances["integrator.rtol_spec"]
    extra_overrides = parse_extra_overrides(args.extra_override)

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
        tolerances = dict(base_tolerances)
        if args.sweep_key in TOLERANCE_KEYS:
            scale = rtol_spec / base_rtol
            tolerances = {
                "integrator.rtol_spec": rtol_spec,
                "integrator.atol_spec": base_tolerances["integrator.atol_spec"] * scale,
                "integrator.rtol_enuc": base_tolerances["integrator.rtol_enuc"] * scale,
                "integrator.atol_enuc": base_tolerances["integrator.atol_enuc"] * scale,
            }
        cases: list[dict[str, object]] = []
        if not args.no_vode:
            vode_overrides: dict[str, float | int] = dict(tolerances)
            vode_overrides.update(extra_overrides)
            if args.sweep_key not in TOLERANCE_KEYS:
                vode_overrides[args.sweep_key] = rtol_spec
            cases.append(
                {
                    "integrator": "VODE",
                    "executable": vode_executable,
                    "rosenbrock_tableau": None,
                    "tableau_name": "",
                    "overrides": vode_overrides,
                }
            )
        for tableau in rosenbrock_tableaux:
            tableau_name = ROSENBROCK_TABLEAUX[tableau]
            overrides: dict[str, float | int] = dict(tolerances)
            overrides.update(extra_overrides)
            if args.sweep_key not in TOLERANCE_KEYS:
                overrides[args.sweep_key] = rtol_spec
            overrides["integrator.rosenbrock_tableau"] = tableau
            cases.append(
                {
                    "integrator": tableau_name,
                    "executable": ros_executable,
                    "rosenbrock_tableau": tableau,
                    "tableau_name": tableau_name,
                    "overrides": overrides,
                }
            )

        for case in cases:
            for repeat in range(1, args.repeats + 1):
                run_dir = (
                    log_root
                    / (
                        f"{case['integrator']}_rtol_{slug_float(rtol_spec)}"
                        f"_run_{repeat:02d}"
                    )
                )
                if args.dry_run:
                    dry_commands.append(
                        command_for_run(
                            case["executable"], inputs, case["overrides"]
                        )
                    )
                    continue
                print(
                    f"running {case['integrator']} rtol_spec={rtol_spec:.3e} "
                    f"repeat {repeat}/{args.repeats}"
                )
                run_rows.append(
                    run_case(
                        integrator=case["integrator"],
                        executable=case["executable"],
                        inputs=inputs,
                        tolerances=tolerances,
                        overrides=case["overrides"],
                        rtol_spec=rtol_spec,
                        rosenbrock_tableau=case["rosenbrock_tableau"],
                        tableau_name=case["tableau_name"],
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
        "rosenbrock_tableau",
        "tableau_name",
        "repeat",
        "rtol_spec",
        "atol_spec",
        "rtol_enuc",
        "atol_enuc",
        "returncode",
        "successful",
        "walltime_s",
        "rhs_evals",
        "jacobian_evals",
        "internal_substeps",
        "rhs_walltime_s",
        "jacobian_walltime_s",
        "lu_walltime_s",
        "linear_solve_walltime_s",
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
        "rosenbrock_tableau",
        "tableau_name",
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
        "rhs_evals_median",
        "rhs_evals_mean",
        "rhs_evals_min",
        "rhs_evals_max",
        "jacobian_evals_median",
        "jacobian_evals_mean",
        "jacobian_evals_min",
        "jacobian_evals_max",
        "internal_substeps_median",
        "internal_substeps_mean",
        "internal_substeps_min",
        "internal_substeps_max",
        "rhs_walltime_s_median",
        "rhs_walltime_s_mean",
        "rhs_walltime_s_min",
        "rhs_walltime_s_max",
        "jacobian_walltime_s_median",
        "jacobian_walltime_s_mean",
        "jacobian_walltime_s_min",
        "jacobian_walltime_s_max",
        "lu_walltime_s_median",
        "lu_walltime_s_mean",
        "lu_walltime_s_min",
        "lu_walltime_s_max",
        "linear_solve_walltime_s_median",
        "linear_solve_walltime_s_mean",
        "linear_solve_walltime_s_min",
        "linear_solve_walltime_s_max",
        "error_l1_rel",
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
