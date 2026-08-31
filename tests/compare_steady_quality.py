#!/usr/bin/env python3
"""Compare the asymptotic sparse quality solve with STACI's time marcher."""

import argparse
import csv
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys


def run(command, cwd):
    completed = subprocess.run(command, cwd=cwd, text=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError("command failed:\n" + " ".join(map(str, command)) +
                           "\n" + completed.stdout)
    return completed.stdout


def rows_by_id(path, key):
    with path.open(newline="", encoding="utf-8") as stream:
        return {row[key]: row for row in csv.DictReader(stream)}


def last_rows(path, key):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    last_time = max(int(row["time_seconds"]) for row in rows)
    return {row[key]: row for row in rows
            if int(row["time_seconds"]) == last_time}


def relative_error(actual, expected, floor):
    return abs(actual - expected) / max(abs(expected), floor)


def replace_pipe_diameter(source, destination, diameter_mm):
    text = source.read_text(encoding="utf-8")
    replacement = f"P1   R1     J1     1000      {diameter_mm:.6f}"
    text, count = re.subn(
        r"^P1\s+R1\s+J1\s+1000\s+[-+0-9.eE]+", replacement,
        text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError("P1 diameter row was not found")
    destination.write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--staci", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()

    args.work_dir.mkdir(parents=True, exist_ok=True)
    chemical_input = args.work_dir / "steady-quality.inp"
    shutil.copyfile(args.input, chemical_input)
    age_input = args.work_dir / "steady-age.inp"
    age_text = chemical_input.read_text(encoding="utf-8").replace(
        "Quality     CHEMICAL Chlorine mg/L", "Quality     AGE")
    age_input.write_text(age_text, encoding="utf-8")

    steady_prefix = args.work_dir / "steady"
    chemical_prefix = args.work_dir / "time-chemical"
    age_prefix = args.work_dir / "time-age"
    log_parts = []
    log_parts.append(run([
        str(args.staci), "-q", str(chemical_input), "-o", str(steady_prefix),
        "--quality-mode", "both", "--quality-sensitivity",
        "-e", "P1", "-p", "diameter"], args.work_dir))
    log_parts.append(run([
        str(args.staci), "-z", str(chemical_input), "-o", str(chemical_prefix)],
        args.work_dir))
    log_parts.append(run([
        str(args.staci), "-z", str(age_input), "-o", str(age_prefix)],
        args.work_dir))

    steady_nodes = rows_by_id(args.work_dir / "steady-steady-nodes.csv", "node_id")
    steady_links = rows_by_id(args.work_dir / "steady-steady-links.csv", "link_id")
    chemical_nodes = last_rows(args.work_dir / "time-chemical-nodes.csv", "node_id")
    chemical_links = last_rows(args.work_dir / "time-chemical-links.csv", "link_id")
    age_nodes = last_rows(args.work_dir / "time-age-nodes.csv", "node_id")
    age_links = last_rows(args.work_dir / "time-age-links.csv", "link_id")

    checks = []
    for identifier in ("J1", "J2"):
        steady = float(steady_nodes[identifier]["water_age_s"])
        marched = float(age_nodes[identifier]["water_age_s"])
        checks.append((f"node {identifier} water age", steady, marched,
                       relative_error(steady, marched, 60.0), 0.03))
        steady = float(steady_nodes[identifier]["chemical_concentration_kgm3"])
        marched = float(chemical_nodes[identifier]["chlorine_kgm3"])
        checks.append((f"node {identifier} chemical", steady, marched,
                       relative_error(steady, marched, 1.0e-8), 0.03))
    for identifier in ("P1", "P2", "P3"):
        steady = float(steady_links[identifier]["water_age_s"])
        marched = float(age_links[identifier]["water_age_s"])
        checks.append((f"link {identifier} water age", steady, marched,
                       relative_error(steady, marched, 60.0), 0.05))
        steady = float(steady_links[identifier]["chemical_concentration_kgm3"])
        marched = float(chemical_links[identifier]["chlorine_kgm3"])
        checks.append((f"link {identifier} chemical", steady, marched,
                       relative_error(steady, marched, 1.0e-8), 0.05))

    diameter_delta_m = 0.0002
    perturbed = []
    for suffix, diameter_mm in (("plus", 200.2), ("minus", 199.8)):
        inp = args.work_dir / f"diameter-{suffix}.inp"
        prefix = args.work_dir / f"diameter-{suffix}"
        replace_pipe_diameter(chemical_input, inp, diameter_mm)
        log_parts.append(run([
            str(args.staci), "-q", str(inp), "-o", str(prefix),
            "--quality-mode", "both"], args.work_dir))
        perturbed.append(rows_by_id(
            args.work_dir / f"diameter-{suffix}-steady-nodes.csv", "node_id"))

    with (args.work_dir / "steady-steady-sensitivity.csv").open(
            newline="", encoding="utf-8") as stream:
        sensitivity_rows = list(csv.DictReader(stream))
    analytic = {(row["result_quantity"], row["result_id"]):
                float(row["derivative_per_si_parameter"])
                for row in sensitivity_rows}
    for quantity, column in (("node_water_age_s", "water_age_s"),
                             ("node_chemical_concentration_kgm3",
                              "chemical_concentration_kgm3")):
        for identifier in ("J1", "J2"):
            finite = (float(perturbed[0][identifier][column]) -
                      float(perturbed[1][identifier][column])) / (2 * diameter_delta_m)
            value = analytic[(quantity, identifier)]
            checks.append((f"{quantity} sensitivity at {identifier}", value,
                           finite, relative_error(value, finite, 1.0e-10), 0.015))

    failures = [item for item in checks if not math.isfinite(item[3]) or
                item[3] > item[4]]
    report = ["STACI steady water-quality benchmark",
              "steady sparse solution vs. 24-hour time-marching solution",
              "all concentrations are kg/m3; ages are seconds",
              ""]
    for label, steady, reference, error, tolerance in checks:
        status = "PASS" if math.isfinite(error) and error <= tolerance else "FAIL"
        report.append(f"{status:4s} {label:48s} steady={steady:.12g} "
                      f"reference={reference:.12g} rel_error={error:.4g} "
                      f"limit={tolerance:.4g}")
    report.append("")
    report.append(f"RESULT: {'FAIL' if failures else 'PASS'} "
                  f"({len(checks) - len(failures)}/{len(checks)} checks passed)")
    report_text = "\n".join(report) + "\n"
    (args.work_dir / "steady-quality-test.log").write_text(
        report_text + "\nPROGRAM OUTPUT\n" + "\n".join(log_parts),
        encoding="utf-8")
    print(report_text, end="")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
