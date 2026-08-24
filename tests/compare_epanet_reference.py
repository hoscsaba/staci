#!/usr/bin/env python3
"""Compare STACI hydraulics with an official EPANET 2.2 report."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Dict, Iterable, Optional, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_DIR = SCRIPT_DIR.parent
SKIP_CODE = 77
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
RESULT_HEADING = re.compile(r"^\s*(Node|Link) Results(?: at (\d+):(\d+):(\d+) hrs)?:")
STATUS_EVENT = re.compile(
    r"^\s*(\d+):(\d+):(\d+):\s+(?:Pipe|Pump|TCV|Valve|CV)\s+(\S+)\s+"
    r"(?:changed from \S+ to )?(open|closed|active)\s*$",
    re.IGNORECASE,
)
NODE_ROW = re.compile(rf"^\s*(\S+)\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})(?:\s+.*)?$")
LINK_ROW = re.compile(rf"^\s*(\S+)\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})(?:\s+.*)?$")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staci", type=Path, required=True)
    parser.add_argument("--epanet", type=Path)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--head-abs", type=float, default=0.002)
    parser.add_argument("--head-rel", type=float, default=0.0)
    parser.add_argument("--flow-abs", type=float, default=1.0e-9)
    parser.add_argument("--flow-rel", type=float, default=0.0)
    parser.add_argument("--demand-abs", type=float, default=5.0e-9)
    parser.add_argument("--velocity-abs", type=float, default=5.0e-6)
    parser.add_argument("--velocity-rel", type=float, default=0.0)
    parser.add_argument("--headloss-abs", type=float)
    parser.add_argument("--headloss-rel", type=float, default=0.0)
    parser.add_argument(
        "--check-status", action="store_true",
        help="Compare STACI enabled states with EPANET Status Full events.",
    )
    return parser.parse_args()


def epanet_candidates(requested: Optional[Path]) -> Iterable[Path]:
    if requested is not None:
        yield requested.expanduser()
    environment = os.environ.get("EPANET_EXECUTABLE")
    if environment:
        yield Path(environment).expanduser()
    found = shutil.which("runepanet")
    if found:
        yield Path(found)
    names = ("runepanet.exe", "runepanet") if os.name == "nt" else ("runepanet", "runepanet.exe")
    roots = (
        REPOSITORY_DIR / "build" / "epanet-reference" / "build" / "bin",
        REPOSITORY_DIR / "build" / "epanet-reference-build" / "bin",
    )
    for root in roots:
        for configuration in (Path(), Path("Release"), Path("RelWithDebInfo"), Path("Debug")):
            for name in names:
                yield root / configuration / name


def resolve_epanet(requested: Optional[Path]) -> Optional[Path]:
    for candidate in epanet_candidates(requested):
        resolved = candidate.resolve()
        if resolved.is_file():
            return resolved
    return None


def run(command: list[str], cwd: Path) -> str:
    completed = subprocess.run(
        command, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", errors="replace", check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"Command failed ({completed.returncode}): {subprocess.list2cmdline(command)}\n"
            f"{completed.stdout}"
        )
    return completed.stdout


def parse_inp(path: Path) -> Tuple[
        Dict[str, float], Dict[str, Tuple[str, str]], Dict[str, bool]]:
    junction_elevations: Dict[str, float] = {}
    links: Dict[str, Tuple[str, str]] = {}
    initial_status: Dict[str, bool] = {}
    flow_units = "LPS"
    section = ""
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        content = raw.split(";", 1)[0].strip()
        if not content:
            continue
        if content.startswith("[") and content.endswith("]"):
            section = content[1:-1].strip().upper()
            continue
        fields = content.split()
        if section == "JUNCTIONS" and len(fields) >= 2:
            junction_elevations[fields[0]] = float(fields[1])
        elif section in {"PIPES", "PUMPS", "VALVES"} and len(fields) >= 3:
            links[fields[0]] = (fields[1], fields[2])
            initial_status[fields[0]] = not (
                section == "PIPES" and len(fields) > 7 and fields[7].upper() == "CLOSED"
            )
        elif section == "STATUS" and len(fields) >= 2:
            value = fields[1].upper()
            initial_status[fields[0]] = value != "CLOSED" and value != "0"
        elif section == "OPTIONS" and len(fields) >= 2 and fields[0].upper() == "UNITS":
            flow_units = fields[1].upper()
    if flow_units in {"CFS", "GPM", "MGD", "IMGD", "AFD"}:
        junction_elevations = {
            identifier: elevation * 0.3048
            for identifier, elevation in junction_elevations.items()
        }
    return junction_elevations, links, initial_status


def parse_epanet_statuses(path: Path, initial_status: Dict[str, bool],
                          report_times: Iterable[int]) -> Dict[Tuple[int, str], bool]:
    events: list[Tuple[int, str, bool]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = STATUS_EVENT.match(line)
        if not match:
            continue
        seconds = int(match.group(1)) * 3600 + int(match.group(2)) * 60 + int(match.group(3))
        events.append((seconds, match.group(4), match.group(5).lower() != "closed"))
    states = dict(initial_status)
    result: Dict[Tuple[int, str], bool] = {}
    event_index = 0
    for report_time in sorted(set(report_times)):
        while event_index < len(events) and events[event_index][0] <= report_time:
            _, identifier, enabled = events[event_index]
            states[identifier] = enabled
            event_index += 1
        for identifier, enabled in states.items():
            result[(report_time, identifier)] = enabled
    return result


def parse_epanet_report(path: Path) -> Tuple[Dict[Tuple[int, str], dict], Dict[Tuple[int, str], dict]]:
    nodes: Dict[Tuple[int, str], dict] = {}
    links: Dict[Tuple[int, str], dict] = {}
    mode = ""
    time_seconds = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        heading = RESULT_HEADING.match(line)
        if heading:
            mode = heading.group(1).lower()
            if heading.group(2) is not None:
                time_seconds = int(heading.group(2)) * 3600 + int(heading.group(3)) * 60 + int(heading.group(4))
            else:
                time_seconds = 0
            continue
        identifier = ""
        values: Tuple[float, float, float]
        try:
            # EPANET's text report uses a 15-character ID column followed by
            # three adjacent 10-character numeric fields. At high precision
            # large values can touch, so whitespace tokenization is unsafe.
            identifier = line[2:17].strip()
            values = tuple(float(line[start:start + 10]) for start in (17, 27, 37))  # type: ignore[assignment]
        except (ValueError, IndexError):
            match = NODE_ROW.match(line) if mode == "node" else LINK_ROW.match(line) if mode == "link" else None
            if not match:
                continue
            identifier = match.group(1)
            values = tuple(float(match.group(index)) for index in range(2, 5))  # type: ignore[assignment]
        if not identifier:
            continue
        if mode == "node":
            nodes[(time_seconds, identifier)] = {
                "demand_m3s": values[0] / 1000.0,
                "head_m": values[1],
                "pressure_m": values[2],
            }
        else:
            links[(time_seconds, identifier)] = {
                "flow_m3s": values[0] / 1000.0,
                "velocity_mps": values[1],
            }
    if not nodes or not links:
        raise RuntimeError(f"Could not parse node/link results from {path}")
    return nodes, links


def read_csv(path: Path, id_column: str) -> Dict[Tuple[int, str], dict]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"No result rows in {path}")
    return {(int(row["time_seconds"]), row[id_column]): row for row in rows}


def main() -> int:
    args = parse_arguments()
    epanet = resolve_epanet(args.epanet)
    if epanet is None:
        print(
            "SKIP: official EPANET executable not found. Run "
            "'python3 tests/setup_epanet_reference.py' or pass --epanet PATH.",
            file=sys.stderr,
        )
        return SKIP_CODE

    staci = args.staci.expanduser().resolve()
    input_path = args.input.expanduser().resolve()
    if not staci.is_file() or not input_path.is_file():
        raise SystemExit("STACI executable or input file does not exist")
    work = args.work_dir.expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)
    local_input = work / input_path.name
    shutil.copy2(input_path, local_input)
    prefix = work / "staci"
    report = work / "epanet.rpt"
    binary_output = work / "epanet.bin"

    staci_output = run([str(staci), "--epanet-eps", str(local_input), "-o", str(prefix)], work)
    epanet_output = run([str(epanet), str(local_input), str(report), str(binary_output)], work)
    staci_nodes = read_csv(Path(str(prefix) + "-nodes.csv"), "node_id")
    staci_links = read_csv(Path(str(prefix) + "-links.csv"), "link_id")
    epanet_nodes, epanet_links = parse_epanet_report(report)
    junction_elevations, topology, initial_status = parse_inp(input_path)
    epanet_statuses = parse_epanet_statuses(
        report, initial_status, (time_seconds for time_seconds, _ in epanet_links)
    ) if args.check_status else {}

    maxima = {"head_m": 0.0, "pressure_m": 0.0, "demand_m3s": 0.0,
              "flow_m3s": 0.0, "velocity_mps": 0.0, "headloss_m": 0.0,
              "status": 0.0}
    element_maxima: Dict[str, Dict[str, float]] = {}
    comparisons = 0
    failures: list[str] = []

    def compare(label: str, actual: float, expected: float,
                absolute_tolerance: float, relative_tolerance: float = 0.0) -> None:
        nonlocal comparisons
        comparisons += 1
        error = abs(actual - expected)
        tolerance = max(
            absolute_tolerance,
            relative_tolerance * max(abs(actual), abs(expected)),
        )
        metric = label.rsplit("/", 1)[-1]
        maxima[metric] = max(maxima[metric], error)
        label_parts = label.split("/")
        if len(label_parts) == 4:
            element_key = f"{label_parts[0]}/{label_parts[2]}"
            element_metrics = element_maxima.setdefault(element_key, {})
            element_metrics[metric] = max(element_metrics.get(metric, 0.0), error)
        if not math.isfinite(actual) or error > tolerance:
            failures.append(
                f"{label}: STACI={actual:.12g}, EPANET={expected:.12g}, "
                f"abs_error={error:.6g}, tolerance={tolerance:.6g}"
            )

    for key, reference in sorted(epanet_nodes.items()):
        if key not in staci_nodes:
            failures.append(f"missing STACI node result {key}")
            continue
        row = staci_nodes[key]
        time_seconds, identifier = key
        compare(f"node/{time_seconds}/{identifier}/head_m", float(row["total_head_m"]), reference["head_m"], args.head_abs, args.head_rel)
        if identifier in junction_elevations:
            reference_pressure_m = reference["head_m"] - junction_elevations[identifier]
            compare(f"node/{time_seconds}/{identifier}/pressure_m", float(row["pressure_head_m"]), reference_pressure_m, args.head_abs, args.head_rel)
            compare(f"node/{time_seconds}/{identifier}/demand_m3s", float(row["demand_m3s"]), reference["demand_m3s"], args.demand_abs)
        if row.get("converged") != "1":
            failures.append(f"STACI did not converge for node result {key}")

    for key, reference in sorted(epanet_links.items()):
        if key not in staci_links:
            failures.append(f"missing STACI link result {key}")
            continue
        row = staci_links[key]
        time_seconds, identifier = key
        compare(f"link/{time_seconds}/{identifier}/flow_m3s", float(row["flow_m3s"]), reference["flow_m3s"], args.flow_abs, args.flow_rel)
        # EPANET reports velocity magnitude even when its signed flow is
        # negative; STACI retains the link-orientation sign in its CSV.
        if row.get("type") != "PUMP":
            compare(f"link/{time_seconds}/{identifier}/velocity_mps", abs(float(row["velocity_mps"])), reference["velocity_mps"], args.velocity_abs, args.velocity_rel)
        if identifier not in topology:
            failures.append(f"missing input topology for EPANET link {identifier}")
        else:
            node_from, node_to = topology[identifier]
            from_head = epanet_nodes[(time_seconds, node_from)]["head_m"]
            to_head = epanet_nodes[(time_seconds, node_to)]["head_m"]
            compare(
                f"link/{time_seconds}/{identifier}/headloss_m",
                float(row["headloss_m"]), from_head - to_head,
                args.headloss_abs if args.headloss_abs is not None else args.head_abs,
                args.headloss_rel,
            )
        if row.get("converged") != "1":
            failures.append(f"STACI did not converge for link result {key}")
        if args.check_status and key in epanet_statuses:
            compare(
                f"link/{time_seconds}/{identifier}/status",
                1.0 if row.get("status") == "1" else 0.0,
                1.0 if epanet_statuses[key] else 0.0,
                0.0,
            )

    summary = {
        "schema": "STACI EPANET reference comparison v1",
        "input": str(input_path),
        "staci_executable": str(staci),
        "epanet_executable": str(epanet),
        "tolerances": {
            "head": {"absolute_m": args.head_abs, "relative": args.head_rel},
            "flow": {"absolute_m3s": args.flow_abs, "relative": args.flow_rel},
            "demand_absolute_m3s": args.demand_abs,
            "velocity": {"absolute_mps": args.velocity_abs, "relative": args.velocity_rel},
            "headloss": {
                "absolute_m": args.headloss_abs if args.headloss_abs is not None else args.head_abs,
                "relative": args.headloss_rel,
            },
            "status_exact": args.check_status,
        },
        "comparison_count": comparisons,
        "max_absolute_error": maxima,
        "element_max_absolute_error": element_maxima,
        "status": "FAIL" if failures else "PASS",
        "failures": failures,
    }
    (work / "comparison.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (work / "staci.stdout.txt").write_text(staci_output, encoding="utf-8")
    (work / "epanet.stdout.txt").write_text(epanet_output, encoding="utf-8")
    print(f"EPANET reference comparison: {summary['status']} ({comparisons} values)")
    for name, value in maxima.items():
        print(f"  max |delta {name}| = {value:.9g}")
    if failures:
        for failure in failures[:50]:
            print("  FAIL:", failure)
        if len(failures) > 50:
            print(f"  ... {len(failures) - 50} additional failures; see comparison.json")
        return 1
    print(f"  details: {work / 'comparison.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
