#!/usr/bin/env python3
"""Run STACI integration tests for every .spr and .inp file in tests/.

The runner uses only Python's standard library. Each network is copied to an
isolated temporary directory because a steady-state STACI run writes results
back to the input .spr file. A human-readable log is always produced and all
networks are tested even if an earlier test fails.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_DIR = SCRIPT_DIR.parent
DEFAULT_LOG = SCRIPT_DIR / "run_tests.log"
HDF5_SIGNATURE = bytes.fromhex("894844460d0a1a0a")
COMMAND_TIMEOUT_SECONDS = 300.0
SPR_HYDRAULICS_MAX_BYTES = 1_000_000
FULL_SPR_HYDRAULICS = False
MAX_LOGGED_OUTPUT_LINES = 200
STACI_LOG_SUFFIXES = (
    ".ros", ".rps", ".rrs", ".rot", ".rpt", ".rrt", ".roc", ".rpc", ".rrc"
)
HDF5_REQUIRED_DATASETS = (
    "/time",
    "/nodes/id",
    "/nodes/head",
    "/nodes/pressure_head",
    "/nodes/demand",
    "/links/id",
    "/links/flow_rate",
    "/links/velocity",
    "/links/headloss",
    "/links/status",
    "/tanks/id",
    "/tanks/level",
    "/tanks/volume",
    "/tanks/inflow",
    "/simulation/converged",
)
HDF5_DYNAMIC_DATASETS = HDF5_REQUIRED_DATASETS[2:5] + HDF5_REQUIRED_DATASETS[6:10] + HDF5_REQUIRED_DATASETS[11:]
HDF5_EXPECTED_UNITS = {
    "/time": "s",
    "/nodes/head": "m",
    "/nodes/pressure_head": "m",
    "/nodes/demand": "m3/s",
    "/links/flow_rate": "m3/s",
    "/links/velocity": "m/s",
    "/links/headloss": "m",
    "/tanks/level": "m",
    "/tanks/volume": "m3",
    "/tanks/inflow": "m3/s",
}


class TestFailure(RuntimeError):
    """An expected integration-test condition was not satisfied."""


class Log:
    def __init__(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self.stream = path.open("w", encoding="utf-8", newline="\n")

    def close(self) -> None:
        self.stream.close()

    def write(self, text: str = "") -> None:
        print(text, file=self.stream, flush=True)

    def section(self, title: str, character: str = "=") -> None:
        self.write()
        self.write(title)
        self.write(character * len(title))


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Test the compiled STACI executable with every tests/*.spr and tests/*.inp file."
    )
    parser.add_argument(
        "--binary",
        type=Path,
        help="Path to staci or staci.exe. If omitted, common build locations are searched.",
    )
    parser.add_argument(
        "--log",
        type=Path,
        default=DEFAULT_LOG,
        help=f"Human-readable output log (default: {DEFAULT_LOG}).",
    )
    parser.add_argument(
        "--tests-dir",
        type=Path,
        default=SCRIPT_DIR,
        help="Directory searched recursively for .spr and .inp files.",
    )
    parser.add_argument(
        "--keep-workdir",
        action="store_true",
        help="Keep generated files in tests/test-results instead of deleting them.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=COMMAND_TIMEOUT_SECONDS,
        help="Maximum seconds allowed for one STACI command (default: 300).",
    )
    parser.add_argument(
        "--full-spr-hydraulics",
        action="store_true",
        help="Run steady-state hydraulics for large SPR files too (can take several minutes).",
    )
    parser.add_argument(
        "--spr-hydraulics-max-mb",
        type=float,
        default=SPR_HYDRAULICS_MAX_BYTES / 1_000_000.0,
        help="Run hydraulics by default for SPR files no larger than this many MB (default: 1).",
    )
    return parser.parse_args()


def candidate_binaries() -> Iterable[Path]:
    environment_binary = os.environ.get("STACI_BINARY")
    if environment_binary:
        yield Path(environment_binary)
    executable_names = ("staci.exe", "staci") if os.name == "nt" else ("staci", "staci.exe")
    build_directories = (
        REPOSITORY_DIR / "build",
        REPOSITORY_DIR / "cmake-build-debug",
        REPOSITORY_DIR / "cmake-build-release",
        Path.cwd() / "build",
        Path.cwd(),
    )
    configurations = (Path(), Path("Release"), Path("Debug"), Path("RelWithDebInfo"))
    for directory in build_directories:
        for configuration in configurations:
            for name in executable_names:
                yield directory / configuration / name


def resolve_binary(requested: Optional[Path]) -> Path:
    if requested is not None:
        candidate = requested.expanduser().resolve()
        if not candidate.is_file():
            raise TestFailure(f"STACI executable does not exist: {candidate}")
        return candidate
    for candidate in candidate_binaries():
        candidate = candidate.expanduser().resolve()
        if candidate.is_file():
            return candidate
    raise TestFailure(
        "STACI executable was not found. Build the project first or pass "
        "--binary /path/to/staci (the STACI_BINARY environment variable is also supported)."
    )


def discover_networks(tests_dir: Path) -> Tuple[List[Path], List[Path]]:
    tests_dir = tests_dir.resolve()
    generated_root = (tests_dir / "test-results").resolve()
    spr_files = sorted(
        path for path in tests_dir.rglob("*")
        if path.is_file() and path.suffix.lower() == ".spr"
        and generated_root not in path.resolve().parents
    )
    inp_files = sorted(
        path for path in tests_dir.rglob("*")
        if path.is_file() and path.suffix.lower() == ".inp"
        and generated_root not in path.resolve().parents
    )
    if not spr_files and not inp_files:
        raise TestFailure(f"No .spr or .inp files were found under {tests_dir}")
    return spr_files, inp_files


def remove_previous_logs(tests_dir: Path, current_log: Path) -> List[Path]:
    """Remove generated STACI sidecar logs; the report itself is truncated by Log."""
    removed: List[Path] = []
    for path in sorted(tests_dir.rglob("*")):
        if not path.is_file() or path.resolve() == current_log:
            continue
        if path.name.lower().endswith(STACI_LOG_SUFFIXES):
            path.unlink()
            removed.append(path)
    return removed


def format_command(command: Sequence[str]) -> str:
    return subprocess.list2cmdline(list(command))


def log_program_output(log: Log, output: str, heading: str = "Program output") -> None:
    if not output:
        return
    lines = output.rstrip().splitlines()
    log.write(f"  {heading} ({len(lines)} lines):")
    if len(lines) > MAX_LOGGED_OUTPUT_LINES:
        half = MAX_LOGGED_OUTPUT_LINES // 2
        selected = lines[:half]
        selected.append(
            f"... {len(lines) - MAX_LOGGED_OUTPUT_LINES} lines omitted from the log ..."
        )
        selected.extend(lines[-half:])
    else:
        selected = lines
    for line in selected:
        cleaned = line.rstrip()
        log.write(f"    | {cleaned}" if cleaned else "    |")


def run_command(log: Log, command: Sequence[str], cwd: Path, label: str) -> str:
    log.write(f"  Command ({label}): {format_command(command)}")
    started = time.monotonic()
    try:
        completed = subprocess.run(
            list(command),
            cwd=str(cwd),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=COMMAND_TIMEOUT_SECONDS,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        elapsed = time.monotonic() - started
        partial = error.stdout or ""
        if isinstance(partial, bytes):
            partial = partial.decode("utf-8", errors="replace")
        log.write(f"  TIMEOUT after {elapsed:.3f} s")
        log_program_output(log, partial, "Partial program output")
        raise TestFailure(
            f"{label} exceeded the {COMMAND_TIMEOUT_SECONDS:g} s command timeout"
        )
    elapsed = time.monotonic() - started
    log.write(f"  Exit code: {completed.returncode}; elapsed: {elapsed:.3f} s")
    log_program_output(log, completed.stdout)
    if completed.returncode != 0:
        raise TestFailure(f"{label} returned exit code {completed.returncode}")
    return completed.stdout


def require_file(path: Path, description: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise TestFailure(f"Missing or empty {description}: {path}")


def require_text(path: Path, expected: str, description: str) -> None:
    require_file(path, description)
    text = path.read_text(encoding="utf-8", errors="replace")
    if expected not in text:
        raise TestFailure(f"{description} does not contain {expected!r}: {path}")


def check_spr(binary: Path, source: Path, case_dir: Path, log: Log) -> None:
    network = case_dir / source.name
    shutil.copy2(source, network)

    run_hydraulics = FULL_SPR_HYDRAULICS or source.stat().st_size <= SPR_HYDRAULICS_MAX_BYTES
    if run_hydraulics:
        run_command(log, (str(binary), "-s", str(network)), case_dir, "steady-state hydraulics")
        convergence_marker = Path(str(network) + ".rrs")
        require_file(convergence_marker, "hydraulic convergence marker")
        marker = convergence_marker.read_text(encoding="utf-8", errors="replace").strip().upper()
        if marker != "OK":
            raise TestFailure(f"Hydraulic solver did not converge; marker contains {marker!r}")
    else:
        size_mb = source.stat().st_size / 1_000_000.0
        limit_mb = SPR_HYDRAULICS_MAX_BYTES / 1_000_000.0
        log.write(
            f"  Hydraulic solve: SKIPPED for {size_mb:.2f} MB SPR "
            f"(automatic limit {limit_mb:.2f} MB; use --full-spr-hydraulics to force it)"
        )

    exported = case_dir / f"{source.stem}-export.inp"
    run_command(
        log,
        (str(binary), "--export-epanet", str(network), "-o", str(exported)),
        case_dir,
        "export SPR to EPANET",
    )
    require_text(exported, "[END]", "exported EPANET file")
    roundtrip = run_command(log, (str(binary), "-l", str(exported)), case_dir, "reload exported INP")
    if "ready." not in roundtrip.lower():
        raise TestFailure("Exported EPANET network could not be loaded")


def read_summary(path: Path) -> Dict[str, str]:
    require_file(path, "EPS summary CSV")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.reader(stream))
    return {row[0]: row[1] for row in rows[1:] if len(row) >= 2}


def validate_si_csv(
    path: Path, expected_columns: Sequence[str], allow_empty: bool = False
) -> None:
    require_file(path, "EPS SI CSV")
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        header = next(reader, [])
        first_data_row = next(reader, [])
    if header != list(expected_columns):
        raise TestFailure(f"Unexpected CSV header in {path.name}: {header}")
    if not first_data_row and not allow_empty:
        raise TestFailure(f"CSV has no result rows: {path}")
    forbidden = {"time_hours", "demand_m3h", "flow_m3h"}
    if forbidden.intersection(header):
        raise TestFailure(f"Non-SI columns found in {path.name}")


def check_inp(binary: Path, source: Path, case_dir: Path, log: Log) -> None:
    network = case_dir / source.name
    shutil.copy2(source, network)

    roundtrip = case_dir / f"{source.stem}-roundtrip.inp"
    run_command(
        log,
        (str(binary), "--export-epanet", str(network), "-o", str(roundtrip)),
        case_dir,
        "lossless EPANET metadata round trip",
    )
    require_file(roundtrip, "round-trip EPANET file")
    if roundtrip.read_bytes() != network.read_bytes():
        raise TestFailure(
            "EPANET export changed metadata, comments, section ordering, or line endings"
        )

    prefix = case_dir / f"{source.stem}-eps"
    eps_output = run_command(
        log,
        (str(binary), "--epanet-eps", str(network), "-o", str(prefix)),
        case_dir,
        "EPANET extended-period simulation",
    )
    if "epanet eps complete:" not in eps_output.lower():
        raise TestFailure("EPS run did not print its completion marker")

    metadata_path = Path(str(prefix) + ".meta.json")
    require_file(metadata_path, "EPS metadata JSON")
    with metadata_path.open(encoding="utf-8") as stream:
        metadata = json.load(stream)
    if metadata.get("units") != "SI":
        raise TestFailure("EPS metadata does not declare SI units")
    simulation = metadata.get("simulation", {})
    if simulation.get("status") != "complete" or int(simulation.get("frames", 0)) < 1:
        raise TestFailure(f"EPS metadata reports an incomplete run: {simulation}")
    if int(simulation.get("failed_frames", 0)) != 0:
        raise TestFailure(f"EPS contains failed report frames: {simulation['failed_frames']}")

    validate_si_csv(
        Path(str(prefix) + "-nodes.csv"),
        ("time_seconds", "node_id", "elevation_m", "pressure_head_m",
         "total_head_m", "demand_m3s", "converged"),
    )
    validate_si_csv(
        Path(str(prefix) + "-links.csv"),
        ("time_seconds", "link_id", "type", "node_from", "node_to",
         "flow_m3s", "velocity_mps", "headloss_m", "status", "converged"),
    )
    validate_si_csv(
        Path(str(prefix) + "-tanks.csv"),
        ("time_seconds", "tank_id", "level_m", "volume_m3", "inflow_m3s",
         "min_level_m", "max_level_m", "converged"),
        allow_empty=True,
    )
    summary = read_summary(Path(str(prefix) + "-summary.csv"))
    if int(summary.get("failed_states", "-1")) != 0:
        raise TestFailure(f"EPS hydraulic states failed: {summary.get('failed_states')}")


def find_h5dump(binary: Path) -> Optional[Path]:
    candidates: List[Path] = []
    configured = os.environ.get("H5DUMP")
    if configured:
        candidates.append(Path(configured))
    on_path = shutil.which("h5dump")
    if on_path:
        candidates.append(Path(on_path))
    candidates.extend((
        binary.parent / ("h5dump.exe" if os.name == "nt" else "h5dump"),
        Path("/opt/homebrew/bin/h5dump"),
        Path("/opt/homebrew/anaconda3/bin/h5dump"),
        Path("/usr/bin/h5dump"),
        Path("/usr/local/bin/h5dump"),
    ))
    vcpkg_root = os.environ.get("VCPKG_ROOT")
    if vcpkg_root:
        candidates.append(Path(vcpkg_root) / "installed" / "x64-windows" / "tools" / "hdf5" / "h5dump.exe")
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    return None


def inspect_hdf5_with_h5py(hdf5_path: Path, log: Log) -> bool:
    try:
        import h5py  # type: ignore
    except ImportError:
        return False

    with h5py.File(hdf5_path, "r") as hdf5_file:
        for name in HDF5_REQUIRED_DATASETS:
            if name not in hdf5_file:
                raise TestFailure(f"Required HDF5 dataset is missing: {name}")
        time_dataset = hdf5_file["/time"]
        if time_dataset.chunks != (16,) or time_dataset.maxshape != (None,):
            raise TestFailure(
                f"Unexpected /time chunk/max shape: {time_dataset.chunks}, {time_dataset.maxshape}"
            )
        frame_count = time_dataset.shape[0]
        for name in HDF5_DYNAMIC_DATASETS:
            dataset = hdf5_file[name]
            if dataset.shape[0] != frame_count:
                raise TestFailure(f"HDF5 time dimension differs for {name}: {dataset.shape}")
            if dataset.chunks is None or dataset.chunks[0] != 16:
                raise TestFailure(f"HDF5 dataset is not chunked by 16 frames: {name} {dataset.chunks}")
            if dataset.maxshape is None or dataset.maxshape[0] is not None:
                raise TestFailure(f"HDF5 time dimension is not extendible: {name} {dataset.maxshape}")
        for name, expected in HDF5_EXPECTED_UNITS.items():
            actual = hdf5_file[name].attrs.get("unit")
            if isinstance(actual, bytes):
                actual = actual.decode("utf-8", errors="replace")
            if actual != expected:
                raise TestFailure(f"Unexpected SI unit for {name}: {actual!r}, expected {expected!r}")
    log.write("  HDF5 inspector: h5py")
    log.write(f"  Verified datasets: {len(HDF5_REQUIRED_DATASETS)}; chunk time dimension: 16")
    return True


def inspect_hdf5_with_h5dump(
    hdf5_path: Path, binary: Path, case_dir: Path, log: Log
) -> bool:
    h5dump = find_h5dump(binary)
    if h5dump is None:
        return False
    names = run_command(log, (str(h5dump), "-n", str(hdf5_path)), case_dir, "list HDF5 datasets")
    for name in HDF5_REQUIRED_DATASETS:
        if name not in names:
            raise TestFailure(f"Required HDF5 dataset is missing: {name}")
    layout = run_command(
        log,
        (str(h5dump), "-pH", "-d", "/time", "-d", "/nodes/head",
         "-d", "/links/flow_rate", "-d", "/tanks/level", str(hdf5_path)),
        case_dir,
        "inspect HDF5 chunking",
    )
    if len(re.findall(r"CHUNKED\s*\(\s*16(?:\s*,|\s*\))", layout)) < 4:
        raise TestFailure("HDF5 datasets are not chunked with a 16-frame time dimension")
    if layout.count("H5S_UNLIMITED") < 4:
        raise TestFailure("HDF5 datasets do not have an extendible time dimension")
    unit_command = [str(h5dump)]
    for name in HDF5_EXPECTED_UNITS:
        unit_command.extend(("-a", f"{name}/unit"))
    unit_command.append(str(hdf5_path))
    unit_output = run_command(log, unit_command, case_dir, "inspect HDF5 SI units")
    actual_units = re.findall(r'\(0\):\s*"([^"]*)"', unit_output)
    expected_units = list(HDF5_EXPECTED_UNITS.values())
    if actual_units != expected_units:
        raise TestFailure(
            f"Unexpected HDF5 unit attributes: {actual_units}; expected {expected_units}"
        )
    log.write(f"  HDF5 inspector: {h5dump}")
    log.write(f"  Verified datasets: {len(HDF5_REQUIRED_DATASETS)}; chunk time dimension: 16")
    return True


def check_hdf5(binary: Path, source: Path, case_dir: Path, log: Log) -> None:
    network = case_dir / source.name
    shutil.copy2(source, network)
    prefix = case_dir / f"{source.stem}-hdf5"
    run_command(
        log,
        (str(binary), "--epanet-eps", str(network), "-o", str(prefix)),
        case_dir,
        "generate HDF5 EPS output",
    )
    metadata_path = Path(str(prefix) + ".meta.json")
    require_file(metadata_path, "EPS metadata JSON")
    with metadata_path.open(encoding="utf-8") as stream:
        metadata = json.load(stream)
    if not bool(metadata.get("format", {}).get("hdf5_available")):
        raise TestFailure("The tested STACI executable was built without HDF5 support")
    hdf5_path = Path(str(prefix) + ".h5")
    require_file(hdf5_path, "EPS HDF5 result")
    with hdf5_path.open("rb") as stream:
        if stream.read(len(HDF5_SIGNATURE)) != HDF5_SIGNATURE:
            raise TestFailure(f"Invalid HDF5 signature: {hdf5_path}")
    if inspect_hdf5_with_h5py(hdf5_path, log):
        return
    if inspect_hdf5_with_h5dump(hdf5_path, binary, case_dir, log):
        return
    raise TestFailure(
        "Cannot verify HDF5 datasets and chunking: install h5py or h5dump, "
        "or set H5DUMP to the h5dump executable"
    )


def run_case(
    kind: str,
    source: Path,
    binary: Path,
    work_root: Path,
    tests_dir: Path,
    log: Log,
) -> Tuple[bool, float, str]:
    relative = source.relative_to(tests_dir)
    safe_name = "_".join(relative.parts).replace(" ", "_")
    case_dir = work_root / f"{kind.lower()}-{safe_name}"
    case_dir.mkdir(parents=True, exist_ok=True)
    log.section(f"{kind}: {relative}", "-")
    print(f"Testing {kind}: {relative} ...", flush=True)
    started = time.monotonic()
    try:
        if kind == "SPR":
            check_spr(binary, source, case_dir, log)
        elif kind == "INP":
            check_inp(binary, source, case_dir, log)
        else:
            check_hdf5(binary, source, case_dir, log)
        elapsed = time.monotonic() - started
        log.write(f"RESULT: PASS ({elapsed:.3f} s)")
        print(f"  PASS ({elapsed:.3f} s)", flush=True)
        return True, elapsed, ""
    except TestFailure as error:
        elapsed = time.monotonic() - started
        details = f"{type(error).__name__}: {error}"
        log.write(f"RESULT: FAIL ({elapsed:.3f} s)")
        log.write(f"REASON: {details}")
        print(f"  FAIL ({elapsed:.3f} s): {error}", flush=True)
        return False, elapsed, details
    except Exception as error:  # Continue so the log contains every failing input.
        elapsed = time.monotonic() - started
        details = f"{type(error).__name__}: {error}"
        log.write(f"RESULT: FAIL ({elapsed:.3f} s)")
        log.write(f"REASON: {details}")
        log.write("TRACEBACK:")
        for line in traceback.format_exc().rstrip().splitlines():
            log.write(f"  {line}")
        print(f"  FAIL ({elapsed:.3f} s): {error}", flush=True)
        return False, elapsed, details


def make_work_root(
    keep: bool, tests_dir: Path, run_stamp: str
) -> Tuple[Path, Optional[tempfile.TemporaryDirectory]]:
    if keep:
        destination = tests_dir / "test-results" / "runs" / run_stamp
        destination.mkdir(parents=True)
        return destination, None
    temporary = tempfile.TemporaryDirectory(prefix="staci-tests-")
    return Path(temporary.name), temporary


def main() -> int:
    global COMMAND_TIMEOUT_SECONDS, FULL_SPR_HYDRAULICS, SPR_HYDRAULICS_MAX_BYTES
    arguments = parse_arguments()
    if arguments.timeout <= 0:
        print("Test runner error: --timeout must be positive", file=sys.stderr)
        return 2
    COMMAND_TIMEOUT_SECONDS = arguments.timeout
    if arguments.spr_hydraulics_max_mb < 0:
        print("Test runner error: --spr-hydraulics-max-mb cannot be negative", file=sys.stderr)
        return 2
    SPR_HYDRAULICS_MAX_BYTES = int(arguments.spr_hydraulics_max_mb * 1_000_000)
    FULL_SPR_HYDRAULICS = arguments.full_spr_hydraulics
    log_path = arguments.log.expanduser().resolve()
    tests_dir = arguments.tests_dir.expanduser().resolve()
    removed_logs = remove_previous_logs(tests_dir, log_path)
    log = Log(log_path)
    temporary: Optional[tempfile.TemporaryDirectory] = None
    try:
        started_utc = dt.datetime.now(dt.timezone.utc)
        run_stamp = started_utc.strftime("%Y%m%dT%H%M%SZ")
        binary = resolve_binary(arguments.binary)
        spr_files, inp_files = discover_networks(tests_dir)
        work_root, temporary = make_work_root(
            arguments.keep_workdir, tests_dir, run_stamp
        )
        hdf5_work_root = tests_dir / "test-results" / "hdf5" / run_stamp

        log.write("STACI integration test report")
        log.write("=============================")
        log.write(f"Started UTC: {started_utc.isoformat()}")
        log.write(f"Platform: {platform.platform()}")
        log.write(f"Python: {platform.python_version()}")
        log.write(f"Executable: {binary}")
        log.write(f"Tests directory: {tests_dir}")
        log.write(f"Work directory: {work_root}")
        log.write(f"Discovered: {len(spr_files)} SPR, {len(inp_files)} INP")
        log.write(f"Removed previous STACI sidecar logs: {len(removed_logs)}")
        log.write(f"Command timeout: {COMMAND_TIMEOUT_SECONDS:g} s")
        log.write(
            "SPR hydraulics: "
            + ("all files" if FULL_SPR_HYDRAULICS else
               f"files <= {SPR_HYDRAULICS_MAX_BYTES / 1_000_000.0:.2f} MB")
        )

        results: List[Tuple[str, Path, bool, float, str]] = []
        for kind, files in (("SPR", spr_files), ("INP", inp_files)):
            for source in files:
                passed, elapsed, reason = run_case(
                    kind, source, binary, work_root, tests_dir, log
                )
                results.append((kind, source.relative_to(tests_dir), passed, elapsed, reason))

        if inp_files:
            smoke_input = tests_dir / "epanet_eps_smoke.inp"
            hdf5_source = smoke_input if smoke_input in inp_files else inp_files[0]
            passed, elapsed, reason = run_case(
                "HDF5", hdf5_source, binary, hdf5_work_root, tests_dir, log
            )
            results.append(("HDF5", hdf5_source.relative_to(tests_dir), passed, elapsed, reason))
            relative = hdf5_source.relative_to(tests_dir)
            safe_name = "_".join(relative.parts).replace(" ", "_")
            hdf5_result = (
                hdf5_work_root / f"hdf5-{safe_name}" /
                f"{hdf5_source.stem}-hdf5.h5"
            )
            if hdf5_result.is_file():
                log.write(f"Persistent HDF5 result: {hdf5_result}")

        passed_count = sum(1 for result in results if result[2])
        failed_count = len(results) - passed_count
        log.section("SUMMARY", "-")
        log.write(f"Total test cases: {len(results)}")
        log.write(f"Passed: {passed_count}")
        log.write(f"Failed: {failed_count}")
        log.write()
        for kind, path, passed, elapsed, reason in results:
            status = "PASS" if passed else "FAIL"
            suffix = "" if passed else f" -- {reason}"
            log.write(f"[{status}] {kind:5s} {path} ({elapsed:.3f} s){suffix}")
        log.write()
        log.write("OVERALL RESULT: " + ("PASS" if failed_count == 0 else "FAIL"))
        log.write(f"Finished UTC: {dt.datetime.now(dt.timezone.utc).isoformat()}")

        print(f"STACI tests: {passed_count} passed, {failed_count} failed")
        print(f"Log: {log_path}")
        if inp_files and hdf5_result.is_file():
            print(f"Persistent HDF5 result: {hdf5_result}")
        if arguments.keep_workdir:
            print(f"Generated files: {work_root}")
        return 0 if failed_count == 0 else 1
    except Exception as error:
        log.section("FATAL ERROR")
        log.write(f"{type(error).__name__}: {error}")
        for line in traceback.format_exc().rstrip().splitlines():
            log.write(f"  {line}")
        print(f"Test runner error: {error}", file=sys.stderr)
        print(f"Log: {log_path}", file=sys.stderr)
        return 2
    finally:
        log.close()
        if temporary is not None:
            temporary.cleanup()


if __name__ == "__main__":
    raise SystemExit(main())
