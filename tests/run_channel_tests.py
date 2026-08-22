#!/usr/bin/env python3
"""Reference tests for the circular STACI channel element.

The test networks are derived from tests/channel.spr.  Reference values are
computed independently with the circular-segment geometry, Manning's equation,
the critical-flow condition and a separate RK4 implementation of the steady
gradually-varied-flow equation.  Only Python's standard library is required.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from typing import Callable, List, Optional, Sequence, Tuple
import xml.etree.ElementTree as ET


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_DIR = SCRIPT_DIR.parent
DEFAULT_TEMPLATE = SCRIPT_DIR / "channel.spr"
DEFAULT_RESULTS = SCRIPT_DIR / "test-results" / "channel"
G = 9.81
RHO = 1000.0


class ReferenceError(RuntimeError):
    pass


@dataclass(frozen=True)
class ChannelGeometry:
    diameter: float
    length: float
    bed_slope: float
    roughness_mm: float

    @property
    def full_area(self) -> float:
        return math.pi * self.diameter * self.diameter / 4.0

    @property
    def full_hydraulic_radius(self) -> float:
        return self.diameter / 4.0

    @property
    def darcy_factor(self) -> float:
        argument = 14.8 * self.full_hydraulic_radius / (self.roughness_mm / 1000.0)
        return 1.0 / (2.0 * math.log10(argument)) ** 2

    @property
    def manning_n(self) -> float:
        return self.full_hydraulic_radius ** (1.0 / 6.0) * math.sqrt(
            self.darcy_factor / (8.0 * G)
        )

    def section(self, depth: float) -> Tuple[float, float, float]:
        """Return area, free-surface top width and hydraulic radius."""
        if depth <= 0.0:
            raise ReferenceError(f"non-positive reference depth: {depth}")
        if depth >= self.diameter:
            return self.full_area, 0.0, self.full_hydraulic_radius
        radius = self.diameter / 2.0
        theta = math.acos(1.0 - depth / radius)
        wetted_perimeter = 2.0 * radius * theta
        area = radius * radius * (theta - math.sin(theta) * math.cos(theta))
        top_width = 2.0 * radius * math.sin(theta)
        return area, top_width, area / wetted_perimeter

    def froude_squared(self, flow: float, depth: float) -> float:
        area, top_width, _ = self.section(depth)
        return flow * flow * top_width / (G * area ** 3)

    def pressure_integral(self, depth: float) -> float:
        """Return I1=int_0^y (y-z) T(z) dz for the circular section."""
        depth = max(0.0, min(depth, self.diameter))
        intervals = 480
        dz = depth / intervals
        total = 0.0
        for index in range(intervals + 1):
            z = index * dz
            width = 2.0 * math.sqrt(max(0.0, z * (self.diameter - z)))
            value = (depth - z) * width
            weight = 1.0 if index in (0, intervals) else (2.0 if index % 2 == 0 else 4.0)
            total += weight * value
        return total * dz / 3.0

    def momentum(self, flow: float, depth: float) -> float:
        area, _, _ = self.section(depth)
        return flow * flow / (G * area) + self.pressure_integral(depth)

    def conveyance_flow(self, depth: float) -> float:
        area, _, hydraulic_radius = self.section(depth)
        return (
            area
            * hydraulic_radius ** (2.0 / 3.0)
            * math.sqrt(self.bed_slope)
            / self.manning_n
        )

    def diffusive_flow(self, start_surface: float, end_surface: float, depth: float) -> float:
        surface_slope = (start_surface - end_surface) / self.length
        area, _, hydraulic_radius = self.section(depth)
        magnitude = (
            area * hydraulic_radius ** (2.0 / 3.0)
            * math.sqrt(abs(surface_slope)) / self.manning_n
        )
        return math.copysign(magnitude, surface_slope)

    def critical_depth(self, flow: float) -> float:
        target = abs(flow)
        if target == 0.0:
            return 0.0
        function = lambda depth: self.froude_squared(target, depth) - 1.0
        return bisect_root(function, self.diameter * 1.0e-7, self.diameter * (1.0 - 1.0e-8))

    def normal_depths(self, flow: float) -> List[float]:
        target = abs(flow)
        if target == 0.0 or self.bed_slope <= 0.0:
            return []
        return scan_roots(
            lambda depth: self.conveyance_flow(depth) - target,
            self.diameter * 1.0e-7,
            self.diameter * (1.0 - 1.0e-8),
        )

    def gvf_rhs(self, depth: float, flow: float) -> float:
        area, top_width, hydraulic_radius = self.section(depth)
        friction_slope = (
            self.manning_n ** 2
            * flow
            * abs(flow)
            / (area ** 2 * hydraulic_radius ** (4.0 / 3.0))
        )
        numerator = self.bed_slope - friction_slope
        if depth >= self.diameter:
            # This is the same steady mixed-flow closure used by Csatorna.cpp.
            # It is a model-consistency reference, not a pressure-wave model.
            return numerator
        denominator = 1.0 - flow * flow * top_width / (G * area ** 3)
        if abs(denominator) < 1.0e-5:
            raise ReferenceError("reference profile reached the critical singularity")
        return numerator / denominator


@dataclass
class Case:
    name: str
    description: str
    length: float
    start_invert: float
    end_invert: float
    start_depth: float
    end_depth: float
    expected_flow: float
    oracle: str
    profile: Optional[List[Tuple[float, float]]] = None
    hydraulic_jump: bool = False


@dataclass
class Result:
    case: Case
    status: str
    actual_flow: Optional[float]
    flow_error: Optional[float]
    profile_error: Optional[float]
    momentum_error: Optional[float]
    branch: str
    elapsed: float
    message: str


def bisect_root(function: Callable[[float], float], lower: float, upper: float) -> float:
    f_lower = function(lower)
    f_upper = function(upper)
    if not math.isfinite(f_lower) or not math.isfinite(f_upper) or f_lower * f_upper > 0.0:
        raise ReferenceError(f"root is not bracketed: f({lower})={f_lower}, f({upper})={f_upper}")
    for _ in range(100):
        middle = (lower + upper) / 2.0
        f_middle = function(middle)
        if f_lower * f_middle <= 0.0:
            upper = middle
            f_upper = f_middle
        else:
            lower = middle
            f_lower = f_middle
    return (lower + upper) / 2.0


def scan_roots(function: Callable[[float], float], lower: float, upper: float) -> List[float]:
    roots: List[float] = []
    previous_x = lower
    previous_f = function(previous_x)
    for index in range(1, 10001):
        current_x = lower + (upper - lower) * index / 10000.0
        current_f = function(current_x)
        if previous_f == 0.0:
            roots.append(previous_x)
        elif previous_f * current_f < 0.0:
            roots.append(bisect_root(function, previous_x, current_x))
        previous_x, previous_f = current_x, current_f
    return roots


def integrate_profile(
    geometry: ChannelGeometry,
    flow: float,
    start_depth: float,
    x_start: float,
    x_end: float,
    step: float = 0.01,
) -> List[Tuple[float, float]]:
    direction = 1.0 if x_end >= x_start else -1.0
    x = x_start
    depth = start_depth
    result = [(x, depth)]
    while direction * (x_end - x) > 1.0e-12:
        h = direction * min(step, abs(x_end - x))
        k1 = geometry.gvf_rhs(depth, flow)
        k2 = geometry.gvf_rhs(depth + h * k1 / 2.0, flow)
        k3 = geometry.gvf_rhs(depth + h * k2 / 2.0, flow)
        k4 = geometry.gvf_rhs(depth + h * k3, flow)
        next_depth = depth + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        if not math.isfinite(next_depth) or next_depth <= 0.0:
            raise ReferenceError(f"invalid reference depth {next_depth} at x={x + h}")
        x += h
        depth = next_depth
        result.append((x, depth))
    return result


def interpolate(profile: Sequence[Tuple[float, float]], x: float) -> float:
    ordered = profile if profile[0][0] <= profile[-1][0] else list(reversed(profile))
    if x <= ordered[0][0]:
        return ordered[0][1]
    if x >= ordered[-1][0]:
        return ordered[-1][1]
    lo, hi = 0, len(ordered) - 1
    while hi - lo > 1:
        middle = (lo + hi) // 2
        if ordered[middle][0] <= x:
            lo = middle
        else:
            hi = middle
    x0, y0 = ordered[lo]
    x1, y1 = ordered[hi]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def build_cases(base: ChannelGeometry) -> List[Case]:
    cases: List[Case] = []

    flow = 0.5
    critical = base.critical_depth(flow)
    normal = base.normal_depths(flow)[0]
    cases.append(Case(
        "positive-normal",
        "Positive uniform flow at the independently computed Manning normal depth.",
        base.length, 1.0, 0.0, normal, normal, flow, "Manning analytic",
        [(0.0, normal), (base.length, normal)],
    ))

    flow = 0.2
    critical = base.critical_depth(flow)
    normal = base.normal_depths(flow)[0]
    for name, depth in (
        ("positive-below-normal", 0.80 * normal),
        ("positive-between-normal-critical", (normal + critical) / 2.0),
        ("positive-above-critical", 1.50 * critical),
    ):
        case_geometry = (
            ChannelGeometry(base.diameter, 20.0, base.bed_slope, base.roughness_mm)
            if name == "positive-above-critical" else base
        )
        start_invert = case_geometry.length * case_geometry.bed_slope
        profile = integrate_profile(case_geometry, flow, depth, 0.0, case_geometry.length)
        cases.append(Case(
            name,
            f"Positive GVF profile; inlet depth={depth:.9f} m, yn={normal:.9f} m, yc={critical:.9f} m.",
            case_geometry.length, start_invert, 0.0, depth, profile[-1][1], flow,
            "independent RK4 GVF", profile,
        ))

    reverse_flow = -0.2
    reverse_critical = base.critical_depth(reverse_flow)
    reverse_start = 1.20 * reverse_critical
    reverse_profile = integrate_profile(base, reverse_flow, reverse_start, 0.0, base.length)
    cases.append(Case(
        "reverse-above-critical",
        f"Reverse-flow GVF profile above critical depth yc={reverse_critical:.9f} m.",
        base.length, 1.0, 0.0, reverse_start, reverse_profile[-1][1], reverse_flow,
        "independent RK4 GVF", reverse_profile,
    ))

    short_geometry = ChannelGeometry(base.diameter, 10.0, base.bed_slope, base.roughness_mm)
    reverse_start = 0.80 * short_geometry.critical_depth(reverse_flow)
    reverse_profile = integrate_profile(short_geometry, reverse_flow, reverse_start, 0.0, short_geometry.length)
    cases.append(Case(
        "reverse-below-critical",
        "Reverse supercritical profile on a shortened 10 m reach, kept away from the dry boundary.",
        short_geometry.length, 0.1, 0.0, reverse_start, reverse_profile[-1][1], reverse_flow,
        "independent RK4 GVF", reverse_profile,
    ))

    jump_start_depth, jump_end_depth = 0.5, 1.3
    jump_flow = base.diffusive_flow(
        1.0 + jump_start_depth,
        jump_end_depth,
        (jump_start_depth + jump_end_depth) / 2.0,
    )
    cases.append(Case(
        "hydraulic-jump-positive",
        "Positive supercritical-to-subcritical transition; discharge closure is diffusive-wave and jump position is selected by momentum matching.",
        base.length, 1.0, 0.0, jump_start_depth, jump_end_depth, jump_flow,
        "signed diffusive-wave Q plus circular-section momentum balance",
        None, True,
    ))

    flow = 5.0
    mixed_start = 1.05 * base.critical_depth(flow)
    mixed_profile = integrate_profile(base, flow, mixed_start, 0.0, base.length)
    cases.append(Case(
        "mixed-positive",
        "Positive flow from an open inlet to a full downstream end.",
        base.length, 1.0, 0.0, mixed_start, mixed_profile[-1][1], flow,
        "independent RK4 mixed-model consistency", mixed_profile,
    ))

    reverse_flow = -2.0
    mixed_start = 0.8
    mixed_profile = integrate_profile(base, reverse_flow, mixed_start, 0.0, base.length)
    cases.append(Case(
        "mixed-reverse",
        "Reverse flow from an open upstream end to a full downstream end.",
        base.length, 1.0, 0.0, mixed_start, mixed_profile[-1][1], reverse_flow,
        "independent RK4 mixed-model consistency", mixed_profile,
    ))

    full_start_depth, full_end_depth = 2.0, 1.8
    head_difference = (1.0 + full_start_depth) - full_end_depth
    full_flow = math.sqrt(
        2.0 * G * base.diameter * base.full_area ** 2 * head_difference
        / (base.darcy_factor * base.length)
    )
    cases.append(Case(
        "pressurised-positive",
        "Fully pressurised positive flow checked by Darcy-Weisbach.",
        base.length, 1.0, 0.0, full_start_depth, full_end_depth, full_flow,
        "Darcy-Weisbach analytic",
    ))

    full_start_depth, full_end_depth = 1.6, 3.0
    head_difference = (1.0 + full_start_depth) - full_end_depth
    full_flow = -math.sqrt(
        2.0 * G * base.diameter * base.full_area ** 2 * abs(head_difference)
        / (base.darcy_factor * base.length)
    )
    cases.append(Case(
        "pressurised-reverse",
        "Fully pressurised reverse flow checked by Darcy-Weisbach.",
        base.length, 1.0, 0.0, full_start_depth, full_end_depth, full_flow,
        "Darcy-Weisbach analytic",
    ))

    cases.append(Case(
        "hydrostatic-zero-flow",
        "Equal absolute water-surface elevations; analytic discharge is zero.",
        base.length, 1.0, 0.0, 0.4, 1.4, 0.0,
        "hydrostatic analytic", [(0.0, 0.4), (base.length, 1.4)],
    ))
    return cases


def decimal(text: Optional[str]) -> float:
    if text is None:
        raise ReferenceError("missing numeric XML value")
    return float(text.strip().replace(",", "."))


def find_by_id(elements: Sequence[ET.Element], identifier: str) -> ET.Element:
    for element in elements:
        if (element.findtext("id") or "").strip() == identifier:
            return element
    raise ReferenceError(f"element {identifier!r} not found in SPR")


def set_text(parent: ET.Element, path: str, value: object) -> None:
    node = parent.find(path)
    if node is None:
        raise ReferenceError(f"missing template XML path: {path}")
    node.text = str(value)


def write_case_network(template: Path, destination: Path, case: Case) -> None:
    tree = ET.parse(template)
    root = tree.getroot()
    nodes = root.findall("./nodes/node")
    edges = root.findall("./edges/edge")
    node8 = find_by_id(nodes, "NODE8")
    node9 = find_by_id(nodes, "NODE9")
    channel = find_by_id(edges, "CHANNEL10")
    pool16 = find_by_id(edges, "POOL16")
    pool19 = find_by_id(edges, "POOL19")

    start_surface = case.start_invert + case.start_depth
    end_surface = case.end_invert + case.end_depth
    for node, surface in ((node8, start_surface), (node9, end_surface)):
        set_text(node, "height", 0.0)
        set_text(node, "head", surface)
        set_text(node, "pressure", surface * RHO * G)
    set_text(pool16, "./edge_spec/pool/water_level", start_surface)
    set_text(pool19, "./edge_spec/pool/water_level", end_surface)

    set_text(channel, "mass_flow_rate", case.expected_flow * RHO)
    set_text(channel, "./edge_spec/channel1/length", case.length)
    set_text(channel, "./edge_spec/channel1/start_height", case.start_invert)
    set_text(channel, "./edge_spec/channel1/end_height", case.end_invert)
    set_text(channel, "./edge_spec/channel1/integral_steps", 200)
    set_text(root, "./settings/iter_max", 2000)
    set_text(root, "./settings/e_p_max", 1.0e-4)
    set_text(root, "./settings/e_mp_max", 1.0e-3)
    initial_mass_flow = case.expected_flow * RHO
    if abs(initial_mass_flow) < 0.1:
        initial_mass_flow = 0.1 if start_surface >= end_surface else -0.1
    set_text(root, "./settings/mp_init", initial_mass_flow)
    # Keep the first Newton evaluation in the intended open-channel regime.
    # The template's generic 5 m pressure initialization makes this 1.5 m
    # channel appear fully pressurised before the fixed-level boundaries have
    # been applied, which needlessly sends the iteration to a remote root.
    set_text(root, "./settings/p_init", max(start_surface, end_surface))

    for curve in channel.findall("./edge_spec/channel1/curve"):
        points = curve.find("points")
        if points is not None:
            points.clear()
    destination.parent.mkdir(parents=True, exist_ok=True)
    tree.write(destination, encoding="utf-8", xml_declaration=True)


def read_channel_result(network: Path) -> Tuple[float, List[Tuple[float, float]]]:
    root = ET.parse(network).getroot()
    channel = find_by_id(root.findall("./edges/edge"), "CHANNEL10")
    flow = decimal(channel.findtext("mass_flow_rate")) / RHO
    curve = None
    for candidate in channel.findall("./edge_spec/channel1/curve"):
        if (candidate.findtext("id") or "").strip() == "curve_p":
            curve = candidate
            break
    profile: List[Tuple[float, float]] = []
    if curve is not None:
        points = curve.find("points")
        if points is not None:
            xs = [decimal(node.text) for node in points.findall("point_x")]
            ys = [decimal(node.text) for node in points.findall("point_y")]
            profile = list(zip(xs, ys))
    return flow, profile


def parse_branch(log_path: Path) -> str:
    if not log_path.is_file():
        return "unknown"
    text = log_path.read_text(encoding="utf-8", errors="replace")
    matches = re.findall(r"(?:case\s+([^:\s]+)|case\s*:\s*([^\s]+)|\s(teltszelveny):)", text)
    values = [(first or second or third) for first, second, third in matches]
    return values[-1] if values else "unknown"


def compare_profile(actual: Sequence[Tuple[float, float]], expected: Sequence[Tuple[float, float]]) -> float:
    if not actual:
        raise ReferenceError("STACI did not write curve_p points")
    return max(abs(depth - interpolate(expected, x)) for x, depth in actual)


def jump_momentum_error(
    geometry: ChannelGeometry, flow: float, profile: Sequence[Tuple[float, float]]
) -> float:
    for index in range(len(profile) - 1):
        x0, depth0 = profile[index]
        x1, depth1 = profile[index + 1]
        if abs(x1 - x0) <= 1.0e-9 and abs(depth1 - depth0) > 1.0e-6:
            return abs(geometry.momentum(flow, depth0) - geometry.momentum(flow, depth1))
    raise ReferenceError("STACI did not write a discontinuous hydraulic-jump pair")


def run_case(binary: Path, template: Path, root: Path, case: Case, timeout: float) -> Result:
    case_dir = root / case.name
    case_dir.mkdir(parents=True, exist_ok=True)
    network = case_dir / "channel.spr"
    write_case_network(template, network, case)
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [str(binary), "-s", str(network)],
            cwd=case_dir,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
            check=False,
        )
        elapsed = time.monotonic() - started
        (case_dir / "console.log").write_text(completed.stdout, encoding="utf-8")
        marker = Path(str(network) + ".rrs")
        if completed.returncode != 0 or not marker.is_file() or marker.read_text().strip().upper() != "OK":
            raise ReferenceError(f"STACI did not converge (exit={completed.returncode})")
        actual_flow, actual_profile = read_channel_result(network)
        flow_error = abs(actual_flow - case.expected_flow)
        # The network Newton solver stops on its own pressure/flow tolerances;
        # 1e-3 m3/s is still below 0.2% for the smallest non-zero reference Q.
        flow_tolerance = max(1.0e-3, 5.0e-4 * max(1.0, abs(case.expected_flow)))
        profile_error = None
        if case.profile is not None and abs(case.expected_flow) > 1.0e-10:
            profile_error = compare_profile(actual_profile, case.profile)
        profile_ok = profile_error is None or profile_error <= 2.0e-3
        momentum_error = None
        if case.hydraulic_jump:
            momentum_error = jump_momentum_error(
                ChannelGeometry(1.5, case.length,
                                (case.start_invert - case.end_invert) / case.length,
                                5.0),
                actual_flow, actual_profile,
            )
        momentum_ok = momentum_error is None or momentum_error <= 5.0e-3
        passed = flow_error <= flow_tolerance and profile_ok and momentum_ok
        if passed:
            status = "PASS"
            message = "independent reference checks agree"
        else:
            status = "FAIL"
            message = "reference tolerance exceeded"
        return Result(
            case, status, actual_flow, flow_error, profile_error, momentum_error,
            parse_branch(case_dir / "CHANNEL10.out"), elapsed, message,
        )
    except (ReferenceError, subprocess.TimeoutExpired, ET.ParseError, OSError, ValueError) as error:
        elapsed = time.monotonic() - started
        return Result(case, "FAIL", None, None, None, None, "unknown", elapsed, str(error))


def resolve_binary(requested: Optional[Path]) -> Path:
    candidates = [requested] if requested else [
        REPOSITORY_DIR / "build" / ("staci.exe" if sys.platform == "win32" else "staci"),
        REPOSITORY_DIR / "cmake-build-debug" / ("staci.exe" if sys.platform == "win32" else "staci"),
    ]
    for candidate in candidates:
        if candidate is not None and candidate.expanduser().resolve().is_file():
            return candidate.expanduser().resolve()
    raise SystemExit("STACI executable not found; build it or pass --binary PATH")


def write_reports(results_root: Path, geometry: ChannelGeometry, results: Sequence[Result]) -> None:
    csv_path = results_root / "results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "case", "status", "oracle", "expected_flow_m3_s", "actual_flow_m3_s",
            "flow_abs_error_m3_s", "profile_max_error_m", "jump_momentum_abs_error_m3",
            "branch", "elapsed_s", "message",
        ])
        for result in results:
            writer.writerow([
                result.case.name, result.status, result.case.oracle,
                f"{result.case.expected_flow:.12g}",
                "" if result.actual_flow is None else f"{result.actual_flow:.12g}",
                "" if result.flow_error is None else f"{result.flow_error:.12g}",
                "" if result.profile_error is None else f"{result.profile_error:.12g}",
                "" if result.momentum_error is None else f"{result.momentum_error:.12g}",
                result.branch, f"{result.elapsed:.6f}", result.message,
            ])

    log_path = results_root / "run_channel_tests.log"
    with log_path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("STACI CHANNEL REFERENCE TESTS\n")
        stream.write("=============================\n\n")
        stream.write("Independent SI reference parameters\n")
        stream.write(f"  diameter: {geometry.diameter:.9g} m\n")
        stream.write(f"  length: {geometry.length:.9g} m\n")
        stream.write(f"  bed slope: {geometry.bed_slope:.9g} m/m\n")
        stream.write(f"  absolute roughness: {geometry.roughness_mm:.9g} mm\n")
        stream.write(f"  Darcy factor: {geometry.darcy_factor:.12g}\n")
        stream.write(f"  derived Manning n: {geometry.manning_n:.12g} s/m^(1/3)\n\n")
        for result in results:
            case = result.case
            stream.write(f"[{result.status}] {case.name}\n")
            stream.write(f"  {case.description}\n")
            stream.write(f"  oracle: {case.oracle}\n")
            stream.write(f"  depths: start={case.start_depth:.12g} m, end={case.end_depth:.12g} m\n")
            stream.write(f"  expected Q: {case.expected_flow:.12g} m3/s\n")
            stream.write(f"  actual Q: {result.actual_flow if result.actual_flow is not None else 'unavailable'} m3/s\n")
            stream.write(f"  flow error: {result.flow_error if result.flow_error is not None else 'unavailable'} m3/s\n")
            stream.write(f"  profile max error: {result.profile_error if result.profile_error is not None else 'not applicable'} m\n")
            stream.write(f"  jump momentum error: {result.momentum_error if result.momentum_error is not None else 'not applicable'} m3\n")
            stream.write(f"  STACI branch: {result.branch}\n")
            stream.write(f"  note: {result.message}\n\n")
        counts = {status: sum(r.status == status for r in results) for status in ("PASS", "FAIL")}
        stream.write("SUMMARY\n-------\n")
        stream.write(f"PASS: {counts['PASS']}\nFAIL: {counts['FAIL']}\n")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run only the circular channel reference tests.")
    parser.add_argument("--binary", type=Path, help="Path to staci or staci.exe")
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE, help="Base channel.spr")
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS, help="Retained result directory")
    parser.add_argument("--timeout", type=float, default=60.0, help="Timeout per hydraulic case")
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    binary = resolve_binary(arguments.binary)
    template = arguments.template.expanduser().resolve()
    results_root = arguments.results.expanduser().resolve()
    if not template.is_file():
        raise SystemExit(f"channel template not found: {template}")
    if results_root.exists():
        shutil.rmtree(results_root)
    results_root.mkdir(parents=True)

    geometry = ChannelGeometry(1.5, 100.0, 0.01, 5.0)
    cases = build_cases(geometry)
    results: List[Result] = []
    for case in cases:
        print(f"Testing channel: {case.name} ...", flush=True)
        result = run_case(binary, template, results_root, case, arguments.timeout)
        results.append(result)
        print(
            f"  {result.status}: expected Q={case.expected_flow:.8g}, "
            f"actual Q={result.actual_flow if result.actual_flow is not None else 'unavailable'}, "
            f"branch={result.branch}",
            flush=True,
        )
    write_reports(results_root, geometry, results)

    failures = sum(result.status == "FAIL" for result in results)
    passes = len(results) - failures
    print(f"Channel reference tests: {passes} passed, {failures} failed")
    print(f"Log: {results_root / 'run_channel_tests.log'}")
    print(f"CSV: {results_root / 'results.csv'}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
