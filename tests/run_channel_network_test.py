#!/usr/bin/env python3
"""Stationary integration test for a branched network of channel1 elements."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

from plot_channel_profiles import (
    PlotError, render_channel_profiles, render_channel_profiles_pdf,
)


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_DIR = SCRIPT_DIR.parent
DEFAULT_NETWORK = SCRIPT_DIR / "channel_network_merge_split.spr"
DEFAULT_RESULTS = SCRIPT_DIR / "test-results" / "channel-network"


class TestFailure(RuntimeError):
    pass


def text(element: ET.Element, path: str) -> str:
    value = element.findtext(path)
    if value is None:
        raise TestFailure(f"Missing XML value: {path}")
    return value.strip()


def number(element: ET.Element, path: str) -> float:
    return float(text(element, path).replace(",", "."))


def resolve_binary(requested: Path | None) -> Path:
    names = ("staci.exe", "staci") if sys.platform == "win32" else ("staci", "staci.exe")
    candidates = [requested] if requested else [
        REPOSITORY_DIR / "build" / name for name in names
    ] + [
        REPOSITORY_DIR / "build" / configuration / name
        for configuration in ("Release", "Debug") for name in names
    ]
    for candidate in candidates:
        if candidate is not None and candidate.expanduser().resolve().is_file():
            return candidate.expanduser().resolve()
    raise TestFailure("STACI executable not found; build it or pass --binary PATH")


def channel_edges(root: ET.Element) -> list[ET.Element]:
    return [
        edge for edge in root.findall("./edges/edge")
        if text(edge, "pipe_type") == "channel1"
    ]


def find_merge_split_junction(channels: list[ET.Element]) -> tuple[str, int, int]:
    nodes = {text(edge, "node_from") for edge in channels}
    nodes.update(text(edge, "node_to") for edge in channels)
    for node in sorted(nodes):
        incoming = sum(text(edge, "node_to") == node for edge in channels)
        outgoing = sum(text(edge, "node_from") == node for edge in channels)
        if incoming >= 2 and outgoing >= 2:
            return node, incoming, outgoing
    raise TestFailure("No channel junction has at least two incoming and two outgoing channels")


def channel_profile(edge: ET.Element) -> list[tuple[float, float]]:
    for curve in edge.findall("./edge_spec/channel1/curve"):
        if (curve.findtext("id") or "").strip() != "curve_p":
            continue
        points = curve.find("points")
        if points is None:
            break
        xs = [number(item, ".") for item in points.findall("point_x")]
        depths = [number(item, ".") for item in points.findall("point_y")]
        if len(xs) == len(depths) and len(xs) >= 2:
            return list(zip(xs, depths))
    raise TestFailure(f"No calculated curve_p profile in {text(edge, 'id')}")


def flow_direction_bed_slope(edge: ET.Element) -> float:
    """Return positive for a falling bed and negative for an adverse bed."""
    stored_slope = (
        number(edge, "./edge_spec/channel1/start_height")
        - number(edge, "./edge_spec/channel1/end_height")
    ) / number(edge, "./edge_spec/channel1/length")
    return stored_slope if number(edge, "mass_flow_rate") >= 0.0 else -stored_slope


def actual_flow_nodes(edge: ET.Element) -> tuple[str, str]:
    if number(edge, "mass_flow_rate") >= 0.0:
        return text(edge, "node_from"), text(edge, "node_to")
    return text(edge, "node_to"), text(edge, "node_from")


def validate_stationary_result(
    network: Path,
) -> tuple[int, str, float, float, float, float, float]:
    root = ET.parse(network).getroot()
    channels = channel_edges(root)
    if len(channels) < 5:
        raise TestFailure(f"Expected at least five channels, found {len(channels)}")
    junction, incoming_count, outgoing_count = find_merge_split_junction(channels)

    node_heads = {
        text(node, "id"): number(node, "head") for node in root.findall("./nodes/node")
    }
    incoming_flow = 0.0
    outgoing_flow = 0.0
    channels_by_id = {text(edge, "id"): edge for edge in channels}
    shallow_inlet = channels_by_id.get("CHANNEL_1")
    steep_inlet = channels_by_id.get("CHANNEL_2")
    if shallow_inlet is None or steep_inlet is None:
        raise TestFailure("Expected source channels CHANNEL_1 and CHANNEL_2")

    shallow_slope = flow_direction_bed_slope(shallow_inlet)
    steep_slope = flow_direction_bed_slope(steep_inlet)
    if not (0.0 < shallow_slope <= 0.006):
        raise TestFailure(f"CHANNEL_1 is not a low-slope inlet: {shallow_slope}")
    if steep_slope < 0.014:
        raise TestFailure(f"CHANNEL_2 is not a steep inlet: {steep_slope}")
    if steep_slope / shallow_slope < 2.5:
        raise TestFailure("Source-channel slopes are not sufficiently different")
    branch_flows: dict[str, float] = {}
    for edge in channels:
        flow = number(edge, "mass_flow_rate")
        if not math.isfinite(flow) or abs(flow) <= 1.0e-6:
            raise TestFailure(f"Invalid stationary flow in {text(edge, 'id')}: {flow}")
        if text(edge, "node_to") == junction:
            incoming_flow += flow
        if text(edge, "node_from") == junction:
            outgoing_flow += flow
            branch_flows[text(edge, "id")] = flow

        diameter = number(edge, "./edge_spec/channel1/diameter")
        start_depth = node_heads[text(edge, "node_from")] - number(
            edge, "./edge_spec/channel1/start_height"
        )
        end_depth = node_heads[text(edge, "node_to")] - number(
            edge, "./edge_spec/channel1/end_height"
        )
        if not (0.0 < start_depth < diameter and 0.0 < end_depth < diameter):
            raise TestFailure(
                f"{text(edge, 'id')} is not open at both ends: "
                f"depths=({start_depth}, {end_depth}), diameter={diameter}"
            )

    residual = abs(incoming_flow - outgoing_flow)
    tolerance = max(1.0e-3, 1.0e-5 * max(abs(incoming_flow), abs(outgoing_flow)))
    if residual > tolerance:
        raise TestFailure(
            f"Mass balance failed at {junction}: incoming={incoming_flow}, "
            f"outgoing={outgoing_flow}, residual={residual} kg/s"
        )
    adverse_edge = channels_by_id.get("CHANNEL_5")
    smooth_edge = channels_by_id.get("CHANNEL_6")
    if adverse_edge is None or smooth_edge is None:
        raise TestFailure("Expected parallel downstream channels CHANNEL_5 and CHANNEL_6")
    adverse_slope = flow_direction_bed_slope(adverse_edge)
    flow_from, flow_to = actual_flow_nodes(adverse_edge)
    if (flow_from, flow_to) != ("BRANCH_A", "SINK_A"):
        raise TestFailure(
            f"CHANNEL_5 flow direction is {flow_from} -> {flow_to}, expected BRANCH_A -> SINK_A"
        )
    if adverse_slope > -4.0e-4:
        raise TestFailure(
            f"CHANNEL_5 is not sufficiently adverse in the flow direction: {adverse_slope}"
        )
    profile = channel_profile(adverse_edge)
    if not all(math.isfinite(x) and math.isfinite(depth) and depth > 0.0 for x, depth in profile):
        raise TestFailure("CHANNEL_5 has an invalid adverse-slope water profile")
    flow_difference = abs(branch_flows["CHANNEL_3"] - branch_flows["CHANNEL_4"])
    if flow_difference <= 0.05 * max(abs(branch_flows["CHANNEL_3"]), abs(branch_flows["CHANNEL_4"])):
        raise TestFailure("Parallel branch flows are still nearly symmetric")
    return (
        len(channels), junction, incoming_flow, residual, adverse_slope,
        shallow_slope, steep_slope,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--network", type=Path, default=DEFAULT_NETWORK)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--timeout", type=float, default=60.0)
    arguments = parser.parse_args()

    try:
        binary = resolve_binary(arguments.binary)
        source = arguments.network.expanduser().resolve()
        results = arguments.results.expanduser().resolve()
        if results.exists():
            shutil.rmtree(results)
        results.mkdir(parents=True)
        network = results / source.name
        shutil.copy2(source, network)
        completed = subprocess.run(
            [str(binary), "-s", str(network)], cwd=results,
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, encoding="utf-8",
            errors="replace", timeout=arguments.timeout, check=False,
        )
        (results / "console.log").write_text(completed.stdout, encoding="utf-8")
        marker = Path(str(network) + ".rrs")
        if completed.returncode != 0 or not marker.is_file() or marker.read_text().strip() != "OK":
            raise TestFailure(f"STACI did not converge (exit={completed.returncode})")

        (
            count, junction, total_flow, residual, adverse_slope,
            shallow_slope, steep_slope,
        ) = validate_stationary_result(network)
        profile_path = results / "channel-network-longitudinal-profile.svg"
        render_channel_profiles(
            network,
            profile_path,
            "STACI channel network with adverse slope — longitudinal profiles",
        )
        pdf_profile_path = profile_path.with_suffix(".pdf")
        render_channel_profiles_pdf(
            network,
            pdf_profile_path,
            "STACI channel network with adverse slope - longitudinal profiles",
        )
        solved_channels = channel_edges(ET.parse(network).getroot())
        expected_topology_labels = {text(edge, "id") for edge in solved_channels}
        expected_topology_labels.update(text(edge, "node_from") for edge in solved_channels)
        expected_topology_labels.update(text(edge, "node_to") for edge in solved_channels)
        svg_document = profile_path.read_text(encoding="utf-8")
        pdf_document = pdf_profile_path.read_bytes().decode("latin-1", errors="replace")
        for label in expected_topology_labels:
            if label not in svg_document or label not in pdf_document:
                raise TestFailure(f"Topology visualization is missing label: {label}")
        report = (
            "STACI STATIONARY MULTI-CHANNEL NETWORK TEST\n"
            "============================================\n"
            "RESULT: PASS\n"
            f"Channels: {count}\n"
            f"Merge/split junction: {junction} (2 incoming, 2 outgoing)\n"
            f"Incoming and outgoing mass flow: {total_flow:.12g} kg/s\n"
            f"Mass-balance residual: {residual:.12g} kg/s\n"
            "All channel end depths are between zero and the diameter.\n"
            f"Source-channel bed slopes: CHANNEL_1={100.0 * shallow_slope:.6g}%, "
            f"CHANNEL_2={100.0 * steep_slope:.6g}%\n"
            f"Flow-direction bed slope in CHANNEL_5: {100.0 * adverse_slope:.6g}%\n"
            "CHANNEL_5 carries open-channel flow from BRANCH_A to SINK_A "
            "against the rising bed.\n"
            f"Longitudinal profile: {profile_path.name}\n"
            f"Longitudinal profile PDF: {pdf_profile_path.name}\n"
            f"Topology labels verified: {len(expected_topology_labels)}\n"
        )
        (results / "run_channel_network_test.log").write_text(report, encoding="utf-8")
        print(f"Stationary multi-channel test PASS: {count} channels, junction={junction}")
        print(f"Adverse slope PASS: CHANNEL_5={100.0 * adverse_slope:.3g}%")
        print(
            f"Source slopes PASS: CHANNEL_1={100.0 * shallow_slope:.3g}%, "
            f"CHANNEL_2={100.0 * steep_slope:.3g}%"
        )
        print(f"Log: {results / 'run_channel_network_test.log'}")
        print(f"Profile: {profile_path}")
        print(f"Profile PDF: {pdf_profile_path}")
        return 0
    except (
        TestFailure, PlotError, OSError, ValueError, ET.ParseError,
        subprocess.TimeoutExpired,
    ) as error:
        print(f"Stationary multi-channel test FAIL: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
