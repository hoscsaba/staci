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


def validate_stationary_result(network: Path) -> tuple[int, str, float, float]:
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
    for edge in channels:
        flow = number(edge, "mass_flow_rate")
        if not math.isfinite(flow) or abs(flow) <= 1.0e-6:
            raise TestFailure(f"Invalid stationary flow in {text(edge, 'id')}: {flow}")
        if text(edge, "node_to") == junction:
            incoming_flow += flow
        if text(edge, "node_from") == junction:
            outgoing_flow += flow

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
    return len(channels), junction, incoming_flow, residual


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

        count, junction, total_flow, residual = validate_stationary_result(network)
        report = (
            "STACI STATIONARY MULTI-CHANNEL NETWORK TEST\n"
            "============================================\n"
            "RESULT: PASS\n"
            f"Channels: {count}\n"
            f"Merge/split junction: {junction} (2 incoming, 2 outgoing)\n"
            f"Incoming and outgoing mass flow: {total_flow:.12g} kg/s\n"
            f"Mass-balance residual: {residual:.12g} kg/s\n"
            "All channel end depths are between zero and the diameter.\n"
        )
        (results / "run_channel_network_test.log").write_text(report, encoding="utf-8")
        print(f"Stationary multi-channel test PASS: {count} channels, junction={junction}")
        print(f"Log: {results / 'run_channel_network_test.log'}")
        return 0
    except (TestFailure, OSError, ValueError, ET.ParseError, subprocess.TimeoutExpired) as error:
        print(f"Stationary multi-channel test FAIL: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
