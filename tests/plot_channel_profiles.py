#!/usr/bin/env python3
"""Render flow-oriented longitudinal profiles from a solved STACI SPR file.

Only the Python standard library is used.  The output is SVG so the elevation
and water-level curves remain sharp when the engineering drawing is enlarged.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from html import escape
import math
from pathlib import Path
from typing import Iterable
import xml.etree.ElementTree as ET


RHO = 1000.0
G = 9.81


class PlotError(RuntimeError):
    pass


def _text(element: ET.Element, path: str) -> str:
    value = element.findtext(path)
    if value is None:
        raise PlotError(f"Missing SPR value: {path}")
    return value.strip()


def _number(element: ET.Element, path: str) -> float:
    return float(_text(element, path).replace(",", "."))


@dataclass(frozen=True)
class Channel:
    identifier: str
    node_from: str
    node_to: str
    length: float
    start_invert: float
    end_invert: float
    diameter: float
    mass_flow: float
    depth_profile: tuple[tuple[float, float], ...]

    @property
    def flow_from(self) -> str:
        return self.node_from if self.mass_flow >= 0.0 else self.node_to

    @property
    def flow_to(self) -> str:
        return self.node_to if self.mass_flow >= 0.0 else self.node_from

    def bed(self, stored_x: float) -> float:
        return self.start_invert + (self.end_invert - self.start_invert) * stored_x / self.length

    def oriented_profile(self) -> list[tuple[float, float, float]]:
        """Return flow distance, invert and absolute water-surface elevation."""
        points = [
            (x, self.bed(x), self.bed(x) + depth)
            for x, depth in self.depth_profile
        ]
        if self.mass_flow < 0.0:
            points = [(self.length - x, bed, surface) for x, bed, surface in reversed(points)]
        return points

    def energy_profile(self) -> list[tuple[float, float]]:
        """Return total head along the actual flow direction."""
        discharge = abs(self.mass_flow) / RHO
        result = []
        for distance, bed, surface in self.oriented_profile():
            depth = min(max(surface - bed, 0.0), self.diameter)
            radius = self.diameter / 2.0
            if depth <= 0.0:
                raise PlotError(f"Cannot calculate energy head for dry channel {self.identifier}")
            if depth >= self.diameter:
                area = math.pi * radius * radius
            else:
                theta = 2.0 * math.acos((radius - depth) / radius)
                area = 0.5 * radius * radius * (theta - math.sin(theta))
            velocity = discharge / area
            result.append((distance, surface + velocity * velocity / (2.0 * G)))
        return result

    def flow_boundary_data(self) -> tuple[
        tuple[str, str, float, float],
        tuple[str, str, float, float],
    ]:
        """Return flow-oriented node id, endpoint symbol, invert and depth."""
        profile = self.oriented_profile()
        upstream = profile[0]
        downstream = profile[-1]
        return (
            (self.flow_from, "e", upstream[1], upstream[2] - upstream[1]),
            (self.flow_to, "v", downstream[1], downstream[2] - downstream[1]),
        )


def _topology_endpoint_data(
    channels: list[Channel],
) -> dict[str, list[tuple[str, float, float]]]:
    """Collect distinct flow-oriented endpoint elevations and depths per node."""
    result: dict[str, list[tuple[str, float, float]]] = {}
    for channel in channels:
        for node_id, symbol, invert, depth in channel.flow_boundary_data():
            entries = result.setdefault(node_id, [])
            if not any(
                existing_symbol == symbol
                and abs(existing_invert - invert) <= 1.0e-6
                and abs(existing_depth - depth) <= 1.0e-6
                for existing_symbol, existing_invert, existing_depth in entries
            ):
                entries.append((symbol, invert, depth))
    for entries in result.values():
        entries.sort(key=lambda item: (item[0] != "e", item[0], item[1], item[2]))
    return result


def _curve_points(edge: ET.Element) -> tuple[tuple[float, float], ...]:
    for curve in edge.findall("./edge_spec/channel1/curve"):
        if (curve.findtext("id") or "").strip() != "curve_p":
            continue
        points = curve.find("points")
        if points is None:
            break
        xs = [_number(item, ".") for item in points.findall("point_x")]
        ys = [_number(item, ".") for item in points.findall("point_y")]
        if len(xs) == len(ys) and len(xs) >= 2:
            return tuple(zip(xs, ys))
    return ()


def _read_channels(network: Path) -> list[Channel]:
    root = ET.parse(network).getroot()
    heads = {
        _text(node, "id"): _number(node, "head")
        for node in root.findall("./nodes/node")
    }
    channels: list[Channel] = []
    for edge in root.findall("./edges/edge"):
        if (edge.findtext("pipe_type") or "").strip() != "channel1":
            continue
        length = _number(edge, "./edge_spec/channel1/length")
        if length <= 0.0:
            raise PlotError(f"Non-positive channel length in {_text(edge, 'id')}")
        start = _number(edge, "./edge_spec/channel1/start_height")
        end = _number(edge, "./edge_spec/channel1/end_height")
        node_from = _text(edge, "node_from")
        node_to = _text(edge, "node_to")
        profile = _curve_points(edge)
        if not profile:
            if node_from not in heads or node_to not in heads:
                raise PlotError(f"No calculated water profile in {_text(edge, 'id')}")
            profile = ((0.0, heads[node_from] - start), (length, heads[node_to] - end))
        else:
            # Some legacy STACI save paths swap a negative-flow edge's nodes
            # and invert elevations but retain curve_p in its original (now
            # opposite) orientation.  Normalize the curve to the stored edge
            # direction by matching its end depths to the solved node heads.
            stored_start_depth = heads[node_from] - start
            stored_end_depth = heads[node_to] - end
            direct_error = (
                abs(profile[0][1] - stored_start_depth)
                + abs(profile[-1][1] - stored_end_depth)
            )
            reverse_error = (
                abs(profile[0][1] - stored_end_depth)
                + abs(profile[-1][1] - stored_start_depth)
            )
            if reverse_error + 1.0e-6 < direct_error:
                profile = tuple(
                    (length - x, depth) for x, depth in reversed(profile)
                )
        channels.append(Channel(
            _text(edge, "id"), node_from, node_to, length, start, end,
            _number(edge, "./edge_spec/channel1/diameter"),
            _number(edge, "mass_flow_rate"), profile,
        ))
    if not channels:
        raise PlotError("The SPR file contains no channel1 elements")
    return channels


def _flow_paths(channels: list[Channel]) -> list[list[Channel]]:
    outgoing: dict[str, list[Channel]] = {}
    incoming_count: dict[str, int] = {}
    nodes: set[str] = set()
    for channel in channels:
        outgoing.setdefault(channel.flow_from, []).append(channel)
        incoming_count[channel.flow_to] = incoming_count.get(channel.flow_to, 0) + 1
        nodes.update((channel.flow_from, channel.flow_to))
    sources = sorted(node for node in nodes if incoming_count.get(node, 0) == 0)
    paths: list[list[Channel]] = []

    def visit(node: str, path: list[Channel], used: set[str]) -> None:
        choices = [edge for edge in outgoing.get(node, ()) if edge.identifier not in used]
        if not choices:
            if path:
                paths.append(path)
            return
        for edge in choices:
            visit(edge.flow_to, path + [edge], used | {edge.identifier})

    for source in sources:
        visit(source, [], set())
    if not paths:  # Cyclic or zero-flow topology: still provide useful edge profiles.
        paths = [[channel] for channel in channels]
    if len(paths) > 64:
        raise PlotError(f"Channel topology expands to too many flow paths ({len(paths)})")
    return paths


def _topology_coordinates(network: Path, channels: list[Channel]) -> dict[str, tuple[float, float]]:
    channel_nodes = {
        node_id
        for channel in channels
        for node_id in (channel.node_from, channel.node_to)
    }
    root = ET.parse(network).getroot()
    coordinates: dict[str, tuple[float, float]] = {}
    for node in root.findall("./nodes/node"):
        node_id = (node.findtext("id") or "").strip()
        if node_id not in channel_nodes:
            continue
        try:
            coordinates[node_id] = (_number(node, "xcoord"), _number(node, "ycoord"))
        except (PlotError, ValueError):
            continue
    if len(coordinates) == len(channel_nodes):
        xs = {point[0] for point in coordinates.values()}
        ys = {point[1] for point in coordinates.values()}
        if len(xs) > 1 or len(ys) > 1:
            return coordinates

    # Fallback for SPR files without drawing coordinates: arrange the nodes in
    # flow-direction layers and spread each layer vertically.
    outgoing: dict[str, list[str]] = {}
    indegree = {node_id: 0 for node_id in channel_nodes}
    for channel in channels:
        outgoing.setdefault(channel.flow_from, []).append(channel.flow_to)
        indegree[channel.flow_to] = indegree.get(channel.flow_to, 0) + 1
    levels = {node_id: 0 for node_id, degree in indegree.items() if degree == 0}
    queue = sorted(levels)
    while queue:
        node_id = queue.pop(0)
        for target in outgoing.get(node_id, ()):
            if target not in levels:
                levels[target] = levels[node_id] + 1
                queue.append(target)
    for node_id in channel_nodes:
        levels.setdefault(node_id, 0)
    by_level: dict[int, list[str]] = {}
    for node_id, level in levels.items():
        by_level.setdefault(level, []).append(node_id)
    return {
        node_id: (float(level), float(index))
        for level, node_ids in by_level.items()
        for index, node_id in enumerate(sorted(node_ids))
    }


def _scaled_topology_positions(
    coordinates: dict[str, tuple[float, float]],
    left: float,
    bottom: float,
    width: float,
    height: float,
) -> dict[str, tuple[float, float]]:
    xs = [point[0] for point in coordinates.values()]
    ys = [point[1] for point in coordinates.values()]
    x_span = max(max(xs) - min(xs), 1.0)
    y_span = max(max(ys) - min(ys), 1.0)
    scale = min(width / x_span, height / y_span)
    used_width = (max(xs) - min(xs)) * scale
    used_height = (max(ys) - min(ys)) * scale
    x_offset = left + (width - used_width) / 2.0
    y_offset = bottom + (height - used_height) / 2.0
    return {
        node_id: (
            x_offset + (x - min(xs)) * scale,
            y_offset + (max(ys) - y) * scale,
        )
        for node_id, (x, y) in coordinates.items()
    }


def _polyline(points: Iterable[tuple[float, float]]) -> str:
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in points)


def _prepared_paths(
    channels: list[Channel],
    include_energy_line: bool = False,
) -> tuple[
    list[list[Channel]],
    list[list[tuple[Channel, float, list[tuple[float, float, float]]]]],
    float,
    float,
]:
    paths = _flow_paths(channels)
    path_profiles: list[list[tuple[Channel, float, list[tuple[float, float, float]]]]] = []
    elevations: list[float] = []
    for path in paths:
        distance = 0.0
        prepared = []
        for channel in path:
            profile = channel.oriented_profile()
            prepared.append((channel, distance, profile))
            for _, bed, surface in profile:
                elevations.extend((bed, surface, bed + channel.diameter))
            if include_energy_line:
                elevations.extend(head for _, head in channel.energy_profile())
            distance += channel.length
        path_profiles.append(prepared)
    z_min, z_max = min(elevations), max(elevations)
    padding = max(0.15, 0.08 * max(z_max - z_min, 1.0))
    return paths, path_profiles, z_min - padding, z_max + padding


def render_channel_profiles(
    network: Path,
    output: Path,
    title: str | None = None,
    show_energy_line: bool = False,
    show_boundary_levels: bool = False,
    show_topology_hydraulics: bool = False,
) -> None:
    network = network.expanduser().resolve()
    channels = _read_channels(network)
    paths, path_profiles, z_min, z_max = _prepared_paths(channels, show_energy_line)
    width, panel_height = 1500, 310
    header, topology_height, footer, left, right = 95, 350, 55, 105, 35
    height = header + topology_height + panel_height * len(paths) + footer
    plot_width = width - left - right
    plot_height = panel_height - 82
    document_title = title or f"STACI channel profiles — {network.name}"
    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<defs><marker id="flow-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="8" markerHeight="8" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#075f91"/></marker></defs>',
        '<style>text{font-family:Arial,sans-serif;fill:#24303a}.axis{stroke:#4d5963;stroke-width:1}.grid{stroke:#dce3e8;stroke-width:1}.bed{stroke:#72502b;stroke-width:4;fill:none}.crown{stroke:#8d969d;stroke-width:1.5;stroke-dasharray:7 5;fill:none}.water{fill:#82c9ee;fill-opacity:.62;stroke:none}.surface{stroke:#087fbd;stroke-width:3;fill:none}.energy{stroke:#c23b22;stroke-width:2.5;stroke-dasharray:10 6;fill:none}.boundary-level{stroke:#d7261e;stroke-width:3;fill:#fff}.boundary-label{font-size:13px;font-weight:bold;fill:#b51f18;paint-order:stroke;stroke:#fff;stroke-width:4px;stroke-linejoin:round}.arrow{fill:#075f91}.small{font-size:13px}.label{font-size:15px;font-weight:bold}.title{font-size:24px;font-weight:bold}.topology-edge{stroke:#087fbd;stroke-width:4;fill:none;marker-end:url(#flow-arrow)}.topology-node{fill:#f4a742;stroke:#59431f;stroke-width:2}.topology-label{font-size:15px;font-weight:bold;paint-order:stroke;stroke:#fff;stroke-width:6px;stroke-linejoin:round}</style>',
        f'<text class="title" x="{left}" y="38">{escape(document_title)}</text>',
        f'<text class="small" x="{left}" y="64">Flow is left to right. Elevations and water levels are read from the solved SPR file in SI units.</text>',
        '<rect x="1080" y="27" width="34" height="14" fill="#82c9ee" fill-opacity=".62"/><text class="small" x="1122" y="39">water</text>',
        '<line x1="1210" y1="35" x2="1250" y2="35" class="bed"/><text class="small" x="1260" y="39">channel invert</text>',
    ]
    if show_energy_line:
        svg.append('<line x1="1080" y1="61" x2="1120" y2="61" class="energy"/><text class="small" x="1130" y="65">energy grade line H = z + y + v²/(2g)</text>')

    topology_top = header
    topology_positions = _scaled_topology_positions(
        _topology_coordinates(network, channels),
        left + 85,
        topology_top + 62,
        plot_width - 170,
        topology_height - 135,
    )
    svg.extend([
        f'<rect x="{left}" y="{topology_top}" width="{plot_width}" height="{topology_height-25}" rx="12" fill="#f7fafc" stroke="#cbd5dc"/>',
        f'<text class="label" x="{left+22}" y="{topology_top+32}">Network topology and calculated flow directions</text>',
    ])
    for channel in channels:
        x0, y0 = topology_positions[channel.flow_from]
        x1, y1 = topology_positions[channel.flow_to]
        svg.append(f'<line class="topology-edge" x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y1:.2f}"/>')
        mx, my = (x0 + x1) / 2.0, (y0 + y1) / 2.0
        topology_edge_label = channel.identifier
        if show_topology_hydraulics:
            flow = abs(channel.mass_flow) / RHO
            topology_edge_label += f" · Q={flow:.4g} m³/s"
        svg.append(f'<text class="topology-label" text-anchor="middle" x="{mx:.2f}" y="{my-10:.2f}">{escape(topology_edge_label)}</text>')
    topology_boundary_data = (
        _topology_endpoint_data(channels) if show_topology_hydraulics else {}
    )
    for node_id, (x, y) in topology_positions.items():
        svg.append(f'<circle class="topology-node" cx="{x:.2f}" cy="{y:.2f}" r="9"/>')
        svg.append(f'<text class="topology-label" text-anchor="middle" x="{x:.2f}" y="{y+29:.2f}">{escape(node_id)}</text>')
        for line_index, (symbol, invert, depth) in enumerate(
            topology_boundary_data.get(node_id, ())
        ):
            svg.append(
                f'<text class="small" text-anchor="middle" x="{x:.2f}" '
                f'y="{y+49+18*line_index:.2f}">'
                f'z_{symbol}={invert:.3f} m · h_{symbol}={depth:.3f} m</text>'
            )

    for panel, (path, prepared) in enumerate(zip(paths, path_profiles)):
        top = header + topology_height + panel * panel_height
        bottom = top + plot_height
        total_length = sum(edge.length for edge in path)
        x_scale = plot_width / max(total_length, 1.0)
        y_scale = plot_height / max(z_max - z_min, 1.0e-9)

        def sx(distance: float) -> float:
            return left + distance * x_scale

        def sy(elevation: float) -> float:
            return bottom - (elevation - z_min) * y_scale

        for tick in range(6):
            elevation = z_min + (z_max - z_min) * tick / 5.0
            y = sy(elevation)
            svg.append(f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{width-right}" y2="{y:.2f}"/>')
            svg.append(f'<text class="small" text-anchor="end" x="{left-10}" y="{y+5:.2f}">{elevation:.2f}</text>')
        svg.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{bottom}"/>')
        svg.append(f'<line class="axis" x1="{left}" y1="{bottom}" x2="{width-right}" y2="{bottom}"/>')
        route = " → ".join([path[0].flow_from] + [edge.flow_to for edge in path])
        svg.append(f'<text class="label" x="{left}" y="{top-12}">Route {panel+1}: {escape(route)}</text>')

        for channel, offset, profile in prepared:
            surface = [(sx(offset + x), sy(level)) for x, _, level in profile]
            bed = [(sx(offset + x), sy(level)) for x, level, _ in profile]
            crown = [(sx(offset + x), sy(level + channel.diameter)) for x, level, _ in profile]
            fill = surface + list(reversed(bed))
            svg.append(f'<polygon class="water" points="{_polyline(fill)}"/>')
            svg.append(f'<polyline class="crown" points="{_polyline(crown)}"/>')
            svg.append(f'<polyline class="bed" points="{_polyline(bed)}"/>')
            svg.append(f'<polyline class="surface" points="{_polyline(surface)}"/>')
            if show_energy_line:
                energy = [
                    (sx(offset + x), sy(head))
                    for x, head in channel.energy_profile()
                ]
                svg.append(f'<polyline class="energy" points="{_polyline(energy)}"/>')
            middle = len(surface) // 2
            ax, ay = surface[middle]
            svg.append(f'<polygon class="arrow" points="{ax-12:.2f},{ay-15:.2f} {ax+12:.2f},{ay-15:.2f} {ax:.2f},{ay-3:.2f}" transform="rotate(90 {ax:.2f} {ay-9:.2f})"/>')
            midpoint = sx(offset + channel.length / 2.0)
            flow = abs(channel.mass_flow) / RHO
            svg.append(f'<text class="small" text-anchor="middle" x="{midpoint:.2f}" y="{bottom+22}">{escape(channel.identifier)} · Q={flow:.4g} m³/s</text>')
            svg.append(f'<line class="grid" x1="{sx(offset):.2f}" y1="{top}" x2="{sx(offset):.2f}" y2="{bottom}"/>')
        if show_boundary_levels:
            upstream_level = prepared[0][0].energy_profile()[0][1]
            downstream_level = prepared[-1][2][-1][2]
            boundary_levels = (
                (sx(0.0), sy(upstream_level), "upstream rest level z_e+h_e+v_e²/(2g)", upstream_level, "start", 14),
                (sx(total_length), sy(downstream_level), "downstream rest level z_v+h_v", downstream_level, "end", -14),
            )
            for x, y, label, level, anchor, dx in boundary_levels:
                svg.extend([
                    f'<g class="boundary-level"><circle cx="{x:.2f}" cy="{y:.2f}" r="8"/>',
                    f'<line x1="{x-5:.2f}" y1="{y:.2f}" x2="{x+5:.2f}" y2="{y:.2f}"/>',
                    f'<line x1="{x:.2f}" y1="{y-5:.2f}" x2="{x:.2f}" y2="{y+5:.2f}"/></g>',
                    f'<text class="boundary-label" text-anchor="{anchor}" x="{x+dx:.2f}" y="{y-14:.2f}">{label}: {level:.3f} m</text>',
                ])
        svg.append(f'<line class="grid" x1="{sx(total_length):.2f}" y1="{top}" x2="{sx(total_length):.2f}" y2="{bottom}"/>')
        svg.append(f'<text class="small" text-anchor="middle" x="{(left+width-right)/2:.2f}" y="{bottom+47}">Distance along flow direction [m] · total {total_length:.3g} m</text>')

    profile_center = header + topology_height + panel_height * len(paths) / 2.0
    svg.append(f'<text class="small" transform="translate(25 {profile_center:.2f}) rotate(-90)" text-anchor="middle">Elevation [m]</text>')
    svg.append('</svg>')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(svg) + "\n", encoding="utf-8")


def _pdf_string(value: str) -> str:
    safe = value.encode("latin-1", errors="replace").decode("latin-1")
    return safe.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def _pdf_text(x: float, y: float, size: float, value: str, bold: bool = False) -> str:
    font = "F2" if bold else "F1"
    return f"BT /{font} {size:.2f} Tf {x:.2f} {y:.2f} Td ({_pdf_string(value)}) Tj ET"


def _pdf_path(points: Iterable[tuple[float, float]], close: bool = False) -> str:
    points = list(points)
    if not points:
        return ""
    commands = [f"{points[0][0]:.2f} {points[0][1]:.2f} m"]
    commands.extend(f"{x:.2f} {y:.2f} l" for x, y in points[1:])
    if close:
        commands.append("h")
    return " ".join(commands)


def _pdf_topology_page(
    network: Path,
    channels: list[Channel],
    document_title: str,
    total_pages: int,
    show_hydraulic_details: bool = False,
) -> str:
    positions = _scaled_topology_positions(
        _topology_coordinates(network, channels), 105.0, 150.0, 632.0, 285.0
    )
    commands = [
        "1 1 1 rg 0 0 842 595 re f",
        "0.14 0.19 0.23 rg",
        _pdf_text(72, 558, 17, document_title, True),
        _pdf_text(72, 530, 13, "Network topology and calculated flow directions", True),
        _pdf_text(72, 512, 8.5, "Only channel1 elements are shown. Arrows follow the solved mass-flow signs."),
        _pdf_text(682, 558, 8.5, f"Source: {network.name}"),
        "0.97 0.98 0.99 rg 72 105 698 375 re f",
        "0.80 0.84 0.87 RG 0.8 w 72 105 698 375 re S",
    ]
    for channel in channels:
        x0, y0 = positions[channel.flow_from]
        x1, y1 = positions[channel.flow_to]
        dx, dy = x1 - x0, y1 - y0
        distance = max(math.hypot(dx, dy), 1.0)
        ux, uy = dx / distance, dy / distance
        start = (x0 + 10.0 * ux, y0 + 10.0 * uy)
        end = (x1 - 12.0 * ux, y1 - 12.0 * uy)
        commands.extend([
            "0.03 0.50 0.74 RG 2.6 w",
            f"{start[0]:.2f} {start[1]:.2f} m {end[0]:.2f} {end[1]:.2f} l S",
        ])
        arrow_x = x0 + 0.68 * dx
        arrow_y = y0 + 0.68 * dy
        px, py = -uy, ux
        arrow = [
            (arrow_x + 9.0 * ux, arrow_y + 9.0 * uy),
            (arrow_x - 7.0 * ux + 5.0 * px, arrow_y - 7.0 * uy + 5.0 * py),
            (arrow_x - 7.0 * ux - 5.0 * px, arrow_y - 7.0 * uy - 5.0 * py),
        ]
        commands.extend(["0.03 0.37 0.57 rg", _pdf_path(arrow, close=True) + " f"])
        mx, my = (x0 + x1) / 2.0, (y0 + y1) / 2.0
        channel_label = channel.identifier
        if show_hydraulic_details:
            flow = abs(channel.mass_flow) / RHO
            channel_label += f", Q={flow:.4g} m3/s"
        label_width = max(43.0, 4.8 * len(channel_label))
        commands.extend([
            f"1 1 1 rg {mx-label_width/2:.2f} {my+5:.2f} {label_width:.2f} 15 re f",
            "0.14 0.19 0.23 rg",
            _pdf_text(mx - 2.4 * len(channel_label), my + 9, 8.5, channel_label, True),
        ])
    topology_boundary_data = (
        _topology_endpoint_data(channels) if show_hydraulic_details else {}
    )
    for node_id, (x, y) in positions.items():
        commands.extend([
            "0.96 0.65 0.26 rg",
            f"{x-7:.2f} {y-7:.2f} 14 14 re f",
            "0.35 0.26 0.12 RG 1.2 w",
            f"{x-7:.2f} {y-7:.2f} 14 14 re S",
            "0.14 0.19 0.23 rg",
            _pdf_text(x - 2.7 * len(node_id), y - 26, 9, node_id, True),
        ])
        for line_index, (symbol, invert, depth) in enumerate(
            topology_boundary_data.get(node_id, ())
        ):
            details = f"z_{symbol}={invert:.3f} m, h_{symbol}={depth:.3f} m"
            commands.append(
                _pdf_text(x - 2.05 * len(details), y - 41 - 13 * line_index, 7.5, details)
            )
    commands.extend([
        "0.03 0.50 0.74 RG 2.6 w 92 72 m 128 72 l S",
        "0.03 0.37 0.57 rg",
        _pdf_path([(128, 72), (118, 77), (118, 67)], close=True) + " f",
        "0.14 0.19 0.23 rg",
        _pdf_text(139, 68, 8.5, "calculated flow direction"),
        "0.96 0.65 0.26 rg 315 65 14 14 re f",
        "0.35 0.26 0.12 RG 1.2 w 315 65 14 14 re S",
        "0.14 0.19 0.23 rg",
        _pdf_text(339, 68, 8.5, "junction / boundary node"),
        _pdf_text(735, 12, 7.5, f"Page 1 of {total_pages}"),
    ])
    return "\n".join(command for command in commands if command)


def _write_pdf(output: Path, pages: list[str], title: str) -> None:
    """Write a compact vector PDF using only built-in PDF drawing operators."""
    page_count = len(pages)
    info_number = 5 + 2 * page_count
    objects: dict[int, bytes] = {
        1: b"<< /Type /Catalog /Pages 2 0 R >>",
        2: (
            f"<< /Type /Pages /Count {page_count} /Kids ["
            + " ".join(f"{5 + 2 * index} 0 R" for index in range(page_count))
            + "] >>"
        ).encode("ascii"),
        3: b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
        4: b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>",
    }
    for index, commands in enumerate(pages):
        page_number = 5 + 2 * index
        content_number = page_number + 1
        stream = commands.encode("latin-1", errors="replace")
        objects[page_number] = (
            f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 842 595] "
            f"/Resources << /Font << /F1 3 0 R /F2 4 0 R >> >> "
            f"/Contents {content_number} 0 R >>"
        ).encode("ascii")
        objects[content_number] = (
            f"<< /Length {len(stream)} >>\nstream\n".encode("ascii")
            + stream
            + b"\nendstream"
        )
    objects[info_number] = (
        f"<< /Title ({_pdf_string(title)}) /Creator (STACI channel profile renderer) >>"
    ).encode("latin-1", errors="replace")

    document = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offsets = [0] * (info_number + 1)
    for number in range(1, info_number + 1):
        offsets[number] = len(document)
        document.extend(f"{number} 0 obj\n".encode("ascii"))
        document.extend(objects[number])
        document.extend(b"\nendobj\n")
    xref = len(document)
    document.extend(f"xref\n0 {info_number + 1}\n".encode("ascii"))
    document.extend(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        document.extend(f"{offset:010d} 00000 n \n".encode("ascii"))
    document.extend(
        f"trailer\n<< /Size {info_number + 1} /Root 1 0 R /Info {info_number} 0 R >>\n"
        f"startxref\n{xref}\n%%EOF\n".encode("ascii")
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(document)


def render_channel_profiles_pdf(
    network: Path,
    output: Path,
    title: str | None = None,
    show_energy_line: bool = False,
    show_boundary_levels: bool = False,
    show_topology_hydraulics: bool = False,
) -> None:
    network = network.expanduser().resolve()
    channels = _read_channels(network)
    paths, path_profiles, z_min, z_max = _prepared_paths(channels, show_energy_line)
    document_title = title or f"STACI channel profiles - {network.name}"
    page_width, page_height = 842.0, 595.0
    left, right, bottom, top = 72.0, 28.0, 105.0, 470.0
    plot_width, plot_height = page_width - left - right, top - bottom
    total_pages = len(paths) + 1
    pages: list[str] = [
        _pdf_topology_page(
            network, channels, document_title, total_pages,
            show_topology_hydraulics,
        )
    ]

    for page_index, (path, prepared) in enumerate(zip(paths, path_profiles)):
        total_length = sum(channel.length for channel in path)
        x_scale = plot_width / max(total_length, 1.0)
        y_scale = plot_height / max(z_max - z_min, 1.0e-9)

        def sx(distance: float) -> float:
            return left + distance * x_scale

        def sy(elevation: float) -> float:
            return bottom + (elevation - z_min) * y_scale

        route = " -> ".join([path[0].flow_from] + [edge.flow_to for edge in path])
        commands = [
            "1 1 1 rg 0 0 842 595 re f",
            "0.14 0.19 0.23 rg",
            _pdf_text(left, 558, 17, document_title, True),
            _pdf_text(left, 535, 11, f"Flow route {page_index + 1}/{len(paths)}: {route}", True),
            _pdf_text(left, 518, 8.5, "Flow is left to right; all distances, elevations and discharges are in SI units."),
            _pdf_text(682, 558, 8.5, f"Source: {network.name}"),
        ]
        for tick in range(6):
            elevation = z_min + (z_max - z_min) * tick / 5.0
            y = sy(elevation)
            commands.extend([
                "0.86 0.89 0.91 RG 0.5 w",
                f"{left:.2f} {y:.2f} m {page_width-right:.2f} {y:.2f} l S",
                "0.14 0.19 0.23 rg",
                _pdf_text(left - 39, y - 3, 8, f"{elevation:.2f}"),
            ])
        commands.extend([
            "0.30 0.35 0.39 RG 0.8 w",
            f"{left:.2f} {bottom:.2f} m {left:.2f} {top:.2f} l S",
            f"{left:.2f} {bottom:.2f} m {page_width-right:.2f} {bottom:.2f} l S",
        ])
        for channel, offset, profile in prepared:
            surface = [(sx(offset + x), sy(level)) for x, _, level in profile]
            bed = [(sx(offset + x), sy(level)) for x, level, _ in profile]
            crown = [(sx(offset + x), sy(level + channel.diameter)) for x, level, _ in profile]
            fill = surface + list(reversed(bed))
            commands.extend([
                "0.51 0.79 0.93 rg",
                _pdf_path(fill, close=True) + " f",
                "0.55 0.59 0.62 RG 0.8 w [5 4] 0 d",
                _pdf_path(crown) + " S",
                "[] 0 d 0.45 0.31 0.17 RG 2.2 w",
                _pdf_path(bed) + " S",
                "0.03 0.50 0.74 RG 1.8 w",
                _pdf_path(surface) + " S",
            ])
            if show_energy_line:
                energy = [
                    (sx(offset + x), sy(head))
                    for x, head in channel.energy_profile()
                ]
                commands.extend([
                    "0.76 0.23 0.13 RG 1.6 w [7 4] 0 d",
                    _pdf_path(energy) + " S",
                    "[] 0 d",
                ])
            middle = len(surface) // 2
            ax, ay = surface[middle]
            arrow = [(ax - 7, ay + 9), (ax + 7, ay + 9), (ax + 7, ay + 14),
                     (ax + 15, ay + 6), (ax + 7, ay - 2), (ax + 7, ay + 3),
                     (ax - 7, ay + 3)]
            commands.extend(["0.03 0.37 0.57 rg", _pdf_path(arrow, close=True) + " f"])
            midpoint = sx(offset + channel.length / 2.0)
            flow = abs(channel.mass_flow) / RHO
            label = f"{channel.identifier}  Q={flow:.4g} m3/s"
            commands.extend([
                "0.14 0.19 0.23 rg",
                _pdf_text(midpoint - min(55, 3.2 * len(label)), 84, 7.5, label),
                "0.86 0.89 0.91 RG 0.5 w",
                f"{sx(offset):.2f} {bottom:.2f} m {sx(offset):.2f} {top:.2f} l S",
            ])
        if show_boundary_levels:
            upstream_level = prepared[0][0].energy_profile()[0][1]
            downstream_level = prepared[-1][2][-1][2]
            upstream_y = sy(upstream_level)
            downstream_y = sy(downstream_level)
            boundary_levels = (
                (sx(0.0), upstream_y, "upstream rest level z_e+h_e+v_e^2/(2g)", upstream_level, 9.0),
                (sx(total_length), downstream_y, "downstream rest level z_v+h_v", downstream_level, -156.0),
            )
            for x, y, label, level, label_dx in boundary_levels:
                commands.extend([
                    "0.84 0.15 0.12 RG 2.2 w",
                    f"{x-6:.2f} {y:.2f} m {x+6:.2f} {y:.2f} l S",
                    f"{x:.2f} {y-6:.2f} m {x:.2f} {y+6:.2f} l S",
                    "0.71 0.12 0.09 rg",
                    _pdf_text(x + label_dx, y + 10, 8, f"{label}: {level:.3f} m", True),
                ])
        commands.extend([
            "0.86 0.89 0.91 RG 0.5 w",
            f"{sx(total_length):.2f} {bottom:.2f} m {sx(total_length):.2f} {top:.2f} l S",
            "0.14 0.19 0.23 rg",
            _pdf_text(330, 54, 9, f"Distance along flow direction [m] - total {total_length:.3g} m"),
            "BT /F1 9 Tf 0 1 -1 0 20 255 Tm (Elevation [m]) Tj ET",
            "0.51 0.79 0.93 rg 595 30 22 10 re f",
            "0.14 0.19 0.23 rg",
            _pdf_text(623, 31, 8, "water"),
            *( [
                "0.76 0.23 0.13 RG 1.6 w [7 4] 0 d 430 35 m 465 35 l S [] 0 d",
                "0.14 0.19 0.23 rg",
                _pdf_text(472, 31, 8, "energy grade line H=z+y+v2/(2g)"),
            ] if show_energy_line else [] ),
            "0.45 0.31 0.17 RG 2.2 w 675 35 m 705 35 l S",
            "0.14 0.19 0.23 rg",
            _pdf_text(712, 31, 8, "channel invert"),
            _pdf_text(735, 12, 7.5, f"Page {page_index + 2} of {total_pages}"),
        ])
        pages.append("\n".join(command for command in commands if command))
    _write_pdf(output, pages, document_title)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", type=Path, required=True, help="Solved SPR file")
    parser.add_argument("--output", type=Path, help="SVG output path")
    parser.add_argument("--pdf-output", type=Path, help="PDF output path")
    parser.add_argument("--title", help="Optional plot title")
    arguments = parser.parse_args()
    output = arguments.output or arguments.network.with_name(arguments.network.stem + "-profile.svg")
    pdf_output = arguments.pdf_output or output.with_suffix(".pdf")
    try:
        render_channel_profiles(arguments.network, output, arguments.title)
        render_channel_profiles_pdf(arguments.network, pdf_output, arguments.title)
        print(f"Channel profile: {output.resolve()}")
        print(f"Channel profile PDF: {pdf_output.resolve()}")
        return 0
    except (PlotError, OSError, ValueError, ET.ParseError) as error:
        print(f"Channel profile generation failed: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
