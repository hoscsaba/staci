#!/usr/bin/env python3
"""Render flow-oriented longitudinal profiles from a solved STACI SPR file.

Only the Python standard library is used.  The output is SVG so the elevation
and water-level curves remain sharp when the engineering drawing is enlarged.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from html import escape
from pathlib import Path
from typing import Iterable
import xml.etree.ElementTree as ET


RHO = 1000.0


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


def _polyline(points: Iterable[tuple[float, float]]) -> str:
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in points)


def _prepared_paths(
    channels: list[Channel],
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
            distance += channel.length
        path_profiles.append(prepared)
    z_min, z_max = min(elevations), max(elevations)
    padding = max(0.15, 0.08 * max(z_max - z_min, 1.0))
    return paths, path_profiles, z_min - padding, z_max + padding


def render_channel_profiles(network: Path, output: Path, title: str | None = None) -> None:
    network = network.expanduser().resolve()
    channels = _read_channels(network)
    paths, path_profiles, z_min, z_max = _prepared_paths(channels)
    width, panel_height = 1500, 310
    header, footer, left, right = 95, 55, 105, 35
    height = header + panel_height * len(paths) + footer
    plot_width = width - left - right
    plot_height = panel_height - 82
    document_title = title or f"STACI channel profiles — {network.name}"
    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#24303a}.axis{stroke:#4d5963;stroke-width:1}.grid{stroke:#dce3e8;stroke-width:1}.bed{stroke:#72502b;stroke-width:4;fill:none}.crown{stroke:#8d969d;stroke-width:1.5;stroke-dasharray:7 5;fill:none}.water{fill:#82c9ee;fill-opacity:.62;stroke:none}.surface{stroke:#087fbd;stroke-width:3;fill:none}.arrow{fill:#075f91}.small{font-size:13px}.label{font-size:15px;font-weight:bold}.title{font-size:24px;font-weight:bold}</style>',
        f'<text class="title" x="{left}" y="38">{escape(document_title)}</text>',
        f'<text class="small" x="{left}" y="64">Flow is left to right. Elevations and water levels are read from the solved SPR file in SI units.</text>',
        '<rect x="1080" y="27" width="34" height="14" fill="#82c9ee" fill-opacity=".62"/><text class="small" x="1122" y="39">water</text>',
        '<line x1="1210" y1="35" x2="1250" y2="35" class="bed"/><text class="small" x="1260" y="39">channel invert</text>',
    ]

    for panel, (path, prepared) in enumerate(zip(paths, path_profiles)):
        top = header + panel * panel_height
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
            middle = len(surface) // 2
            ax, ay = surface[middle]
            svg.append(f'<polygon class="arrow" points="{ax-12:.2f},{ay-15:.2f} {ax+12:.2f},{ay-15:.2f} {ax:.2f},{ay-3:.2f}" transform="rotate(90 {ax:.2f} {ay-9:.2f})"/>')
            midpoint = sx(offset + channel.length / 2.0)
            flow = abs(channel.mass_flow) / RHO
            svg.append(f'<text class="small" text-anchor="middle" x="{midpoint:.2f}" y="{bottom+22}">{escape(channel.identifier)} · Q={flow:.4g} m³/s</text>')
            svg.append(f'<line class="grid" x1="{sx(offset):.2f}" y1="{top}" x2="{sx(offset):.2f}" y2="{bottom}"/>')
        svg.append(f'<line class="grid" x1="{sx(total_length):.2f}" y1="{top}" x2="{sx(total_length):.2f}" y2="{bottom}"/>')
        svg.append(f'<text class="small" text-anchor="middle" x="{(left+width-right)/2:.2f}" y="{bottom+47}">Distance along flow direction [m] · total {total_length:.3g} m</text>')

    svg.append(f'<text class="small" transform="translate(25 {height/2:.2f}) rotate(-90)" text-anchor="middle">Elevation [m]</text>')
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


def render_channel_profiles_pdf(network: Path, output: Path, title: str | None = None) -> None:
    network = network.expanduser().resolve()
    channels = _read_channels(network)
    paths, path_profiles, z_min, z_max = _prepared_paths(channels)
    document_title = title or f"STACI channel profiles - {network.name}"
    page_width, page_height = 842.0, 595.0
    left, right, bottom, top = 72.0, 28.0, 105.0, 470.0
    plot_width, plot_height = page_width - left - right, top - bottom
    pages: list[str] = []

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
        commands.extend([
            f"{sx(total_length):.2f} {bottom:.2f} m {sx(total_length):.2f} {top:.2f} l S",
            "0.14 0.19 0.23 rg",
            _pdf_text(330, 54, 9, f"Distance along flow direction [m] - total {total_length:.3g} m"),
            "BT /F1 9 Tf 0 1 -1 0 20 255 Tm (Elevation [m]) Tj ET",
            "0.51 0.79 0.93 rg 595 30 22 10 re f",
            "0.14 0.19 0.23 rg",
            _pdf_text(623, 31, 8, "water"),
            "0.45 0.31 0.17 RG 2.2 w 675 35 m 705 35 l S",
            "0.14 0.19 0.23 rg",
            _pdf_text(712, 31, 8, "channel invert"),
            _pdf_text(735, 12, 7.5, f"Page {page_index + 1} of {len(paths)}"),
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
