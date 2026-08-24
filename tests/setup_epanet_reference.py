#!/usr/bin/env python3
"""Download and build the pinned official OWA EPANET 2.2 reference solver."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request


EPANET_COMMIT = "4d8d82ddc260fce216af9321fc3d9a4646ac6827"
EPANET_ARCHIVE = (
    "https://github.com/OpenWaterAnalytics/EPANET/archive/"
    f"{EPANET_COMMIT}.tar.gz"
)
SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_DIR = SCRIPT_DIR.parent


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--build-root", type=Path,
        default=REPOSITORY_DIR / "build" / "epanet-reference",
        help="Source and build directory (default: build/epanet-reference).",
    )
    parser.add_argument(
        "--jobs", type=int, default=max(1, os.cpu_count() or 1),
        help="Maximum parallel build jobs.",
    )
    return parser.parse_args()


def run(command: list[str]) -> None:
    print("+", subprocess.list2cmdline(command), flush=True)
    subprocess.run(command, check=True)


def executable_candidates(build_dir: Path) -> list[Path]:
    names = ("runepanet.exe", "runepanet") if os.name == "nt" else ("runepanet", "runepanet.exe")
    return [
        build_dir / "bin" / name
        for name in names
    ] + [
        build_dir / "bin" / configuration / name
        for configuration in ("Release", "RelWithDebInfo", "Debug")
        for name in names
    ]


def safely_extract(archive: Path, destination: Path) -> Path:
    with tarfile.open(archive, "r:gz") as bundle:
        members = bundle.getmembers()
        roots = {Path(member.name).parts[0] for member in members if Path(member.name).parts}
        if len(roots) != 1:
            raise RuntimeError("Unexpected EPANET source archive layout")
        destination_resolved = destination.resolve()
        for member in members:
            target = (destination / member.name).resolve()
            if target != destination_resolved and destination_resolved not in target.parents:
                raise RuntimeError(f"Unsafe archive member: {member.name}")
            if member.issym() or member.islnk():
                raise RuntimeError(f"Archive links are not accepted: {member.name}")
        if sys.version_info >= (3, 12):
            bundle.extractall(destination, members=members, filter="data")
        else:
            bundle.extractall(destination, members=members)
    return destination / roots.pop()


def main() -> int:
    args = parse_arguments()
    if args.jobs <= 0:
        raise SystemExit("--jobs must be positive")
    if shutil.which("cmake") is None:
        raise SystemExit("CMake was not found on PATH")

    root = args.build_root.expanduser().resolve()
    source_dir = root / "src"
    build_dir = root / "build"
    marker = source_dir / ".staci-epanet-commit"
    if not marker.is_file() or marker.read_text(encoding="ascii").strip() != EPANET_COMMIT:
        if source_dir.exists():
            shutil.rmtree(source_dir)
        root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="staci-epanet-download-") as temporary:
            temporary_dir = Path(temporary)
            archive = temporary_dir / "epanet.tar.gz"
            print(f"Downloading official OWA EPANET commit {EPANET_COMMIT} ...", flush=True)
            urllib.request.urlretrieve(EPANET_ARCHIVE, archive)
            extracted = safely_extract(archive, temporary_dir / "source")
            shutil.move(str(extracted), str(source_dir))
        marker.write_text(EPANET_COMMIT + "\n", encoding="ascii")

    run([
        "cmake", "-S", str(source_dir), "-B", str(build_dir),
        "-DCMAKE_BUILD_TYPE=Release", "-DBUILD_TESTS=OFF", "-DBUILD_PY_LIB=OFF",
        "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
    ])
    run(["cmake", "--build", str(build_dir), "--config", "Release", "--parallel", str(args.jobs)])

    executable = next((path for path in executable_candidates(build_dir) if path.is_file()), None)
    if executable is None:
        raise SystemExit(f"Build succeeded but runepanet was not found below {build_dir}")
    print(f"EPANET_EXECUTABLE={executable}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
