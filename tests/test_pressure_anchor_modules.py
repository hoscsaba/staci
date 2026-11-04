#!/usr/bin/env python3
"""Exercise an EPANET-compatible passive pressure anchor in every STACI app.

The anchor is a zero-demand terminal junction, joined by a fixed connector.
Its target head is an external design/calibration constraint, not a new
fixed-head source.  This test catches accidental exchange of water or age.
"""

import argparse
import csv
import json
from pathlib import Path
import subprocess
import tempfile


def run(binary, directory, *arguments):
    result = subprocess.run([str(binary.resolve()), *map(str, arguments)],
                            cwd=directory, capture_output=True, text=True,
                            timeout=60)
    assert result.returncode == 0, (binary.name, arguments,
                                    result.stdout[-1500:], result.stderr[-1500:])
    return result


def read_csv(path, key):
    with path.open(newline="") as stream:
        return {row[key]: row for row in csv.DictReader(stream)}


def main():
    parser = argparse.ArgumentParser()
    for app in ("staci", "split", "calibrate"):
        parser.add_argument("--" + app, type=Path, required=True)
    args = parser.parse_args()
    network = """[JUNCTIONS]
J1 0 0
D1 0 1
ANCHOR 0 0
[RESERVOIRS]
R1 50
[PIPES]
P0 R1 J1 100 100 120 0 Open
P1 J1 D1 100 100 120 0 Open
P_ANCHOR J1 ANCHOR 100 100 120 0 Open
[QUALITY]
ANCHOR 0.25
[SOURCES]
R1 CONCEN 1
[REACTIONS]
GLOBAL BULK -1
[OPTIONS]
UNITS LPS
HEADLOSS H-W
QUALITY CHEMICAL Chlorine mg/L
[TIMES]
DURATION 24:00
QUALITY TIMESTEP 0:05
[END]
"""
    with tempfile.TemporaryDirectory(prefix="staci pressure anchor modules ") as temp:
        work = Path(temp)
        inp = work / "anchor.inp"
        inp.write_text(network)
        base = work / "base.inp"
        base.write_text(network.replace("ANCHOR 0 0\n", "")
                        .replace("P_ANCHOR J1 ANCHOR 100 100 120 0 Open\n", "")
                        .replace("[QUALITY]\nANCHOR 0.25\n", ""))
        for name, path in (("anchor", inp), ("base", base)):
            run(args.staci, work, "-s", path)
            run(args.staci, work, "-q", path, "-o", work / name,
                "--quality-mode", "both")
        anchored = json.loads(Path(str(inp) + ".hydraulics.json").read_text())
        original = json.loads(Path(str(base) + ".hydraulics.json").read_text())
        for node in ("J1", "D1"):
            assert abs(anchored["nodes"][node]["head_m"] -
                       original["nodes"][node]["head_m"]) < 1e-7
        assert abs(anchored["links"]["P_ANCHOR"]["flow_m3s"]) < 1e-9
        anchor_nodes = read_csv(work / "anchor-steady-nodes.csv", "node_id")
        base_nodes = read_csv(work / "base-steady-nodes.csv", "node_id")
        for node in ("J1", "D1"):
            for field in ("water_age_s", "chemical_concentration_kgm3"):
                assert abs(float(anchor_nodes[node][field]) -
                           float(base_nodes[node][field])) < 1e-8
        assert float(anchor_nodes["ANCHOR"]["water_age_s"]) == 0
        assert abs(float(anchor_nodes["ANCHOR"]["chemical_concentration_kgm3"])
                   - 0.00025) < 1e-12

        # The legacy residence-time command also loads and solves this INP.
        run(args.staci, work, "-t", inp)

        split = {"global_debug_level": 0, "n_comm": 2,
                 "weight_type": "topology", "weight_type_mod": "diameter",
                 "fname": str(inp), "logfilename": "split.log",
                 "obj_type": "modularity", "popsize": 4, "ngen": 1,
                 "pmut": 0.25, "pcross": 0.8}
        split_file = work / "split.json"
        split_file.write_text(json.dumps(split))
        run(args.split, work, "--settings", split_file, "--seed", 12345)
        assert (work / "membership.txt").exists()

        calibrate = {"global_debug_level": 0, "Staci_debug_level": 0,
                     "network_file": str(inp), "dir_name": str(work) + "/",
                     "logfilename": "calibrate.log",
                     "best_logfilename": "calibrate-best.log",
                     "sollwert_dfile": "targets.json", "Start_of_Periods": 0,
                     "Num_of_Periods": 1, "Spoil_Active_Pipes": False,
                     "dt": 1.0, "weight_p_err": 1.0,
                     "type_of_pipe_selection": "list", "pipe_ids": ["P1"],
                     "popsize": 4, "ngen": 1, "pmut": 0.2, "pcross": 0.8}
        (work / "calibrate.json").write_text(json.dumps(calibrate))
        target_head = anchored["nodes"]["ANCHOR"]["pressure_m"]
        (work / "targets.json").write_text(json.dumps({"measurements": [
            {"id": "ANCHOR", "type": "node", "values": [target_head]}]}))
        run(args.calibrate, work, "--settings", work / "calibrate.json",
            "--seed", 12345)
        assert (work / "calibrate-best.log").exists()
    print("PASS zero-flow pressure anchor in hydraulics, residence time, "
          "steady quality, split and calibration")


if __name__ == "__main__":
    main()
