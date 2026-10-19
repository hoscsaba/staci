#!/usr/bin/env python3
"""Fail hosted CI when an EPANET comparison is missing/skipped or a test fails."""
import argparse
from pathlib import Path
import xml.etree.ElementTree as ET

# The 18 formerly skipped comparisons plus two tests enabled by the shared library.
REQUIRED_REFERENCE_TESTS = {
    'epanet_reference_pipe', 'epanet_reference_tcv', 'epanet_reference_pattern',
    'epanet_reference_anytown', 'epanet_reference_anytown_chlorine',
    'epanet_node_metadata', 'epanet_pumps', 'epanet_eps_smoke', 'epanet_controls',
    'epanet_rules', 'epanet_volume_curve', 'epanet_reference_tcv_controls',
    'epanet_reference_tcv_rules', 'public_epanet_reference', 'epanet_physical_models',
    'epanet_reference_valves', 'epanet_pressure_elements', 'epanet_reference_units',
    'epanet_eps_chemical_reference', 'epanet_eps_check_valve_reference',
}


def check_report(path):
    cases = ET.parse(path).getroot().findall('.//testcase')
    missing = REQUIRED_REFERENCE_TESTS - {case.get('name') for case in cases}
    unsuccessful = [case.get('name', '<unnamed>') for case in cases
                    if any(case.find(tag) is not None for tag in ('skipped', 'failure', 'error'))
                    or case.get('status', 'run') not in ('run', 'passed')]
    if missing or unsuccessful:
        raise ValueError(f"Incomplete CI test coverage: missing reference tests={sorted(missing)}; "
                         f"skipped/failed tests={unsuccessful}")
    print(f"PASS: {len(cases)} tests, no skips/failures; all {len(REQUIRED_REFERENCE_TESTS)} "
          "required EPANET reference comparisons executed")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    args = parser.parse_args()
    check_report(args.report)
