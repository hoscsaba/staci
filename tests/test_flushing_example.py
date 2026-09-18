"""Run the documented example in isolation and check its reference and hydraulics."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

exe = str(Path(sys.argv[1]).resolve())
example = Path(sys.argv[2]).resolve()
original = {name: (example/name).read_bytes() for name in ('network.inp','flushing_config.json')}
with tempfile.TemporaryDirectory(prefix='staci worked example ') as temp:
    root = Path(temp)
    for name, data in original.items():
        (root/name).write_bytes(data)
    result = subprocess.run([exe, '--inp', str(root/'network.inp'),
                             '--config', str(root/'flushing_config.json')],
                            cwd=root, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, (result.stdout, result.stderr)
    output = root/'results'
    def rows(path):
        with path.open() as f:
            return list(csv.DictReader(f))
    actual = rows(output/'flushing_plan.csv')
    expected = rows(example/'expected_plan.csv')
    assert len(actual) == len(expected) == 3
    for row, reference in zip(actual,expected):
        assert {key:row[key] for key in reference} == reference, (row,reference)
    assert '1. C: additional' in (output/'flushing_plan.txt').read_text()
    # Solve h_supply = h_MAIN(Q) + h_branch(Q) + K*(Q/A)^2/(2g)
    # independently by bisection, using STACI's documented H-W resistance.
    def resistance(length, diameter):
        return length/120**1.85/diameter**4.87*7.88/.85**1.85
    for row in rows(output/'scenarios.csv'):
        node = row['node_id']
        pipe_resistance = resistance(100,.15)
        if node != 'A':
            pipe_resistance += resistance(100 if node=='B' else 150,.1)
        lo, hi = 0., .1
        for _ in range(100):
            q = (lo+hi)/2
            head = pipe_resistance*q**1.85 + 2*(q/.002)**2/(2*9.81)
            if head > 40: hi=q
            else: lo=q
        assert row['status']=='ok'
        assert math.isclose(float(row['flow_m3s']), (lo+hi)/2, rel_tol=1e-6)
    detail = rows(output/'scenario_pipes.csv')
    reverse = next(r for r in detail if r['node_id']=='C' and r['pipe_id']=='BRANCH_C')
    assert float(reverse['velocity_mps']) < -.5 and reverse['above_threshold']=='1'
    exported = (output/'networks/network_hydrant_C.inp').read_text()
    assert 'LINK MAIN FLUSHED' in exported and 'LINK BRANCH_C FLUSHED' in exported
    assert 'LINK BRANCH_B BELOW_THRESHOLD' in exported
    assert '[EMITTERS]\nC ' in exported and '[COORDINATES]' in exported
    assert json.loads((output/'run.json').read_text())['status']=='complete'
    for name,data in original.items():
        assert (example/name).read_bytes()==data and (root/name).read_bytes()==data
    assert not (root/'network.inp.ros').exists()
print('Worked flushing example: reference plan, independent hydraulics and exports passed')
