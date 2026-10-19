"""Run and regenerate the documented Anytown calibration example in isolation."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--staci', type=Path, required=True)
p.add_argument('--example', type=Path, required=True)
a = p.parse_args()
example = a.example.resolve()
original = {f.name: f.read_bytes() for f in example.iterdir() if f.is_file()}
with tempfile.TemporaryDirectory(prefix='staci anytown calibration ') as tmp:
    job = Path(tmp)/'job'
    shutil.copytree(example, job)
    provenance = json.loads((job/'provenance.json').read_text())
    assert hashlib.sha256((job/'anytown.inp').read_bytes()).hexdigest() == provenance['source_sha256']
    # Reproducibility of all snapshot inputs and synthetic measurements.
    regenerated = Path(tmp)/'regenerated'
    import sys
    run = subprocess.run([sys.executable, str(job/'prepare.py'), '--staci', str(a.staci.resolve()),
                          '--output', str(regenerated)], cwd=job, capture_output=True, text=True, timeout=60)
    assert run.returncode == 0, run.stdout+run.stderr
    for name in ('network_0.inp', 'network_1.inp', 'network_2.inp'):
        generated_bytes = (regenerated/name).read_bytes()
        assert b'\r' not in generated_bytes, f'{name}: generated INP must use LF line endings'
        assert (job/name).read_bytes() == generated_bytes, name
    for name in ('settings.json', 'provenance.json'):
        assert json.loads((job/name).read_text()) == json.loads((regenerated/name).read_text()), name
    observed = json.loads((job/'measurements.json').read_text())['measurements']
    fresh = json.loads((regenerated/'measurements.json').read_text())['measurements']
    assert len(observed) == len(fresh) == 8
    for reference, actual in zip(observed, fresh):
        assert reference['id'] == actual['id'] and reference['type'] == actual['type']
        assert len(actual['values']) == 3
        assert all(math.isclose(x,y,abs_tol=1e-5) for x,y in zip(reference['values'],actual['values']))
    run = subprocess.run([str(a.binary.resolve()), '--settings', 'settings.json', '--seed', '12345'],
                         cwd=job, capture_output=True, text=True, timeout=120)
    assert run.returncode == 0, run.stdout[-2000:]+run.stderr
    best = [[s.strip() for s in line.split(';') if s.strip()] for line in (job/'best.dat').read_text().splitlines()]
    assert best[0] == ['1','2']
    for i, identifier in enumerate(best[0]):
        truth = provenance['true_diameters_m'][identifier]
        fitted = float(best[1][i])
        assert math.isclose(fitted, truth, rel_tol=.01), (identifier, fitted, truth)
        assert math.isclose(float(best[2][i]), fitted/(truth*provenance['input_diameter_factors'][identifier]), rel_tol=1e-4)
    assert len(best[4:]) == 3
    for row in best[4:]:
        assert len(row) == 16
        assert all(abs(float(row[i])-float(row[i+1])) < .1 for i in range(0,len(row),2))
    assert (job/'sensitivity_Pressures.dat').is_file()
    for name,data in original.items():
        assert (example/name).read_bytes() == data and (job/name).read_bytes() == data
print('PASS: Anytown regeneration, two selected diameters, three periods and source preservation')
