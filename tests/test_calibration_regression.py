"""Multi-period objective invariance and strict legacy CSV measurement parsing."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET

p = argparse.ArgumentParser()
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--state-test', type=Path, required=True)
a = p.parse_args()
root = Path(__file__).resolve().parents[1]
source = root / 'tests/anytown_1med.spr'
tree = ET.parse(source)
pool = next(e for e in tree.findall('edges/edge') if e.find('edge_spec/pool') is not None)
pool_id = pool.findtext('id').strip()
level = float(pool.findtext('edge_spec/pool/water_level'))

with tempfile.TemporaryDirectory(prefix='staci calibration regression ') as tmp:
    work = Path(tmp)
    base = json.loads((root / 'examples/config/staci_calibrate_settings.json').read_text())
    base.update(dir_name=str(work) + '/', global_debug_level=0)
    settings = work / 'settings.json'
    targets = work / 'calibrate_targets.json'
    for offset in (0, 3):
        config = dict(base, Start_of_Periods=offset, Num_of_Periods=2)
        settings.write_text(json.dumps(config))
        for i in range(offset, offset + 2):
            shutil.copy2(source, work / f'calibrate_network_{i}.spr')
        targets.write_text(json.dumps({'measurements': [
            {'id': 'NODE64', 'type': 'node', 'values': [0] * (offset + 2)},
            {'id': pool_id, 'type': 'pool', 'values': [level] * (offset + 2)}]}))
        result = subprocess.run([str(a.state_test.resolve()), '--settings', str(settings)],
                                cwd=work, capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, result.stdout[-2000:] + result.stderr

    settings.write_text(json.dumps(dict(base, sollwert_dfile='targets.csv')))
    csv = work / 'targets.csv'
    for kind, identifier in [('node', 'NODE64'), ('pool', pool_id)]:
        for value in ('not-a-number', '', '12oops', 'nan', 'inf', '1e999'):
            csv.write_text(f'{identifier};{kind};{value};\n')
            log = work / 'diagnostics.jsonl'
            if log.exists():
                log.unlink()
            run = subprocess.run([str(a.binary.resolve()), '--settings', str(settings),
                                  '--diagnostics-file', str(log)], cwd=work,
                                 capture_output=True, text=True, timeout=30)
            assert run.returncode == 2, (kind, value, run.stderr)
            records = [json.loads(line) for line in log.read_text().splitlines()]
            messages = [r['message'] for r in records if r.get('code') == 'INPUT.CONFIG']
            assert any('targets.csv' in m and identifier in m and 'period 0' in m and
                       'finite number' in m for m in messages), records
            assert records[-1]['exit_code'] == 2
    csv.write_text('NODE64;node; 0.0e+0 ;\n')
    run = subprocess.run([str(a.binary.resolve()), '--settings', str(settings)],
                         cwd=work, capture_output=True, text=True, timeout=30)
    assert run.returncode == 0, run.stderr
print('PASS: multi-period calibration and numeric CSV validation')
