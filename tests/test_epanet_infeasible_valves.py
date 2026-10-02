"""Physical valve/pump conflicts must not pass by relaxing numerical tolerances."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--binary', type=Path, required=True)
a = p.parse_args()
root = Path(__file__).parent / 'public_networks'
cases = [
 ('wntr/wntr/tests/networks_for_testing/fcv_open_no_downstream_sources.inp', 'HYDRAULICS.VALVE_CONSTRAINT', 'VALVE'),
 ('wntr/wntr/tests/networks_for_testing/prv_closed_no_upstream_sources.inp', 'HYDRAULICS.VALVE_CONSTRAINT', 'VALVE'),
 ('owa/epanet-tests/valves/2fcvs.inp', 'HYDRAULICS.VALVE_CONSTRAINT', '44'),
 ('wntr/wntr/tests/networks_for_testing/io.inp', 'EPANET.POWER_PUMP_DEAD_END', 'pump2'),
]
with tempfile.TemporaryDirectory(prefix='staci physical constraints ') as tmp:
 for path, code, element in cases:
  for tolerance in ('0.0001', '1e-12'):
   work = Path(tmp) / (Path(path).stem + tolerance)
   work.mkdir()
   net = work / Path(path).name
   shutil.copy2(root / path, net)
   digest = hashlib.sha256(net.read_bytes()).hexdigest()
   run = subprocess.run([str(a.binary.resolve()), '-s', str(net),
       '--head-tolerance-m', tolerance, '--max-iterations', '1000'],
       cwd=work, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, timeout=30)
   records = [json.loads(line) for line in (work / 'staci-diagnostics.jsonl').read_text().splitlines()]
   assert run.returncode == 1, (path, tolerance, run.stderr)
   assert any(r.get('code') == code and element in r.get('message', '') for r in records), (path, records)
   assert records[-1]['event'] == 'run_end' and records[-1]['exit_code'] == 1
   output = Path(str(net) + '.hydraulics.json')
   assert not output.exists() or json.loads(output.read_text())['converged'] is False
   assert hashlib.sha256(net.read_bytes()).hexdigest() == digest
print('PASS: four physical conflicts rejected with element diagnostics at default and strict tolerances.')
