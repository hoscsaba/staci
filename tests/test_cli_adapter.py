"""Exercise the documented Python adapter against a real solver process."""
import argparse
import importlib.util
from pathlib import Path
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--binary', type=Path, required=True)
a = p.parse_args()
root = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('run_staci', root / 'examples/integration/run_staci.py')
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
network = root / 'tests/epanet_eps_smoke.inp'
with tempfile.TemporaryDirectory(prefix='staci adapter spaces ') as tmp:
 for eps in (False, True):
  result = adapter.run_network(a.binary, network, tmp, eps=eps)
  assert result['state']=='success' and result['exit_code']==0, result
  assert Path(result['result_file']).is_file()
  assert result['run_id']
 bad = Path(tmp) / 'bad.inp'
 bad.write_text('[JUNCTIONS]\nJ bad-number\n[END]\n')
 result = adapter.run_network(a.binary, bad, tmp)
 assert result['state']=='failed' and result['exit_code']!=0 and result['result_file'] is None
 result = adapter.run_network(a.binary, network, tmp, eps=True, timeout=1e-6)
 assert result['state']=='timeout' and result['exit_code'] is None and result['result_file'] is None
print('PASS: Python adapter steady/EPS success, bad input, Unicode-safe argument arrays and timeout.')
