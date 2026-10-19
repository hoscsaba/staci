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

    # INP periods: real EPANET tank IDs, different demands, selected pipe not
    # first in the importer order, and nonzero measurement-window offsets.
    inp_template = """[JUNCTIONS]
J1 0 {demand}
[RESERVOIRS]
R1 40
[TANKS]
T1 0 30 0 100 20 0
[PIPES]
Small R1 J1 100 100 120 0 Open
Large J1 T1 100 200 120 0 Open
[OPTIONS]
UNITS LPS
HEADLOSS H-W
[END]
"""
    for offset in (0, 3):
        config = dict(base, network_extension='.inp', Start_of_Periods=offset,
                      Num_of_Periods=2, dt=.01)
        settings.write_text(json.dumps(config))
        for i in range(offset, offset+2):
            (work/f'calibrate_network_{i}.inp').write_text(inp_template.format(demand=i+1))
        targets.write_text(json.dumps({'measurements': [
            {'id': 'J1', 'type': 'node', 'values': [30] * (offset+2)},
            {'id': 'T1', 'type': 'pool', 'values': [30] * (offset+2)}]}))
        run = subprocess.run([str(a.state_test.resolve()), '--settings', str(settings), '--seed', '123'],
                             cwd=work, capture_output=True, text=True, timeout=60)
        assert run.returncode == 0, run.stdout[-2000:]+run.stderr
        best = (work/'best.dat').read_text().splitlines()
        assert best[0].strip() == 'Large ;', best
        assert abs(float(best[2].split(';')[0]) - float(best[1].split(';')[0])/.2) < 1e-4

    # Direct filenames work with JSON and XML settings; no prefix required.
    direct = dict(base, network_file='calibrate_network_0.inp', dt=.01)
    del direct['fname_prefix']
    targets.write_text(json.dumps({'measurements': [
        {'id': 'J1', 'type': 'node', 'values': [30]},
        {'id': 'T1', 'type': 'pool', 'values': [30]}]}))
    for suffix in ('json', 'xml'):
        path = work / ('direct.'+suffix)
        if suffix == 'json':
            path.write_text(json.dumps(direct))
        else:
            xml = ET.Element('settings')
            for key,value in direct.items():
                ET.SubElement(xml,key).text = ('yes' if value else 'no') if isinstance(value,bool) else str(value)
            ET.ElementTree(xml).write(path, encoding='unicode')
        run = subprocess.run([str(a.binary.resolve()), '--settings', str(path)],
                             cwd=work, capture_output=True, text=True, timeout=30)
        assert run.returncode == 0, run.stdout[-2000:]+run.stderr
        assert 'CALIBRATION.UNINFORMATIVE_PRESSURE' in run.stderr
    # Explicit lists and all: no count or automatic-selection fields required.
    for selection, ids, expected_ids in [('list', ['Large'], ['Large']),
                                         ('list', ['Large','Small'], ['Small','Large']),
                                         ('all', None, ['Small','Large'])]:
        selected = dict(direct, type_of_pipe_selection=selection)
        selected.pop('num_of_active_pipes')
        if ids is not None: selected['pipe_ids'] = ids
        for suffix in ('json', 'xml'):
            path = work / ('selection.'+suffix)
            if suffix == 'json': path.write_text(json.dumps(selected))
            else:
                xml = ET.Element('settings')
                for key,value in selected.items():
                    ET.SubElement(xml,key).text = (json.dumps(value) if isinstance(value,list) else
                        ('yes' if value else 'no') if isinstance(value,bool) else str(value))
                ET.ElementTree(xml).write(path,encoding='unicode')
            run = subprocess.run([str(a.binary.resolve()), '--settings', str(path)],
                                 cwd=work, capture_output=True, text=True, timeout=30)
            assert run.returncode == 0, run.stderr
            best = (work/'best.dat').read_text().splitlines()
            assert [x.strip() for x in best[0].split(';') if x.strip()] == expected_ids
            assert len([x for x in best[1].split(';') if x.strip()]) == len(expected_ids)
    for overrides, message in [
        ({'pipe_ids': []}, 'nonempty array'),
        ({'pipe_ids': ['Large','Large']}, 'duplicate pipe ID'),
        ({'pipe_ids': ['missing']}, 'not a calibratable pipe ID'),
        ({'pipe_ids': ['J1']}, 'not a calibratable pipe ID'),
        ({'pipe_ids': ['EPANET_TANK_T1']}, 'not a calibratable pipe ID'),
        ({'pipe_ids': [42]}, 'nonempty string IDs'),
        ({'pipe_ids': ['']}, 'nonempty string IDs'),
        ({'type_of_pipe_selection': 'all', 'pipe_ids': ['Large']}, 'requires type_of_pipe_selection=list')]:
        settings.write_text(json.dumps(dict(direct, **{'type_of_pipe_selection': 'list', **overrides})))
        run = subprocess.run([str(a.binary.resolve()), '--settings', str(settings)],
                             cwd=work, capture_output=True, text=True, timeout=30)
        assert run.returncode == 2 and message in run.stderr, run.stderr

    for overrides, message in [
        ({'network_file': 'missing.inp'}, 'does not exist'),
        ({'network_file': 'network.csv'}, '.spr or .inp'),
        ({'Num_of_Periods': 2}, 'requires Num_of_Periods=1'),
        ({'network_extension': '.inp'}, 'cannot be combined'),
        ({'dt': 0}, 'dt must be positive')]:
        settings.write_text(json.dumps(dict(direct, **overrides)))
        run = subprocess.run([str(a.binary.resolve()), '--settings', str(settings)],
                             cwd=work, capture_output=True, text=True, timeout=30)
        assert run.returncode == 2 and message in run.stderr, run.stderr
    # Absolute bounds override (rather than intersect) the old relative bounds.
    for bounds, lower, upper in [({}, .1, .4),
                                 ({'min_diameter_m': .45, 'max_diameter_m': .5}, .45, .5),
                                 ({'min_diameter_m': .25}, .25, .4),
                                 ({'max_diameter_m': .15}, .1, .15)]:
        for suffix in ('json', 'xml'):
            data = dict(direct, type_of_pipe_selection='list', pipe_ids=['Large'], **bounds)
            path = work/('bounds.'+suffix)
            if suffix == 'json': path.write_text(json.dumps(data))
            else:
                xml = ET.Element('settings')
                for key,value in data.items():
                    ET.SubElement(xml,key).text = (json.dumps(value) if isinstance(value,list) else
                        ('yes' if value else 'no') if isinstance(value,bool) else str(value))
                ET.ElementTree(xml).write(path,encoding='unicode')
            run = subprocess.run([str(a.binary.resolve()), '--settings', str(path), '--seed', '123'],
                                 cwd=work, capture_output=True, text=True, timeout=30)
            assert run.returncode == 0, run.stderr
            diameter = float((work/'best.dat').read_text().splitlines()[1].split(';')[0])
            assert lower <= diameter <= upper, (bounds,diameter)
            history = (work/base['best_logfilename']).read_text().splitlines()
            assert history and all(lower <= float(row.split(';')[2]) <= upper for row in history)
    for bounds, message in [({'min_diameter_m': 0}, 'positive finite'),
                             ({'max_diameter_m': -1}, 'positive finite'),
                             ({'min_diameter_m': 'nan'}, 'finite number'),
                             ({'min_diameter_m': .3, 'max_diameter_m': .3}, 'strictly less'),
                             ({'min_diameter_m': .4, 'max_diameter_m': .2}, 'strictly less'),
                             ({'min_diameter_m': .5}, "bounds for pipe 'Large'"),
                             ({'max_diameter_m': .05}, "bounds for pipe 'Large'")]:
        data = dict(direct, type_of_pipe_selection='list', pipe_ids=['Large'], **bounds)
        settings.write_text(json.dumps(data))
        run = subprocess.run([str(a.binary.resolve()), '--settings', str(settings)],
                             cwd=work, capture_output=True, text=True, timeout=30)
        assert run.returncode == 2 and message in run.stderr, run.stderr

    # Prevent silently applying measurements/tank history to different elements.
    settings.write_text(json.dumps(dict(base, network_extension='.inp', Num_of_Periods=2)))
    (work/'calibrate_network_1.inp').write_text(inp_template.format(demand=2).replace('Large ', 'Changed '))
    run = subprocess.run([str(a.binary.resolve()), '--settings', str(settings)],
                         cwd=work, capture_output=True, text=True, timeout=30)
    assert run.returncode == 2 and 'identical node/link IDs' in run.stderr, run.stderr

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
