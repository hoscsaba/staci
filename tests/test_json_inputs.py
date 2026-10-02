"""Check legacy/JSON equivalence and diagnosed malformed auxiliary inputs."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET

p = argparse.ArgumentParser()
for name in ('staci', 'split', 'calibrate'):
    p.add_argument('--' + name, type=Path, required=True)
a = p.parse_args()
root = Path(__file__).resolve().parents[1]


def run(binary, work, args, expected=0):
    log = work / 'diagnostics.jsonl'
    if log.exists():
        log.unlink()
    result = subprocess.run([str(binary.resolve()), '--diagnostics-file', str(log), *args],
                            cwd=work, capture_output=True, text=True, timeout=120)
    assert result.returncode == expected, (result.returncode, result.stdout[-3000:], result.stderr)
    records = [json.loads(line) for line in log.read_text().splitlines()]
    assert records[-1]['event'] == 'run_end' and records[-1]['exit_code'] == expected
    if expected == 2:
        assert any(r.get('code') == 'INPUT.CONFIG' for r in records), records
    return result


with tempfile.TemporaryDirectory(prefix='staci JSON inputs ') as temporary:
    base = Path(temporary)
    for name in ('split', 'calibrate'):
        outputs = []
        for form in ('xml', 'json'):
            work = base / (name + form)
            work.mkdir()
            template = ET.parse(root / 'tests' / ('staci_' + name + '_settings.xml.in'))
            settings = template.getroot()
            if name == 'split':
                shutil.copy2(root / 'tests/LOV-LOVOTV-2-input_mod.spr', work / 'network.spr')
                settings.find('fname').text = 'network.spr'
            else:
                shutil.copy2(root / 'tests/anytown_1med.spr', work / 'calibrate_network_0.spr')
                settings.find('dir_name').text = './'
                if form == 'xml':
                    shutil.copy2(root / 'tests/calibrate_targets.csv', work / 'calibrate_targets.csv')
                else:
                    shutil.copy2(root / 'examples/config/calibrate_targets.json', work / 'calibrate_targets.json')
                    settings.find('sollwert_dfile').text = 'calibrate_targets.json'
            path = work / ('custom settings.' + form.upper())
            if form == 'xml':
                template.write(path)
            else:
                data = {}
                for child in settings:
                    try:
                        data[child.tag] = json.loads(child.text)
                    except ValueError:
                        data[child.tag] = child.text
                if name == 'calibrate':
                    data['Spoil_Active_Pipes'] = False
                path.write_text(json.dumps(data))
            binary = getattr(a, name)
            run(binary, work, ['--settings', str(path), '--seed', '12345'])
            output = work / ('membership.txt' if name == 'split' else 'calibrate-best.log')
            outputs.append(output.read_text())
            if form == 'json':
                # JSON-only default discovery, without an explicit CLI path.
                default = work / ('staci_' + name + '_settings.json')
                shutil.copy2(path, default)
                run(binary, work, ['--seed', '12345'])
                assert output.read_text() == outputs[-1]
                bad = work / 'bad.json'
                for value in ('{', '[]', json.dumps({**data, 'popsize': 'bad'}),
                              json.dumps({**data, 'ngen': 1.5})):
                    bad.write_text(value)
                    run(binary, work, ['--settings', str(bad)], 2)
                bad.write_text('<settings/>')
                run(binary, work, ['--settings', str(bad)], 2)
                run(binary, work, ['--settings'], 2)
                wrong_extension = work / 'settings.txt'
                wrong_extension.write_text(json.dumps(data))
                run(binary, work, ['--settings', str(wrong_extension)], 2)
                if name == 'calibrate':
                    (work / 'calibrate_targets.json').write_text('{"measurements":[{"id":"NODE64","type":"node","values":["bad"]}]}')
                    run(binary, work, ['--settings', str(path)], 2)
                    (work / 'calibrate_targets.json').write_text('{"measurements":[{"id":"NODE64","type":"node","values":[]}]}')
                    run(binary, work, ['--settings', str(path)], 2)
        assert outputs[0] == outputs[1], (name, 'XML/JSON optimization results differ')

    work = base / 'initial values'
    work.mkdir()
    source = root / 'tests/anytown_1med.spr'
    initial = ET.parse(source).getroot()
    values = {'nodes': [], 'edges': []}
    for tag, key, property_name, json_key in (
            ('nodes/node', 'nodes', 'pressure', 'pressure_pa'),
            ('edges/edge', 'edges', 'mass_flow_rate', 'mass_flow_rate_kg_s')):
        for item in initial.findall(tag):
            value = item.findtext(property_name)
            if value is not None:
                values[key].append({'id': item.findtext('id'), json_key: float(value)})
    # End nodes are boundary metadata, not hydraulic unknowns: use an existing junction.
    values['nodes'] = [n for n in values['nodes'] if n['id'] == 'NODE64']
    values['edges'] = []
    path = work / 'initial.JSON'
    path.write_text(json.dumps(values))
    results = []
    for suffix in ('xml', 'json'):
        net = work / (suffix + '.spr')
        shutil.copy2(source, net)
        run(a.staci, work, ['-s', str(net), '-i', str(source if suffix == 'xml' else path)])
        result = ET.parse(net).getroot()
        results.append({n.findtext('id'): float(n.findtext('pressure'))
                        for n in result.findall('nodes/node') if n.findtext('pressure') is not None})
    assert results[0].keys() == results[1].keys()
    assert all(abs(value - results[1][key]) < 100 for key, value in results[0].items())  # < 1 cm head
    path.write_text('{"nodes":[{"id":"missing","pressure_pa":10}]}')
    run(a.staci, work, ['-s', str(net), '-i', str(path)], 2)
print('PASS: XML/JSON optimizer equivalence, JSON measurements/defaults/initial values and input diagnostics.')
