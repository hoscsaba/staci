"""End-to-end hydraulic assertions; no third-party Python packages required."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

exe, fixture = sys.argv[1:]
with tempfile.TemporaryDirectory(prefix='staci flushing ') as d:
    root = Path(d)
    inp = root / 'network.inp'
    original = Path(fixture).read_bytes()
    inp.write_bytes(original)
    nodes = root / 'nodes.txt'
    nodes.write_text('J1\nJ2\n')

    def run(label, *, model=inp, hydrants=nodes, extra=(), expected=0, threshold='.5'):
        out = root / label
        args = [exe, '--inp', str(model), '--hydrants', str(hydrants),
                '--hydrant-area-m2', '.002', '--loss-coefficient', '2',
                '--velocity-threshold-mps', threshold, '--output-dir', str(out), *extra]
        p = subprocess.run(args, capture_output=True, text=True, timeout=120)
        assert p.returncode == expected, (args, p.returncode, p.stdout, p.stderr)
        return out

    def rows(out, name):
        with (out / name).open() as f:
            return list(csv.DictReader(f))

    out = run('forward')
    scenarios = rows(out, 'scenarios.csv')
    assert len(scenarios) == 2 and all(r['status'] == 'ok' for r in scenarios)
    # Independently solve H-W pipe + outlet head balance by bisection.
    # Use the published resistance form implemented by STACI, independently of Newton.
    for r in scenarios:
        q = float(r['flow_m3s']); h = float(r['hydrant_pressure_head_m'])
        assert math.isclose(h, 2*(q/.002)**2/(2*9.81), abs_tol=1e-5)
        assert float(r['max_mass_residual_kgs']) < 1e-6
        pipe_count = 1 if r['hydrant_id'] == 'J1' else 2
        resistance = pipe_count * 100 / 120**1.85 / .1**4.87 * 7.88 / .85**1.85
        lo, hi = 0., .1
        for _ in range(100):
            mid = (lo + hi)/2
            loss = resistance*mid**1.85 + 2*(mid/.002)**2/(2*9.81)
            if loss > 40: hi = mid
            else: lo = mid
        assert math.isclose(q, (lo+hi)/2, rel_tol=1e-6)
    pipes = rows(out, 'scenario_pipes.csv')
    assert len(pipes) == 4
    j2p2 = next(r for r in pipes if r['hydrant_id']=='J2' and r['pipe_id']=='P2')
    assert float(j2p2['velocity_mps']) < 0 and j2p2['above_threshold']=='1'
    j1p2 = next(r for r in pipes if r['hydrant_id']=='J1' and r['pipe_id']=='P2')
    assert abs(float(j1p2['flow_m3s'])) < 1e-9  # other hydrant stayed closed
    for r in pipes:
        q = float(next(s for s in scenarios if s['hydrant_id']==r['hydrant_id'])['flow_m3s'])
        if r is not j1p2:
            assert math.isclose(abs(float(r['flow_m3s'])), q, abs_tol=1e-8)
    nodes.write_text('J2\nJ1\n')
    reverse = run('reverse')
    indexed = lambda rr: {(r['hydrant_id'],r['pipe_id']):r for r in rr}
    reversed_pipes = indexed(rows(reverse,'scenario_pipes.csv'))
    for key, a in indexed(pipes).items():
        b = reversed_pipes[key]
        for field in ('above_threshold', 'newly_above_threshold'):
            assert a[field] == b[field]
        for field in ('flow_m3s', 'velocity_mps'):
            assert math.isclose(float(a[field]), float(b[field]), rel_tol=1e-6, abs_tol=1e-9)
    pump = root/'closed-pump.inp'
    pump.write_text(original.decode().replace('R1 40', 'R1 40\nR2 100').replace(
        '[OPTIONS]', '[PUMPS]\nPump R2 J1 POWER 10\n[STATUS]\nPump Closed\n[OPTIONS]'))
    pump_out = run('closed-pump', model=pump)
    assert json.loads((pump_out/'run.json').read_text())['initial_pump_states'] == {'Pump': False}
    for r in rows(pump_out, 'scenarios.csv'):
        ref = next(s for s in scenarios if s['hydrant_id']==r['hydrant_id'])
        assert math.isclose(float(r['flow_m3s']), float(ref['flow_m3s']), rel_tol=1e-6)
    assert inp.read_bytes()==original and not Path(str(inp)+'.ros').exists()
    # JSON adapter, including two separate assets sharing a node.
    js = root/'hydrants.json'
    js.write_text(json.dumps({'hydrants':[{'status':'matched','node_id':'J1','dxf_handle':'A'},
                                         {'status':'matched','node_id':'J1','dxf_handle':'B'}]}))
    shared = run('shared', hydrants=js)
    sr = rows(shared,'scenarios.csv')
    assert len(sr)==2 and math.isclose(float(sr[0]['flow_m3s']), float(sr[1]['flow_m3s']), rel_tol=1e-6)
    low = run('low-pressure',extra=('--min-pressure-head-m','35'),expected=2)
    assert all(r['status']=='below_min_pressure' for r in rows(low,'scenarios.csv'))
    assert rows(low,'pipes_above_threshold.csv') == []
    assert rows(low,'scenario_pipes.csv') == []
    # Ranking is adaptive, and opening times include the full reverse-flow route.
    plan = rows(out, 'flushing_plan.csv')
    assert [r['node_id'] for r in plan] == ['J2', 'J1']
    fields = list(plan[0])
    assert fields[fields.index('cumulative_volume_m3')+1] == 'cumulative_volume_percent'
    assert plan[0]['cumulative_volume_percent'] == '100.00'
    for row in plan:
        for key in ('qualifying_volume_m3','additional_volume_m3','cumulative_volume_m3'):
            assert row[key] == f"{float(row[key]):.2f}"
    single = root/'single-hydrant.txt'
    single.write_text('J1\n')
    partial = run('partial-coverage', hydrants=single)
    partial_plan = rows(partial, 'flushing_plan.csv')
    # Denominator includes the uncovered P2, so covering P1 is 50%, not 100%.
    assert len(partial_plan) == 1
    assert partial_plan[0]['cumulative_volume_percent'] == '50.00'
    assert float(plan[1]['additional_volume_m3']) == 0 and plan[1]['redundant'] == '1'
    assert all(float(r['opening_time_min']) > 0 for r in plan)
    for step in plan:
        active = [p for p in pipes if p['node_id']==step['node_id'] and p['above_threshold']=='1']
        expected = sum(100/abs(float(p['velocity_mps'])) for p in active)
        assert step['opening_time_min'] == f"{expected/60:.1f}"
        scenario = next(s for s in scenarios if s['node_id']==step['node_id'])
        flow = float(scenario['flow_m3s'])
        expected_vq = len(active)*math.pi*.1**2/4*100/flow
        assert math.isclose(float(step['hydrant_flow_m3s']), flow, rel_tol=1e-10)
        assert step['volume_over_flow_time_min'] == f"{expected_vq/60:.1f}"
        assert step['volume_over_flow_status']=='volume_over_flow_estimate'
    # The second step is redundant, but its full-volume exchange time is positive.
    assert float(plan[1]['volume_over_flow_time_min']) > 0
    no_coverage = run('no-coverage', threshold='100')
    for step in rows(no_coverage, 'flushing_plan.csv'):
        assert step['volume_over_flow_time_min'] == '0.0'
        assert step['volume_over_flow_status'] == 'no_qualifying_pipes'

    assert [r['node_id'] for r in rows(reverse,'flushing_plan.csv')] == ['J2','J1']
    assert rows(low,'flushing_plan.csv') == []
    assert float(plan[-1]['cumulative_volume_m3']) == float(plan[0]['qualifying_volume_m3'])
    # Unknown node, duplicate node/asset, unsupported hydraulic inputs.
    nodes.write_text('unknown\n');run('unknown',expected=1)
    nodes.write_text('J1\nJ1\n');run('duplicate',expected=1)
    nodes.write_text('J1\n')
    for label, before, after in [
        ('closed','120 0 Open','120 0 Closed'),
        ('minor','120 0 Open','120 3 Open'),
        ('emitter','[END]','[EMITTERS]\nJ1 2\n[END]'),
        ('negative','R1 40','R1 -10'),
    ]:
        bad=root/(label+'.inp');bad.write_text(original.decode().replace(before,after))
        run(label,model=bad,expected=1)
    run('zero-threshold',threshold='0',expected=1)
    # Prevent accidental overwrite of old output.
    run('forward',expected=1)
    manifest=json.loads((out/'run.json').read_text())
    assert manifest['status']=='complete' and manifest['scenario_count']==2

    # Config paths are relative to the config, even when launched elsewhere.
    nodes.write_text('J1\nJ2\n')
    config_dir = root / 'config with spaces'
    config_dir.mkdir()
    config = dict(hydrant_area_m2=.002, total_loss_coefficient=2,
                  velocity_threshold_mps=.5, output_dir='results')
    config_path = config_dir / 'flushing.json'
    def config_run(data, expected=0, extra=()):
        config_path.write_text(json.dumps(data))
        p = subprocess.run([exe, '--inp', str(inp), '--hydrants', str(nodes),
                            '--config', str(config_path), *extra], cwd=root,
                           capture_output=True, text=True, timeout=120)
        assert p.returncode == expected, (p.stdout, p.stderr)
    config_run(config)
    result = json.loads((config_dir / 'results' / 'run.json').read_text())
    assert result['status'] == 'complete' and result['min_pressure_head_m'] == 0
    assert len(rows(config_dir / 'results', 'scenarios.csv')) == 2
    for changes in [dict(hydrant_area_m2='0.002'), dict(total_loss_coefficient=True),
                    dict(typo=1), dict(output_dir=''), dict(min_pressure_head_m=-1)]:
        config_run({**config, 'output_dir': 'invalid', **changes}, expected=1)
    missing = dict(config); del missing['hydrant_area_m2']
    config_run(missing, expected=1)
    config_run(dict(config, output_dir='duplicate'), expected=1,
               extra=('--loss-coefficient', '2'))
    config_run(dict(config, output_dir='filtered', min_pressure_head_m=35), expected=2)

    # Node IDs now live in config; export geometry/marking is independent of hydraulics.
    inp.write_text(original.decode().replace('[END]',
        '[COORDINATES]\nR1 0 0\nJ1 100 0\nJ2 200 0\n[VERTICES]\nP2 150 5\n[END]'))
    inline = dict(config, output_dir='exports', hydrant_node_ids=['J1', 'J2'],
                  write_network_files=True)
    def inline_run(data, expected=0):
        config_path.write_text(json.dumps(data))
        p = subprocess.run([exe, '--inp', str(inp), '--config', str(config_path)],
                           capture_output=True, text=True, timeout=120)
        assert p.returncode == expected, (p.stdout, p.stderr)
        return p
    exported = inline_run(inline)
    export_dir = config_dir / 'exports'
    assert 'qualifying pipe volume:' in exported.stdout
    assert 'Unique pipe volume' in (export_dir/'summary.txt').read_text()
    index = rows(export_dir, 'networks.csv')
    assert len(index) == 2
    def sections(path):
        result = {}; section = ''
        for line in path.read_text().splitlines():
            line = line.split(';')[0].strip()
            if not line: continue
            if line.startswith('['):
                section=line.strip('[]'); result.setdefault(section, [])
            else: result.setdefault(section, []).append(line.split())
        return result
    for entry in index:
        assert entry['network_file'] == f"networks/network_hydrant_{entry['node_id']}.inp"
        path = export_dir/entry['network_file']
        text = path.read_text()
        assert text.index('[JUNCTIONS]') < text.index('[PIPES]')
        doc = sections(path)
        assert doc['COORDINATES'] == [['R1','0','0'], ['J1','100','0'], ['J2','200','0']]
        assert doc['VERTICES'] == [['P2','150','5']]
        assert len(doc['EMITTERS']) == 1 and doc['EMITTERS'][0][0] == entry['node_id']
        assert math.isclose(float(doc['EMITTERS'][0][1]), .002*math.sqrt(2*9.81/2)*1000, rel_tol=1e-10)
        tags = {row[1]: row[2] for row in doc['TAGS'] if row[0]=='LINK'}
        qualifying = [row for row in rows(export_dir,'pipes_above_threshold.csv') if row['node_id']==entry['node_id']]
        assert {pipe for pipe,tag in tags.items() if tag=='FLUSHED'} == {row['pipe_id'] for row in qualifying}
        scenario = next(row for row in rows(export_dir,'scenarios.csv') if row['node_id']==entry['node_id'])
        assert scenario['qualifying_volume_m3'] == f"{len(qualifying)*math.pi*.1**2/4*100:.2f}"
        assert doc['TIMES'] == [['DURATION','0']]
        assert 'CONTROLS' not in doc and 'RULES' not in doc
    inline_run(dict(inline, output_dir='exports-off', write_network_files=False))
    assert not (config_dir/'exports-off'/'networks').exists()
    for change in [dict(hydrant_node_ids=[]), dict(hydrant_node_ids=['J1','J1']),
                   dict(hydrant_node_ids=['missing']), dict(hydrant_node_ids=[1]),
                   dict(write_network_files='true')]:
        inline_run({**inline, 'output_dir':'bad-inline', **change}, expected=1)
    inline_run(dict(inline, output_dir='invalid-exports', min_pressure_head_m=35), expected=2)
    for entry in rows(config_dir/'invalid-exports', 'networks.csv'):
        doc = sections(config_dir/'invalid-exports'/entry['network_file'])
        assert all(row[2]=='INVALID_SCENARIO' for row in doc['TAGS'] if row[0]=='LINK')
    # Supplying both node IDs and an external list is ambiguous and rejected.
    config_run(inline, expected=1)

print('Flushing end-to-end checks passed')
