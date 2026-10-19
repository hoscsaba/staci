"""Create Anytown snapshot calibration inputs and synthetic STACI measurements."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

FACTORS = (0.75, 1.0, 1.15)
DT_HOURS = 0.01
PIPE_FACTORS = {'1': 1.10, '2': 0.90}
MEASURED_NODES = ('1', '4', '8', '12', '17', '19')


def records(source):
    result = {}
    section = ''
    for line in source.read_text().splitlines():
        line = line.split(';', 1)[0].strip()
        if line.startswith('['):
            section = line.strip('[]').upper()
            result.setdefault(section, [])
        elif line:
            result.setdefault(section, []).append(line.split())
    return result


def snapshot(data, factor, levels, perturb):
    lines = ['[TITLE]', 'Anytown calibration snapshot; synthetic experiment, not field data']
    for section in ('JUNCTIONS', 'RESERVOIRS', 'TANKS', 'PIPES', 'PUMPS', 'CURVES', 'OPTIONS', 'COORDINATES', 'VERTICES'):
        lines.append('\n['+section+']')
        for original in data.get(section, []):
            row = original.copy()
            if section == 'JUNCTIONS':
                row = row[:3]
                row[2] = str(float(row[2])*factor)
            elif section == 'RESERVOIRS':
                row = row[:2]
            elif section == 'TANKS':
                row[2] = str(levels[row[0]] / 0.3048)  # Anytown uses feet.
            elif section == 'PIPES' and perturb and row[0] in PIPE_FACTORS:
                row[4] = str(float(row[4])*PIPE_FACTORS[row[0]])
            elif section == 'PUMPS':
                row = row[:5] + ['SPEED', '1.0']
            elif section == 'OPTIONS' and row[0].upper() in ('QUALITY', 'PATTERN', 'DEMAND'):
                continue
            lines.append(' '.join(row))
    lines.extend(['\n[TIMES]', 'DURATION 0', '\n[END]', ''])
    return '\n'.join(lines)


def prepare(binary, source, output):
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output directory must be new or empty')
    output.mkdir(parents=True, exist_ok=True)
    data = records(source)
    levels = {row[0]: 20*0.3048 for row in data['TANKS']}
    initial = levels.copy()
    areas = {row[0]: math.pi*(float(row[5])*0.3048)**2/4 for row in data['TANKS']}
    measurements = [{'id': n, 'type': 'node', 'values': []} for n in MEASURED_NODES]
    measurements += [{'id': n, 'type': 'pool', 'values': []} for n in levels]
    with tempfile.TemporaryDirectory(prefix='anytown truth ') as tmp:
        for period, factor in enumerate(FACTORS):
            truth = Path(tmp)/f'truth_{period}.inp'
            truth.write_text(snapshot(data, factor, levels, False))
            run = subprocess.run([str(binary), '-s', str(truth)], cwd=tmp,
                                 capture_output=True, text=True, timeout=60)
            if run.returncode:
                raise RuntimeError(run.stdout[-2000:]+run.stderr)
            result = json.loads(Path(str(truth)+'.hydraulics.json').read_text())
            if not result['converged']:
                raise RuntimeError('Synthetic reference did not converge')
            for row in measurements:
                row['values'].append(result['nodes'][row['id']]['pressure_m']
                                     if row['type']=='node' else levels[row['id']])
            (output/f'network_{period}.inp').write_text(snapshot(data, factor, initial, True))
            for node in levels:
                levels[node] += result['links']['EPANET_TANK_'+node]['flow_m3s'] * DT_HOURS*3600/areas[node]
    (output/'measurements.json').write_text(json.dumps({'measurements': measurements}, indent=2)+'\n')
    settings = dict(global_debug_level=1, Staci_debug_level=0, dir_name='./',
                    fname_prefix='network_', network_extension='.inp',
                    logfilename='calibrate.log', best_logfilename='calibrate-best.log',
                    sollwert_dfile='measurements.json', Start_of_Periods=0, Num_of_Periods=3,
                    Spoil_Active_Pipes=False, dt=DT_HOURS, weight_p_err=.5,
                    type_of_pipe_selection='list', pipe_ids=list(PIPE_FACTORS),
                    popsize=20, ngen=25, pmut=.2, pcross=.8)
    (output/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    truth_d = {row[0]: float(row[4])*.0254 for row in data['PIPES'] if row[0] in PIPE_FACTORS}
    (output/'provenance.json').write_text(json.dumps(dict(
        source='tests/public_networks/wntr/wntr/tests/networks_for_testing/Anytown.inp',
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        measurements='Synthetic STACI steady solutions with sequential constant-area tank balance; not independent validation or field data',
        demand_multipliers=FACTORS, initial_tank_level_m=initial,
        true_diameters_m=truth_d, input_diameter_factors=PIPE_FACTORS,
        counts={k:len(data[k]) for k in ('JUNCTIONS','RESERVOIRS','TANKS','PIPES','PUMPS')}
    ),indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--staci',type=Path,required=True)
    parser.add_argument('--source',type=Path,default=Path(__file__).with_name('anytown.inp'))
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    prepare(args.staci.resolve(), args.source.resolve(), args.output.resolve())
