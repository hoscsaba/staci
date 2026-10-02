#!/usr/bin/env python3
"""Check CLI diagnostics for unsupported physics and malformed hydraulic inputs."""
import argparse
from pathlib import Path
import subprocess
import tempfile

BASE = '''[JUNCTIONS]
J1 0 1
[RESERVOIRS]
R1 30
[PIPES]
P1 R1 J1 100 150 120 0 OPEN
[OPTIONS]
UNITS LPS
HEADLOSS H-W
[END]
'''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve()
    cases = [
        ('unknown_node', BASE.replace('P1 R1 J1', 'P1 R1 MISSING'), ('[PIPES]', "'P1'", "'MISSING'")),
        ('bad_number', BASE.replace('100 150 120', '100 NaN 120'), ('[PIPES]', 'finite number', 'NaN')),
        ('zero_diameter', BASE.replace('100 150 120', '100 0 120'), ('Pipe diameter', 'positive')),
        ('negative_length', BASE.replace('100 150 120', '-100 150 120'), ('Pipe length', 'non-negative')),
        ('duplicate_node', BASE.replace('J1 0 1', 'J1 0 1\nJ1 1 2'), ('Duplicate node ID', "'J1'")),
        ('duplicate_link', BASE.replace('[OPTIONS]', 'P1 R1 J1 100 150 120\n[OPTIONS]'), ('Duplicate link ID', "'P1'")),
        ('unknown_units', BASE.replace('UNITS LPS', 'UNITS UNKNOWN'), ('Unknown flow units', "'UNKNOWN'")),
        ('unsupported_headloss', BASE.replace('HEADLOSS H-W', 'HEADLOSS C-M'), ('Head-loss formula', "'C-M'", 'not implemented')),
        ('pda', BASE.replace('[END]', 'DEMAND MODEL UNKNOWN\n[END]'), ('Pressure-dependent demand', 'UNKNOWN')),
        ('emitter', BASE.replace('[END]', '[EMITTERS]\nJ1 -0.2\n[END]'), ('[EMITTERS]', "'J1'", 'non-negative')),
        ('emitter_unknown', BASE.replace('[END]', '[EMITTERS]\nABSENT 1\n[END]'), ('[EMITTERS]', 'ABSENT', 'unknown junction')),
        ('emitter_nan', BASE.replace('[END]', '[EMITTERS]\nJ1 NaN\n[END]'), ('[EMITTERS]', 'finite number')),
        ('emitter_exponent', BASE.replace('[END]', 'EMITTER EXPONENT 0\n[END]'), ('[OPTIONS]', 'exponent', 'positive')),
        ('pda_range', BASE.replace('[END]', 'DEMAND MODEL PDA\nMINIMUM PRESSURE 20\nREQUIRED PRESSURE 10\n[END]'), ('[OPTIONS]', 'Required pressure', 'minimum pressure')),
        ('unknown_section', BASE.replace('[END]', '[FOO]\nJ1 0.2\n[END]'), ('[FOO]', 'Unknown EPANET section')),
        ('disconnected', BASE.replace('120 0 OPEN', '120 0 CLOSED'), ('[TOPOLOGY]', "'J1'", 'no reservoir or tank')),
        ('missing_curve', BASE.replace('[OPTIONS]', '[PUMPS]\nPU1 R1 J1 HEAD ABSENT\n[OPTIONS]'), ('[PUMPS]', "'PU1'", "'ABSENT'", 'not found')),
    ]
    cases += [
        ('gpv_nan', BASE.replace('[OPTIONS]', '[VALVES]\nV1 R1 J1 150 GPV C 0\n[CURVES]\nC 0 NaN\nC 1 1\n[OPTIONS]'), ('[CURVES]', 'GPV curve headloss', 'finite number')),
        ('gpv_duplicate', BASE.replace('[OPTIONS]', '[VALVES]\nV1 R1 J1 150 GPV C 0\n[CURVES]\nC 1 1\nC 1 2\n[OPTIONS]'), ('[VALVES]', 'strictly increasing flow')),
        ('gpv_decreasing', BASE.replace('[OPTIONS]', '[VALVES]\nV1 R1 J1 150 GPV C 0\n[CURVES]\nC 0 2\nC 1 1\n[OPTIONS]'), ('[VALVES]', 'nondecreasing headloss')),
        ('gpv_missing', BASE.replace('[OPTIONS]', '[VALVES]\nV1 R1 J1 150 GPV ABSENT 0\n[OPTIONS]'), ('GPV headloss curve', 'ABSENT', 'two points')),
        ('fcv_negative', BASE.replace('[OPTIONS]', '[VALVES]\nV1 R1 J1 150 FCV -1 0\n[OPTIONS]'), ('Valve setting', 'non-negative')),
        ('fcv_nan', BASE.replace('[OPTIONS]', '[VALVES]\nV1 R1 J1 150 FCV NaN 0\n[OPTIONS]'), ('Valve setting', 'finite number')),
    ]
    with tempfile.TemporaryDirectory(prefix='staci-diagnostics-') as temporary:
        work = Path(temporary)
        for name, text, tokens in cases:
            path = work / (name + '.inp')
            path.write_text(text)
            Path(str(path) + ".rrs").write_text("OK\n")
            Path(str(path) + ".hydraulics.json").write_text('{"converged":true}\n')
            process = subprocess.run([str(binary), '-s', str(path)], cwd=work,
                                     capture_output=True, text=True, timeout=15)
            output = process.stdout + process.stderr
            assert process.returncode == 2, (name, process.returncode, output)
            assert path.name in process.stderr and 'ERROR [EPANET][COMPATIBILITY]' in process.stderr, output
            assert all(token in process.stderr for token in tokens), (name, tokens, output)
            assert 'line ' in process.stderr or '[TOPOLOGY]' in process.stderr, output
            assert not Path(str(path) + '.rrs').exists(), name
            assert not Path(str(path) + '.hydraulics.json').exists(), name
            assert path.read_text() == text, name
            print('PASS', name)
        # EPANET allows a node and a link to use the same ID.
        path = work / 'independent_namespaces.inp'
        path.write_text(BASE.replace('P1 R1 J1', 'J1 R1 J1'))
        process = subprocess.run([str(binary), '-s', str(path)], cwd=work,
                                 capture_output=True, text=True, timeout=15)
        assert process.returncode == 0, process.stdout + process.stderr
        assert Path(str(path) + '.rrs').read_text().strip() == 'OK'
        print('PASS independent_namespaces')
        for value in ('0', '0.0', '0e0'):
            path = work / ('zero-emitter-' + value + '.inp')
            path.write_text(BASE.replace('[END]', '[EMITTERS]\nJ1 ' + value + '\n[END]'))
            process = subprocess.run([str(binary), '-s', str(path)], cwd=work,
                                     capture_output=True, text=True, timeout=15)
            assert process.returncode == 0, process.stdout + process.stderr
            assert Path(str(path) + '.rrs').read_text().strip() == 'OK'
            print('PASS zero emitter', value)

        path = work / 'stopped-pump-empty-tank.inp'
        path.write_text("""[JUNCTIONS]
J1 0 1
J2 0 1
[RESERVOIRS]
R1 30
[TANKS]
T1 20 0 0 10 10 0
[PIPES]
P1 J1 J2 100 150 120 0 OPEN
P2 T1 J2 100 150 120 0 OPEN
[PUMPS]
PU R1 J1 HEAD C1 PATTERN OFF
[CURVES]
C1 5 20
[PATTERNS]
OFF 0
[OPTIONS]
UNITS LPS
HEADLOSS H-W
[END]
""")
        process = subprocess.run([str(binary), '-s', str(path)], cwd=work,
                                 capture_output=True, text=True, timeout=15)
        text = process.stdout + process.stderr
        assert process.returncode == 1, text
        assert 'EPANET.NO_AVAILABLE_SUPPLY' in text and 'zero-speed pumps' in text, text
        assert 'minimum level' in text and '2 nodes with positive demand' in text, text
        assert 'Initial guesses or relaxation cannot replace the missing supply' in text, text
        print('PASS stopped pump / empty tank supply diagnostic')

    print(f'{len(cases)+5} diagnostic/namespace/emitter/supply checks passed')


if __name__ == '__main__':
    main()
