#!/usr/bin/env python3
"""Test pinned public EPANET networks through the standalone staci -s CLI.

Success, compatibility rejection and known nonconvergence are distinct outcomes.
Unexpected failures, signals, timeouts, changed diagnostics and changed source
files always fail the suite. No third-party Python packages or network needed.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time
from epanet_reference import snapshot, compare, ReferenceError, resolve_library

CORPUS = Path(__file__).resolve().parent / 'public_networks'


def run_case(binary, entry, output, timeout, reference_library=None, head_tolerance=None):
    source = CORPUS / entry['path']
    case = output / entry['path'][:-4]
    case.mkdir(parents=True, exist_ok=True)
    result = dict(network=entry['path'], expected=entry['expected'], counts=entry['counts'])
    if hashlib.sha256(source.read_bytes()).hexdigest() != entry['sha256']:
        return dict(result, passed=False, outcome='integrity_error', diagnostic='Upstream input SHA256 changed.')
    network = case / source.name
    shutil.copyfile(source, network)
    marker = Path(str(network) + '.rrs')
    marker.unlink(missing_ok=True)
    started = time.monotonic()
    try:
        with (case / 'staci-output.log').open('w', encoding='utf-8') as log:
            solver_options = ['--head-tolerance-m', str(head_tolerance if head_tolerance is not None else 1e-12), '--mass-tolerance-kg-s', '1e-8'] if reference_library or head_tolerance is not None else []
            if reference_library: solver_options += ['--max-iterations','1000']
            process = subprocess.run([str(binary), '-s', str(network)] + solver_options, cwd=case,
                                     stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
        text = (case / 'staci-output.log').read_text(encoding='utf-8', errors='replace')
        diagnostic_lines = []
        lines = text.splitlines()
        for i, line in enumerate(lines):
            if 'ERROR [' in line:
                diagnostic_lines.append(line[line.index('ERROR ['):])
                for detail in lines[i+1:i+24]:
                    if not detail.startswith('  ') or not detail.strip():
                        break
                    diagnostic_lines.append(detail)
        diagnostic = '\n'.join(diagnostic_lines)
        if process.returncode == 0 and marker.is_file() and marker.read_text().strip() == 'OK':
            outcome = 'solved'
        elif process.returncode == 2 and 'ERROR [EPANET][COMPATIBILITY]' in text:
            outcome = 'compatibility'
        elif process.returncode == 1 and 'ERROR [HYDRAULICS][NONCONVERGENCE]' in text:
            outcome = 'nonconvergence'
        else:
            outcome = 'unexpected_failure'
            diagnostic = diagnostic or '\n'.join(lines[-25:])
        expected = entry['expected']
        passed = outcome == expected or (expected == 'nonconvergence' and entry.get('accept_solved_improvement', False) and outcome == 'solved')
        if outcome in ('compatibility', 'nonconvergence'):
            passed = passed and source.name in diagnostic and all(
                token in diagnostic for token in entry.get('diagnostic_contains', []))
        if outcome == 'nonconvergence':
            passed = passed and all(token in diagnostic for token in
                                    ('RMS head residual=', 'Worst link', 'Worst node'))
        if outcome == 'compatibility' and marker.exists():
            passed = False
        if outcome == 'nonconvergence' and (not marker.exists() or marker.read_text().strip() != 'ERROR!'):
            passed = False
        result['warnings'] = [line[line.index('WARNING ['):] for line in lines if 'WARNING [' in line]
        result.update(returncode=process.returncode, outcome=outcome, passed=passed,
                      diagnostic=diagnostic, log=str(case / 'staci-output.log'))
        if hashlib.sha256(network.read_bytes()).hexdigest() != entry['sha256']:
            result.update(passed=False, diagnostic='STACI modified the EPANET input file.')
    except subprocess.TimeoutExpired:
        result.update(passed=False, outcome='timeout', diagnostic=f'Calculation exceeded {timeout:g} seconds; inspect staci-output.log. A timeout is not an accepted compatibility result.')
    if reference_library:
        reference_path = case / 'epanet-hydraulics.json'
        reference_path.unlink(missing_ok=True)
        try:
            reference = snapshot(reference_library, network, case)
            reference_path.write_text(json.dumps(reference, indent=2) + '\n')
            result['reference_warnings'] = reference['warnings']
            if entry.get('reference_expected') == 'invalid_input': result['passed'] = False
            if result['outcome'] == 'solved':
                actual = json.loads(Path(str(network) + '.hydraulics.json').read_text())
                result['comparison'] = compare(actual, reference)
                if not actual.get('converged'): result['comparison']['passed'] = False
                result['reference_outcome'] = 'match' if result['comparison']['passed'] else 'mismatch'
                result['passed'] = result['passed'] and result['comparison']['passed']
            else:
                result['reference_outcome'] = 'staci_unsupported' if result['outcome'] == 'compatibility' else 'staci_numerical_failure'
        except ReferenceError as error:
            result['reference_outcome'] = 'invalid_input' if error.code >= 200 else 'reference_numerical_failure'
            result['reference_error'] = dict(code=error.code, message=str(error))
            if result['outcome'] == 'solved': result['passed'] = entry.get('reference_expected') == 'invalid_input' and result['reference_outcome'] == 'invalid_input' and error.code == entry.get('reference_error_code')
            if entry.get('reference_expected') == 'invalid_input' and result['reference_outcome'] != 'invalid_input': result['passed'] = False
        except Exception as error:
            result['reference_outcome'] = 'comparison_error'
            result['reference_error'] = str(error)
            result['passed'] = False
        reference_expected = entry.get('reference_expected')
        improved = entry.get('accept_solved_improvement') and result['reference_outcome'] == 'match'
        if reference_expected and result['reference_outcome'] != reference_expected and not improved:
            result['passed'] = False
    else:
        result['reference_outcome'] = 'not_checked'
    result['seconds'] = round(time.monotonic() - started, 3)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--timeout', type=float, default=120)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--head-tolerance-m', type=float, help='Explicit STACI residual limit for testing a selected convergence profile.')
    parser.add_argument('--reference-library', type=Path, help='Official EPANET 2.2 shared library; compare initial-time hydraulics.')
    parser.add_argument('--require-reference', action='store_true')
    parser.add_argument('--reference-only', action='store_true', help='Return skip code 77 when the independent reference is unavailable.')
    parser.add_argument('--network', help='Run only this manifest path.')
    args = parser.parse_args()
    entries = json.loads((CORPUS / 'manifest.json').read_text())['networks']
    if args.network:
        entries = [e for e in entries if e['path'] == args.network]
    if not entries:
        parser.error('No manifest networks selected.')
    if args.jobs < 1 or args.timeout <= 0:
        parser.error('Jobs and timeout must be positive.')
    if args.reference_library and not args.reference_library.is_file():
        parser.error('Explicit reference library does not exist.')
    args.reference_library = resolve_library(args.reference_library)
    if args.reference_only and not args.reference_library:
        print('SKIP: official EPANET library unavailable; numerical equivalence was NOT checked.')
        return 77
    if args.require_reference and not args.reference_library:
        parser.error('--require-reference needs --reference-library.')
    if args.reference_library and not args.reference_library.is_file():
        parser.error('Reference library does not exist.')
    if args.head_tolerance_m is not None and (not math.isfinite(args.head_tolerance_m) or args.head_tolerance_m <= 0):
        parser.error('--head-tolerance-m must be positive and finite.')
    binary = args.binary.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda e: run_case(binary, e, output, args.timeout, args.reference_library, args.head_tolerance_m), entries))
    counts = {k:sum(r['outcome'] == k for r in results) for k in
              ('solved', 'compatibility', 'nonconvergence', 'unexpected_failure', 'timeout', 'integrity_error')}
    failed = sum(not r['passed'] for r in results)
    reference_counts = {k:sum(r['reference_outcome'] == k for r in results) for k in sorted(set(r['reference_outcome'] for r in results))}
    report = dict(head_tolerance_m=args.head_tolerance_m if args.head_tolerance_m is not None else (1e-12 if args.reference_library else 0.0001), reference_outcomes=reference_counts, reference_library=str(args.reference_library) if args.reference_library else None, total=len(results), failed=failed, outcomes=counts, results=results)
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    lines = ['STACI public EPANET steady-snapshot report',
             f'Total: {len(results)}; regression failures: {failed}',
             'Outcomes: ' + ', '.join(f'{k}={v}' for k,v in counts.items()),
             'Reference outcomes: ' + ', '.join(f'{k}={v}' for k,v in reference_counts.items()),
             'Compatibility rejections and known numerical failures are diagnostic tests, not solved networks.', '']
    for r in results:
        lines.append(f"{'PASS' if r['passed'] else 'FAIL'} {r['network']}: {r['outcome']} ({r.get('seconds',0):.3f}s)")
        lines.append('Reference: ' + r['reference_outcome'])
        if r.get('comparison'):
            lines.append(json.dumps(r['comparison'], indent=2))
        if r.get('reference_error'):
            lines.append(str(r['reference_error']))
        if r.get('diagnostic'):
            lines.append(r['diagnostic'])
        lines.extend(r.get('warnings', []))
        if not r['passed']:
            lines.append(f"Expected: {r['expected']}; log: {r.get('log','')}")
    (output / 'report.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines[:4]), flush=True)
    for r in results:
        print(f"{'PASS' if r['passed'] else 'FAIL'} {r['network']}: {r['outcome']}", flush=True)
    print(f'Report: {output / "report.txt"}', flush=True)
    return int(failed != 0)


if __name__ == '__main__':
    sys.exit(main())
