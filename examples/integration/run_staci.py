"""Isolated CLI adapter for a Python service or desktop GUI (Python 3.9+)."""
import json
import shutil
import subprocess
import uuid
from pathlib import Path


def run_network(binary, network, job_root, *, eps=False, timeout=120):
    binary = Path(binary).resolve(strict=True)
    network = Path(network).resolve(strict=True)
    job = Path(job_root).resolve() / uuid.uuid4().hex
    job.mkdir(parents=True)
    copied = job / ('network' + network.suffix.lower())
    shutil.copy2(network, copied)
    diagnostics = job / 'diagnostics.jsonl'
    args = [str(binary), '--diagnostics-file', str(diagnostics)]
    if eps:
        args += ['--epanet-eps', str(copied), '-o', str(job / 'result')]
    else:
        args += ['-s', str(copied)]
    timed_out = False
    with (job / 'console.log').open('wb') as log:
        try:
            process = subprocess.run(args, cwd=job, stdout=log,
                                     stderr=subprocess.STDOUT, timeout=timeout,
                                     check=False, shell=False)
            exit_code = process.returncode
        except subprocess.TimeoutExpired:
            timed_out = True
            exit_code = None
    records = []
    if diagnostics.exists():
        for line in diagnostics.read_text(encoding='utf-8').splitlines():
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                # An interrupted writer can leave a final incomplete record.
                records.append({'severity': 'warning', 'code': 'ADAPTER.INCOMPLETE_LOG',
                                'message': 'Incomplete diagnostic record.'})
    starts = [r for r in records if r.get('event') == 'run_start']
    run_id = starts[-1].get('run_id') if starts else None
    end = next((r for r in reversed(records)
                if r.get('event') == 'run_end' and r.get('run_id') == run_id), None)
    complete = end is not None and end.get('exit_code') == exit_code
    state = ('timeout' if timed_out else 'partial' if complete and exit_code == 3
             else 'success' if complete and exit_code == 0 else 'failed')
    result = job / ('result.meta.json' if eps else copied.name + '.hydraulics.json')
    return {'state': state, 'exit_code': exit_code, 'run_id': run_id,
            'job_dir': str(job), 'diagnostics': records,
            'result_file': str(result) if state in ('success', 'partial') and result.exists() else None}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', required=True)
    parser.add_argument('--network', required=True)
    parser.add_argument('--job-root', required=True)
    parser.add_argument('--eps', action='store_true')
    parser.add_argument('--timeout', type=float, default=120)
    options = parser.parse_args()
    response = run_network(options.binary, options.network, options.job_root,
                           eps=options.eps, timeout=options.timeout)
    print(json.dumps(response, ensure_ascii=False, indent=2))
    raise SystemExit(0 if response['state'] == 'success' else 1)
