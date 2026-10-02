#!/usr/bin/env python3
"""Verify the four CLI programs' shared GUI diagnostics protocol."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    parser=argparse.ArgumentParser()
    for name in ('staci','split','calibrate','flush'):
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    apps={'staci':args.staci.resolve(),'staci_split':args.split.resolve(),
          'staci_calibrate':args.calibrate.resolve(),'staci_flush':args.flush.resolve()}
    with tempfile.TemporaryDirectory(prefix='staci gui diagnostics ') as directory:
        work=Path(directory)
        log=work/'közös diagnostics.jsonl'
        expected=[]
        def run(program, options, expected_code):
            process=subprocess.run([str(apps[program]), *options, '--diagnostics-file', str(log)],
                                   cwd=work, capture_output=True, text=True, timeout=30)
            assert process.returncode==expected_code,(program,options,process.returncode,process.stdout,process.stderr)
            expected.append((program,expected_code))
            return process
        for name in apps:run(name,['--help'],0)
        run('staci',['-s',str(work/'missing.inp')],2)
        run('staci',['-s'],2)
        run('staci',['--unknown-option'],2)
        run('staci_split',[],2)
        run('staci_calibrate',[],2)
        run('staci_flush',[],2)
        for name in ('staci_split','staci_calibrate'):
            (work/(name+'_settings.xml')).write_text('<settings><global_debug_level>0</global_debug_level></settings>')
            process=run(name,[],2)
            assert 'missing or empty field' in process.stderr,process.stderr
            assert ('n_comm' if name == 'staci_split' else 'Staci_debug_level') in process.stderr
        network=work/'network.inp'
        network.write_text('[JUNCTIONS]\nJ1 0 1\n[RESERVOIRS]\nR1 40\n[PIPES]\nP1 R1 J1 100 100 120\n[OPTIONS]\nUNITS LPS\nHEADLOSS H-W\n[END]\n')
        run('staci',['-s',str(network)],0)
        (work/'hydrants.txt').write_text('J1\n')
        run('staci_flush',['--inp',str(network),'--hydrants',str(work/'hydrants.txt'),
                          '--hydrant-area-m2','.002','--loss-coefficient','2',
                          '--velocity-threshold-mps','.5','--min-pressure-head-m','35',
                          '--output-dir',str(work/'flushing')],3)
        # Parallel invocations must append complete records without truncating
        # the earlier runs. Include all four programs in the same physical file.
        with ThreadPoolExecutor(max_workers=8) as pool:
            list(pool.map(lambda name:run(name,['--help'],0),list(apps)*8))
        records=[json.loads(line) for line in log.read_text().splitlines()]
        groups={}
        for record in records:
            assert record['schema_version']==1 and record['severity'] in ('info','warning','error')
            assert record['program'] in apps and record['code'] and record['message'] and record['timestamp'].endswith('Z')
            groups.setdefault(record['run_id'],[]).append(record)
        assert len(groups)==len(expected),(len(groups),len(expected))
        for events in groups.values():
            assert [e['sequence'] for e in events]==list(range(1,len(events)+1)),events
            assert events[0]['event']=='run_start' and events[-1]['event']=='run_end',events
            assert events[-1]['error_count']==sum(e['severity']=='error' for e in events)
            assert events[-1]['warning_count']==sum(e['severity']=='warning' for e in events)
        ends=[events[-1] for events in groups.values()]
        assert sorted((e['program'],e['exit_code']) for e in ends)==sorted(expected)
        assert any(e['severity']=='warning' for e in records)
        assert any(e['code']=='FLUSH.SCENARIO' for e in records)
        assert any(e['code']=='FLUSH.PARTIAL_FAILURE' for e in records)
        assert any(e['code']=='INPUT.XML' for e in records)
        # The environment option is also shared, without a CLI override.
        import os
        environment=os.environ.copy(); environment['STACI_DIAGNOSTICS_FILE']=str(work/'environment.jsonl')
        p=subprocess.run([str(apps['staci']),'--help'],cwd=work,env=environment,capture_output=True,text=True,timeout=15)
        assert p.returncode==0 and (work/'environment.jsonl').is_file()
        bad=work/'absent-parent'/'diagnostics.jsonl'
        p=subprocess.run([str(apps['staci']),'--help','--diagnostics-file',str(bad)],cwd=work,capture_output=True,text=True,timeout=15)
        assert p.returncode==1 and 'DIAGNOSTICS_OPEN' in p.stderr
        print(f'{len(expected)} shared-file runs verified, including concurrent writers, warnings, errors and partial results.')

if __name__=='__main__':main()
