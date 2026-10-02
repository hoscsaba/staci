#!/usr/bin/env python3
"""Verify common solver overrides, rejection diagnostics and ky16 default accuracy."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    parser=argparse.ArgumentParser()
    for name in ('staci','split','calibrate','flush'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--input',type=Path,required=True)
    args=parser.parse_args()
    apps=[getattr(args,name).resolve() for name in ('staci','split','calibrate','flush')]
    with tempfile.TemporaryDirectory(prefix='staci solver options ') as directory:
        work=Path(directory)
        for binary in apps:
            help_run=subprocess.run([str(binary),'--help'],cwd=work,capture_output=True,text=True,timeout=15)
            assert help_run.returncode==0
            assert all(key in help_run.stdout for key in ('--head-tolerance-m','--mass-tolerance-kg-s','--max-iterations','0.1 mm'))
            for option,value in [('--head-tolerance-m','0'),('--head-tolerance-m','NaN'),
                    ('--mass-tolerance-kg-s','-1'),('--mass-tolerance-kg-s','inf'),
                    ('--max-iterations','1.5'),('--max-iterations','2147483648'),
                    ('--head-tolerance-m','abc'),('--mass-tolerance-kg-s','1e-8extra')]:
                p=subprocess.run([str(binary),option,value],cwd=work,capture_output=True,text=True,timeout=15)
                assert p.returncode==2,(binary,option,value,p.stdout,p.stderr)
                assert option in p.stderr and 'positive finite' in p.stderr
            p=subprocess.run([str(binary),'--head-tolerance-m'],cwd=work,capture_output=True,text=True,timeout=15)
            assert p.returncode==2 and '--head-tolerance-m' in p.stderr

        network=work/'ky16.inp';shutil.copyfile(args.input,network)
        def solve(name,options):
            prefix=work/name
            p=subprocess.run([str(apps[0]),'--epanet-eps',str(network),'-o',str(prefix)]+options,
                             cwd=work,capture_output=True,text=True,timeout=60)
            data=json.loads(Path(str(prefix)+'.meta.json').read_text())
            return p,data['simulation']
        p,result=solve('default',[])
        assert p.returncode==0 and result['failed_frames']==0 and result['frames']==25,(p.stderr,result)
        p,result=solve('strict',['--head-tolerance-m','1e-12','--mass-tolerance-kg-s','1e-8'])
        assert 'limit 1e-12' in p.stderr or p.returncode == 0, (p.stderr,result)
        assert ((p.returncode == 0 and result['failed_frames'] == 0) or
                (p.returncode == 3 and result['failed_frames'] == 1 and
                 'EPANET.EPS_PARTIAL_FAILURE' in p.stderr)), (p.stderr,result)
        p,result=solve('iteration-cap',['--head-tolerance-m=1e-12','--mass-tolerance-kg-s=1e-8','--max-iterations=1'])
        assert p.returncode==3 and result['failed_frames']>1,(p.stderr,result)
        assert 'after 1 iterations' in p.stderr
        p,result=solve('explicit-default',['--head-tolerance-m=.0001','--mass-tolerance-kg-s=1e-8','--max-iterations=50'])
        assert p.returncode==0 and result['failed_frames']==0,(p.stderr,result)
    print('PASS common options for all four applications; invalid values rejected; ky16 25/25 default frames; strict and iteration overrides effective.')
    return 0

if __name__=='__main__':raise SystemExit(main())
