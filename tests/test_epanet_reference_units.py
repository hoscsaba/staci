#!/usr/bin/env python3
"""Independent EPANET checks for all flow units and hydraulic import semantics."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile
from epanet_reference import snapshot, compare, resolve_library, validate_reference_balance, ReferenceError


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--reference-library',type=Path)
    args=parser.parse_args();library=resolve_library(args.reference_library)
    if not library:
        print('SKIP: independent EPANET reference unavailable.');return 77
    # Conservation audit must reject nominal convergence with unserved demand,
    # including signed/reverse flow; valid source-boundary exchange is exempt.
    nodes = {'R': {'type': 1, 'head_m': 20., 'pressure_m': 0., 'demand_m3s': -.001},
             'J': {'type': 0, 'head_m': 10., 'pressure_m': 10., 'demand_m3s': .001}}
    links = {'P': {'from': 'J', 'to': 'R', 'flow_m3s': -.001}}
    validate_reference_balance(nodes, links)
    for flow in (0., float('nan')):
        broken = copy.deepcopy(links); broken['P']['flow_m3s'] = flow
        try:
            validate_reference_balance(nodes, broken)
        except ReferenceError as error:
            assert error.code == 110
        else:
            raise AssertionError('Unreliable reference balance was accepted.')
    # Unit factors are specified here independently of the reference reader.
    units={'CFS':.028316846592,'GPM':.0000630901964,'MGD':.0438126363889,
        'IMGD':.0526167824074,'AFD':.0142764101852,'LPS':.001,'LPM':1/60000,
        'MLD':1000/86400,'CMH':1/3600,'CMD':1/86400,'SI':.001}
    cases=[(unit,'H-W',.001,False,False) for unit in units]
    cases += [('LPS','D-W',flow,False,False) for flow in (1e-5,2e-4,1e-3)]
    cases += [('LPS','H-W',.001,True,False),('SI','D-W',.001,False,True)]
    binary=args.binary.resolve();checks=0
    with tempfile.TemporaryDirectory(prefix='staci reference units ') as directory:
        for number,(unit,formula,flow,pattern,legacy) in enumerate(cases):
            work=Path(directory)/str(number);work.mkdir()
            us=unit in {'CFS','GPM','MGD','IMGD','AFD'}
            length=.3048 if us else 1.0
            diameter=.1/.0254 if us else 100.0
            roughness=120 if formula=='H-W' else .01
            source_section='TANKS' if legacy else 'RESERVOIRS'
            pattern_text='[PATTERNS]\n1 0.5 2\n' if pattern else ''
            text=f'''[JUNCTIONS]
J1 {5/length:.15g} 0
J2 {5/length:.15g} {flow/units[unit]:.15g}
[{source_section}]
R {20/length:.15g}
[PIPES]
P1 R J1 {100/length:.15g} {diameter:.15g} {roughness} 2 OPEN
P2 J1 J2 {100/length:.15g} {diameter:.15g} {roughness} 0 OPEN
{pattern_text}[OPTIONS]
UNITS {unit}
HEADLOSS {formula}
[END]
'''
            network=work/'network.inp';network.write_text(text)
            digest=hashlib.sha256(network.read_bytes()).hexdigest()
            run=subprocess.run([str(binary),'-s',str(network),'--head-tolerance-m','1e-12','--mass-tolerance-kg-s','1e-8'],cwd=work,capture_output=True,text=True,timeout=30)
            assert run.returncode==0,(unit,formula,run.stdout,run.stderr)
            assert hashlib.sha256(network.read_bytes()).hexdigest()==digest
            actual=json.loads(Path(str(network)+'.hydraulics.json').read_text())
            reference=snapshot(library,network,work)
            expected_flow=flow*(.5 if pattern else 1)
            assert math.isclose(reference['links']['P1']['flow_m3s'],expected_flow,rel_tol=1e-6,abs_tol=1e-10)
            assert math.isclose(reference['nodes']['R']['head_m'],20,abs_tol=1e-8)
            result=compare(actual,reference);assert result['passed'],(unit,formula,result['failures'])
            # Demonstrate the numerical test cannot pass an incorrect solution.
            for kind,name,field,amount in [('nodes','J2','head_m',1),('links','P1','flow_m3s',.01),
                                           ('links','P1','velocity_mps',1)]:
                changed=copy.deepcopy(actual);changed[kind][name][field]+=amount
                assert not compare(changed,reference)['passed']
            changed=copy.deepcopy(actual);changed['links']['P1']['enabled']=False
            assert not compare(changed,reference)['passed']
            checks+=1
    print(f'PASS {checks} EPANET unit/pipe/pattern/legacy cases and deliberate wrong-result rejection checks.')
    return 0

if __name__=='__main__':raise SystemExit(main())
