#!/usr/bin/env python3
"""Verify physical adapted Anytown and zero-speed HEAD/POWER pump behavior."""
import argparse,csv,hashlib,json,math,shutil,subprocess,tempfile
from pathlib import Path
from epanet_reference import snapshot,compare,resolve_library


def stopped_pump(kind,mode,speed=None):
    definition='HEAD HC' if kind=='HEAD' else 'POWER 1'
    if speed is None:speed=0 if mode=='direct' else 1
    pattern=' PATTERN PS' if mode in ('pattern','eps') else ''
    status='[STATUS]\nPU 0\n' if mode=='status' else ''
    values='0 1 0' if mode=='eps' else '0'
    return f'''[JUNCTIONS]
J1 0 1
J2 0 1
[RESERVOIRS]
R1 10
R2 50
[PIPES]
BYPASS R2 J2 100 150 120 0 OPEN
P2 J2 J1 100 150 120 0 OPEN
[PUMPS]
PU R1 J1 {definition} SPEED {speed}{pattern}
{status}[PATTERNS]
PS {values}
[CURVES]
HC 0 60
HC 5 40
HC 10 0
[TIMES]
DURATION 2:00
HYDRAULIC TIMESTEP 1:00
PATTERN TIMESTEP 1:00
REPORT TIMESTEP 1:00
[OPTIONS]
UNITS LPS
HEADLOSS H-W
[END]
'''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--anytown',type=Path,required=True)
    parser.add_argument('--reference-library',type=Path)
    args=parser.parse_args();library=resolve_library(args.reference_library)
    if not library:print('SKIP: independent EPANET reference unavailable.');return 77
    binary=args.binary.resolve();checks=0
    with tempfile.TemporaryDirectory(prefix='staci physical models ') as directory:
        work=Path(directory);model=work/'Anytown.inp';shutil.copyfile(args.anytown,model)
        digest=hashlib.sha256(model.read_bytes()).hexdigest()
        def run(command):
            p=subprocess.run([str(binary)]+command,cwd=work,capture_output=True,text=True,timeout=120)
            assert p.returncode==0,(command,p.stdout,p.stderr)
            return p
        # Independent snapshot checks at strict and normal tolerances.
        reference=snapshot(library,model,work)
        assert not reference['warnings'],reference['warnings']
        for options in ([],['--head-tolerance-m','1e-12']):
            run(['-s',str(model)]+options)
            actual=json.loads(Path(str(model)+'.hydraulics.json').read_text())
            result=compare(actual,reference);assert result['passed'],result['failures']
            pressure=[n['pressure_m'] for n in actual['nodes'].values() if n['demand_m3s']>0]
            assert min(pressure)>20 and max(pressure)<100,pressure
            assert all(actual['links'][name]['enabled'] and actual['links'][name]['flow_m3s']>0 for name in ('78','79','80'))
        prefix=work/'anytown-eps';run(['--epanet-eps',str(model),'-o',str(prefix),'--head-tolerance-m','1e-9'])
        meta=json.loads(Path(str(prefix)+'.meta.json').read_text())['simulation']
        assert meta['frames']==25 and meta['hydraulic_states']>=1441 and meta['failed_hydraulic_states']==0,meta
        rows=list(csv.DictReader(Path(str(prefix)+'-nodes.csv').open()))
        junction_rows=[r for r in rows if r['node_id'] not in ('40','41','42')]
        assert min(float(r['pressure_head_m']) for r in junction_rows)>20
        assert max(float(r['pressure_head_m']) for r in junction_rows)<100
        links=list(csv.DictReader(Path(str(prefix)+'-links.csv').open()))
        assert max(abs(float(r['velocity_mps'])) for r in links if r['type']=='PIPE')<3
        for row in csv.DictReader(Path(str(prefix)+'-tanks.csv').open()):
            assert float(row['min_level_m'])-1e-8<=float(row['level_m'])<=float(row['max_level_m'])+1e-8,row
        assert hashlib.sha256(model.read_bytes()).hexdigest()==digest
        print('PASS adapted Anytown: default/strict EPANET agreement; at least 1441 hydraulic states, positive pressure and tank limits.')
        for kind in ('HEAD','POWER'):
            for mode in ('direct','pattern','status'):
                network=work/f'{kind}-{mode}.inp';network.write_text(stopped_pump(kind,mode))
                run(['-s',str(network),'--head-tolerance-m','1e-12'])
                actual=json.loads(Path(str(network)+'.hydraulics.json').read_text());ref=snapshot(library,network,work)
                comparison=compare(actual,ref);assert comparison['passed'],(kind,mode,comparison['failures'])
                assert not actual['links']['PU']['enabled'] and abs(actual['links']['PU']['flow_m3s'])<1e-10
                assert actual['nodes']['J1']['pressure_m']>0
                checks+=1
            network=work/f'{kind}-eps.inp';network.write_text(stopped_pump(kind,'eps'));prefix=work/f'{kind}-eps'
            run(['--epanet-eps',str(network),'-o',str(prefix),'--head-tolerance-m','1e-12'])
            meta=json.loads(Path(str(prefix)+'.meta.json').read_text())['simulation'];assert meta['failed_frames']==0 and meta['frames']==3,meta
            nodes=list(csv.DictReader(Path(str(prefix)+'-nodes.csv').open()));links=list(csv.DictReader(Path(str(prefix)+'-links.csv').open()))
            for t,speed in ((0,0),(3600,1),(7200,0)):
                refnetwork=work/f'{kind}-static-{t}.inp';refnetwork.write_text(stopped_pump(kind,'direct',speed))
                ref=snapshot(library,refnetwork,work)
                actual={'nodes':{},'links':{}}
                for row in nodes:
                    if int(row['time_seconds'])==t:
                        actual['nodes'][row['node_id']]={'head_m':float(row['total_head_m']),'pressure_m':float(row['pressure_head_m']),'demand_m3s':float(row['demand_m3s'])}
                for row in links:
                    if int(row['time_seconds'])==t:
                        actual['links'][row['link_id']]={'flow_m3s':float(row['flow_m3s']),'velocity_mps':abs(float(row['velocity_mps'])),'enabled':int(row['status'])!=0}
                for name in ('R1','R2'):
                    actual['nodes'][name]['demand_m3s']=math.fsum((1 if link['to']==name else -1)*actual['links'][identifier]['flow_m3s'] for identifier,link in ref['links'].items() if name in (link['from'],link['to']))
                result=compare(actual,ref);assert result['passed'],(kind,t,result['failures'])
                assert actual['links']['PU']['enabled']==bool(speed)
                if not speed:assert abs(actual['links']['PU']['flow_m3s'])<1e-10
                else:assert actual['links']['PU']['flow_m3s']>0
                checks+=1
        print(f'PASS {checks} zero-speed HEAD/POWER independent cases, including restart/shutdown, status and patterns.')
    return 0

if __name__=='__main__':raise SystemExit(main())
