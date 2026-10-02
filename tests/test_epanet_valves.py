#!/usr/bin/env python3
"""Independent hydraulic checks for PRV, FCV and GPV operating states and units."""
import argparse, csv, hashlib, json, math, subprocess, tempfile
from pathlib import Path
from epanet_reference import snapshot,compare,resolve_library

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--reference-library',type=Path)
    args=parser.parse_args();library=resolve_library(args.reference_library)
    if not library: print('SKIP: EPANET reference unavailable.');return 77
    units={'CFS':.028316846592,'GPM':.0000630901964,'MGD':.0438126363889,
      'IMGD':.0526167824074,'AFD':.0142764101852,'LPS':.001,'LPM':1/60000,
      'MLD':1000/86400,'CMH':1/3600,'CMD':1/86400}
    scenarios=[('PRV','active'),('PRV','open'),('PRV','forced_open'),('PRV','reverse'),('PRV','closed'),
      ('FCV','active'),('FCV','unattainable'),('FCV','reverse'),('FCV','closed'),
      ('GPV','normal'),('GPV','extrapolation'),('GPV','reverse'),('GPV','closed')]
    checks=0
    with tempfile.TemporaryDirectory(prefix='staci valves ') as directory:
      for unit,factor in units.items():
       for kind,mode in scenarios:
        if unit not in ('LPS','GPM','CFS') and mode not in ('active','normal'): continue
        work=Path(directory)/f'{unit}-{kind}-{mode}';work.mkdir()
        us=unit in {'CFS','GPM','MGD','IMGD','AFD'};length=.3048 if us else 1
        diameter=.1/.0254 if us else 100
        pressurefactor=.3048/.4333 if us else 1
        head=15 if mode=='open' else 60
        flow=.03 if mode=='extrapolation' else .005
        # Two source paths keep reverse/closed cases physically determined.
        two_sources=kind=='FCV' or mode in ('reverse','closed')
        r2=80 if mode=='reverse' else 0
        reservoirs=f'R1 {head/length:.15g}' + (f'\nR2 {r2/length:.15g}' if two_sources else '')
        demand=0 if two_sources else flow/factor
        pipes=f'P1 R1 J1 {100/length:.15g} {diameter:.15g} 120 0'
        if two_sources: pipes+=f'\nP2 J2 R2 {100/length:.15g} {diameter:.15g} 120 0'
        setting=20/pressurefactor if kind=='PRV' else (.5 if mode=='unattainable' else .005)/factor
        setting='C' if kind=='GPV' else f'{setting:.15g}'
        status='V OPEN' if mode=='forced_open' else 'V CLOSED' if mode=='closed' else ''
        curve='\n'.join(f'C {q/factor:.15g} {h/length:.15g}' for q,h in [(0,0),(.01,10),(.02,40)])
        text=f'[JUNCTIONS]\nJ1 0 0\nJ2 0 {demand:.15g}\n[RESERVOIRS]\n{reservoirs}\n[PIPES]\n{pipes}\n[VALVES]\nV J1 J2 {diameter:.15g} {kind} {setting} 0\n[CURVES]\n{curve}\n[STATUS]\n{status}\n[OPTIONS]\nUNITS {unit}\nHEADLOSS H-W\nTRIALS 100\n[END]\n'
        if unit=='LPS' and kind=='GPV' and mode=='normal':
            # A pump curve deliberately collides with the GPV export prefix.
            text=text.replace('[VALVES]','[PUMPS]\nPU R1 J1 HEAD STACI_GPV_V\n[VALVES]').replace('[CURVES]','[CURVES]\nSTACI_GPV_V 0 20\nSTACI_GPV_V 10 10\nSTACI_GPV_V 20 0').replace('[STATUS]','[STATUS]\nPU CLOSED')
        network=work/'network.inp';network.write_text(text);digest=hashlib.sha256(network.read_bytes()).hexdigest()
        run=subprocess.run([str(args.binary.resolve()),'-s',str(network),'--head-tolerance-m','1e-9'],cwd=work,capture_output=True,text=True,timeout=30)
        assert run.returncode==0,(unit,kind,mode,run.stdout,run.stderr)
        assert hashlib.sha256(network.read_bytes()).hexdigest()==digest
        actual=json.loads(Path(str(network)+'.hydraulics.json').read_text())
        reference=snapshot(library,network,work);result=compare(actual,reference)
        assert result['passed'],(unit,kind,mode,result['failures'])
        if kind=='PRV' and mode=='active': assert math.isclose(actual['nodes']['J2']['pressure_m'],20,abs_tol=1e-6)
        if kind=='FCV' and mode=='active': assert math.isclose(actual['links']['V']['flow_m3s'],.005,abs_tol=1e-8)
        if unit=='LPS' and mode in ('active','normal','forced_open','closed'):
            exported=work/'exported.inp'
            run=subprocess.run([str(args.binary.resolve()),'-y',str(network),'-o',str(exported)],cwd=work,capture_output=True,text=True,timeout=30)
            assert run.returncode==0,(kind,mode,'export',run.stdout,run.stderr)
            export_work=work/'export-reference';export_work.mkdir()
            exported_reference=snapshot(library,exported,export_work)
            result=compare(exported_reference,reference)
            assert result['passed'],(kind,mode,'export round trip',result['failures'])
        if unit=='LPS' and mode in ('active','normal'):
            modified=work/'modified.inp'
            run=subprocess.run([str(args.binary.resolve()),'-m',str(network),'-e','V','-p','tcv_minor_loss','-n','5','-o',str(modified)],cwd=work,capture_output=True,text=True,timeout=30)
            assert run.returncode==0,(kind,'modify',run.stdout,run.stderr)
            run=subprocess.run([str(args.binary.resolve()),'-s',str(modified),'--head-tolerance-m','1e-9'],cwd=work,capture_output=True,text=True,timeout=30)
            assert run.returncode==0,(kind,'modified solve',run.stdout,run.stderr)
            modified_work=work/'modified-reference';modified_work.mkdir()
            reference_modified=snapshot(library,modified,modified_work)
            actual_modified=json.loads(Path(str(modified)+'.hydraulics.json').read_text())
            result=compare(actual_modified,reference_modified)
            assert result['passed'],(kind,'modified solve',result['failures'])
        checks+=1
      for pressure,gravity,raw,target in [('METERS',2,20,10),('KPA',1,196.1334,20)]:
        work=Path(directory)/('pressure-'+pressure);work.mkdir()
        text=(Path(directory)/'LPS-PRV-active/network.inp').read_text().replace('UNITS LPS',f'UNITS LPS\nPRESSURE {pressure}\nSPECIFIC GRAVITY {gravity}').replace('100 PRV 20 0',f'100 PRV {raw} 0')
        network=work/'network.inp';network.write_text(text)
        run=subprocess.run([str(args.binary.resolve()),'-s',str(network),'--head-tolerance-m','1e-9'],cwd=work,capture_output=True,text=True,timeout=30)
        assert run.returncode==0,(pressure,run.stdout,run.stderr)
        actual=json.loads(Path(str(network)+'.hydraulics.json').read_text());reference=snapshot(library,network,work)
        result=compare(actual,reference);assert result['passed'],(pressure,result['failures'])
        assert math.isclose(actual['nodes']['J2']['pressure_m'],target,abs_tol=1e-4)
        checks+=1
      eps_checks=0
      base=(Path(directory)/'LPS-FCV-active/network.inp').read_text()
      for kind in ('PRV','FCV','GPV'):
        work=Path(directory)/('eps-'+kind);work.mkdir()
        lines=base.splitlines()
        for i,line in enumerate(lines):
            if line.startswith('V J1 J2 '): lines[i]='V J1 J2 100 '+kind+' '+({'PRV':'20','FCV':'5','GPV':'C'}[kind])+' 0'
        static='\n'.join(lines)+'\n'
        action='CLOSED' if kind=='GPV' else '15' if kind=='PRV' else '10'
        text=static.replace('[END]',f'[CONTROLS]\nLINK V {action} AT TIME 1\nLINK V OPEN AT TIME 2\n[TIMES]\nDURATION 2:00\nHYDRAULIC TIMESTEP 1:00\nREPORT TIMESTEP 1:00\n[END]')
        network=work/'network.inp';network.write_text(text);prefix=work/'eps'
        run=subprocess.run([str(args.binary.resolve()),'--epanet-eps',str(network),'-o',str(prefix),'--head-tolerance-m','1e-9'],cwd=work,capture_output=True,text=True,timeout=30)
        assert run.returncode==0,(kind,'EPS',run.stdout,run.stderr)
        nodes=list(csv.DictReader(Path(str(prefix)+'-nodes.csv').open()))
        links=list(csv.DictReader(Path(str(prefix)+'-links.csv').open()))
        for time in (0,3600,7200):
            reference_text=static
            if time==3600:
                if kind=='GPV': reference_text=reference_text.replace('[STATUS]','[STATUS]\nV CLOSED')
                else:
                    old='20' if kind=='PRV' else '5'
                    reference_text=reference_text.replace('100 '+kind+' '+old+' 0','100 '+kind+' '+action+' 0')
            if time==7200: reference_text=reference_text.replace('[STATUS]','[STATUS]\nV OPEN')
            ref_work=work/str(time);ref_work.mkdir();ref_input=ref_work/'reference.inp';ref_input.write_text(reference_text)
            reference=snapshot(library,ref_input,ref_work)
            actual={'nodes':{},'links':{}}
            for row in nodes:
                if int(row['time_seconds'])!=time: continue
                assert row['converged']=='1'
                actual['nodes'][row['node_id']]={'head_m':float(row['total_head_m']),'pressure_m':float(row['pressure_head_m']),'demand_m3s':float(row['demand_m3s'])}
            for row in links:
                if int(row['time_seconds'])!=time: continue
                assert row['converged']=='1'
                actual['links'][row['link_id']]={'flow_m3s':float(row['flow_m3s']),'velocity_mps':abs(float(row['velocity_mps'])),'enabled':row['status']=='1'}
            # EPS node demand is prescribed demand; the API returns net
            # boundary exchange for reservoirs. Reconstruct it from link flows.
            for name,node in reference['nodes'].items():
                if node['type']==0: continue
                actual['nodes'][name]['demand_m3s']=sum(
                    actual['links'][key]['flow_m3s']*((1 if link['to']==name else 0)-(1 if link['from']==name else 0))
                    for key,link in reference['links'].items())
            result=compare(actual,reference)
            assert result['passed'],(kind,'EPS',time,result['failures'])
            eps_checks+=1
    print(f'PASS {checks} valve/reference cases across all ten EPANET flow units, export round trips and {eps_checks} EPS controlled snapshots.')
    return 0
if __name__=='__main__': raise SystemExit(main())
