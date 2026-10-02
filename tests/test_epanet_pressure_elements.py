#!/usr/bin/env python3
import sys,json,subprocess,math,shutil,argparse,tempfile,hashlib,csv
from pathlib import Path
from epanet_reference import snapshot,compare,resolve_library
parser=argparse.ArgumentParser(description="PSV/PBV/emitter/PDA initial hydraulic comparisons against EPANET.")
parser.add_argument('--binary',type=Path,required=True)
parser.add_argument('--reference-library',type=Path)
args=parser.parse_args();lib=resolve_library(args.reference_library)
if not lib:
 print('SKIP: EPANET reference unavailable.');sys.exit(77)
temporary=tempfile.TemporaryDirectory(prefix='staci pressure elements ')
root=Path(temporary.name);results=[]

for unit,flowfactor,length,diam,pressure in [('LPS',.001,1,100,1),('GPM',.0000630901964,.3048,.1/.0254,.3048/.4333),('CFS',.028316846592,.3048,.1/.0254,.3048/.4333)]:
 for kind in ['PRV','PSV','PBV','EMITTER','PDA']:
  for mode in (['active','open','closed','reverse','minor','negative'] if kind in ['PRV','PSV'] else ['active','open','closed','reverse','minor'] if kind in ['PRV','PSV','PBV'] else ['negative','zero','partial','full']):
   w=root/(unit+'-'+kind+'-'+mode);w.mkdir(exist_ok=True)
   if kind in ['PRV','PSV','PBV']:
    head=60 if mode!='reverse' else -60
    setting=(-20 if mode=='negative' else 20)/pressure;pipe_length=1000 if kind=='PSV' and mode=='active' else 100
    text=f'[JUNCTIONS]\nJ1 0\nJ2 0\n[RESERVOIRS]\nR1 {head/length}\nR2 0\n[PIPES]\nP1 R1 J1 {pipe_length/length} {diam} 120 0\nP2 J2 R2 {100/length} {diam} 120 0\n[VALVES]\nV J1 J2 {diam} {kind} {setting} {100 if mode=="minor" else 0}\n[STATUS]\n'+('V CLOSED\n' if mode=='closed' else 'V OPEN\n' if mode=='open' else '')
   else:
    head={'negative':-10,'zero':0,'partial':5,'full':40}[mode]
    text=f'[JUNCTIONS]\nJ 0 {(.005 if kind=="PDA" else .001)/flowfactor}\n[RESERVOIRS]\nR {head/length}\n[PIPES]\nP R J {100/length} {diam} 120 0\n'
    if kind=='EMITTER':text+=f'[EMITTERS]\nJ {.0002/flowfactor*pressure**.5}\n'
   text+=f'[OPTIONS]\nUNITS {unit}\nHEADLOSS H-W\nTRIALS 1000\n'
   if kind=='PDA':text+=f'DEMAND MODEL PDA\nMINIMUM PRESSURE 0\nREQUIRED PRESSURE {20/pressure}\nPRESSURE EXPONENT 0.5\n'
   network=w/'network.inp';network.write_text(text+'[END]\n')
   digest=hashlib.sha256(network.read_bytes()).hexdigest()
   run=subprocess.run([str(args.binary.resolve()),'-s',str(network),'--head-tolerance-m','1e-9','--max-iterations','1000'],cwd=w,capture_output=True,text=True,timeout=30);(w/'staci.log').write_text(run.stdout+run.stderr)
   result={'case':w.name,'exit':run.returncode}
   if not run.returncode:
    try:
     ref=snapshot(lib,network,w);comp=compare(json.loads(Path(str(network)+'.hydraulics.json').read_text()),ref);result.update(passed=comp['passed'],failures=comp['failures'][:3]);(w/'comparison.json').write_text(json.dumps(comp,indent=2))
    except Exception as e:result['reference_error']=str(e)
   else:result['error']=(run.stdout+run.stderr)[-500:]
   assert hashlib.sha256(network.read_bytes()).hexdigest()==digest
   assert result.get('passed'),result
   actual=json.loads(Path(str(network)+'.hydraulics.json').read_text())
   if kind=='PSV' and mode=='active': assert math.isclose(actual['nodes']['J1']['pressure_m'],20,abs_tol=1e-6)
   if kind=='PBV' and mode=='active': assert math.isclose(actual['nodes']['J1']['head_m']-actual['nodes']['J2']['head_m'],20,abs_tol=1e-6)
   if kind=='PDA':
    p=actual['nodes']['J']['pressure_m'];delivered=.005*min(1,max(0,p/20))**.5
    assert math.isclose(actual['nodes']['J']['demand_m3s'],delivered,abs_tol=1e-9)
   exported=w/'exported.inp'
   export=subprocess.run([str(args.binary.resolve()),'-y',str(network),'-o',str(exported)],cwd=w,capture_output=True,text=True,timeout=30)
   assert export.returncode==0,(w.name,export.stdout,export.stderr)
   reference_work=w/'export-reference';reference_work.mkdir()
   roundtrip=compare(snapshot(lib,exported,reference_work),ref)
   assert roundtrip['passed'],(w.name,'export',roundtrip['failures'])
   results.append(result)
# Check EPS demand output after source pressure changes; compare each report
# to an independently solved EPANET steady model at that physical state.
for kind in ('EMITTER','PDA'):
 w=root/('eps-'+kind);w.mkdir()
 original=(root/('LPS-'+kind+'-partial')/'network.inp').read_text()
 network=w/'network.inp'
 network.write_text(original.replace('R 5.0','R 5.0 H').replace('[OPTIONS]','[PATTERNS]\nH 1 8 0.2\n[OPTIONS]').replace('[END]','[TIMES]\nDURATION 2:00\nHYDRAULIC TIMESTEP 1:00\nPATTERN TIMESTEP 1:00\nREPORT TIMESTEP 1:00\n[END]'))
 prefix=w/'eps'
 run=subprocess.run([str(args.binary.resolve()),'--epanet-eps',str(network),'-o',str(prefix),'--head-tolerance-m','1e-9'],cwd=w,capture_output=True,text=True,timeout=30)
 assert run.returncode==0,(kind,'EPS',run.stdout,run.stderr)
 rows=list(csv.DictReader(Path(str(prefix)+'-nodes.csv').open()))
 for time,head in ((0,5),(3600,40),(7200,1)):
  rw=w/str(time);rw.mkdir();static=rw/'network.inp';static.write_text(original.replace('R 5.0',f'R {head}'))
  reference=snapshot(lib,static,rw)
  row=next(x for x in rows if x['node_id']=='J' and int(x['time_seconds'])==time)
  assert math.isclose(float(row['demand_m3s']),reference['nodes']['J']['demand_m3s'],rel_tol=1e-4,abs_tol=1e-8),(kind,time,row,reference)
  assert math.isclose(float(row['pressure_head_m']),reference['nodes']['J']['pressure_m'],abs_tol=1e-5)
print(f'PASS {len(results)} PSV/PBV/emitter/PDA scenarios and {len(results)} INP export round trips.')
print('PASS 6 emitter/PDA EPS report frames against EPANET.')
temporary.cleanup()
