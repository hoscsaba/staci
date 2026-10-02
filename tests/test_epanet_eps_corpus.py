"""Run original public INP periods and compare hydraulic/quality report frames."""
import argparse,csv,json,shutil,subprocess,time,math,sys,hashlib
from pathlib import Path
from epanet_reference import compare,ReferenceError
from epanet_eps_reference import simulate
parser=argparse.ArgumentParser();parser.add_argument('--binary',type=Path,required=True);parser.add_argument('--library',type=Path,required=True);parser.add_argument('--manifest',type=Path,default=Path(__file__).parent/'public_networks/manifest.json');parser.add_argument('--output',type=Path,required=True);parser.add_argument('--only');parser.add_argument('--timeout',type=float,default=180);parser.add_argument('--require-equivalence',action='store_true');args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
results=[]
def load_frames(prefix):
 frames={}
 for typ in ['nodes','links']:
  p=Path(str(prefix)+'-'+typ+'.csv')
  if not p.exists():return frames
  with p.open() as f:
   for row in csv.DictReader(f):
    t=int(row['time_seconds']);frame=frames.setdefault(t,{'nodes':{},'links':{},'converged':True,'quality_nodes':{}});frame['converged']=frame['converged'] and row['converged']=='1'
    if typ=='nodes':
     frame['nodes'][row['node_id']]={'head_m':float(row['total_head_m']),'pressure_m':float(row['pressure_head_m']),'demand_m3s':float(row['demand_m3s'])}
     frame['quality_nodes'][row['node_id']]={'age':float(row['water_age_s']),'chemical':float(row['chlorine_kgm3'])}
    else:frame['links'][row['link_id']]={'flow_m3s':float(row['flow_m3s']),'velocity_mps':abs(float(row['velocity_mps'])),'headloss_m':float(row['headloss_m']),'enabled':int(row['status'])>0}
 return frames
for entry in json.load(args.manifest.open())['networks']:
 if args.only and Path(entry['path']).stem not in args.only.split(','):continue
 source=args.manifest.parent/entry['path'];work=args.output/Path(entry['path']).with_suffix('');work.mkdir(parents=True,exist_ok=True);network=work/source.name;shutil.copy2(source,network);prefix=work/'staci';started=time.monotonic();r={'network':entry['path']};digest=hashlib.sha256(source.read_bytes()).hexdigest();assert digest==entry['sha256'],entry['path']
 try:
  with (work/'staci-output.log').open('w') as log:
   run=subprocess.run([str(args.binary),'-z',str(network),'-o',str(prefix),'--head-tolerance-m','1e-9','--max-iterations','1000'],cwd=work,stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout)
  r['staci_exit_code']=run.returncode
 except subprocess.TimeoutExpired:r['staci_exit_code']='timeout'
 actual=load_frames(prefix);reference=None;attempts=[]
 for accuracy,damping in [(1e-6,0),(1e-5,0),(1e-6,.01),(1e-5,.01)]:
  try:reference=simulate(args.library,network,work,damping,accuracy);break
  except ReferenceError as e:
   attempts.append({'accuracy':accuracy,'damping':damping,'code':e.code,'error':str(e)})
   if e.code!=110:break
  except Exception as e:attempts.append({'error':repr(e)});break
 r['reference_attempts']=attempts
 if reference is None:r['outcome']='reference_unavailable';r['actual_report_frames']=len(actual)
 elif r['staci_exit_code'] not in (0,3):r['outcome']='staci_failure';r.update(duration_s=reference['duration_s'],reference_report_frames=len(reference['frames']))
 else:
  r.update(duration_s=reference['duration_s'],actual_report_frames=len(actual),reference_report_frames=len(reference['frames']),reference_settings=reference['settings'],quality_type=reference['quality_type'],quality_reference_status=reference['quality_status'],reference_warnings=reference['warnings'])
  failures=[];nfail=0;quantities=0;frames=0;quality_count=0;quality_bad=0;quality_max=0;quality_examples=[];partial=[]
  for t,ref in reference['frames'].items():
   if t not in actual:failures.append({'time_s':t,'error':'missing STACI report frame'});nfail+=1;continue
   if not actual[t]['converged']:partial.append(t);continue
   comp=compare(actual[t],ref);frames+=1;nfail+=comp['failure_count'];quantities+=comp['comparisons']
   failures.extend(dict(f,time_s=t) for f in comp['failures'][:max(0,30-len(failures))])
   for name,expected in (ref['quality_nodes'].items() if reference['quality_status']=='available' else []):
    field='age' if reference['quality_type']==2 else 'chemical';value=actual[t]['quality_nodes'][name][field];quality_count+=1
    # AGE: 1 minute / 3%; chemical: 0.01 mg/L / 3%.
    error=abs(value-expected);limit=max(60 if field=='age' else 1e-5,abs(expected)*.03)
    if math.isfinite(error):quality_max=max(quality_max,error)
    if not math.isfinite(error) or error>limit:
     quality_bad+=1
     if len(quality_examples)<10:quality_examples.append({'time_s':t,'node':name,'field':field,'actual':value if math.isfinite(value) else None,'reference':expected,'tolerance':limit})
  r.update(compared_frames=frames,compared_quantities=quantities,hydraulic_failure_count=nfail,hydraulic_examples=failures[:30],partial_times_s=partial,quality_compared_values=quality_count,quality_failure_count=quality_bad,quality_max_absolute_error=quality_max,quality_examples=quality_examples)
  if partial or r['staci_exit_code']==3:r['outcome']='staci_partial'
  elif not frames:r['outcome']='no_report_comparison'
  elif nfail:r['outcome']='hydraulic_mismatch'
  else:r['outcome']='hydraulic_match'
  r['quality_outcome']='not_configured' if reference['quality_type']==0 else 'reference_unavailable' if reference['quality_status']!='available' else 'not_compared' if not quality_count else 'mismatch' if quality_bad else 'match'
 assert hashlib.sha256(network.read_bytes()).hexdigest()==digest,entry['path']
 r['seconds']=round(time.monotonic()-started,2);results.append(r);(work/'comparison-summary.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k in ['network','outcome','quality_outcome','seconds','compared_frames','hydraulic_failure_count']}),flush=True)
 report={'scope':'Original full configured periods, unchanged physical inputs; STACI 1e-9 m RMS/1000 iterations; EPANET independent reliability checks','results':results,'total':len(results),'outcomes':{k:sum(x['outcome']==k for x in results) for k in sorted(set(x['outcome'] for x in results))}}
 (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')

if args.require_equivalence:
 sys.exit(0 if results and all(r['outcome']=='hydraulic_match' and r.get('quality_outcome') in ('match','not_configured') for r in results) else 1)
