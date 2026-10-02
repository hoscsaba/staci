"""Independent initial-time hydraulic reference using pinned EPANET 2.2 C API."""
import ctypes as c
import os
import math
import threading
_REFERENCE_LOCK = threading.Lock()
from pathlib import Path

def resolve_library(requested=None):
    roots = [Path(__file__).resolve().parent.parent / 'build' / 'epanet-reference' / 'build']
    candidates = [Path(requested)] if requested else []
    if os.environ.get('EPANET_LIBRARY'): candidates.append(Path(os.environ['EPANET_LIBRARY']))
    for root in roots:
        for directory in ('lib','bin','lib/Release','bin/Release'):
            for name in ('libepanet2.dylib','libepanet2.so','epanet2.dll'):
                candidates.append(root / directory / name)
    return next((p.resolve() for p in candidates if p.is_file()),None)

class ReferenceError(RuntimeError):
    def __init__(self, code, message):
        super().__init__(message); self.code=code

def validate_reference_balance(nodes, links):
    """Audit SI flows independently of EPANET's relative-flow stopping flag."""
    import math
    balance = {name: [] for name in nodes}
    throughput = {name: [] for name in nodes}
    for name, link in links.items():
        flow = link['flow_m3s']
        if not math.isfinite(flow):
            raise ReferenceError(110, "EPANET returned nonfinite flow for link '" + name + "'.")
        balance[link['from']].append(-flow)
        balance[link['to']].append(flow)
        throughput[link['from']].append(abs(flow))
        throughput[link['to']].append(abs(flow))
    failures = []
    for name, node in nodes.items():
        if not all(math.isfinite(node[field]) for field in ('head_m', 'pressure_m', 'demand_m3s')):
            raise ReferenceError(110, "EPANET returned nonfinite hydraulic values for node '" + name + "'.")
        if node['type'] != 0:
            continue
        residual = math.fsum(balance[name] + [-node['demand_m3s']])
        limit = max(1e-6, 1e-6 * max(abs(node['demand_m3s']), math.fsum(throughput[name])))
        if abs(residual) > limit:
            failures.append((abs(residual), name, residual, limit))
    if failures:
        _, name, residual, limit = max(failures)
        raise ReferenceError(110, "EPANET reports convergence but returned flows violate junction continuity: "
                             + "node '" + name + "', net inflow minus demand=" + format(residual, '.9g')
                             + " m3/s (limit " + format(limit, '.9g') + "). No reliable reference exists.")

def snapshot(library, network, work):
    # EPANET input parsing contains shared tokenizer state. Serialize C API calls.
    with _REFERENCE_LOCK:
        attempts = []
        for accuracy, damping in ((1e-6, 0.0), (1e-5, 0.0), (1e-6, 0.01), (1e-5, 0.01)):
            try:
                result = _snapshot(library, network, work, accuracy, damping)
                result['reference_settings'] = {'accuracy':accuracy,'trials':1000,'damping_limit':damping,'head_error_limit_m':0.0001,'flow_change_limit_m3s':1e-8 if damping else None}
                result['preceding_failed_attempts'] = attempts
                return result
            except ReferenceError as error:
                if error.code != 110: raise
                attempts.append({'accuracy':accuracy,'damping_limit':damping,'error':str(error)})
        raise ReferenceError(110, 'No reliable EPANET initial hydraulic reference: ' + str(attempts))

def _snapshot(library, network, work, accuracy, damping=0.0):
    loader = c.WinDLL if os.name == "nt" else c.CDLL
    api=loader(str(library))
    project=c.c_void_p()
    signatures={'EN_getversion':[c.POINTER(c.c_int)], 'EN_createproject':[c.POINTER(c.c_void_p)], 'EN_deleteproject':[c.c_void_p],
        'EN_open':[c.c_void_p,c.c_char_p,c.c_char_p,c.c_char_p], 'EN_close':[c.c_void_p],
        'EN_openH':[c.c_void_p], 'EN_initH':[c.c_void_p,c.c_int],
        'EN_runH':[c.c_void_p,c.POINTER(c.c_long)], 'EN_closeH':[c.c_void_p],
        'EN_getcount':[c.c_void_p,c.c_int,c.POINTER(c.c_int)],
        'EN_getflowunits':[c.c_void_p,c.POINTER(c.c_int)],
        'EN_getnodeid':[c.c_void_p,c.c_int,c.c_char_p], 'EN_getlinkid':[c.c_void_p,c.c_int,c.c_char_p],
        'EN_getlinktype':[c.c_void_p,c.c_int,c.POINTER(c.c_int)],
        'EN_getnodetype':[c.c_void_p,c.c_int,c.POINTER(c.c_int)],
        'EN_getnodevalue':[c.c_void_p,c.c_int,c.c_int,c.POINTER(c.c_double)],
        'EN_getlinkvalue':[c.c_void_p,c.c_int,c.c_int,c.POINTER(c.c_double)],
        'EN_getlinknodes':[c.c_void_p,c.c_int,c.POINTER(c.c_int),c.POINTER(c.c_int)],
        'EN_setoption':[c.c_void_p,c.c_int,c.c_double],
        'EN_getstatistic':[c.c_void_p,c.c_int,c.POINTER(c.c_double)],
        'EN_setlinkvalue':[c.c_void_p,c.c_int,c.c_int,c.c_double],
        'EN_geterror':[c.c_int,c.c_char_p,c.c_int]}
    for name,args in signatures.items():
        getattr(api,name).argtypes=args;getattr(api,name).restype=c.c_int
    warnings=[]
    def call(name,*args):
        code=getattr(api,name)(*args)
        if code:
            text=c.create_string_buffer(256);api.EN_geterror(code,text,256)
            if code>=100: raise ReferenceError(code,text.value.decode(errors='replace'))
            warnings.append({'code':code,'message':text.value.decode(errors='replace')})
    version=c.c_int();call('EN_getversion',c.byref(version))
    if version.value // 100 != 202:
        raise ReferenceError(250, 'Expected EPANET 2.2 reference, got version code ' + str(version.value))
    call('EN_createproject',c.byref(project))
    opened=hydraulics=False
    try:
        call('EN_open',project,str(network).encode(),str(work/'epanet.rpt').encode(),b'');opened=True
        # Tighten numerical accuracy, retaining physical parameters and controls.
        call('EN_setoption',project,1,accuracy)
        call('EN_setoption',project,0,1000.0)
        units=c.c_int();call('EN_getflowunits',project,c.byref(units))
        flow=[.028316846592,.0000630901964,.0438126363889,.0526167824074,.0142764101852,.001,1/60000,1000/86400,1/3600,1/86400][units.value]
        length=.3048 if units.value<5 else 1.0
        if damping:
            call('EN_setoption',project,5,0.0001/length)
            call('EN_setoption',project,6,1e-8/flow)
            call('EN_setoption',project,17,damping)
        call('EN_openH',project);hydraulics=True
        call('EN_initH',project,0);t=c.c_long();call('EN_runH',project,c.byref(t))
        disabled_pumps=[];count=c.c_int();call('EN_getcount',project,2,c.byref(count))
        for index in range(1,count.value+1):
            typ=c.c_int();call('EN_getlinktype',project,index,c.byref(typ))
            if typ.value!=2:continue
            speed=c.c_double();status=c.c_double()
            call('EN_getlinkvalue',project,index,12,c.byref(speed));call('EN_getlinkvalue',project,index,11,c.byref(status))
            if speed.value==0 and status.value>0:
                call('EN_setlinkvalue',project,index,11,0.0);disabled_pumps.append(index)
        if disabled_pumps:
            # Zero speed is a physical off state. Correct EPANET's inconsistent
            # current OPEN flag before auditing energy; do not alter initial
            # pump settings or any input bytes.
            warnings.clear();call('EN_runH',project,c.byref(t))
        units=c.c_int();call('EN_getflowunits',project,c.byref(units))
        flow=[.028316846592,.0000630901964,.0438126363889,.0526167824074,.0142764101852,.001,1/60000,1000/86400,1/3600,1/86400][units.value]
        length=.3048 if units.value<5 else 1.0
        def value(kind,index,prop):
            v=c.c_double();call('EN_get'+kind+'value',project,index,prop,c.byref(v));return v.value
        if any(w['code'] in (1,2) for w in warnings):
            raise ReferenceError(110, 'EPANET initial hydraulics did not converge reliably: ' + str(warnings))
        maximum=c.c_double();call('EN_getstatistic',project,2,c.byref(maximum))
        maximum_head_error_m=maximum.value*length
        if not math.isfinite(maximum_head_error_m) or maximum_head_error_m>0.0001:
            raise ReferenceError(110,'EPANET reports convergence but maximum head-equation error is '+format(maximum_head_error_m,'.9g')+' m (limit 0.0001 m).')
        nodes={}; indices={};links={}
        for kind,countprop in [('node',0),('link',2)]:
            count=c.c_int();call('EN_getcount',project,countprop,c.byref(count))
            for i in range(1,count.value+1):
                name=c.create_string_buffer(256);call('EN_get'+kind+'id',project,i,name);name=name.value.decode()
                if kind=='node':
                    typ=c.c_int();call('EN_getnodetype',project,i,c.byref(typ))
                    head=value(kind,i,10)*length;elev=value(kind,i,0)*length
                    nodes[name]={'head_m':head,'pressure_m':head-elev,'demand_m3s':value(kind,i,9)*flow,'type':typ.value};indices[i]=name
                else:
                    a=c.c_int();b=c.c_int();call('EN_getlinknodes',project,i,c.byref(a),c.byref(b))
                    typ=c.c_int();call('EN_getlinktype',project,i,c.byref(typ))
                    raw_status = value(kind,i,11)
                    pump_speed = value(kind,i,12) if typ.value == 2 else None
                    # SPEED 0 in the input can retain EPANET's raw OPEN/XHEAD
                    # flag despite zero hydraulic operation. Preserve the raw
                    # flag; normalized enabled means an operating link.
                    power=value(kind,i,18) if typ.value==2 else 0.0
                    links[name]={'constant_power_input':power,'type':typ.value,'velocity_mps':value(kind,i,9)*length,'flow_m3s':value(kind,i,8)*flow,
                        'raw_status':raw_status,'pump_speed':pump_speed,
                        'enabled':raw_status>0 and (pump_speed is None or pump_speed>0),
                        'area_m2':3.141592653589793/4*(value(kind,i,0)*(.0254 if units.value<5 else .001))**2,
                        'headloss_m':nodes[indices[a.value]]['head_m']-nodes[indices[b.value]]['head_m'],
                        'from':indices[a.value],'to':indices[b.value]}
        for name,link in links.items():
            if not link['enabled'] or link['constant_power_input']<=0: continue
            coefficient=link['constant_power_input']*(8.814*.3048*.028316846592 if units.value<5 else 8.814*.3048*.028316846592/.7457)*link['pump_speed']**3
            q=link['flow_m3s'];gain=-link['headloss_m']
            if q<=0 or abs(q*gain-coefficient)>max(1e-8,coefficient*1e-5):
                raise ReferenceError(110,"Constant-power pump '"+name+"' has no reliable physical operating point: Q="+format(q,'.9g')+' m3/s, head gain='+format(gain,'.9g')+' m; Q*head must equal '+format(coefficient,'.9g')+' m4/s.')
        validate_reference_balance(nodes, links)
        return {'api_version':version.value,'nodes':nodes,'links':links,'warnings':warnings,'time_seconds':t.value,'flow_units':units.value,'maximum_head_error_m':maximum_head_error_m,'zero_speed_status_normalizations':disabled_pumps}
    finally:
        if hydraulics: api.EN_closeH(project)
        if opened: api.EN_close(project)
        api.EN_deleteproject(project)

def compare(actual, reference):
    # Fixed engineering tolerances; do not adapt to individual discrepancies.
    tolerances={'head_m':(0.01,1e-6),'pressure_m':(0.01,1e-6),
                'demand_m3s':(1e-8,1e-4),'flow_m3s':(1e-6,1e-4),'headloss_m':(.02,1e-5),'velocity_mps':(1e-5,1e-4)}
    errors=[];maxima={};count=0
    def check(asset,field,a,b,closed=False,area=None,boundary=False):
        import math
        nonlocal count
        count+=1
        absolute,relative=tolerances[field]
        if field == 'velocity_mps' and area and area > 0: absolute = max(absolute, tolerances['flow_m3s'][0]/area)
        if boundary and field == 'demand_m3s': absolute = tolerances['flow_m3s'][0]
        if closed and field == 'flow_m3s': absolute = 1e-6  # EPANET closed-link penalty permits tiny leakage.
        error=abs(a-b);limit=max(absolute,relative*max(abs(a),abs(b)))
        if error>maxima.get(field,{}).get('absolute_error',-1):maxima[field]={'asset':asset,'absolute_error':error,'staci':a,'epanet':b}
        if not math.isfinite(a) or not math.isfinite(b) or error>limit:
            errors.append({'asset':asset,'field':field,'staci':a,'epanet':b,'absolute_error':error,'tolerance':limit})
    for kind in ['nodes','links']:
        for name,ref in reference[kind].items():
            asset=kind+'/'+name
            if name not in actual[kind]:errors.append({'asset':asset,'field':'missing'});continue
            act=actual[kind][name]
            if kind=='nodes':
                check(asset,'head_m',act['head_m'],ref['head_m'])
                check(asset,'demand_m3s',act['demand_m3s'],ref['demand_m3s'],boundary=ref['type']!=0)
                if ref['type']==0:
                    check(asset,'pressure_m',act['pressure_m'],ref['pressure_m'])
            else:
                check(asset,'flow_m3s',act['flow_m3s'],ref['flow_m3s'],closed=not ref['enabled'])
                if ref['type'] != 2:
                    check(asset,'velocity_mps',act['velocity_mps'],ref['velocity_mps'],area=ref['area_m2'])
                check(asset,'headloss_m',actual['nodes'][ref['from']]['head_m']-actual['nodes'][ref['to']]['head_m'],ref['headloss_m'])
                count+=1
                if act['enabled']!=ref['enabled']:errors.append({'asset':asset,'field':'enabled','staci':act['enabled'],'epanet':ref['enabled']})
    return {'passed':not errors,'comparisons':count,'failure_count':len(errors),'failures':errors[:100],
        'maxima':maxima,'tolerances':tolerances,'closed_link_flow_absolute_m3s':1e-6}
