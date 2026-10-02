"""Independent EPANET 2.2 EPS reference: hydraulic events and quality reports."""
import ctypes as c
import math
import os
from epanet_reference import ReferenceError, validate_reference_balance
FLOW=(.028316846592,.0000630901964,.0438126363889,.0526167824074,.0142764101852,.001,1/60000,1000/86400,1/3600,1/86400)

def _simulate(library,network,work,damping=0.0,accuracy=1e-6):
    api=(c.WinDLL if os.name=='nt' else c.CDLL)(str(library))
    p=c.c_void_p();I=c.POINTER(c.c_int);D=c.POINTER(c.c_double);L=c.POINTER(c.c_long)
    signatures={'EN_getversion':[I],'EN_createproject':[c.POINTER(c.c_void_p)],'EN_deleteproject':[c.c_void_p],
      'EN_open':[c.c_void_p,c.c_char_p,c.c_char_p,c.c_char_p],'EN_close':[c.c_void_p],
      'EN_openH':[c.c_void_p],'EN_initH':[c.c_void_p,c.c_int],'EN_runH':[c.c_void_p,L],'EN_nextH':[c.c_void_p,L],'EN_closeH':[c.c_void_p],
      'EN_openQ':[c.c_void_p],'EN_initQ':[c.c_void_p,c.c_int],'EN_runQ':[c.c_void_p,L],'EN_nextQ':[c.c_void_p,L],'EN_closeQ':[c.c_void_p],
      'EN_gettimeparam':[c.c_void_p,c.c_int,L],'EN_getcount':[c.c_void_p,c.c_int,I],'EN_getflowunits':[c.c_void_p,I],
      'EN_getnodeid':[c.c_void_p,c.c_int,c.c_char_p],'EN_getlinkid':[c.c_void_p,c.c_int,c.c_char_p],
      'EN_getnodetype':[c.c_void_p,c.c_int,I],'EN_getlinktype':[c.c_void_p,c.c_int,I],'EN_getlinknodes':[c.c_void_p,c.c_int,I,I],
      'EN_getnodevalue':[c.c_void_p,c.c_int,c.c_int,D],'EN_getlinkvalue':[c.c_void_p,c.c_int,c.c_int,D],
      'EN_getstatistic':[c.c_void_p,c.c_int,D],'EN_setoption':[c.c_void_p,c.c_int,c.c_double],
      'EN_getqualinfo':[c.c_void_p,I,c.c_char_p,c.c_char_p,I],'EN_geterror':[c.c_int,c.c_char_p,c.c_int]}
    for name,args in signatures.items():getattr(api,name).argtypes=args;getattr(api,name).restype=c.c_int
    warnings=[]
    def call(name,*args):
        code=getattr(api,name)(*args)
        if code:
            msg=c.create_string_buffer(256);api.EN_geterror(code,msg,256);text=msg.value.decode(errors='replace')
            if code>=100 or code in (1,2):raise ReferenceError(code if code>=100 else 110,name+': '+text)
            warnings.append({'code':code,'message':text})
    def val(kind,index,prop):
        value=c.c_double();call('EN_get'+kind+'value',p,index,prop,c.byref(value));return value.value
    def count(prop):
        value=c.c_int();call('EN_getcount',p,prop,c.byref(value));return value.value
    def time(prop):
        value=c.c_long();call('EN_gettimeparam',p,prop,c.byref(value));return value.value
    version=c.c_int();call('EN_getversion',c.byref(version))
    if version.value//100!=202:raise ReferenceError(250,'Expected EPANET 2.2')
    call('EN_createproject',c.byref(p));opened=hopen=qopen=False
    try:
        call('EN_open',p,str(network).encode(),str(work/'epanet-eps.rpt').encode(),b'');opened=True
        unit=c.c_int();call('EN_getflowunits',p,c.byref(unit));flow=FLOW[unit.value];length=.3048 if unit.value<5 else 1
        call('EN_setoption',p,0,1000.0);call('EN_setoption',p,1,accuracy)
        if damping:
            call('EN_setoption',p,5,.0001/length);call('EN_setoption',p,6,1e-8/flow);call('EN_setoption',p,17,damping)
        duration,report_step,report_start=time(0),time(5),time(6)
        nmeta=[];lmeta=[]
        for i in range(1,count(0)+1):
            name=c.create_string_buffer(256);typ=c.c_int();call('EN_getnodeid',p,i,name);call('EN_getnodetype',p,i,c.byref(typ))
            nmeta.append((i,name.value.decode(),typ.value,val('node',i,0)*length))
        ids={i:name for i,name,typ,elev in nmeta};node_index={name:i for i,name,typ,elev in nmeta}
        for i in range(1,count(2)+1):
            name=c.create_string_buffer(256);typ=c.c_int();a=c.c_int();b=c.c_int();call('EN_getlinkid',p,i,name);call('EN_getlinktype',p,i,c.byref(typ));call('EN_getlinknodes',p,i,c.byref(a),c.byref(b))
            area=math.pi/4*(val('link',i,0)*(.0254 if unit.value<5 else .001))**2
            power=val('link',i,18) if typ.value==2 else 0
            lmeta.append((i,name.value.decode(),typ.value,ids[a.value],ids[b.value],area,power))
        frames={};events=0;max_energy=0
        call('EN_openH',p);hopen=True;call('EN_initH',p,1)
        while True:
            t=c.c_long();call('EN_runH',p,c.byref(t));events+=1
            error=c.c_double();call('EN_getstatistic',p,2,c.byref(error));head_error=error.value*length;max_energy=max(max_energy,head_error)
            if not math.isfinite(head_error) or head_error>.0001:raise ReferenceError(110,f'At {t.value}s head-equation error={head_error:.9g}m exceeds 0.0001m')
            # Check the physical POWER law at every event, including events
            # that fall between report times.
            for i,name,typ,a,b,area,power in lmeta:
                if not power or val('link',i,11)<=0:continue
                speed=val('link',i,12)
                if speed<=0:continue
                q=val('link',i,8)*flow;gain=(val('node',node_index[b],10)-val('node',node_index[a],10))*length
                coefficient=power*(8.814*.3048*.028316846592 if unit.value<5 else 8.814*.3048*.028316846592/.7457)*speed**3
                if q<=0 or abs(q*gain-coefficient)>max(1e-8,coefficient*1e-5):raise ReferenceError(110,f"At {t.value}s POWER pump '{name}' has no physical operating point")
            if t.value>=report_start and (t.value-report_start)%report_step==0:
                nodes={name:{'type':typ,'head_m':val('node',i,10)*length,'pressure_m':val('node',i,10)*length-elev,'demand_m3s':val('node',i,9)*flow} for i,name,typ,elev in nmeta}
                links={}
                for i,name,typ,a,b,area,power in lmeta:
                    speed=val('link',i,12) if typ==2 else None
                    raw=val('link',i,11)
                    links[name]={'type':typ,'from':a,'to':b,'area_m2':area,'flow_m3s':val('link',i,8)*flow,'velocity_mps':val('link',i,9)*length,'headloss_m':nodes[a]['head_m']-nodes[b]['head_m'],'enabled':raw>0 and (speed is None or speed>0),'raw_status':raw}
                validate_reference_balance(nodes,links)
                frames[t.value]={'nodes':nodes,'links':links,'maximum_head_error_m':head_error,'quality_nodes':{}}
            step=c.c_long();call('EN_nextH',p,c.byref(step))
            if step.value==0:break
        call('EN_closeH',p);hopen=False
        qtype=c.c_int();trace=c.c_int();chemical=c.create_string_buffer(256);qunits=c.create_string_buffer(256)
        call('EN_getqualinfo',p,c.byref(qtype),chemical,qunits,c.byref(trace))
        quality_status='not_configured'
        if qtype.value in (1,2):
            unit_text=qunits.value.decode().lower();factor=3600 if qtype.value==2 else .001 if unit_text in ('mg/l','mg/liter') else 1e-6 if unit_text=='ug/l' else None
            if factor is None:quality_status='unsupported_units:'+unit_text
            else:
                quality_status='available';call('EN_openQ',p);qopen=True;call('EN_initQ',p,0)
                while True:
                    t=c.c_long();call('EN_runQ',p,c.byref(t))
                    if t.value in frames:
                        frames[t.value]['quality_nodes']={name:val('node',i,12)*factor for i,name,typ,elev in nmeta}
                    step=c.c_long();call('EN_nextQ',p,c.byref(step))
                    if step.value==0:break
                # Independent sanity audit: the pinned reference can leave
                # SortedNodes uninitialized when every flow is zero, exposing
                # unchanged node quality despite reacting tank segments.
                if qtype.value==1 and len(frames)>1 and all(abs(link['flow_m3s'])<0.005*FLOW[1] for f in frames.values() for link in f['links'].values()):
                    first,last=min(frames),max(frames)
                    for i,name,typ,elev in nmeta:
                        if typ!=2:continue
                        kb=val('node',i,23)/86400
                        initial=frames[first]['quality_nodes'].get(name,0)
                        final=frames[last]['quality_nodes'].get(name,0)
                        if initial>0 and kb<0 and -kb*(last-first)>.03 and abs(final-initial)<initial*1e-10:
                            quality_status='numerical_failure_stagnant_tank'
                            warnings.append({'code':'QUALITY.STAGNANT_TANK','message':f"Tank {name}: reference concentration stays unchanged despite nonzero bulk decay and zero exchange; chemical equivalence cannot be established."})
                call('EN_closeQ',p);qopen=False
        elif qtype.value==3:quality_status='unsupported_trace'
        return {'frames':frames,'duration_s':duration,'hydraulic_events':events,'quality_type':qtype.value,'quality_status':quality_status,'warnings':warnings,'maximum_head_error_m':max_energy,'settings':{'accuracy':accuracy,'damping':damping,'trials':1000}}
    finally:
        if qopen:api.EN_closeQ(p)
        if hopen:api.EN_closeH(p)
        if opened:api.EN_close(p)
        api.EN_deleteproject(p)

def simulate(library,network,work,damping=0.0,accuracy=1e-6):
    # EPANET 2.2 creates its hydraulic scratch files in the current directory.
    # Serialize tokenizer/cwd use and keep all scratch files in the case folder.
    from epanet_reference import _REFERENCE_LOCK
    from pathlib import Path
    library=Path(library).resolve();network=Path(network).resolve();work=Path(work).resolve()
    with _REFERENCE_LOCK:
        previous=os.getcwd()
        try:
            os.chdir(work)
            return _simulate(library,network,work,damping,accuracy)
        finally:os.chdir(previous)
