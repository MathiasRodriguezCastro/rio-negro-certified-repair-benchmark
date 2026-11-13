#!/usr/bin/env python3
"""Audit this self-contained package, optionally replaying all small LP solves."""
import argparse,csv,json,math,platform,subprocess,sys,tempfile
from fractions import Fraction as R
from pathlib import Path
import numpy as np
import benchmark as b

ROOT=Path(__file__).resolve().parent


def production_expanded(coeff,x,y,v):
    a,c1,c2,d,e=coeff
    return a*x+c1*x*v+c2*x*v*v+d*x*x+d*x*y+e*x*x*x+2*e*x*x*y+e*x*y*y


def independent_bound(row):
    """Complete the square for the exported row, not its pre-rounded shift."""
    a,c1,c2,d,e=map(R,b.CONFIG['coefficients'][row['plant']])
    p,q,k=map(R,(row['px'],row['py'],row['intercept']))
    l,u,s,t,v0,v1=map(R,row['box'])
    assert l>=0 and s>=0 and e>=0 and -2*d-e*(6*u+4*t)>=0
    signs=[c1+2*c2*v for v in (v0,v1)]
    assert min(signs)>=0 or max(signs)<=0
    v=v0 if min(signs)>=0 else v1
    for x in (l,u):
        alpha=k+(p-a-c1*v-c2*v*v)*x-d*x*x-e*x*x*x
        beta=q-d*x-2*e*x*x;gamma=-e*x
        values=[alpha+beta*y+gamma*y*y for y in (s,t)]
        if gamma<0:
            center=-beta/(2*gamma)
            if s<=center<=t:values.append(alpha-beta*beta/(4*gamma))
        assert max(values)<=0, ('Exported affine row is not a lower bound',row)


def independent_point(instance,point,rows=None,require_service=True):
    assert point is not None
    T=b.CONFIG['periods'];dt=b.CONFIG['time_step'];tol=b.CONFIG['point_tolerance']
    x,y,v,g=(np.array(point[k],dtype=float) for k in ('x','y','v','g'));q=np.array(point['q'])
    assert all(z.shape==(3,T) and np.isfinite(z).all() for z in (x,y,v,g)) and q.shape==(T,)
    assert min(x.min(),y.min(),q.min())>=-tol
    assert x.max()<=b.CONFIG['xmax']+tol and y.max()<=b.CONFIG['ymax']+tol
    assert v.min()>=b.CONFIG['vmin']-tol and v.max()<=b.CONFIG['vmax']+tol
    assert np.max(q-instance['demand'])<=tol and np.max(q-g.sum(axis=0))<=tol
    for plant in range(3):
        for hour in range(T):
            before=instance['initial'][plant] if hour==0 else v[plant,hour-1]
            upstream=0.
            if plant>0 and hour>=plant:upstream=x[plant-1,hour-plant]+y[plant-1,hour-plant]
            residual=v[plant,hour]-before+dt*(x[plant,hour]+y[plant,hour]-instance['inflows'][plant][hour]-upstream)
            assert abs(residual)<=tol
    power=np.array([production_expanded(c,x[i],y[i],v[i]) for i,c in enumerate(b.CONFIG['coefficients'])])
    assert power.min()>=-tol
    if require_service:assert np.max(q-power.sum(axis=0))<=tol
    if rows:
        for row in rows:
            i,j=row['plant'],row['time'];vals=(x[i,j],y[i,j],v[i,j])
            for z,lo,hi in zip(vals,row['box'][::2],row['box'][1::2]):assert lo-tol<=z<=hi+tol
            assert g[i,j]-row['px']*x[i,j]-row['py']*y[i,j]-row['intercept']<=tol
    water=float(np.dot(b.CONFIG['water_values'],np.array(instance['initial'])-v[:,-1]))
    return dict(scheduled_cost=water+b.CONFIG['failure_cost']*float(np.sum(np.array(instance['demand'])-q)),
                final_cost=water+b.CONFIG['failure_cost']*float(np.maximum(0,np.array(instance['demand'])-power.sum(axis=0)).sum()),
                unsupported=float(np.maximum(0,q-power.sum(axis=0)).max()))


def check_milp(instance,initial):
    point=initial['point'];witness=initial['solver'];tol=b.CONFIG['point_tolerance']
    independent_point(instance,point,require_service=False)
    z=np.asarray(witness['head_selection']);split=np.asarray(witness['split_flow'])
    assert np.max(np.abs(z-np.rint(z)))<=1e-6 and np.max(np.abs(z.sum(2)-1))<=tol
    assert min(z.min(),split.min())>=-tol and np.max(split-z)<=tol
    assert np.max(np.abs(split.sum(2)-np.asarray(point['x'])))<=tol
    v=np.asarray(point['v']);assert np.min(v-2*z[:,:,0]-4*z[:,:,1])>=-tol
    assert np.max(v-4*z[:,:,0]-6*z[:,:,1])<=tol
    for i,(a,c1,c2,_,_) in enumerate(b.CONFIG['coefficients']):
        head=np.array([a+c1*value+c2*value*value for value in (3.,5.)])
        assert np.max(np.asarray(point['g'])[i]-split[i]@head)<=tol


def coordination_audit(point,baseline,variant):
    tol=b.CONFIG['point_tolerance'];T=b.CONFIG['periods'];q=np.asarray(point['q'])
    if variant==4:assert sum(point['v'][i][-1] for i in range(3))>=sum(baseline['v'][i][-1] for i in range(3))-.005-tol
    if variant==5:
        for start in range(0,T,8):assert sum(point['q'][start:start+8])>=.995*sum(baseline['q'][start:start+8])-tol
    if variant==6:
        original=np.array(baseline['x']).sum(0);limit=np.max(np.abs(np.diff(original)))+.01
        assert np.max(np.abs(np.diff(np.array(point['x']).sum(0))))<=limit+tol
        for t in range(T):
            if t%8<2:assert q[t]<=.98*baseline['q'][t]+tol
    multiplier=1000 if variant in (3,5) else 1
    prices=np.ones(T)
    if variant in (1,2,4,6):prices+=.2*np.sin(np.arange(T)*2*np.pi/T+(np.pi if variant==2 else 0))
    return -multiplier*float(np.dot(b.CONFIG['water_values'],np.array(point['v'])[:,-1]))-b.CONFIG['failure_cost']*float(np.dot(prices,q))


def close(a,c):return abs(a-c)<=max(1e-6,1e-9*max(abs(a),abs(c)))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--audit-only',action='store_true');args=parser.parse_args()
    protocol=json.loads((ROOT/'protocol.json').read_text())
    assert protocol['benchmark_sha256']==b.digest(ROOT/'benchmark.py') and protocol['configuration_sha256']==b.digest(b.CONFIG)
    subprocess.run([sys.executable,str(ROOT/'test_certificates.py'),'-v'],check=True)
    exact=points=coord=initials=unsupported=0
    for path in sorted((ROOT/'data/test').glob('*/scenario.json')):
        instance=json.loads(path.read_text());init=json.loads((path.parent/'initial.json').read_text());case=instance['id']
        regenerated=b.generate_instance('test',int(case[-3:]))
        for field in ('inflows','demand','initial'):assert np.allclose(instance[field],regenerated[field],rtol=0,atol=1e-13)
        check_milp(instance,init);initials+=1
        for method in ('taylor','derivative','interval','candidate'):
            directory=ROOT/'results/test'/case/method;rows=json.loads((directory/'certificates.json').read_text())
            computed,_=b.certificates(instance,init['point'],method);assert b.digest(rows)==b.digest(computed)
            if method!='taylor':
                for row in rows:independent_bound(row);exact+=1
            point=json.loads((directory/'point.json').read_text());record=json.loads((directory/'record.json').read_text())
            if point is not None:
                audit=independent_point(instance,point,rows,method!='taylor');points+=1
                assert close(audit['scheduled_cost'],record['point_audit']['scheduled_cost'])
                assert close(audit['final_cost'],record['point_audit']['final_cost'])
                unsupported+=method=='taylor' and audit['unsupported']>b.CONFIG['point_tolerance']
            if method=='candidate' and point is not None:
                for variant in range(1,7):
                    pair=[]
                    for label in ('cached','rebuilt'):
                        saved=json.loads((directory/f'coordination_{variant}_{label}.json').read_text());candidate=saved['point']
                        if candidate is not None:
                            independent_point(instance,candidate,rows);objective=coordination_audit(candidate,point,variant)
                            assert close(objective,saved['execution']['objective']);coord+=1
                        assert saved['execution']['production_rows_sha256']==b.digest(rows)
                        pair.append(saved)
                    assert pair[0]['execution']['success']==pair[1]['execution']['success']
                    if pair[0]['point'] is not None:assert close(pair[0]['execution']['objective'],pair[1]['execution']['objective'])
    solves=0
    if not args.audit_only:
        with tempfile.TemporaryDirectory(prefix='synthetic-replay-') as directory:
            fresh=Path(directory);b.run('test',fresh,reuse=True)
            old=list(csv.DictReader((ROOT/'results/test/metrics.csv').open()));new=list(csv.DictReader((fresh/'metrics.csv').open()))
            old={(r['instance'],r['method']):r for r in old};new={(r['instance'],r['method']):r for r in new};assert old.keys()==new.keys()
            for key,r in old.items():
                n=new[key];assert r['solver_success']==n['solver_success'] and r['admissible']==n['admissible']
                for column in ('scheduled_cost','final_cost'):assert close(float(r[column]),float(n[column]))
                reference=ROOT/'results/test'/key[0]/key[1]/'certificates.json';replica=fresh/key[0]/key[1]/'certificates.json'
                assert reference.read_bytes()==replica.read_bytes()
            solves=len(new)+len(list(csv.DictReader((fresh/'reuse.csv').open())))
    assert initials==12 and exact==2592 and points==48 and coord==144
    result=dict(status='passed',initial_MILP_witnesses_verified=initials,independent_full_box_rational_checks=exact,
                raw_quality_points_audited=points,raw_coordination_points_audited=coord,taylor_outputs_with_unsupported_service=unsupported,
                replayed_LP_solves=solves,regenerated_MILPs=0,proprietary_solver_required=False,
                proprietary_modules_loaded=any(n.startswith(('cplex','casadi','gurobipy')) for n in sys.modules),
                inherited_data_or_code_required=False,network_access_required=False,python=platform.python_version(),numpy=np.__version__,scipy=b.scipy.__version__)
    assert not result['proprietary_modules_loaded']
    b.write_json(ROOT/'REPRODUCTION_CHECKS.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
