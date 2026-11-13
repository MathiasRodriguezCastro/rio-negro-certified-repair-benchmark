#!/usr/bin/env python3
"""Independent synthetic cascade: rational certificates and open-source LP/MILP.

Only newly generated normalized data are used. No Rio Negro coefficients,
historical scenarios, dispatches, native implementation, or proprietary solvers.
"""
from __future__ import annotations
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse, csv, hashlib, json, math, platform, time, warnings
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, linprog, milp
from scipy.sparse import coo_matrix, csr_matrix, hstack, vstack

ROOT = Path(__file__).resolve().parent
CONFIG = dict(seed=20261113, periods=24, plants=3, development=3, test=12,
              radius=.08, margin=.04, interval_subdivisions=4,
              time_step=.2, xmax=1., ymax=.8, vmin=2., vmax=6.,
              coefficients=[[1.25,.12,-.006,-.08,0.],
                            [1.9,-.08,.003,-.07,0.],
                            [1.45,.10,-.004,-.08,.006]],
              failure_cost=100., water_values=[1.2,1.6,2.],
              mip_time_limit_s=30., mip_relative_gap=1e-7,
              lp_tolerance=1e-9, point_tolerance=1e-7,
              coordination_variants=6)


def digest(value):
    if isinstance(value, Path): value=value.read_bytes()
    if not isinstance(value, bytes): value=json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
    return hashlib.sha256(value).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def write_csv(path, rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=keys);w.writeheader();w.writerows(rows)


def upward(value):
    f=float(value)
    return math.nextafter(f,math.inf) if Q(f)<value else f


def downward(value):
    f=float(value)
    return math.nextafter(f,-math.inf) if Q(f)>value else f


def polynomial(coeff,x,y,v):
    a,b,c,d,e=coeff
    return x*(a+b*v+c*v*v+d*(x+y)+e*(x+y)**2)


def fixed_tangent(coeff,x,y,v):
    a,b,c,d,e=coeff
    px=a+b*v+c*v*v+d*(2*x+y)+e*(3*x*x+4*x*y+y*y)
    py=x*(d+2*e*(x+y))
    return px,py,polynomial(coeff,x,y,v)-px*x-py*y


def control_hessian(coeff,x,y):
    a,b,c,d,e=coeff
    return 2*d+e*(6*x+4*y),d+e*(4*x+2*y),2*e*x


def geometry(coeff,box):
    a,b,c,d,e=map(Q,coeff);xl,xh,yl,yh,vl,vh=map(Q,box)
    assert 0<=xl<=xh and 0<=yl<=yh and vl<=vh and e>=0
    assert -2*d-e*(6*xh+4*yh)>=0, 'Control convexity is required'
    dv=[b+2*c*v for v in (vl,vh)]
    if min(dv)>=0:return vl
    if max(dv)<=0:return vh
    raise ValueError('State monotonicity is required')


def exact_candidates(coeff,tangent,box):
    coeff=list(map(Q,coeff));px,py,k=map(Q,tangent)
    xl,xh,yl,yh,_,_=map(Q,box);v=geometry(coeff,box);e=coeff[4];d=coeff[3]
    values=[]
    for x in {xl,xh}:
        ys={yl,yh}
        if e*x>0:
            vertex=(py-d*x)/(2*e*x)-x
            ys.add(min(yh,max(yl,vertex)))
        for y in ys:values.append(px*x+py*y+k-polynomial(coeff,x,y,v))
    return max(Q(0),max(values))


def square_oracle(coeff,tangent,box):
    """Independent completed-square check, without candidate point evaluation."""
    a,b,c,d,e=map(Q,coeff);px,py,k=map(Q,tangent)
    xl,xh,yl,yh,_,_=map(Q,box);v=geometry(coeff,box);values=[]
    for x in (xl,xh):
        A=k+x*(px-a-b*v-c*v*v)-d*x*x-e*x**3
        B=py-d*x-2*e*x*x;C=-e*x
        values.extend(A+B*y+C*y*y for y in (yl,yh))
        if C<0 and yl<=-B/(2*C)<=yh:values.append(A-B*B/(4*C))
    return max(values)


def derivative_shift(coeff,tangent,anchor,box):
    a,b,c,d,e=map(Q,coeff);px,py,k=map(Q,tangent)
    x0,y0,v0=map(Q,anchor);xl,xh,yl,yh,vl,vh=map(Q,box)
    ideal=fixed_tangent((a,b,c,d,e),x0,y0,v0)
    rounding=sum(abs(u-v)*s for u,v,s in zip((px,py,k),ideal,(xh,yh,Q(1))))
    dx=max(x0-xl,xh-x0);dy=max(y0-yl,yh-y0)
    hx=(min(xl,x0),max(xh,x0));hy=(min(yl,y0),max(yh,y0));hv=(min(vl,v0),max(vh,v0))
    hxx=max(abs(control_hessian((a,b,c,d,e),x,y)[0]) for x in hx for y in hy)
    hxy=max(abs(control_hessian((a,b,c,d,e),x,y)[1]) for x in hx for y in hy)
    hyy=2*e*max(abs(x) for x in hx);lv=xh*max(abs(b+2*c*v) for v in hv)
    return hxx*dx*dx/2+hxy*dx*dy+hyy*dy*dy/2+lv*max(abs(v0-vl),abs(vh-v0))+rounding


def interval_shift(coeff,tangent,box,subdivisions):
    """Centered polynomial intervals, using exact rational arithmetic."""
    a,b,c,d,e=map(Q,coeff);px,py,k=map(Q,tangent)
    xl,xh,yl,yh,_,_=map(Q,box);v=geometry(coeff,box)
    terms={(0,0):k,(1,0):px-a-b*v-c*v*v,(0,1):py,
           (2,0):-d,(1,1):-d,(3,0):-e,(2,1):-2*e,(1,2):-e}
    largest=Q(0)
    for ix in range(subdivisions):
        lx=xl+(xh-xl)*ix/subdivisions;ux=xl+(xh-xl)*(ix+1)/subdivisions
        for iy in range(subdivisions):
            ly=yl+(yh-yl)*iy/subdivisions;uy=yl+(yh-yl)*(iy+1)/subdivisions
            cx,cy=(lx+ux)/2,(ly+uy)/2;rx,ry=(ux-lx)/2,(uy-ly)/2
            centered={}
            for (i,j),coef in terms.items():
                for p in range(i+1):
                    for q in range(j+1):
                        key=(p,q)
                        centered[key]=centered.get(key,Q(0))+coef*math.comb(i,p)*math.comb(j,q)*cx**(i-p)*cy**(j-q)
            upper=Q(0)
            for (p,q),coef in centered.items():
                radius=rx**p*ry**q
                lo=-radius if p%2 or q%2 else (Q(1) if p==q==0 else Q(0))
                hi=radius
                upper+=coef*(hi if coef>=0 else lo)
            largest=max(largest,upper)
    return largest


def certificates(instance,anchor,method):
    start=time.perf_counter();x,y,v=(np.asarray(anchor[k]) for k in ('x','y','v'))
    rows=[]
    for i,coeff in enumerate(CONFIG['coefficients']):
        for t in range(CONFIG['periods']):
            ref=(x[i,t],y[i,t],v[i,t]);tan=fixed_tangent(coeff,*ref)
            r,m=CONFIG['radius'],CONFIG['margin']
            box=(max(0.,ref[0]-r*CONFIG['xmax']),min(CONFIG['xmax'],ref[0]+r*CONFIG['xmax']),
                 max(0.,ref[1]-r*CONFIG['ymax']),min(CONFIG['ymax'],ref[1]+r*CONFIG['ymax']),
                 max(CONFIG['vmin'],(1-m)*ref[2]),min(CONFIG['vmax'],(1+m)*ref[2]))
            if method=='candidate':shift=exact_candidates(coeff,tan,box)
            elif method=='derivative':shift=derivative_shift(coeff,tan,ref,box)
            elif method=='interval':shift=interval_shift(coeff,tan,box,CONFIG['interval_subdivisions'])
            elif method=='taylor':shift=Q(0)
            else:raise ValueError(method)
            delta=upward(shift);intercept=downward(Q(tan[2])-Q(delta))
            rows.append(dict(plant=i,time=t,box=list(box),reference=list(ref),px=tan[0],py=tan[1],
                             tangent_intercept=tan[2],shift=delta,intercept=intercept,
                             rational_shift=str(shift)))
    return rows,time.perf_counter()-start


class Cascade:
    def __init__(self,instance):
        self.instance=instance;self.T=CONFIG['periods'];self.n=13*self.T
        self.rows=[];self.lower=[];self.upper=[]
        self.lb=np.full(self.n,-np.inf);self.ub=np.full(self.n,np.inf)
        self.c=np.zeros(self.n)
        for i in range(3):
            for k,lo,hi in [('x',0,CONFIG['xmax']),('y',0,CONFIG['ymax']),('v',CONFIG['vmin'],CONFIG['vmax'])]:
                self.lb[self.block(k,i)]=lo;self.ub[self.block(k,i)]=hi
            self.c[self.index('v',i,self.T-1)]=-CONFIG['water_values'][i]
        self.lb[self.block('q')]=0;self.ub[self.block('q')]=instance['demand']
        self.c[self.block('q')]=-CONFIG['failure_cost']
        dt=CONFIG['time_step'];inflows=np.asarray(instance['inflows'])
        for t in range(self.T):
            for i in range(3):
                row={self.index('v',i,t):1.,self.index('x',i,t):dt,self.index('y',i,t):dt}
                rhs=dt*inflows[i,t]
                if t:row[self.index('v',i,t-1)]=-1.
                else:rhs+=instance['initial'][i]
                delay=i
                if i and t>=delay:
                    row[self.index('x',i-1,t-delay)]=-dt;row[self.index('y',i-1,t-delay)]=-dt
                self.add(row,rhs,rhs)
            self.add({self.index('q',0,t):1.,**{self.index('g',i,t):-1. for i in range(3)}},-np.inf,0.)

    def index(self,k,i=0,t=0):
        return ({'x':0,'y':3,'v':6,'g':9,'q':12}[k]+(i if k!='q' else 0))*self.T+t

    def block(self,k,i=0):return slice(self.index(k,i),self.index(k,i)+self.T)

    def add(self,row,lo,hi):self.rows.append(row);self.lower.append(lo);self.upper.append(hi)

    def sparse(self,n=None):
        r=[];c=[];d=[]
        for j,row in enumerate(self.rows):
            for k,val in row.items():r.append(j);c.append(k);d.append(val)
        return coo_matrix((d,(r,c)),shape=(len(self.rows),n or self.n)).tocsr()

    def decode(self,values):
        return {k:np.array([values[self.block(k,i)] for i in range(3)]).tolist() for k in ('x','y','v','g')}|{'q':values[self.block('q')].tolist()}


def generate_instance(split,index):
    seed=CONFIG['seed']+(100000 if split=='test' else 0)+index
    rng=np.random.default_rng(seed);t=np.arange(CONFIG['periods'])
    inflows=np.array([.25+.07*np.sin(2*np.pi*t/24+rng.uniform(0,6))+.02*rng.uniform(-1,1,len(t)) for _ in range(3)])
    demand=2.6+.5*np.sin(2*np.pi*t/24+rng.uniform(0,6))+.1*rng.uniform(-1,1,len(t))
    return dict(id=f'synthetic-{split}-{index:03d}',split=split,seed=seed,inflows=inflows.tolist(),demand=demand.tolist(),initial=rng.uniform(3.8,4.2,3).tolist())


def solve_initial_milp(instance):
    model=Cascade(instance);T=model.T;offset=model.n;nb=3*T*2;N=model.n+2*nb
    lb=np.r_[model.lb,np.zeros(2*nb)];ub=np.r_[model.ub,np.ones(2*nb)]
    for i in range(3):lb[model.block('g',i)]=0
    integrality=np.zeros(N,dtype=int);integrality[offset+nb:]=1
    for i,coeff in enumerate(CONFIG['coefficients']):
        for t in range(T):
            x=model.index('x',i,t);v=model.index('v',i,t);g=model.index('g',i,t)
            qs=[offset+(i*T+t)*2+j for j in range(2)];zs=[q+nb for q in qs]
            model.add({x:1.,qs[0]:-1.,qs[1]:-1.},0,0)
            model.add({z:1. for z in zs},1,1)
            model.add({v:1.,zs[0]:-2.,zs[1]:-4.},0,np.inf)
            model.add({v:1.,zs[0]:-4.,zs[1]:-6.},-np.inf,0)
            for q,z in zip(qs,zs):model.add({q:1.,z:-1.},-np.inf,0.)
            head=[coeff[0]+coeff[1]*u+coeff[2]*u*u for u in (3.,5.)]
            model.add({g:1.,qs[0]:-head[0],qs[1]:-head[1]},-np.inf,0)
    start=time.perf_counter()
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Unrecognized options detected:.*')
        result=milp(np.r_[model.c,np.zeros(2*nb)],integrality=integrality,bounds=Bounds(lb,ub),
                    constraints=LinearConstraint(model.sparse(N),model.lower,model.upper),
                    options={'time_limit':CONFIG['mip_time_limit_s'],'mip_rel_gap':CONFIG['mip_relative_gap'],'threads':1})
    record=dict(status=int(result.status),success=bool(result.success),message=result.message,wall_s=time.perf_counter()-start,
                mip_gap=float(result.mip_gap) if result.x is not None else None)
    if result.x is None:raise RuntimeError('Synthetic initial MILP has no incumbent: '+str(record))
    record['head_selection']=result.x[offset+nb:].reshape((3,T,2)).tolist()
    record['split_flow']=result.x[offset:offset+nb].reshape((3,T,2)).tolist()
    return model.decode(result.x),record


def coordination(model,baseline,variant):
    c=model.c.copy();T=model.T;phase=np.arange(T)*2*np.pi/T
    extra=[];rhs=[]
    if variant in (1,2,4,6):c[model.block('q')]*=1+.2*np.sin(phase+(np.pi if variant==2 else 0))
    if variant in (3,5):
        for i in range(3):c[model.index('v',i,T-1)]*=1000
    if variant==4:
        extra.append({model.index('v',i,T-1):-1. for i in range(3)})
        rhs.append(-sum(baseline['v'][i][-1] for i in range(3))+.005)
    if variant==5:
        for start in range(0,T,8):
            extra.append({model.index('q',0,t):-1. for t in range(start,min(T,start+8))})
            rhs.append(-.995*sum(baseline['q'][start:start+8]))
    if variant==6:
        flow=np.sum(np.asarray(baseline['x']),axis=0);limit=float(np.max(np.abs(np.diff(flow))))+.01
        for t in range(1,T):
            for sign in (-1.,1.):
                extra.append({model.index('x',i,t):sign for i in range(3)}|{model.index('x',i,t-1):-sign for i in range(3)})
                rhs.append(limit)
        for t in range(T):
            if t%8<2:
                extra.append({model.index('q',0,t):1.});rhs.append(.98*baseline['q'][t])
    return c,extra,rhs


def solve_lp(instance,rows,baseline=None,variant=0):
    started=time.perf_counter();model=Cascade(instance)
    for r in rows:
        i,t=r['plant'],r['time'];xl,xh,yl,yh,vl,vh=r['box']
        for k,lo,hi in [('x',xl,xh),('y',yl,yh),('v',vl,vh)]:model.lb[model.index(k,i,t)]=lo;model.ub[model.index(k,i,t)]=hi
        model.add({model.index('g',i,t):1.,model.index('x',i,t):-r['px'],model.index('y',i,t):-r['py']},-np.inf,r['intercept'])
    c,extra,rhs=coordination(model,baseline,variant) if variant else (model.c,[],[])
    for row,b in zip(extra,rhs):model.add(row,-np.inf,b)
    A=model.sparse();lo=np.array(model.lower);hi=np.array(model.upper);eq=lo==hi
    assembly=time.perf_counter()-started;begin=time.perf_counter()
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Unrecognized options detected:.*')
        result=linprog(c,A_ub=A[~eq],b_ub=hi[~eq],A_eq=A[eq],b_eq=hi[eq],bounds=np.c_[model.lb,model.ub],method='highs-ds',
                       options={'primal_feasibility_tolerance':CONFIG['lp_tolerance'],'dual_feasibility_tolerance':CONFIG['lp_tolerance'],'threads':1})
    record=dict(success=bool(result.success),status=int(result.status),message=result.message,assembly_s=assembly,
                solve_s=time.perf_counter()-begin,procedure_s=time.perf_counter()-started,
                objective=float(result.fun) if result.x is not None else None,extra_rows=len(extra),production_rows_sha256=digest(rows))
    if result.x is not None:
        activity=A@result.x
        record['linear_residual']=max(0.,float(np.max(lo-activity)),float(np.max(activity-hi)),float(np.max(model.lb-result.x)),float(np.max(result.x-model.ub)))
    return model.decode(result.x) if result.x is not None else None,record


def point_audit(instance,point,rows=None):
    if point is None:return dict(passed=False,reason='no numerical point')
    x,y,v,g=(np.asarray(point[k]) for k in ('x','y','v','g'));q=np.asarray(point['q']);T=CONFIG['periods']
    powers=np.array([polynomial(p,x[i],y[i],v[i]) for i,p in enumerate(CONFIG['coefficients'])]);available=powers.sum(0)
    mass=[]
    for i in range(3):
        previous=np.r_[instance['initial'][i],v[i,:-1]];up=np.zeros(T)
        if i:up[i:]=x[i-1,:T-i]+y[i-1,:T-i]
        mass.extend(v[i]-previous-CONFIG['time_step']*(np.asarray(instance['inflows'][i])+up-x[i]-y[i]))
    technical=max(0.,float(np.max(-x)),float(np.max(x-CONFIG['xmax'])),float(np.max(-y)),float(np.max(y-CONFIG['ymax'])),float(np.max(CONFIG['vmin']-v)),float(np.max(v-CONFIG['vmax'])),float(np.max(-q)),float(np.max(q-np.array(instance['demand']))))
    affine=domain=0.
    if rows:
        for r in rows:
            i,t=r['plant'],r['time'];affine=max(affine,g[i,t]-r['px']*x[i,t]-r['py']*y[i,t]-r['intercept'])
            for value,lo,hi in zip((x[i,t],y[i,t],v[i,t]),r['box'][::2],r['box'][1::2]):domain=max(domain,lo-value,value-hi)
    unsupported=max(0.,float(np.max(q-available)));balance=max(abs(z) for z in mass)
    credit=max(0.,float(np.max(q-g.sum(0))));tol=CONFIG['point_tolerance']
    water=float(np.dot(CONFIG['water_values'],np.array(instance['initial'])-v[:,-1]))
    scheduled=water+CONFIG['failure_cost']*float(np.sum(np.array(instance['demand'])-q))
    final=water+CONFIG['failure_cost']*float(np.sum(np.maximum(0,np.array(instance['demand'])-available)))
    hidden=float(np.maximum(0,np.minimum(g.sum(0),instance['demand'])-np.minimum(available,instance['demand'])).sum())
    passed=max(balance,technical,affine,domain,credit,unsupported)<=tol
    return dict(passed=bool(passed),mass_residual=balance,technical_residual=technical,affine_residual=max(0.,affine),domain_residual=max(0.,domain),credit_residual=credit,unsupported_service=unsupported,hidden_energy=hidden,scheduled_cost=scheduled,final_cost=final,final_credit_energy=float(np.sum(np.minimum(instance['demand'],available)-q)))


def audit_certificates(rows):
    worst=Q(0)
    for r in rows:
        coeff=CONFIG['coefficients'][r['plant']]
        excess=square_oracle(coeff,(r['px'],r['py'],r['intercept']),r['box'])
        assert excess<=0, ('Invalid mathematical production bound',r['plant'],r['time'],excess)
        worst=max(worst,excess)
    return dict(passed=True,exact_completed_square_checks=len(rows),maximum_positive_excess=str(worst))


def run(split,output,reuse=False):
    if split=='test':
        protocol=json.loads((ROOT/'protocol.json').read_text())
        assert protocol['benchmark_sha256']==digest(Path(__file__)) and protocol['configuration_sha256']==digest(CONFIG)
    output.mkdir(parents=True,exist_ok=True);metrics=[];reuse_records=[]
    paths=sorted((ROOT/'data'/split).glob('*/scenario.json'))
    assert len(paths)==CONFIG[split]
    for position,path in enumerate(paths):
        instance=json.loads(path.read_text());anchor=json.loads((path.parent/'initial.json').read_text())['point'];case=instance['id']
        methods=['taylor','derivative','interval','candidate'];methods=methods[position%4:]+methods[:position%4]
        for method in methods:
            folder=output/case/method;folder.mkdir(parents=True,exist_ok=True)
            rows,build=certificates(instance,anchor,method)
            math_audit=audit_certificates(rows) if method!='taylor' else dict(passed=None)
            point,record=solve_lp(instance,rows);check=point_audit(instance,point,rows)
            write_json(folder/'certificates.json',rows);write_json(folder/'point.json',point)
            write_json(folder/'record.json',dict(construction_s=build,solver=record,point_audit=check,mathematical_audit=math_audit))
            metrics.append(dict(instance=case,method=method,solver_success=record['success'],admissible=check['passed'],construction_s=build,assembly_s=record['assembly_s'],solve_s=record['solve_s'],**{k:check.get(k) for k in ('hidden_energy','unsupported_service','scheduled_cost','final_cost','final_credit_energy')}))
            if reuse and method=='candidate' and record['success'] and check['passed']:
                frozen=digest(rows)
                for variant in range(1,7):
                    pair={};cache_first=(int(case[-3:])+variant)%2==0
                    for label in (('cached','rebuilt') if cache_first else ('rebuilt','cached')):
                        if label=='rebuilt':current,seconds=certificates(instance,anchor,'candidate')
                        else:current,seconds=rows,0.
                        assert digest(current)==frozen
                        new,execution=solve_lp(instance,current,point,variant);new_audit=point_audit(instance,new,current)
                        if new is not None:
                            new_audit['coordination_linear_residual']=execution['linear_residual']
                            new_audit['passed']=new_audit['passed'] and execution['linear_residual']<=CONFIG['point_tolerance']
                        assert digest(rows)==frozen
                        item=dict(instance=case,variant=variant,path=label,order='cached-first' if cache_first else 'rebuilt-first',certificate_s=seconds,**execution,admissible=new_audit['passed'],changed_from_base=max(float(np.max(np.abs(np.asarray(new[k])-np.asarray(point[k])))) for k in ('x','y','v','q'))>1e-7 if new else False)
                        reuse_records.append(item);pair[label]=item
                        write_json(folder/f'coordination_{variant}_{label}.json',dict(point=new,execution=item,point_audit=new_audit))
                    assert pair['cached']['success']==pair['rebuilt']['success']
                    if pair['cached']['success']:assert abs(pair['cached']['objective']-pair['rebuilt']['objective'])<=1e-6
        print('Audited synthetic case',case,flush=True)
    write_csv(output/'metrics.csv',metrics)
    if reuse:write_csv(output/'reuse.csv',reuse_records)
    summary=dict(status='passed',split=split,attempts=len(metrics),methods={m:dict(attempted=sum(r['method']==m for r in metrics),converged=sum(r['method']==m and r['solver_success'] for r in metrics),admissible=sum(r['method']==m and r['admissible'] for r in metrics)) for m in ('taylor','derivative','interval','candidate')},reuse_attempts=len(reuse_records),reuse_changed=sum(r['path']=='cached' and r['changed_from_base'] for r in reuse_records),configuration_sha256=digest(CONFIG),python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    write_json(output/'summary.json',summary);print(json.dumps(summary,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('generate','run'))
    p.add_argument('--split',choices=('development','test'),required=True);p.add_argument('--output',type=Path);p.add_argument('--reuse',action='store_true');args=p.parse_args()
    if args.action=='generate':
        for i in range(CONFIG[args.split]):
            instance=generate_instance(args.split,i);folder=ROOT/'data'/args.split/instance['id']
            if (folder/'initial.json').exists():continue
            point,record=solve_initial_milp(instance)
            write_json(folder/'scenario.json',instance);write_json(folder/'initial.json',dict(point=point,solver=record))
            print('Saved',instance['id'],record,flush=True)
        write_json(ROOT/'configuration.json',CONFIG)
    else:run(args.split,args.output or ROOT/'results'/args.split,args.reuse)


if __name__=='__main__':main()
