#!/usr/bin/env python3
"""v1.1.0 extension: full Taylor, including storage, on the original fixed boxes.

New independent rational implementation; the v1.0.1 scientific files are retained.
No real Rio Negro data, coefficients, native classes or licensed solvers are used.
"""
import argparse, itertools, json, math, time, warnings
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
import benchmark as b

CONFIG = dict(version='1.1.0', methods=['derivative','interval','candidate'],
              subdivisions=[4,4,2], radius=b.CONFIG['radius'], margin=b.CONFIG['margin'],
              boxes='exactly those of the supplied fixed-storage candidate rows',
              statistical_scope='retrospective reproducibility examples, not new independent test data')
ROOT=Path(__file__).resolve().parent


def full_tangent(coeff,x,y,v):
    px,py,k=b.fixed_tangent(coeff,x,y,v)
    pv=x*(coeff[1]+2*coeff[2]*v)
    return px,py,pv,k-pv*v


def candidates(coeff,affine,box):
    a,bb,c,d,e=map(Q,coeff);px,py,pv,k=map(Q,affine)
    xl,xh,yl,yh,vl,vh=map(Q,box)
    assert xl>=0 and e>=0 and -2*d-e*(6*xh+4*yh)>=0
    values=[];interiors=0
    for x in (xl,xh):
        ys=[yl,yh];vs=[vl,vh]
        if e*x>0:ys.append(min(yh,max(yl,(py-d*x)/(2*e*x)-x)))
        if c*x>0:
            z=min(vh,max(vl,(pv-bb*x)/(2*c*x)));vs.append(z)
            interiors+=vl<z<vh
        values.extend(px*x+py*y+pv*v+k-b.polynomial((a,bb,c,d,e),x,y,v) for y in ys for v in vs)
    return max(Q(0),max(values)),len(values),interiors


def oracle(coeff,affine,box):
    """Independent separated completed squares; return an uncapped maximum."""
    a,bb,c,d,e=map(Q,coeff);px,py,pv,k=map(Q,affine)
    xl,xh,yl,yh,vl,vh=map(Q,box)
    assert xl>=0 and e>=0 and min(-2*d-e*(6*x+4*y) for x in (xl,xh) for y in (yl,yh))>=0
    def quadratic(A,B,lo,hi):
        choices=[(A*z*z+B*z,z) for z in (lo,hi)]
        if A<0 and lo < -B/(2*A) < hi:choices.append((-B*B/(4*A),-B/(2*A)))
        return max(choices,key=lambda r:r[0])
    out=[]
    for x in (xl,xh):
        yy,y=quadratic(-e*x,py-d*x-2*e*x*x,yl,yh)
        vv,v=quadratic(-c*x,pv-bb*x,vl,vh)
        out.append((k+px*x-a*x-d*x*x-e*x*x*x+yy+vv,x,y,v))
    return max(out,key=lambda r:r[0])


def derivative(coeff,affine,anchor,box):
    a,bb,c,d,e=map(Q,coeff);px,py,pv,k=map(Q,affine);x,y,v=map(Q,anchor)
    xl,xh,yl,yh,vl,vh=map(Q,box)
    rx=max(abs(xl-x),abs(xh-x));ry=max(abs(yl-y),abs(yh-y));rv=max(abs(vl-v),abs(vh-v))
    hx=(min(xl,x),max(xh,x));hy=(min(yl,y),max(yh,y));hv=(min(vl,v),max(vh,v))
    hxx=max(abs(2*d+e*(6*xx+4*yy)) for xx in hx for yy in hy)
    hxy=max(abs(d+e*(4*xx+2*yy)) for xx in hx for yy in hy)
    hyy=max(abs(2*e*xx) for xx in hx);hxv=max(abs(bb+2*c*vv) for vv in hv)
    hvv=max(abs(2*c*xx) for xx in hx)
    gx,gy,_=b.fixed_tangent((a,bb,c,d,e),x,y,v);gv=x*(bb+2*c*v)
    rounding=abs(px*x+py*y+pv*v+k-b.polynomial((a,bb,c,d,e),x,y,v))
    rounding+=abs(px-gx)*rx+abs(py-gy)*ry+abs(pv-gv)*rv
    return rounding+(hxx*rx*rx+hyy*ry*ry+hvv*rv*rv)/2+hxy*rx*ry+hxv*rx*rv


def subdivision(coeff,affine,box):
    a,bb,c,d,e=map(Q,coeff);px,py,pv,k=map(Q,affine)
    lo=list(map(Q,box[::2]));hi=list(map(Q,box[1::2]));ns=CONFIG['subdivisions'];largest=Q(0)
    for ix in itertools.product(*(range(n) for n in ns)):
        ls=[l+(h-l)*i/n for l,h,i,n in zip(lo,hi,ix,ns)]
        us=[l+(h-l)*(i+1)/n for l,h,i,n in zip(lo,hi,ix,ns)]
        x,y,v=[(l+u)/2 for l,u in zip(ls,us)];rx,ry,rv=[(u-l)/2 for l,u in zip(ls,us)]
        gx,gy,_=b.fixed_tangent((a,bb,c,d,e),x,y,v)
        terms={(0,0,0):px*x+py*y+pv*v+k-b.polynomial((a,bb,c,d,e),x,y,v),
               (1,0,0):px-gx,(0,1,0):py-gy,(0,0,1):pv-x*(bb+2*c*v),
               (2,0,0):-d-e*(3*x+2*y),(1,1,0):-d-2*e*(2*x+y),
               (0,2,0):-e*x,(1,0,1):-bb-2*c*v,(0,0,2):-c*x,
               (3,0,0):-e,(2,1,0):-2*e,(1,2,0):-e,(1,0,2):-c}
        upper=Q(0)
        for powers,coef in terms.items():
            radius=rx**powers[0]*ry**powers[1]*rv**powers[2]
            low=-radius if any(n%2 for n in powers) else (Q(1) if powers==(0,0,0) else Q(0))
            upper+=coef*(radius if coef>=0 else low)
        largest=max(largest,upper)
    return largest


def certificates(fixed_rows,method):
    begin=time.perf_counter();rows=[]
    for old in fixed_rows:
        coeff=b.CONFIG['coefficients'][old['plant']];ref=old['reference'];box=old['box'];affine=full_tangent(coeff,*ref)
        if method=='candidate':shift,count,interior=candidates(coeff,affine,box)
        elif method=='derivative':shift=derivative(coeff,affine,ref,box);count=interior=0
        elif method=='interval':shift=subdivision(coeff,affine,box);count=interior=0
        else:raise ValueError(method)
        delta=b.upward(shift)
        rows.append(dict(plant=old['plant'],time=old['time'],box=box,reference=ref,px=affine[0],py=affine[1],pv=affine[2],
                         tangent_intercept=affine[3],shift=delta,intercept=b.downward(Q(affine[3])-Q(delta)),
                         rational_shift=str(shift),candidate_count=count,interior_volume_candidates=interior))
    return rows,time.perf_counter()-begin


def solve(instance,rows):
    begin=time.perf_counter();model=b.Cascade(instance)
    for r in rows:
        i,t=r['plant'],r['time'];xl,xh,yl,yh,vl,vh=r['box']
        for name,l,u in [('x',xl,xh),('y',yl,yh),('v',vl,vh)]:model.lb[model.index(name,i,t)]=l;model.ub[model.index(name,i,t)]=u
        model.add({model.index('g',i,t):1.,model.index('x',i,t):-r['px'],model.index('y',i,t):-r['py'],model.index('v',i,t):-r['pv']},-np.inf,r['intercept'])
    A=model.sparse();lo=np.array(model.lower);hi=np.array(model.upper);eq=lo==hi
    assembly=time.perf_counter()-begin;start=time.perf_counter()
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Unrecognized options detected:.*')
        result=linprog(model.c,A_ub=A[~eq],b_ub=hi[~eq],A_eq=A[eq],b_eq=hi[eq],bounds=np.c_[model.lb,model.ub],method='highs-ds',
                       options={'primal_feasibility_tolerance':b.CONFIG['lp_tolerance'],'dual_feasibility_tolerance':b.CONFIG['lp_tolerance'],'threads':1})
    record=dict(success=bool(result.success),status=int(result.status),message=result.message,assembly_s=assembly,solve_s=time.perf_counter()-start,procedure_s=time.perf_counter()-begin,objective=float(result.fun) if result.x is not None else None)
    return model.decode(result.x) if result.x is not None else None,record


def audit(instance,point,rows):
    # Translate the storage term at the returned point solely for reusing the
    # independent original physical/linear auditor. This is not a certificate.
    translated=[dict(r,intercept=r['intercept']+r['pv']*point['v'][r['plant']][r['time']]) for r in rows]
    return b.point_audit(instance,point,translated)


def run(output):
    metrics=[];output.mkdir(parents=True,exist_ok=True)
    for n,path in enumerate(sorted((ROOT/'data/test').glob('*/scenario.json'))):
        instance=json.loads(path.read_text());case=instance['id']
        fixed=json.loads((ROOT/'results/test'/case/'candidate/certificates.json').read_text())
        old=json.loads((ROOT/'results/test'/case/'candidate/record.json').read_text())['point_audit']
        order=CONFIG['methods'][n%3:]+CONFIG['methods'][:n%3]
        for method in order:
            folder=output/case/method;folder.mkdir(parents=True,exist_ok=True)
            rows,seconds=certificates(fixed,method)
            for r in rows:
                maximum=oracle(b.CONFIG['coefficients'][r['plant']],(r['px'],r['py'],r['pv'],r['intercept']),r['box'])
                assert maximum[0]<=0
            point,rec=solve(instance,rows);check=audit(instance,point,rows) if point else dict(passed=False)
            b.write_json(folder/'certificates.json',rows);b.write_json(folder/'point.json',point)
            b.write_json(folder/'record.json',dict(construction_s=seconds,solver=rec,point_audit=check))
            metrics.append(dict(instance=case,method=method,construction_s=seconds,solve_s=rec['solve_s'],admissible=check['passed'],
                                **{k:check.get(k) for k in ('scheduled_cost','final_cost','hidden_energy')},
                                scheduled_advantage_over_fixed=old['scheduled_cost']-check['scheduled_cost'] if point else None,
                                final_advantage_over_fixed=old['final_cost']-check['final_cost'] if point else None))
        print('Full-state synthetic',case,flush=True)
    b.write_csv(output/'metrics.csv',metrics)
    b.write_json(output/'summary.json',dict(status='passed',cases=12,LP_solves=len(metrics),admissible=sum(r['admissible'] for r in metrics),
                  exact_exported_rows=len(metrics)*72,candidate_counts_by_plant=[8,12,12],
                  paired_with_original_fixed_storage=True,configuration=CONFIG))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=ROOT/'results/full_state');args=parser.parse_args()
    run(args.output)
