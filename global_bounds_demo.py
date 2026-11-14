#!/usr/bin/env python3
"""Complete rational outer-bound example on a new synthetic delayed cascade.

The input coefficients and profiles are independently distributable. This module
builds/solves models; verify_global_bounds_demo.py checks its evidence separately.
Costs are synthetic monetary units, not real-bank USD observations.
"""
import argparse,copy,json,math,os,warnings
from fractions import Fraction as Q
from pathlib import Path
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parent


def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def polynomial(p,x,y,v):
    a,b,c,d,e=p
    return x*(a+b*v+c*v*v+d*(x+y)+e*(x+y)**2)


def tangent(p,x,y,v):
    a,b,c,d,e=p
    px=a+b*v+c*v*v+d*(2*x+y)+e*(3*x*x+4*x*y+y*y)
    py=x*(d+2*e*(x+y));pv=x*(b+2*c*v)
    return px,py,pv,polynomial(p,x,y,v)-px*x-py*y-pv*v


def quadmax(A,B,D,E,C,ru,rv):
    points=[(u,v) for u in (Q(0),ru) for v in (Q(0),rv)]
    for u in (Q(0),ru):
        if B<0:points.append((u,min(rv,max(Q(0),-D*u/B))))
    for v in (Q(0),rv):
        if A<0:points.append((min(ru,max(Q(0),-(D*v+E)/A)),v))
    determinant=A*B-D*D
    if determinant:
        u=-B*E/determinant;v=D*E/determinant
        if 0<=u<=ru and 0<=v<=rv:points.append((u,v))
    return max(A*u*u/2+D*u*v+B*v*v/2+E*u+C for u,v in points)


def upper(p,center,box):
    x,y,v=center;a,b,c,d,e=p;xl,xh,yl,yh,vl,vh=box
    rx=max(abs(x-xl),abs(x-xh));ry=max(abs(y-yl),abs(y-yh));rv=max(abs(v-vl),abs(v-vh))
    X=(min(x,xl),max(x,xh));Y=(min(y,yl),max(y,yh));V=(min(v,vl),max(v,vh))
    A=max(2*d+e*(6*xx+4*yy) for xx in X for yy in Y)
    B=max(2*c*xx for xx in X);D=max(abs(b+2*c*vv) for vv in V)
    E=max(abs(d+e*(4*xx+2*yy)) for xx in X for yy in Y)*ry
    C=max(Q(0),max(2*e*xx for xx in X))*ry*ry/2
    shift=max(Q(0),quadmax(A,B,D,E,C,rx,rv));affine=tangent(p,x,y,v)
    exact_intercept=affine[3]+shift;rounded=float(exact_intercept)
    if Q(rounded)<exact_intercept:rounded=math.nextafter(rounded,math.inf)
    assert Q(rounded)>=exact_intercept
    return (*affine[:3],Q(rounded)),shift


def lower(p,center,box):
    a,b,c,d,e=p;px,py,pv,k=tangent(p,*center);xl,xh,yl,yh,vl,vh=box
    assert xl>=0 and e>=0 and min(-2*d-e*(6*x+4*y) for x in (xl,xh) for y in (yl,yh))>=0
    def quadratic(A,B,l,u):
        points=[l,u]
        if A<0:points.append(min(u,max(l,-B/(2*A))))
        return max(A*z*z+B*z for z in points)
    delta=max(Q(0),max(k+px*x-a*x-d*x*x-e*x*x*x+
        quadratic(-e*x,py-d*x-2*e*x*x,yl,yh)+quadratic(-c*x,pv-b*x,vl,vh) for x in (xl,xh)))
    return px,py,pv,k-delta


def read(path):
    raw=json.loads(path.read_text());data=dict(raw)
    for key in ('coefficients','inflows','anchor'):
        data[key]=[list(map(Q,row)) for row in raw[key]]
    for key in ('demand','initial','water_values'):data[key]=list(map(Q,raw[key]))
    for key in ('hydraulic_factor','failure_cost','radius','margin','xmax','ymax','vmin','vmax'):data[key]=Q(raw[key])
    return data


def base(data):
    T=data['periods'];kappa=data['hydraulic_factor'];bounds=[];eq=[];c=[Q(0)]*(9*T)
    for t in range(T):
        anchor=data['anchor'][t]
        for j in range(8):
            if j<6:
                cap=data['xmax' if j%2==0 else 'ymax'];width=data['radius']*cap
                bounds.append((max(Q(0),anchor[j]-width),min(cap,anchor[j]+width)))
            else:bounds.append((max(data['vmin'],anchor[j]*(1-data['margin'])),min(data['vmax'],anchor[j]*(1+data['margin']))))
        bounds.append((Q(0),data['demand'][t]));c[9*t+8]=-data['failure_cost']
        for i,flow,state in ((0,0,6),(2,4,7)):
            row={9*t+state:Q(1),9*t+flow:kappa,9*t+flow+1:kappa}
            rhs=kappa*data['inflows'][i][t]
            if t:row[9*(t-1)+state]=Q(-1)
            else:rhs+=data['initial'][0 if i==0 else 1]
            if i==2 and t>=data['delays'][1]:
                for j in (2,3):row[9*(t-data['delays'][1])+j]=-kappa
            eq.append((row,rhs))
        row={9*t+2:Q(1),9*t+3:Q(1)}
        if t>=data['delays'][0]:
            for j in (0,1):row[9*(t-data['delays'][0])+j]=Q(-1)
        eq.append((row,data['inflows'][1][t]))
    for i in (0,1):c[9*(T-1)+6+i]=-data['water_values'][i]
    constant=data['failure_cost']*sum(data['demand'])+sum(a*b for a,b in zip(data['initial'],data['water_values']))
    return dict(bounds=bounds,equalities=eq,inequalities=[],cost=c,constant=constant)


def box(model,i,t):
    return (*model['bounds'][9*t+2*i],*model['bounds'][9*t+2*i+1],*model['bounds'][9*t+(6 if i==0 else 7)])


def add(model,planes):
    for t,row in enumerate(planes):
        lhs={9*t+8:Q(1)};rhs=Q(0)
        for i,(px,py,pv,k) in enumerate(row):
            lhs[9*t+2*i]=-px;lhs[9*t+2*i+1]=-py
            j=9*t+(6 if i==0 else 7);lhs[j]=lhs.get(j,Q(0))-pv;rhs+=k
        model['inequalities'].append((lhs,rhs))


def solve(model):
    n=len(model['cost'])
    def matrix(rows):return [[float(row.get(j,Q(0))) for j in range(n)] for row,_ in rows]
    eq=model['equalities'];ub=model['inequalities']
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Unrecognized options detected:.*')
        result=linprog(list(map(float,model['cost'])),A_eq=matrix(eq),b_eq=[float(b) for _,b in eq],
            A_ub=matrix(ub),b_ub=[float(b) for _,b in ub],bounds=[tuple(map(float,b)) for b in model['bounds']],
            method='highs-ds',options={'threads':1,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    if not result.success:raise RuntimeError(result.message)
    return result,[[Q(float(result.x[9*t+j])) for j in range(9)] for t in range(n//9)]


def encode(model):
    rows=lambda rs:[[{str(k):str(v) for k,v in row.items()},str(rhs)] for row,rhs in rs]
    return dict(bounds=[[str(l),str(u)] for l,u in model['bounds']],cost=list(map(str,model['cost'])),constant=str(model['constant']),
                equalities=rows(model['equalities']),inequalities=rows(model['inequalities']))


def dual(model,result):
    lam=[max(Q(0),-Q(float(v))) for v in result.ineqlin.marginals]
    nu=[-Q(float(v)) for v in result.eqlin.marginals];r=model['cost'][:];offset=model['constant']
    for multiplier,(row,rhs) in zip(lam+nu,model['inequalities']+model['equalities']):
        offset-=multiplier*rhs
        for j,v in row.items():r[j]+=multiplier*v
    correction=sum(min(v*l,v*u) for v,(l,u) in zip(r,model['bounds']))
    return dict(inequality_multipliers=list(map(str,lam)),equality_multipliers=list(map(str,nu)),residual=list(map(str,r)),
                lagrangian_constant=str(offset),box_correction=str(correction),lower=str(offset+correction))


def realize(data,model,target):
    T=data['periods'];dt=data['hydraulic_factor'];out=[[Q(0)]*9 for _ in range(T)];bounds=model['bounds']
    def clip(v,l,u):
        if l>u:raise ValueError('Empty exact reachability interval; no witness claimed')
        return min(u,max(l,v))
    def reservoir(flow,state,initial,inflows,extra=None):
        releases=[];lo=[];hi=[]
        for t in range(T):
            xl,xh=bounds[9*t+flow];yl,yh=bounds[9*t+flow+1];l=xl+yl;u=xh+yh
            if extra is not None:l=max(l,extra[t][0]);u=min(u,extra[t][1])
            if l>u:raise ValueError('Empty exact release interval')
            releases.append((l,u));l,u=bounds[9*t+state];lo.append(l);hi.append(u)
        for t in range(T-2,-1,-1):
            lo[t]=max(lo[t],lo[t+1]-dt*(inflows[t+1]-releases[t+1][0]))
            hi[t]=min(hi[t],hi[t+1]-dt*(inflows[t+1]-releases[t+1][1]))
        previous=initial
        for t in range(T):
            l=max(lo[t],previous+dt*(inflows[t]-releases[t][1]));u=min(hi[t],previous+dt*(inflows[t]-releases[t][0]))
            v=clip(target[t][state],l,u);R=inflows[t]+(previous-v)/dt
            xl,xh=bounds[9*t+flow];yl,yh=bounds[9*t+flow+1]
            x=clip(target[t][flow],max(xl,R-yh),min(xh,R-yl))
            out[t][flow]=x;out[t][flow+1]=R-x;out[t][state]=v;previous=v
    extra=[];delay=data['delays'][0]
    for t in range(T):
        if t+delay<T:
            j=t+delay;xl,xh=bounds[9*j+2];yl,yh=bounds[9*j+3];a=data['inflows'][1][j]
            extra.append((xl+yl-a,xh+yh-a))
        else:extra.append((Q(0),data['xmax']+data['ymax']))
    reservoir(0,6,data['initial'][0],data['inflows'][0],extra)
    for t in range(T):
        R=data['inflows'][1][t]+(sum(out[t-delay][:2]) if t>=delay else Q(0))
        xl,xh=bounds[9*t+2];yl,yh=bounds[9*t+3]
        out[t][2]=clip(target[t][2],max(xl,R-yh),min(xh,R-yl));out[t][3]=R-out[t][2]
    delay=data['delays'][1]
    inflows=[data['inflows'][2][t]+(out[t-delay][2]+out[t-delay][3] if t>=delay else Q(0)) for t in range(T)]
    reservoir(4,7,data['initial'][1],inflows)
    final=copy.deepcopy(out)
    for t in range(T):
        P=sum(polynomial(data['coefficients'][i],out[t][2*i],out[t][2*i+1],out[t][6 if i==0 else 7]) for i in range(3))
        assert P>=0
        out[t][8]=min(max(Q(0),target[t][8]),data['demand'][t],P)
        final[t]=out[t][:];final[t][8]=min(data['demand'][t],P)
    cost=lambda point:model['constant']+sum(c*point[j//9][j%9] for j,c in enumerate(model['cost']))
    return out,final,cost(out),cost(final)


def generate(path,output):
    data=read(path);model=base(data);rounds=[];centers=data['anchor'];lower_values=[]
    for k in range(3):
        planes=[];records=[]
        for t in range(data['periods']):
            row=[]
            for i in range(3):
                center=(centers[t][2*i],centers[t][2*i+1],centers[t][6 if i==0 else 7])
                aff,shift=upper(data['coefficients'][i],center,box(model,i,t));row.append(aff)
                records.append(dict(plant=i,time=t,center=list(map(str,center)),affine=list(map(str,aff)),majorant_shift=str(shift)))
            planes.append(row)
        add(model,planes);result,centers=solve(model);bound=dual(model,result);lower_values.append(Q(bound['lower']))
        rounds.append(dict(index=k,production_upper_planes=records,model=encode(model),dual=bound,numerical_LP_success=True))
    inner=base(data);lower_rows=[[lower(data['coefficients'][i],(data['anchor'][t][2*i],data['anchor'][t][2*i+1],data['anchor'][t][6 if i==0 else 7]),box(inner,i,t)) for i in range(3)] for t in range(data['periods'])]
    add(inner,lower_rows);result,target=solve(inner);scheduled,final,Js,Jf=realize(data,inner,target)
    witness=dict(scheduled_point=[[str(v) for v in row] for row in scheduled],final_point=[[str(v) for v in row] for row in final],
                 scheduled_cost=str(Js),final_cost=str(Jf),maximum_physical_change=str(max(abs(target[t][j]-scheduled[t][j]) for t in range(data['periods']) for j in range(8))))
    best=max(lower_values);assert best<=Jf<=Js
    norm=data['failure_cost']*sum(data['demand'])
    evidence=dict(format='global-bound-example-v1',input_file=path.name,rounds=rounds,witness=witness,
        lower=str(best),scheduled_upper=str(Js),final_upper=str(Jf),scheduled_gap_cap=str(Js-best),final_gap_cap=str(Jf-best),
        normalization_scale=str(norm),normalized_final_cap_percent=str(100*(Jf-best)/norm),LP_solves=4,
        scope='Continuous nonlinear reference optimum in the declared local boxes; synthetic monetary units; no original MILP global claim')
    dump(output,evidence)
    return evidence


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,default=ROOT/'global_bounds_example/input.json');ap.add_argument('--output',type=Path,default=ROOT/'global_bounds_example/certificate.json');args=ap.parse_args()
    r=generate(args.input,args.output);print(json.dumps({k:v for k,v in r.items() if k not in ('rounds','witness')},indent=2))
