#!/usr/bin/env python3
"""Independent stdlib-only verifier: whole-box Bernstein proof, dual, witnesses.

Does not import the generator, SciPy, NumPy or a solver. Bernstein coefficients
certify an entire box; checking point samples is never accepted as a proof.
"""
import argparse,itertools,json,math
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent


def eval_poly(terms,point):
    return sum(c*math.prod(z**p for z,p in zip(point,powers)) for powers,c in terms.items())


def envelope_difference(coeff,affine):
    a,b,c,d,e=map(F,coeff);px,py,pv,k=map(F,affine)
    terms={(0,0,0):k,(1,0,0):px-a,(0,1,0):py,(0,0,1):pv,(1,0,1):-b,(1,0,2):-c,
           (2,0,0):-d,(1,1,0):-d,(3,0,0):-e,(2,1,0):-2*e,(1,2,0):-e}
    return {p:c for p,c in terms.items() if c}


def bernstein_lower(terms,box):
    """Convex-hull lower enclosure on a physical three-dimensional box."""
    degrees=tuple(max(p[j] for p in terms) for j in range(3)) if terms else (0,0,0)
    coefficients={}
    for powers,value in terms.items():
        for exponents in itertools.product(*(range(n+1) for n in powers)):
            term=value
            for j,(n,k) in enumerate(zip(powers,exponents)):
                lo,hi=box[j];term*=math.comb(n,k)*lo**(n-k)*(hi-lo)**k
            coefficients[exponents]=coefficients.get(exponents,F(0))+term
    values=[]
    for index in itertools.product(*(range(n+1) for n in degrees)):
        value=F(0)
        for powers,c in coefficients.items():
            if all(k<=i for k,i in zip(powers,index)):
                value+=c*math.prod(F(math.comb(i,k),math.comb(n,k)) for k,i,n in zip(powers,index,degrees))
        values.append(value)
    return min(values)


def prove_nonnegative(terms,box,depth=0):
    assert all(l<=u for l,u in box),'Empty box'
    # A negative corner rejects a bad envelope immediately; nonnegative corners
    # alone are never sufficient to accept it.
    assert all(eval_poly(terms,p)>=0 for p in itertools.product(*[(l,u) for l,u in box])),'Envelope fails at a corner'
    if bernstein_lower(terms,box)>=0:return 1
    assert depth<12,'Whole-box Bernstein verification inconclusive within its fixed depth'
    j=depth%3
    if box[j][0]==box[j][1]:
        choices=[i for i in range(3) if box[i][0]<box[i][1]]
        assert choices,'Negative constant on a singleton box';j=choices[0]
    left=list(box);right=list(box);l,u=box[j];mid=(l+u)/2
    left[j]=(l,mid);right[j]=(mid,u)
    return prove_nonnegative(terms,tuple(left),depth+1)+prove_nonnegative(terms,tuple(right),depth+1)


def normalized_rows(raw):
    return [({int(j):F(v) for j,v in row.items() if F(v)},F(rhs)) for row,rhs in raw]


def canonical(data):
    """Physical unscaled constraints independently reconstructed from the input."""
    T=data['periods'];K=F(data['hydraulic_factor']);bounds=[];eq=[]
    for t in range(T):
        for j in range(8):
            v=F(data['anchor'][t][j])
            if j>=6:
                bounds.append((max(F(data['vmin']),v*(1-F(data['margin']))),min(F(data['vmax']),v*(1+F(data['margin'])))))
            else:
                cap=F(data['xmax'] if j%2==0 else data['ymax']);r=F(data['radius'])*cap
                bounds.append((max(F(0),v-r),min(cap,v+r)))
        bounds.append((F(0),F(data['demand'][t])))
        for reservoir,flow,state,inflow in ((0,0,6,0),(1,4,7,2)):
            lhs={t*9+state:F(1),t*9+flow:K,t*9+flow+1:K}
            rhs=K*F(data['inflows'][inflow][t])
            if t==0:rhs+=F(data['initial'][reservoir])
            else:lhs[(t-1)*9+state]=F(-1)
            if reservoir==1 and t>=data['delays'][1]:
                lag=t-data['delays'][1];lhs[lag*9+2]=-K;lhs[lag*9+3]=-K
            eq.append((lhs,rhs))
        lhs={t*9+2:F(1),t*9+3:F(1)}
        if t>=data['delays'][0]:
            lag=t-data['delays'][0];lhs[lag*9]=F(-1);lhs[lag*9+1]=F(-1)
        eq.append((lhs,F(data['inflows'][1][t])))
    cost=[F(0)]*(9*T)
    for t in range(T):cost[t*9+8]=-F(data['failure_cost'])
    cost[(T-1)*9+6]=-F(data['water_values'][0]);cost[(T-1)*9+7]=-F(data['water_values'][1])
    constant=F(data['failure_cost'])*sum(map(F,data['demand']))+sum(F(x)*F(y) for x,y in zip(data['initial'],data['water_values']))
    return bounds,eq,cost,constant


def check_dual(model,dual):
    cost=list(map(F,model['cost']));bounds=[tuple(map(F,b)) for b in model['bounds']]
    ub=normalized_rows(model['inequalities']);eq=normalized_rows(model['equalities'])
    lam=list(map(F,dual['inequality_multipliers']));nu=list(map(F,dual['equality_multipliers']))
    assert len(lam)==len(ub) and len(nu)==len(eq) and len(cost)==len(bounds)
    assert all(v>=0 for v in lam) and all(l<=u for l,u in bounds)
    residual=cost[:];offset=F(model['constant'])
    for multiplier,(row,rhs) in zip(lam+nu,ub+eq):
        offset-=multiplier*rhs
        for j,v in row.items():assert 0<=j<len(cost);residual[j]+=multiplier*v
    correction=sum(min(r*l,r*u) for r,(l,u) in zip(residual,bounds))
    assert list(map(F,dual['residual']))==residual
    assert F(dual['lagrangian_constant'])==offset and F(dual['box_correction'])==correction
    value=offset+correction;assert value==F(dual['lower'])
    return value


def check_point(data,point,bounds):
    T=data['periods'];assert len(point)==T and all(len(row)==9 for row in point)
    available=[]
    for t,row in enumerate(point):
        for j,z in enumerate(row):assert bounds[t*9+j][0]<=z<=bounds[t*9+j][1],'Point outside physical box'
        previous=(list(map(F,data['initial'])) if t==0 else [point[t-1][6],point[t-1][7]])
        up1=sum(point[t-data['delays'][0]][:2]) if t>=data['delays'][0] else F(0)
        up2=sum(point[t-data['delays'][1]][2:4]) if t>=data['delays'][1] else F(0)
        K=F(data['hydraulic_factor'])
        assert row[6]-previous[0]==K*(F(data['inflows'][0][t])-row[0]-row[1]),'First-reservoir balance'
        assert row[2]+row[3]==F(data['inflows'][1][t])+up1,'Run-of-river balance'
        assert row[7]-previous[1]==K*(F(data['inflows'][2][t])+up2-row[4]-row[5]),'Last-reservoir balance'
        powers=[]
        for i,p in enumerate(data['coefficients']):
            a,b,c,d,e=map(F,p);x,y,v=row[2*i],row[2*i+1],row[6 if i==0 else 7]
            powers.append(a*x+b*x*v+c*x*v*v+d*x*x+d*x*y+e*x*x*x+2*e*x*x*y+e*x*y*y)
        assert min(powers)>=0 and row[8]<=sum(powers),'Service has no nonlinear production support'
        available.append(sum(powers))
    water=sum(F(w)*(F(v0)-point[-1][j]) for w,v0,j in zip(data['water_values'],data['initial'],(6,7)))
    cost=water+F(data['failure_cost'])*sum(F(d)-row[8] for d,row in zip(data['demand'],point))
    return cost,available


def verify(data,evidence):
    assert evidence['format']=='global-bound-example-v1' and len(evidence['rounds'])==3
    bounds,eq,cost,constant=canonical(data);T=data['periods'];ubs=[];lower=[];planes=leaves=0
    check_point(data,[list(map(F,row)) for row in data['anchor']],bounds)
    for index,record in enumerate(evidence['rounds']):
        assert record['index']==index
        rs=record['production_upper_planes'];assert [(r['time'],r['plant']) for r in rs]==[(t,i) for t in range(T) for i in range(3)]
        for r in rs:
            i,t=r['plant'],r['time'];offset=t*9;dims=(offset+2*i,offset+2*i+1,offset+(6 if i==0 else 7))
            physical_box=tuple(bounds[j] for j in dims)
            # Algebraically independent whole-box verification of the exported
            # affine plane: no call to the signed-Hessian construction.
            terms=envelope_difference(data['coefficients'][i],r['affine'])
            leaves+=prove_nonnegative(terms,physical_box);planes+=1
        for t in range(T):
            lhs={9*t+8:F(1)};rhs=F(0)
            for r in rs[3*t:3*t+3]:
                i=r['plant'];px,py,pv,k=map(F,r['affine']);lhs[9*t+2*i]=-px;lhs[9*t+2*i+1]=-py
                j=9*t+(6 if i==0 else 7);lhs[j]=lhs.get(j,F(0))-pv;rhs+=k
            ubs.append(({j:v for j,v in lhs.items() if v},rhs))
        model=record['model']
        assert [tuple(map(F,b)) for b in model['bounds']]==bounds
        assert list(map(F,model['cost']))==cost and F(model['constant'])==constant
        assert normalized_rows(model['equalities'])==eq and normalized_rows(model['inequalities'])==ubs
        lower.append(check_dual(model,record['dual']))
    witness=evidence['witness'];points=[[[F(v) for v in row] for row in witness[k]] for k in ('scheduled_point','final_point')]
    Js,_=check_point(data,points[0],bounds);Jf,available=check_point(data,points[1],bounds)
    assert all(a[:8]==b[:8] for a,b in zip(*points))
    assert all(row[8]==min(F(d),p) for row,d,p in zip(points[1],data['demand'],available))
    assert Js==F(witness['scheduled_cost'])==F(evidence['scheduled_upper'])
    assert Jf==F(witness['final_cost'])==F(evidence['final_upper'])
    for point,expected in zip(points,(Js,Jf)):
        assert constant+sum(c*point[j//9][j%9] for j,c in enumerate(cost))==expected
    assert F(evidence['lower'])==max(lower)<=Jf<=Js
    assert F(evidence['scheduled_gap_cap'])==Js-max(lower) and F(evidence['final_gap_cap'])==Jf-max(lower)
    scale=F(data['failure_cost'])*sum(map(F,data['demand']));assert scale>0 and scale==F(evidence['normalization_scale'])
    assert F(evidence['normalized_final_cap_percent'])==100*(Jf-max(lower))/scale
    return dict(status='passed',whole_box_upper_planes=planes,exact_Bernstein_leaf_certificates=leaves,
                exact_dual_bounds=3,exact_feasible_reference_points=2,stored_LP_solves=evidence['LP_solves'],optimization_solves_in_verifier=0,
                nonzero_residual_rounds=sum(any(F(v) for v in r['dual']['residual']) for r in evidence['rounds']),
                lower=str(max(lower)),scheduled_upper=str(Js),final_upper=str(Jf),final_gap_cap=str(Jf-max(lower)),
                scale='Synthetic monetary units; normalized cap relative to all-demand failure penalty, not optimal cost',
                exact_hydraulic_anchor_checked=True,generator_imported=False,proprietary_solver_required=False,network_required=False)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,default=ROOT/'global_bounds_example/input.json');ap.add_argument('--certificate',type=Path,default=ROOT/'global_bounds_example/certificate.json');args=ap.parse_args()
    result=verify(json.loads(args.input.read_text()),json.loads(args.certificate.read_text()));print(json.dumps(result,indent=2))
