#!/usr/bin/env python3
"""Independent full-volume maxima, Hessian and essential interior-state cases."""
import random, unittest
from fractions import Fraction as Q
import benchmark as b
import full_state as f


def interior_example():
    coeff=b.CONFIG['coefficients'][1];x,y,v=1.,0.,4.
    pv=coeff[1]+2*coeff[2]*v
    affine=(0.,0.,pv,b.polynomial(coeff,x,y,v)-pv*v+.01)
    box=(1.,1.,0.,0.,2.,6.)
    maximum=f.oracle(coeff,affine,box)
    ends=[Q(affine[2])*Q(z)+Q(affine[3])-b.polynomial(list(map(Q,coeff)),Q(x),Q(y),Q(z)) for z in (2.,6.)]
    return dict(plant_family='Baygorria-type curvature; independent synthetic coefficients',coefficients=coeff,affine=list(affine),box=list(box),
                exact_maximum=str(maximum[0]),witness=[str(z) for z in maximum[1:]],endpoint_errors=[str(z) for z in ends])


class FullStateChecks(unittest.TestCase):
    def test_interior_volume_is_essential(self):
        example=interior_example()
        self.assertGreater(Q(example['exact_maximum']),0)
        self.assertTrue(all(Q(v)<0 for v in example['endpoint_errors']))
        self.assertGreater(Q(example['witness'][2]),2)
        self.assertLess(Q(example['witness'][2]),6)
        shift,count,interiors=f.candidates(example['coefficients'],example['affine'],example['box'])
        self.assertEqual(shift,Q(example['exact_maximum']))
        self.assertEqual(count,12);self.assertEqual(interiors,2)

    def test_all_three_bounds_against_exact_whole_box_maximum(self):
        rng=random.Random(2026111305)
        for n in range(36):
            coeff=b.CONFIG['coefficients'][n%3]
            x,y,v=rng.uniform(.1,.9),rng.uniform(.1,.7),rng.uniform(2.5,5.5)
            box=(x-.08,x+.08,y-.04,y+.04,v-.12,v+.12)
            affine=f.full_tangent(coeff,x,y,v)
            if n%6==0:box=(0.,0.,box[2],box[3],box[4],box[5])
            if n%7==0:box=box[:2]+(y,y)+box[4:]
            exact=max(Q(0),f.oracle(coeff,affine,box)[0])
            candidate,count,_=f.candidates(coeff,affine,box)
            self.assertEqual(candidate,exact)
            self.assertLessEqual(count,[8,12,12][n%3])
            for shift in (candidate,f.derivative(coeff,affine,(x,y,v),box),f.subdivision(coeff,affine,box)):
                self.assertGreaterEqual(shift,exact)
                intercept=b.downward(Q(affine[3])-Q(b.upward(shift)))
                self.assertLessEqual(f.oracle(coeff,(*affine[:3],intercept),box)[0],0)

    def test_cross_storage_hessian_by_exact_finite_differences(self):
        h=Q(1,16);x,y,v=Q(2,3),Q(1,4),Q(4)
        for coeff in b.CONFIG['coefficients']:
            a,bb,c,d,e=map(Q,coeff)
            p=lambda xx,vv:b.polynomial((a,bb,c,d,e),xx,y,vv)
            mixed=(p(x+h,v+h)-p(x+h,v-h)-p(x-h,v+h)+p(x-h,v-h))/(4*h*h)
            vv=(p(x,v+h)-2*p(x,v)+p(x,v-h))/(h*h)
            self.assertEqual(mixed,bb+2*c*v);self.assertEqual(vv,2*c*x)


if __name__=='__main__':unittest.main()
