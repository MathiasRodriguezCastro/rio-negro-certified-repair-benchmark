#!/usr/bin/env python3
"""Independent maxima, weak-concavity cases and missing-hypothesis examples."""
import math,random,unittest
from fractions import Fraction as Q
from benchmark import (CONFIG,exact_candidates,square_oracle,certificates,
                       audit_certificates,upward,downward,geometry,polynomial,control_hessian)


def concave_slice_candidate(derivative,lo,hi):
    if derivative(lo)<=0:return lo
    if derivative(hi)>=0:return hi
    for _ in range(100):
        middle=(lo+hi)/2
        if derivative(middle)>0:lo=middle
        else:hi=middle
    return hi


class GeneralReductionTest(unittest.TestCase):
    def test_control_hessian_by_exact_finite_differences(self):
        h=Q(1,16)
        for values in CONFIG['coefficients']:
            coeff=list(map(Q,values))
            for x,y,v in [(Q(1,3),Q(1,5),Q(3)),(Q(3,4),Q(2,3),Q(5))]:
                P=lambda u,w:polynomial(coeff,u,w,v)
                xx=(P(x+h,y)-2*P(x,y)+P(x-h,y))/(h*h)
                yy=(P(x,y+h)-2*P(x,y)+P(x,y-h))/(h*h)
                xy=(P(x+h,y+h)-P(x+h,y-h)-P(x-h,y+h)+P(x-h,y-h))/(4*h*h)
                self.assertEqual(control_hessian(coeff,x,y),(xx,xy,yy))

    def test_nonquadratic_and_nonpolynomial_slices(self):
        for center in (-.7,.2,.8):
            y=concave_slice_candidate(lambda z:-4*(z-center)**3,-1.,1.)
            self.assertAlmostEqual(y,center,places=12)
            candidates=[x*x-(y-center)**4+2 for x in (-.5,1.)]
            self.assertAlmostEqual(max(candidates),3.)
        y=concave_slice_candidate(lambda z:2-math.exp(z),-2.,2.)
        self.assertAlmostEqual(y,math.log(2),places=12)
        self.assertAlmostEqual(1-math.exp(y)+2*y, -1+2*math.log(2),places=12)

    def test_flat_stationary_interval(self):
        # C2 concave function with a whole interval of maximizers.
        def h(z):return -max(0.,abs(z)-.25)**4
        def derivative(z):return -4*max(0.,abs(z)-.25)**3*(1 if z>=0 else -1)
        y=concave_slice_candidate(derivative,-1.,1.)
        self.assertAlmostEqual(y,-.25,places=12)
        self.assertEqual(h(y),0.)
        self.assertEqual(h(0.),0.)
        self.assertEqual(h(.25),0.)

    def test_affine_and_degenerate_slices(self):
        self.assertEqual(concave_slice_candidate(lambda z:0.,2.,4.),2.)
        self.assertEqual(concave_slice_candidate(lambda z:1.,2.,4.),4.)
        self.assertEqual(concave_slice_candidate(lambda z:-1.,2.,4.),2.)
        self.assertEqual(concave_slice_candidate(lambda z:0.,3.,3.),3.)

    def test_omitted_hypotheses_can_fail(self):
        # Missing x convexity: interior x beats every x endpoint.
        F=lambda x,y,v:-(x-.5)**2-y*y+v
        self.assertGreater(F(.5,0,1),max(F(x,0,1) for x in (0,1)))
        # Missing slice concavity: picking an arbitrary stationary point fails.
        h=lambda y:-(y*y-1)**2
        self.assertGreater(h(1),max(h(-2),h(0),h(2)))
        # Missing state monotonicity: neither state endpoint maximizes F.
        self.assertGreater(-(.5-.5)**2,max(-(v-.5)**2 for v in (0,1)))

    def test_multiple_monotone_states_and_convex_controls(self):
        # Two convex coordinates: four vertices, rather than only two endpoints.
        values=[x*x+2*u*u+v-3*w for x in (-1,2) for u in (-2,1) for v in (0,1) for w in (-1,1)]
        self.assertEqual(max(values),4+8+1+3)

    def test_exact_polynomial_oracle_and_float_export(self):
        rng=random.Random(20261113)
        interior=0
        for n in range(120):
            coeff=CONFIG['coefficients'][n%3]
            x0=rng.uniform(.1,.9);y0=rng.uniform(.1,.7);v0=rng.uniform(2.5,5.5)
            box=(max(0.,x0-.1),min(1.,x0+.1),max(0.,y0-.12),min(.8,y0+.12),v0-.1,v0+.1)
            tangent=(rng.uniform(1,2),rng.uniform(-.1,.1),rng.uniform(-.2,.2))
            if n%3==2:
                x=box[1];target=(box[2]+box[3])/2
                tangent=(tangent[0],x*(coeff[3]+2*coeff[4]*(x+target)),tangent[2])
                interior+=1
            exact=exact_candidates(coeff,tangent,box)
            self.assertEqual(exact,max(Q(0),square_oracle(coeff,tangent,box)))
            shift=upward(exact);intercept=downward(Q(tangent[2])-Q(shift))
            self.assertGreaterEqual(Q(shift),exact)
            self.assertLessEqual(square_oracle(coeff,(tangent[0],tangent[1],intercept),box),0)
        self.assertEqual(interior,40)

    def test_invalid_geometry_is_rejected(self):
        with self.assertRaises(AssertionError):geometry([1,.1,0,.1,.1],(0,1,0,1,2,6))
        with self.assertRaises(ValueError):geometry([1,-.6,.1,-1,0],(0,1,0,1,2,6))


if __name__=='__main__':unittest.main()
