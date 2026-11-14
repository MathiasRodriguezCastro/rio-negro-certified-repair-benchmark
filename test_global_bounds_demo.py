#!/usr/bin/env python3
"""Whole-box and adversarial checks for the public global-bound certificate."""
import copy,json,unittest
from fractions import Fraction as F
from pathlib import Path
import verify_global_bounds_demo as v
ROOT=Path(__file__).resolve().parent


class GlobalBoundDemoChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((ROOT/'global_bounds_example/input.json').read_text())
        cls.evidence=json.loads((ROOT/'global_bounds_example/certificate.json').read_text())

    def test_whole_box_subdivision_and_hidden_interior_failure(self):
        box=((F(0),F(1)),(F(0),F(1)),(F(0),F(1)))
        good={(2,0,0):F(1),(1,0,0):F(-1),(0,0,0):F(26,100)}
        self.assertLess(v.bernstein_lower(good,box),0)
        self.assertGreater(v.prove_nonnegative(good,box),1)
        bad=dict(good);bad[(0,0,0)]=F(24,100)
        self.assertGreater(v.eval_poly(bad,(F(0),F(0),F(0))),0)
        self.assertGreater(v.eval_poly(bad,(F(1),F(1),F(1))),0)
        with self.assertRaises(AssertionError):v.prove_nonnegative(bad,box)

    def test_arbitrary_dual_residual_and_required_multiplier_sign(self):
        model=dict(cost=['-1'],constant='0',bounds=[['0','1']],equalities=[],inequalities=[[{'0':'1'},'1/3']])
        for lam,expected in ((F(0),F(-1)),(F(1),F(-1,3)),(F(2),F(-2,3))):
            residual=lam-1;offset=-lam/3;correction=min(F(0),residual)
            dual=dict(inequality_multipliers=[str(lam)],equality_multipliers=[],residual=[str(residual)],
                      lagrangian_constant=str(offset),box_correction=str(correction),lower=str(expected))
            self.assertEqual(v.check_dual(model,dual),expected)
        dual['inequality_multipliers']=['-1']
        with self.assertRaises(AssertionError):v.check_dual(model,dual)

    def test_complete_stored_certificate(self):
        r=v.verify(self.data,self.evidence)
        self.assertEqual(r['whole_box_upper_planes'],216)
        self.assertEqual(r['exact_feasible_reference_points'],2)
        self.assertEqual(r['nonzero_residual_rounds'],3)

    def test_corrupt_upper_plane_and_dual_are_rejected(self):
        bad=copy.deepcopy(self.evidence);bad['rounds'][0]['production_upper_planes'][0]['affine'][3]='-100'
        with self.assertRaises(AssertionError):v.verify(self.data,bad)
        bad=copy.deepcopy(self.evidence);bad['rounds'][0]['dual']['lower']=str(F(bad['rounds'][0]['dual']['lower'])+1)
        with self.assertRaises(AssertionError):v.verify(self.data,bad)

    def test_inexact_hydraulics_and_unsupported_service_are_rejected(self):
        bad=copy.deepcopy(self.evidence);p=bad['witness']['scheduled_point'];p[0][6]=str(F(p[0][6])+F(1,10000))
        with self.assertRaises(AssertionError):v.verify(self.data,bad)
        bad=copy.deepcopy(self.evidence);bad['witness']['final_point'][0][8]='2'
        with self.assertRaises(AssertionError):v.verify(self.data,bad)

    def test_objective_constant_and_finite_box_are_required(self):
        bad=copy.deepcopy(self.evidence);bad['rounds'][0]['model']['constant']='0'
        with self.assertRaises(AssertionError):v.verify(self.data,bad)
        bad=copy.deepcopy(self.evidence);bad['rounds'][0]['model']['bounds'][0][1]='inf'
        with self.assertRaises((ValueError,AssertionError)):v.verify(self.data,bad)


if __name__=='__main__':unittest.main()
