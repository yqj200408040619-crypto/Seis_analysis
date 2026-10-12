"""Optional solver cross-check: units, root completeness and mode indexing."""
import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dispersion.inversion import forward,independent_phase

class ForwardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try: import disba,pysurf96
        except ImportError: raise unittest.SkipTest("Optional disba/pysurf96 not installed")

    def test_homogeneous_analytic_rayleigh_root(self):
        # β=Vs²/Vp²=1/4; Rayleigh secular polynomial in x=(c/Vs)².
        roots=np.roots([1.,-8.,20.,-12.]);x=[r.real for r in roots if abs(r.imag)<1e-10 and 0<r.real<1][0]
        model=np.array([[0.,1.,.5,2.]])
        f=np.array([2.,3.,4.,5.]);expected=np.sqrt(x)*500
        p=forward(model,f,0);q=independent_phase(model,f,0)
        self.assertLess(abs(p-expected).max(),.5);self.assertLess(abs(p-q).max(),.5)

    def test_layered_fundamental_and_first_higher_mode(self):
        vs=np.array([.2806,.566,.789]);model=np.column_stack(([.0204,.1003,0.],2*vs,vs,[2.,2.,2.]))
        f=np.linspace(4.3,5.,15)
        for mode in (0,1):
            p=forward(model,f,mode);q=independent_phase(model,f,mode)
            self.assertLess(abs(p-q).max(),.5);self.assertTrue((p<789.).all())

    def test_independent_solver_more_than_sixty_frequencies(self):
        model=np.array([[.06,1.1,.55,2.],[0.,1.6,.8,2.]])
        f=np.linspace(2.,5.,131)
        p=forward(model,f,0);q=independent_phase(model,f,0)
        self.assertEqual(len(q),131);self.assertLess(abs(p-q).max(),.5)

if __name__=="__main__": unittest.main()
