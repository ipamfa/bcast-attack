import unittest
from L3 import *

# Test Gram-schmidt, ambil dari Bremner      
class GramSTest(unittest.TestCase) :
    
    def test_v1equal_bstar1(self):
        X = np.array([
            [3, -1, 5],
            [-5, 2, -1],
            [-3, 9, 2]
        ])
        gram(X)
        #
        y = getBstar(0)
        x = X[0]
        self.assertEqual( all( x==y ), True)
        
    def test_bstar_before_run_isNone(self) :
        # test for 3 dim
        self.assertIsNone( getBstar(0) )
        self.assertIsNone( getBstar(1) )
        self.assertIsNone( getBstar(2) )
        
    def test_lll_fullrank_3dim(self) :
        # apply reduction algorithm
        b = np.array([
            [-2, 7, 7, -5],
            [3, -2, 6, -1],
            [2, -8, -9, -7],
            [8, -9, 6, -4],
        ])
        LLL(b, 1.0)

        
Y = np.array([
    [-11,  12,  -4],
    [3, -1, 5],
    [-5, 2, -1],
    [-3, 9, 2]
])
Z = np.array([
    [13, 15, -5,  8],
    [8, 6, 2, 6],
    [2, 4, 8, 4],
    [1, -1, -6, -8],
    [2, 6, -9, 6],
])

if __name__ == "__main__" :
    unittest.main()
