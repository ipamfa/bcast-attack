import unittest
from L3 import *
from kannan import orthoproject

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
        
    def test_orthogonal_project_to_b0(self) :
        print "orthogonal projection to 1st vector basis"
        print "========================================="
        b1 = np.array([
            [4,5,1],
            [4,8,2],
            [6,2,6]
        ])
        b1_ = orthoproject(b1)
        print b1_
        
        b2 = np.array([
            [-2, 7, 7, -5],
            [3, -2, 6, -1],
            [2, -8, -9, -7],
            [8, -9, 6, -4],
        ])
        b2_ = orthoproject(b2)
        print b2_

if __name__ == "__main__" :
    unittest.main()
