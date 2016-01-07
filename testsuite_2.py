import unittest
from L3 import *
from kannan import orthoproject, shortest

# Test Gram-schmidt, ambil dari Bremner      
class GramSTest(unittest.TestCase) :
    
    def test_v1equal_bstar1(self):
        X = np.array([
            [3, -1, 5],
            [-5, 2, -1],
            [-3, 9, 2]
        ])
        g = Gram(X)
        #
        y = g.getBstar(0)
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
        b1_ = orthoproject(b1[0], b1[1:])
        print b1_
        
        b2 = np.array([
            [-2, 7, 7, -5],
            [3, -2, 6, -1],
            [2, -8, -9, -7],
            [8, -9, 6, -4],
        ])
        b2_ = orthoproject(b2[0], b2[1:])
        print b2_

    def test_orthoproject_2(self) :
        x = np.array([[4.0, 0.2, -0.6],
                      [0.9268292682926829, -1.853658536585366, 5.560975609756098]] )
        x_ = orthoproject( x[0], x[1:] )
        print x_

if __name__ == "__main__" :
    unittest.main()
