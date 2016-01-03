import unittest
from L3 import *
from kannan import shortest, enumerate, select_basis

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
        print "Test 3 dimensi"
        print "=============="
        x = np.array([
            [4,5,1],
            [4,8,2],
            [6,2,6]
        ])
        b = LLL(x)
        print b
        print "eval norm |b1|=", norm(b[0]), " |b2|=", norm(b[1]), " |b3|=", norm(b[2])
        shortest(x)
        
    def test_lll_fullrank_4dim(self) :
        # apply reduction algorithm
        x = np.array([
            [-2, 7, 7, -5],
            [3, -2, 6, -1],
            [2, -8, -9, -7],
            [8, -9, 6, -4],
        ])
        b = LLL(x, 1.0)
        print b
        print "eval norm |b1|=", norm(b[0]), " |b2|=", norm(b[1]), " |b3|=", norm(b[2]), " |b4|=", norm(b[3])

    def test_lll_8dim(self) :
        x = np.array([
            [ 8, -3, -3, -9, 1, 9, -3, -9],
            [ -7, 5, 1, 1, 9, -3, -4, -2],
            [ -5, 2, 1, -3, -4, 5, 5, 4],
            [ -4, 9, -6, -5, -7, 2, -1, 5],
            [ -5, 0, 2, 2, 0, 5, 6, -5],
            [ -8, -2, 3, 5, -1, 7, 7, 4],
            [ 3, -9, 3, -7, 3, 2, -3, 2],
            [ -4, -2, -8, 6, 0, 4, -9, 7]
        ])
        b = LLL(x)
        print b
        print "eval norm |b1|=", norm(b[0]), " |b2|=", norm(b[1]), " |b3|=", norm(b[2]), " |b4|=", norm(b[3]), "|b5|", norm(b[4]), " |b6|=", norm(b[5]), " |b7|=", norm(b[6]), " |b8|=", norm(b[7])

if __name__ == "__main__" :
    unittest.main()
