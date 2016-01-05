import unittest
from L3 import *
from util import solve_linear
from kannan import shortest, orthoproject

class GramSTest(unittest.TestCase) :

    def test_lll_fullrank_3dim(self) :
        print "Test Kannan"
        print "=============="
        x = np.array([
            [4,5,1],
            [4,8,2],
            [6,2,6]
        ])
        shortest(x)

    def test_solvelin_again(self) :
        x = np.array( [
            [0.0, 0.2, -0.6]
        ])
        y = np.array([
            [0.0, 0.2, -0.6]
        ])
        print "kasus Kannan 2"
        print "=============="
        print x
        print y
        print solve_linear( x.T, y )
        
if __name__ == "__main__" :
    unittest.main()
