import unittest
from L3 import *
from util import solve_linear
from kannan import shortest, init_global
from numpy.linalg import det

class KannnanTest(unittest.TestCase) :

    def test_kannan_fullrank_3dim(self) :
        print "Test Kannan"
        print "=============="
        x = np.array([
            [4,5,1],
            [4,8,2],
            [6,2,6]
        ])
        init_global()
        print shortest(x)   

    def test_kannan_4dim(self) :
        x = np.array([
            [-2, 7, 7, -5],
            [3, -2, 6, -1],
            [2, -8, -9, -7],
            [8, -9, 6, -4],
        ])
        init_global()
        print shortest(x)

    def test_kannan_8dim(self) :
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
        init_global()
        print shortest(x)
        
if __name__ == "__main__" :
    unittest.main()
