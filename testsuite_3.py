import unittest
from numpy import array
# test algoritma eliminasi Gauss
from util import solve_linear

class GaussElimTest(unittest.TestCase) :
    
    def test_gaussel0(self) :
        # matrix dari vektor baris
        X = array([            # kasus j = k : solusi unik
            [1,2,3],
            [4,5,6],
            [7,8,9]            
        ])
        Y = array([
            [10,11,12]
        ])
        print "kasus unik normal"
        print solve_linear( (X.astype(float)).T, Y)

    def test_gaussel1(self):
        X = array([            # kasus j = k : solusi unik
            [1,2,3],
            [4,5,6]            
        ])
        Y = array([
            [7,8,9],
            [10,11,12]
        ])
        print "kasus unik not full rank"
        T = solve_linear( (X.astype(float)).T, Y) # input kolom perlu transp
        print T.T
        
    def test_gaussel2(self):
        X = array([            # kasus j = k : solusi unik
            [1,2,3],
            [7,8,9]            
        ])
        Y = array([
            [4,5,6],
            [10, 11, 12]
        ])
        print "kasus unik not full rank"
        T = solve_linear( (X.astype(float)).T, Y) # input kolom perlu transp
        print T
if __name__ == "__main__" :
    unittest.main()

