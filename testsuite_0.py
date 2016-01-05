import unittest
from numpy import array
# test algoritma eliminasi Gauss
from util import solve_linear

class GaussElimTest(unittest.TestCase) :

    def test_gaussel0(self) :
        X = array( [[1,2,3]] )
        Y = array( [[1,2,3]] )
        print "kasus paling dasar"
        print "=================="
        T = solve_linear( (X.astype(float)).T, Y.astype(float) ) # input kolom perlu transp
        print T
        
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
        T = solve_linear( (X.astype(float)).T, Y.astype(float) ) # input kolom perlu transp
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
        T = solve_linear( (X.astype(float)).T, Y.astype(float) ) # input kolom perlu transp
        print T

    def test_gaussel3(self) :
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
        print solve_linear( (X.astype(float)).T, Y.astype(float) )

    def test_gaussel4(self) :
        X = array([            # kasus j = k : solusi unik
            [0,-2,-2],
            [1,4,5],
            [-1,-1,-4]
        ])
        Y = array([
            [3, 1, -2]
        ])
        print "kasus unik full rank"
        T = solve_linear( (X.astype(float)).T, Y.astype(float) ) # input kolom perlu transp
        print T

    def test_gaussel5(self) :
        X = array([            # kasus j = k : solusi unik
            [1,1,1],
            [1,2,2],
            [1,2,3]
        ])
        Y = array([
            [1, 1, 1]
        ])
        print "kasus unik full rank"
        T = solve_linear( (X.astype(float)).T, Y.astype(float) ) # input kolom perlu transp
        print T
        
    def test_gaussel6(self) :
        X = array([            # kasus j = k : solusi unik
            [2,-1,0],
            [-1,2,-1],
            [0,-1,1]
        ])
        Y = array([
            [0, 0, 1]
        ])
        print "kasus unik full rank"
        T = solve_linear( (X.astype(float)).T, Y.astype(float) ) # input kolom perlu transp
        print T

    def test_gaussel7(self) :
        X = array([            # kasus j = k : solusi unik
            [0,-1,-1],
            [4,7,8],
            [-3,-5,-6]
        ])
        Y = array([
            [3, 4, 5]
        ])
        print "kasus unik full rank"
        T = solve_linear( (X.astype(float)).T, Y.astype(float) ) # input kolom perlu transp
        print T
 
        
if __name__ == "__main__" :
    unittest.main()

