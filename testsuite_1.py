import unittest
from numpy import array
# test algoritma eliminasi Gauss
from util import solve_linear

class GaussElimTest(unittest.TestCase) :    # uji kasus bukan normal
    def test_gaussel0(self) :
        # matrix dari vektor baris
        X = array([            # kasus j = k : solusi unik
            [1,2,3],
            [4,5,6]
        ])
        Y = array([
            [5,7,10]
        ])
        print "linear independen (1)"
        print "================="
        print solve_linear( (X.astype(float)).T, Y.astype(float) )

    def test_gaussel1(self) :
        # matrix dari vektor baris
        X = array([            # kasus j = k : solusi unik
            [1,2,3],
            [4,5,6],
            [5,7,9]
        ])
        Y = array([
            [5,7,10]
        ])
        print "linear independen (2)"
        print "================="
        print solve_linear( (X.astype(float)).T, Y.astype(float) )
        
    def test_gaussel2(self) :
        # matrix dari vektor baris
        X = array([            # kasus j = k : solusi unik
            [1,2,3],
            [4,5,6],
            [5,7,10]
        ])
        Y = array([
            [0,0,0]
        ])
        print "sistem homogen"
        print "=============="
        print solve_linear( (X.astype(float)).T, Y.astype(float) )   

    def test_gaussel3(self) :
        # matrix dari vektor baris
        X = array([            # kasus j = k : solusi unik
            [1,2,3],
            [4,5,6],
            [5,7,9],
            [5,7,10]
        ])
        Y = array([
            [0,0,0]
        ])
        print "kasus Kannan"
        print "============"
        print solve_linear( (X.astype(float)).T, Y.astype(float) )
            
if __name__ == "__main__" :
    unittest.main()
