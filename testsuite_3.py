import unittest
from numpy import array
from numpy.matlib import dot, matrix
# test algoritma Kannan
from kannan import lifting, shortest, gauss_elim
from kannan import getTransform        

class GaussElimTest(unittest.TestCase) :
    
    def test_gaussel1(self) :        
        
        print "=*=*=*=*=*=*=*tes 1=*=*=*=*=*=*=*"
        X = array([                     # kasus j < k
            [1,2,3],
            [4,5,6]
        ])
        Y = matrix([[1,2], [3,4]])*X
        t = getTransform( X,Y )
        print t

        # membalik eliminasi gauss
        X = array([               # kasus linear dependen
            [1,2,3],
            [4,5,6],
            [7,8,9],
            [10,11,12]
        ])
        # tes bagaimana b didapatkan
        (a,b) = gauss_elim(X,True)
        T = getTransform(X,a)
        print T

    def test_gaussel0(self) :
        # cek
        X = array([            # kasus j > k : linear dependen
            [1,2,3],
            [4,5,6],
            [7,8,9],
            [10,11,12]
        ])
        Y = array([
            [1,2,3],
            [4,5,6],
            [7,8,9],
            [0,0,0]
        ])
        print "tes 0"
        print getTransform(X,Y)

    def test_gaussel2(self):
        X = array([            # kasus linear dependen v3 = 3v2-2v1
            [1,2,3],
            [4,5,6],            
            [10,11,12]
        ])
        print "tes 2"
        Y = array([
            [1,2,3],
            [4,5,6], 
            [0,0,0]
        ])
        print getTransform(X,Y)

    def test_transform0(self) :
        print "tes transform 0"
        X = array([
            [1,2,3],
            [4,5,6],
            [7,8,9]
        ])
        t = array([
            [1,2,3],
            [4,5,6],
            [7,8,9]
        ])
        print gauss_elim(X.T)
        Y = matrix(X)*matrix(t)
        print getTransform(X,Y)
        
if __name__ == "__main__" :
    unittest.main()

