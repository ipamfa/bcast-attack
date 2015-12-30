import unittest
from numpy import array
from numpy.matlib import dot, matrix
# test algoritma Kannan
from L3 import LLL
from kannan import lifting, shortest, gauss_elim

class KannanTest(unittest.TestCase) :
    # lifting from 1 dim vector to 2 dim
    def test_lifting_1to2(self) :
        v1 = array([2,3])
        v2 = array([4,5])
        r  = lifting(v1, v2)
        # nilai proyeksi nya hrus minimal
        print "=*=*=*=*=*=*=*proyeksi:=*=*=*=*=*=*=*=*"
        print float( dot(r,v1) )/dot(v1,v1)

    def test_lifting_2to3(self) :
        pass

    def test_gausselim(self) :        
        X = array([
            [1,2,3],
            [4,5,6],
            [7,8,9],
            [10,11,12]
        ])
        # eliminasi gauss sbg OPERASI BARIS T*X
        (a,b) = gauss_elim(X,True)
        print "eliminasi gauss"
        print a
        print b
        print matrix(b)*X # harusnya a

    def test_transform(self) :        
        from kannan import getTransform        
        print "=*=*=*=*=*=*=*transform=*=*=*=*=*=*=*"
        X = array([                     # kasus j < k
            [1,2,3],
            [4,5,6]
        ])
        Y = matrix([[1,2], [3,4]])*X
        t = getTransform( X,Y )
        print t

        # membalik eliminasi gauss
        X = array([                     # kasus j > k
            [1,2,3],
            [4,5,6],
            [7,8,9],
            [10,11,12]
        ])
        # tes bagaimana b didapatkan
        (a,b) = gauss_elim(X,True)
        T = getTransform(X,a)
        print T
        
if __name__ == "__main__" :
    unittest.main()
