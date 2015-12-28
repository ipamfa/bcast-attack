import unittest
from numpy import array
from numpy.matlib import dot
# test algoritma Kannan
from L3 import LLL
from kannan import lifting, shortest

class KannanTest(unittest.TestCase) :
    # lifting from 1 dim vector to 2 dim
    def test_lifting_1to2(self) :
        v1 = array([2,3])
        v2 = array([4,5])
        r  = lifting(v1, v2)
        # nilai proyeksi nya hrus minimal
        print float( dot(r,v1) )/dot(v1,v1)

    def test_lifting_2to3(self) :
        pass

    def test_kannan1(self) :
        x = array([
            [4,5,1],
            [4,8,2],
            [6,2,6]
        ])
        b = shortest(x)
        
if __name__ == "__main__" :
    unittest.main()
