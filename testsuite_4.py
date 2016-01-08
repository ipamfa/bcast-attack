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
       # init_global()
       # print shortest(x)   

    def test_kannan_4dim(self) :
        x = np.array([
            [-2, 7, 7, -5],
            [3, -2, 6, -1],
            [2, -8, -9, -7],
            [8, -9, 6, -4],
        ])
        init_global()
        print shortest(x)

    def test_step_4dim(self) :
        b_ = np.array([
            [2.46632124352332, 0.25906735751295, -2.11917098445596, -3.59067357512953], [3.4786437692646395, -3.6177895200352275, 4.313518273888157, -0.41743725231175377]], dtype=object)
        b = np.array([[ 2.46632124,  0.25906736, -2.11917098, -3.59067358],
                      [ 3.24352332, -3.64248705,  4.51554404, -0.07512953]])
        T = solve_linear(b.T,b_)
        
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
      #  init_global()
      #  print shortest(x)
        
if __name__ == "__main__" :
    unittest.main()
