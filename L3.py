import numpy as np
from numpy.linalg import inv, norm
from numpy.matlib import matrix

from util import swap, gcd, lcm
#
# ingat kata The Practice of Programming : Designing Interface :D
#
class Gram :
    def __init__(self,b) :       # b adalah numpy array object
        (m,n) = b.shape
        # star x star
        self.D = np.array( range(0,m), float )
        # basis hasil reduksi, di buku Bremner notasinya y
        self.B = np.copy(b)
        # basis orthogonal Gram-schmidt, y*
        self.B_ = B.astype(float,copy=True)  # konversi ke float
        # konstanta mu
        self.U = {}
        for i in range(0,m) :                # proses Gram-schmidt thd basis b
            for j in range(0,i) :
                self.U[(i,j)] = np.dot( self.B[i], self.B_[j] )/self.D[j]
                self.B_[i] -= self.U[(i,j)]*self.B_[j]
            self.D[i] = np.dot( self.B_[i], self.B_[i] )

    # param  = i, j : index integer 
    # return = konstanta mu
    def getU(self, i, j) :                   # aksesor ke konstanta mu
        return self.U[ (i,j) ]

    def getBstar(self,i) :                   # aksesor ke basis orthogonal
        return self.B_[i]
        

# in our program, k=1..(n-1) and l=k-1, k-2..0 -> setup l < k
def reduce(k,l) :
    global U, B
    if abs( U[(k,l)] ) > 0.5 :
        ui = int( round( U[(k,l)] ) )
        B[k] -= ui*B[l]
        for j in range(0,l) :
            U[(k,j)] -= ui*U[(l,j)]
        U[(k,l)] -= ui

# mengubah urutan basis dan hitung ulang nilai gamma*
# nilai mu, dll
def exchange(k) :       # in our program, k=1..(n-1)
    global B, U, D
    # swap b_k-1 dan b_k
    z = np.copy( B[k] )
    B[k] = B[k-1]
    B[k-1] = z

    v = U[(k,k-1)]
    delta = D[k] + v**2*D[k-1]
    U[(k,k-1)] = v*D[k-1]/delta
    D[k] *= (D[k-1]/delta)
    D[k-1] = delta

    # exchange U_k-1 dan U_k
    for j in range(0, k-1) :  # only run if k > 1 s/d k-2
        t = U[(k-1,j)]
        U[(k-1,j)] = U[(k,j)]
        U[(k,j)] = t
    # update U
    (n,n) = B.shape           # asumsi full-rank
    for i in range(k+1, n) :
        e = U[(i,k)]
        U[(i,k)] = U[(i,k-1)] - v*U[(i,k)]
        U[(i,k-1)] = U[(k,k-1)]*U[(i,k)] + e

def pdebug( i, k, l=None) :   # print for debug
    global U
    if l is None :        
        print "DEBUG: iterasi", i, "exchange", "k=",k
    else :                
        print "DEBUG: iterasi", i, "reduce", "k=",k, "l=",l, "ukl=", U[(k,l)]
    
def LLL(b,alpha=0.75) :       # main procedure
    global B, U, D
    (m,n) = b.shape           # asumsi b numpy array
    gram(b)
    k = 1
    while k < m :        
        reduce(k,k-1)
        if D[k] >= (alpha - U[(k,k-1)]**2)*D[k-1] :
            for l in range(k-2, -1, -1) :                                
                reduce(k,l)
                
            k += 1
        else :            
            exchange(k)
            if k > 1 :
                k -= 1
    # hasil akhir reduksi adalah basis B
    return B

# TODO : pindah is_reduced ke sini
