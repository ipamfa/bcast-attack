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
        
        # basis orthogonal Gram-schmidt, y*
        self.B_ = b.astype(float,copy=True)  # konversi ke float
        # konstanta mu
        self.U = {}
        for i in range(0,m) : # proses Gram-schmidt thd basis b
            for j in range(0,i) :
                self.U[(i,j)] = np.dot( b[i], self.B_[j] )/self.D[j]
                self.B_[i] -= self.U[(i,j)]*self.B_[j]
            self.D[i] = np.dot( self.B_[i], self.B_[i] )

    # param  = i, j : index integer 
    # return = konstanta mu
    def getU(self, i, j) :           # aksesor ke konstanta mu
        return self.U[ (i,j) ]

    # param = i, j : index integer
    # update konstanta mu
    def setU(self, i, j, val) :      # mutator mu
        self.U[(i,j)] = val

    def swapU(self, k, j, k_, j_ ) :
        t = self.U[(k,j)]
        self.U[(k,j)] = self.U[(k_,j_)]
        self.U[(k_,j_)] = t      
      
    # param = i : index integer
    # return B* (basis orthogonal)
    def getBstar(self,i) :         # aksesor ke basis orthogonal
        return self.B_[i]

    def getD(self, k) :            # aksesor
        return self.D[k]

    def setD(self, k, v) :
        self.D[k] = v

        
# enkapsulasi data (pembungkusan desain)
class LLL :
    def __init__(self, b, alpha=0.75) :
        # basis hasil reduksi, di buku Bremner notasinya y
        self.B = np.copy(b)
        self.gram_data = Gram(b)        # enkapsulasi data gram
        (m,n) = b.shape
        k = 1
        while k < m :
            # OOP desain dasar
            self.reduce(k,k-1)
            u = self.gram_data.getU(k,k-1)
            if self.gram_data.getD(k) >= (alpha - u**2)*self.gram_data.getD(k-1) :
                for l in range(k-2, -1, -1):
                    self.reduce(k,l)
                k += 1
            else :
                self.exchange(k)
                if k > 1 :
                    k -= 1

    # sharing dg algoritma Kannan
    def getGram(self) :              # aksesor
        return self.gram_data
    
    # param = i : index, val = vektor
    # update vektor pada basis
    def setB(self, i, val) :         # mutator B
        self.B[i] = val

    # param = i : index integer
    # return individual vektor basis
    def getB(self, i) :              # aksesor B
        return self.B[i]

    # param = i, j : index integer
    # tukar urutan basis
    def swapB(self, i, j) :          # penyingkat
        z = np.copy( self.B[i] )
        self.B[i] = self.B[j]
        self.B[j] = z

    # dibutuhkan oleh algoritma Kannan
    def getBasis(self) :             # aksesor B entire
        return self.B
    
    # param = k, l : index integer
    def reduce(self, k, l) :
        u = self.gram_data.getU(k,l)
        if abs( u ) > 0.5 :
            ui = int( round(u) )
            v = self.getB(k)
            v -= ui* self.getB(l)
            self.setB( k,v )        # update B
            for j in range(0,l) :
                # update U
                v = self.gram_data.getU(k,j)
                v -= ui*self.gram_data.getU(l,j)
                self.gram_data.setU(k,j,v)

            v = self.gram_data.getU(k,l)
            v -= ui
            self.gram_data.setU(k,l,v)

    def exchange(self, k) :
        # tukar b_k dan b_k-1
        self.swapB(k,k-1)
        # hitung nilai
        v = self.gram_data.getU(k,k-1)
        dk = self.gram_data.getD(k)
        dk_ = self.gram_data.getD(k-1)
        delta = dk + v**2*dk_
        # update
        self.gram_data.setU( k,k-1,v*dk_/delta )
        self.gram_data.setD( k, dk*(dk_/delta))
        self.gram_data.setD( k-1, delta)
        
        # tukar u_k-1 dan u_k
        for j in range(0,k-1) : # run only if k>1 s/d k-2
            self.gram_data.swapU(k-1,j,k,j)
        # update U
        (n,n) = self.getBasis().shape
        for i in range(k+1, n) :
            e = self.gram_data.getU( i,k )
            u_ = self.gram_data.getU( i,k-1 )
            u = u_ - v*self.gram_data.getU( i,k )
            t = u*self.gram_data.getU( k,k-1 ) + e
            self.gram_data.setU( i,k,u )
            self.gram_data.setU( i,k-1,t )
        
    def pdebug( i, k, l=None) :   # print for debug    
        if l is None :        
            print "DEBUG: iterasi", i, "exchange", "k=",k
        else :                
            print "DEBUG: iterasi", i, "reduce", "k=",k, "l=",l, "ukl=", 

# check whether basis reduced            
def is_reduced(b) :              # reduced in term of Lenstra
    return 0
