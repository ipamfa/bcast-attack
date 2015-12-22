
import numpy as np
from numpy.linalg import inv, norm
from numpy.matlib import matrix


from util import swap, gcd, lcm

B = None   # basis hasil reduksi, di buku goblok notasinya y
B_ = None  # hasil Gram-schmidt
D = None   # star x star
U = {}     # konstanta mu

def gram(b) :    # proses Gram-schmidt thd basis b
    (m,n) = b.shape
    global D, B, B_
    D = np.array( range(0,m), float )
    B = np.copy(b)
    B_ = B.astype(float,copy=True)  # konversi ke float
    for i in range(0,m) :
        for j in range(0,i) :
            #
            U[(i,j)] = np.dot( B[i],B_[j] )/D[j]
            B_[i] -= U[(i,j)]*B_[j]
            
        D[i] = np.dot( B_[i],B_[i] )
        
# ingat kata the practice of programming : Designing Interface :D
# aksesor ke b* basis orthogonal
def getBstar(i) :
    return B_[i]

def getU(i,j) :      # indeks dimulai dari 0 (nol)
    return [ (i,j) ]

# in our program, k=1..(n-1) and l=k-1, k-2..0 -> setup l < k
def reduce(k,l) :
    global U, B
    if abs( U[(k,l)] ) > 0.5 :
        ui = int( round( U[(k,l)] ) )
        B[k] -= ui*B[l]
        for j in range(0,l) :
            U[(k,j)] -= ui*U[(l,j)]
        U[(k,l)] -= ui

def exchange(k) :       # in our program, k=1..(n-1)
    global B, U, D
    # swap
    z = np.copy( B[k] )
    B[k] = B[k-1]
    B[k-1] = z

    v = U[(k,k-1)]
    delta = D[k] + v**2*D[k-1]
    U[(k,k-1)] = v*D[k-1]/delta
    D[k] *= (D[k-1]/delta)

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
        
def LLL(b,alpha=0.75) :       # main procedure
    (n,n) = b.shape           # asumsi b numpy array
    gram(b)
    k = 1
    while k < n :
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
