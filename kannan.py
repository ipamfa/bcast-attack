from fractions import Fraction
from numpy import empty, copy, zeros, array, append, ravel, identity
from numpy.linalg import norm
from numpy.matlib import dot, matrix
from util import gauss_elim

# asumsi desain tanpa thread dan lock, gunakan glob variable
Z = None                                # vektor zero (nol)

def lifting(a,b) :                      # proyeksi thd a
    x = float( dot(a, b) )/dot(a, a)
    while True :        
        b = b - a
        y = float( dot(a, b) )/dot(a, a)
        if abs(x) < abs(y) :
            break
        x = y
    # undo subtraction of b
    return b + a

# solve sistem linier T.a = b dengan mencari nilai T (Kannan p.)
# input : a, b matrix baris asumsi keduanya konsisten
def getTransform(a,b) :
    # cek jika b adalah bentuk eliminasi gauss    
    (a_, t) = gauss_elim(a.T, True)
    (j,k) = a_.shape     
    T = zeros((k,k))        # transformasi selalu square
    a = a_
    b = matrix(t)*b.T    
    # hitung batas y    
    Ymax = min(j-1, k-1)
    # loop per baris T
    for x in range(0,k) :
        for y in range(Ymax, -1, -1) :
            if a[y,y]== 0 :
                print "DEBUG: error!!"
            s = 0
            for z in range(Ymax, y, -1) :
                # jml s
                s += a[y,z]*T[z,x]
            T[y,x] = ( b[y,x] - s )/a[y,y]
    return T.T

def orthoproject(b) :             # proyeksi orthogonal thd b0
    (j,k) = b.shape
    b_ = empty( [j-1,k], )
    for i in range( 0,j-1 ) :
        b_[i] = b[i+1] - (float( dot( b[i+1],b[0] ))/dot( b[0],b[0] ))*b[0]
    
    return b_

# ref: Kannan p.
def select_basis(b) :     # b is numpy array    
    (j,k) = b.shape       # j = k+1
    global Z
    if Z is None :
        Z = zeros(k)      # TODO : cek apakah rekursif berpengaruh    
    if all( b[0]==Z ) :
        if len(b) > 1 :
            basis = b[1:]
        else :
            basis = []
    else :
        # tes apakah b[0]= kombinasi linear dari b[1]..b[j-1]
        # menggunakan eliminasi Gauss
        (c,xs) = gauss_elim(b, True)    
        # TODO : numerical precission issue
        c[j-1] = map(lambda x: round(x,14), c[j-1])
        
        a = empty([j-1,k], dtype=object)
        if not all( c[j-1]==Z ) :          # b[0] bebas linear thd b lain
            a[0] = copy(b[0])
        else :
            # xs masih matrix, hanya perlu baris terakhir
            xs = map(lambda x: (-1/xs[j-1,0])*x, xs[j-1, 1:])
            f = Fraction(str(xs[0]))
            m = f.denominator
            for x in xs[1:] :
                f = Fraction(str(x))
                m = lcm(f.denominator, m)
            d = long( round( m*xs[0] ) )
            for x in xs[1:] :
                x = long( round(m*x) )
                d = gcd(x, d)
            if d < 0 :
                m *= -1
                d *= -1
            f = Fraction( '%d/%d' % (m,d) )
            if m < 0 :
                a[0] = ( -1.0/f.denominator )*b[0]
            else :
                a[0] = ( 1.0/f.denominator )*b[0]

        # proyeksi b ke a[0]. karena b[0] sdh diproses, maka
        # b yg diproyeksi berikutnya adalah 1,2..j-1
        b_ = orthoproject( append([ a[0]], b[1:], axis=0 ) )
        # spek dr Kannan menyatakan set b adalah linear dependent
        # TODO : cek nilai 0
        if len(b_)==1 :
            b_[0] =  map(lambda x: round(x,14), b_[0])  # harusnya close to 0
        # rekursif
        c = select_basis(b_)          
        # hasil dari rekursif digunakan utk lifting
        if len(c) > 0 :
            # unique minimal lifting dari c ke a
            j = len(c)
            for i in range(0, j):
                # lihat buku catatan 17/12
                a[i+1] = lifting(a[0], b[i+1] + c[i] - b_[i] )
        
        basis = a
    return basis

from L3 import gram, getBstar, getU

alpha = []    # list of array of possible integer combination

# ref: Kannan p.
# input : k = basis ke, m = jml basis (indexing diawali dari 1/ aturan umum)
# return list of integer
def List(k,m) :
    if k == 0 :
        # list telah komplit, tambahkan ke alpha
        global alpha
        alpha = append( alpha, [ copy(alpha[0]) ], axis=0 )
    # 
    # hitung beta0_k proposisi 2.13 ref Kannan p.
    x = norm( getBstar(0) )/norm( getBstar(k-1) )
    t = 0
    for j in range( k+1, m+1 ) :     # lihat catatan 18/12
        t += getU( j-1,k-1 )         # lihat spek interface
    beta0_k = long( round(- x - t) )
    d = long( round(2*x) )
    # tes utk tiap kombinasi nilai alpha
    global alpha
    for i in range( beta0_k , beta0_k + d + 1 ) :
        alpha[0, k-1] = i
        List(k-1,m) # rekursif

# ref: Kannan p.
# find shortest vector in lattice L(b1, b2, ... b_m)
def enumerate(b) :
    (m,n) = b.shape
    if m == 1 :
        return b
    global alpha
    alpha = array( [ zeros(m) ] , dtype=object)
    # hitung bm(m) = b*_m
    gram(b)    
    lim = norm(b[0])/norm( getBstar( m-1 ) )
    for j in range( long( round(-lim )), long( round(lim )) + 1 ) :
        alpha[0, m-1] = j
        List(m-1,m)
    # alpha sdh selesai, tinggal di iterasi cari SVP :)
    B = matrix(b)
    svp = ravel( alpha[1]*B )
    sn = norm(svp)
    k = len(alpha)
    for i in range(2, k) :
        x = ravel( alpha[i]*B )
        xn = norm(x)
        if sn > xn :
            sn = xn
            svp = x
    return svp

from math import sqrt
from numpy.linalg import inv
from util import toLong
from L3 import LLL

# ref: Kannan p.
# main procedure of Kannan, find v1 and construct basis
def shortest(B) :    
    (j,k) = B.shape
    if j==1 :
        return B
    
    b = LLL(B[n:])    
    # aproximasi reduced basis
    b_ = orthoproject(b)
    
    # rekursif
    b__ = shortest(x1,)
    
    T = 0 # 
    toLong(T)
    # lifting dg transformasi linear, catatan 13/12
    v = T*b[1:]               # karena matrix baris, T di sebelah kiri
    b[1:] = numpy.asarray(v)  # override
    b0n = norm(b[0])
    b1n = norm(b[1])
    if b1n > b0n*sqrt(3)/2 :
        t = copy(b[0]) # swap
        b[0] = b[1]
        b[1] = t
        b0n = b1n
    # temukan nilai batas j0
    gram(b)
    # cari vektor yg norm-nya lebih besar dari b0
    i = 1
    while i < j :
        v = getBstar(i)
        if norm(v) >= b0n :
            break
        i += 1
    # enumerate
    v1 = enumerate(b[:i-1])
    # select-basis
    B = select_basis( append([v1], b, axis=0) )
    # ulangi lagi step
    B_ = orthoproject(B)
    B__ = shortest(B_,n-1)                # rekursif
    T = 0
    toLong(T)
    V = T*B[1:]
    B[1:] = numpy.asarray(V)
    return B                              
