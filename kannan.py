from fractions import Fraction
import numpy as np
from numpy.linalg import norm
from numpy.matlib import dot, matrix
from util import solve_linear, lcm, gcd, findNZ1

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

def orthoproject(b) :             # proyeksi orthogonal thd b0
    (j,k) = b.shape
    b_ = np.empty( [j-1,k], )
    for i in range( 0,j-1 ) :
        b_[i] = b[i+1] - (float( dot( b[i+1],b[0] ))/dot( b[0],b[0] ))*b[0]
    
    return b_

# ref: Kannan p.
def select_basis(b) :     # b is numpy array    
    (j,k) = b.shape       # j = k+1
    global Z
    if Z is None :
        Z = np.zeros(k)      # TODO : cek apakah rekursif berpengaruh    
    if all( b[0]==Z ) :
        if len(b) > 1 :
            basis = b[1:]
        else :
            basis = []
    else :
        # tes apakah b[0]= kombinasi linear dari b[1]..b[j-1]
        # menggunakan eliminasi Gauss
        c = solve_linear( (b[1:j].astype(float)).T , b[:1] )        
        
        a = np.empty([j-1,k], dtype=object)
        if c==[] :                      # b[0] bebas linear thd b lain
            a[0] = np.copy(b[0])
        else :
            xs = c.T[0]                 # transpos c msh dalam represnt kolom
            #
            # dont you look like an idiot who write codes ? even you yourself
            # cant grasp and understand a single piece of code you write by
            # your own hand a couple weeks ago and now you perplex by your own
            # thoughts
            
            # khusus bagian ini lihat buku Bremner p.185
            # kita tidak perlu y dan z karena sdh ditangani oleh Fraction
            # tidak perlu menangani x_i 0
            i = findNZ1( xs )
            f = Fraction( str(xs[i]) )
            m = f.denominator
            i2 = i+1
            while i2 <  len(xs) :
                if xs[i2] :
                    f = Fraction( str(xs[i2]) )
                    m = lcm( f.denominator, m )
                i2 += 1
                
            d = long( round( m*xs[i] ) )
            i += 1
            while i <  len(xs) :
                if xs[i] :
                    x = long( round(m*xs[i]) )
                    d = gcd( x,d )
                i += 1
            
            # python perlu denominator tetap positif
            if d < 0 :
                m *= -1
                d *= -1
            # f adalah bentuk rasional m/d yg jika relative prima = p/q
            # yg berarti faktor bersama nya dikeluarkan dan diambil q nya
            f = Fraction( '%d/%d' % (m,d) )
            if m < 0 :
                a[0] = ( -1.0/f.denominator )*b[0]
            else :
                a[0] = ( 1.0/f.denominator )*b[0]

        # proyeksi b ke a[0]. karena b[0] sdh diproses, maka
        # b yg diproyeksi berikutnya adalah 1,2..j-1
        b_ = orthoproject( np.append([ a[0]], b[1:], axis=0 ) )
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

from L3 import Gram

alpha = []    # list of array of possible integer combination

# ref: Kannan p.
# input : k = basis ke, m = jml basis (indexing diawali dari 1/ aturan umum), g = Gram schmidt data
# return list of integer
def List(k,m,g) :
    if k == 0 :
        # list telah komplit, tambahkan ke alpha
        global alpha
        alpha = np.append( alpha, [ np.copy(alpha[0]) ], axis=0 )
        return
    # 
    # hitung beta0_k proposisi 2.13 ref Kannan p.
    x = norm( g.getBstar(0) )/norm( g.getBstar(k-1) )
    t = 0
    for j in range( k+1, m+1 ) :     # lihat catatan 18/12
        t += g.getU( j-1,k-1 )       # lihat spek interface
    beta0_k = long( round(- x - t) )
    d = long( round(2*x) )
    # tes utk tiap kombinasi nilai alpha
    global alpha
    for i in range( beta0_k , beta0_k + d + 1 ) :
        alpha[0, k-1] = i
        List(k-1,m,g) # rekursif

# ref: Kannan p.
# find shortest vector in lattice L(b1, b2, ... b_m)
def enumerate(b) :
    (m,n) = b.shape
    if m == 1 :
        return b[0]
    global alpha
    alpha = np.array( [ np.zeros(m) ] , dtype=object)
    # hitung bm(m) = b*_m
    
    g = Gram(b)
    lim = norm(b[0])/norm( g.getBstar( m-1 ) )
    for j in range( long( round(-lim )), long( round(lim )) + 1 ) :
        alpha[0, m-1] = j
        List(m-1,m,g)
    # alpha sdh selesai, tinggal di iterasi cari SVP :)
    B = matrix(b)
    svp = np.ravel( alpha[1]*B )
    sn = norm(svp)
    k = len(alpha)
    for i in range(2, k) :
        x = np.ravel( alpha[i]*B )
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
    b = LLL(B).getBasis()              # representasi basis dalam array 2d

    while True :
        # aproximasi reduced basis
        b_ = orthoproject(b)
    
        # rekursif
        b__ = shortest(b_)
        # lifting dg transformasi linear, catatan 13/12       
        T =  solve_linear( b_.T,b__ )  # b_ perlu ditransform, liat spek
        
        # transpose T karena msh kolom 
        v = matrix( T.T )*matrix( b[1:] )            # ubah ke matrix dulu
        b[1:] = np.asarray(v)                        # override array lagi
        b0n = norm(b[0])
        b1n = norm(b[1])
        if b1n**2 > (3.0/4)*b0n**2 :
            break                   
        # swap
        t = np.copy(b[0])
        b[0] = b[1]
        b[1] = t
        
    # temukan nilai batas j0
    g = Gram(b)        
    # cari vektor yg norm-nya lebih besar dari b0
    i = 1
    while i < j :
        v = g.getBstar(i)
        if norm(v) >= b0n :
            break
        i += 1
    # enumerate   
    v1 = enumerate(b[:i])
    # select-basis
    B = select_basis( np.append([v1], b, axis=0) )
    # ulangi lagi step spt yg didalam loop
    B_ = orthoproject(B)
    B__ = shortest(B_)                               # rekursif
    T = solve_linear( B_.T, B__ )
    
    V = matrix( T.T )*matrix( B[1:] )                # perlu bentuk matrix
    B[1:] = np.asarray(V)
    return B                              
