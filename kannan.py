from fractions import Fraction
from numpy import zeros, array, append, ravel
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
        
        a = np.empty([j-1,k], dtype=object)
        if not all( c[j-1]==Z ) :           # b[0] bebas linear thd b lain
            a[0] = np.copy(b[0])
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
        b_ = np.empty( [j-1,k], )
        for i in range( 0,j-1 ) :
            b_[i] = b[i+1] - (np.dot( b[i+1],a[0] )/np.dot( a[0],a[0] ))*a[0]
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
    m = len(b)
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

# ref: Kannan p.
# main procedure of Kannan, find v1 and construct basis
def shortest(b) :
    n = len(b)
    if n==1 :
        return b
    (j,k) = b.shape
    # aproximasi reduced basis
    b_ = np.empty( [j-1,k], )
    for i in range( 0,j-1 ) :      # proyeksi perpendicular thd b1
        b_[i] = b[i+1] - (np.dot( b[i+1],b[0] )/np.dot( b[0],b[0] ))*b[0]
    # enumerate

    # select-basis
    return 0
