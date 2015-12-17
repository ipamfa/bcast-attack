from fractions import Fraction
from numpy.matlib import dot
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
        Z = np.zeros(k)   # TODO : cek apakah rekursif berpengaruh    
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

# ref: Kannan p.
# find shortest vector in lattice L(b1, b2, ... b_n)
def enumerate() :
    return 0

# ref: Kannan p.
# main procedure of Kannan, find v1 and construct basis
def shortest(b) :        
    return 0
