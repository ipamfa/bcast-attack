# 29 Sep 2015
# library utk representasi Matrix
import numpy as np

# semua menggunakan asumsi matrix baris

def swap(b,i,j) :
        t = np.copy(b[i])
        b[i]  = b[j]
        b[j] = t

#
# menghitung Hermite Normal Form dari a
#
def hnf(a) : # a is array n dimensi
	(m,n) = a.shape
	h = np.zeros((m,n))
	i = j = 0
	while i < m and j < n :
		# cek nilai utk semua baris di bawah i (dalam kolom j) adalah 0
		if all(map(lambda x: not x ,a[i:m,j])) :
			j += 1
			continue
		# kondisi normal maka lakukan operasi baris
		while not (a[i,j] > 0 and all(map(lambda x: not x, a[i+1:m,j]) )) :
			s = min(filter(lambda x: x, map(abs, a[i:m,j])))
			rows = range(i,m)
			tabs = np.array([rows, abs(a[rows,j])]) .transpose()
			k = min( map(lambda r: r[0], filter(lambda r: r[1]==s, tabs)) )
			if i != k :
				t = np.copy(a[i])# swap baris
				a[i] = a[k]
				a[k] = t
			if a[i,j] < 0 :
				a[i]*=-1
			for k in range(i+1,m) :
				q = a[k,j] / a[i,j]
				if q != 0 :
					a[k] -= q*a[i]

		# lakukan operasi ke elemen yg di baris atas (msh satu kolom j)
		for k in range(0,i) :
			q = a[k,j] / a[i,j]
			if q != 0 :
				a[k] -= q*a[i]
		i += 1
		j += 1

# Gauss elimination function on matrix
# params :
# return :
def gauss_elim(b, withq=False) :            # b matrix baris, numpy type
        a = b.astype(float, copy=True)
        i = j = 0
        (m,n) = a.shape
        if withq :
                # matrix representasi operasi baris dari b
                qs = np.identity(m)
        while i < m and j < n :
                if all(map(lambda x: not x, a[i:m, j])) :
                        j += 1
                        continue
                while not (a[i,j] > 0 and all(map(lambda x: not x, a[i+1:m,j]) )) :
                        # s nilai minimum dalam satu kolom j
                        s = min(filter(lambda x: x, map(abs, a[i:m,j])))
			rows = range(i,m)
			tabs = np.array([rows, abs(a[rows,j])]) .transpose()
                        # index nilai minimum
			k = min( map(lambda r: r[0], filter(lambda r: r[1]==s, tabs)) )
                        if i != k :
				t = np.copy(a[i])# swap baris
				a[i] = a[k]
				a[k] = t
                                # koefisien
                                if withq :
                                        swap(qs, i, k)
			if a[i,j] < 0 :
				a[i]*=-1
                                if withq :
                                        qs[i] *= -1
			for k in range(i+1,m) :
				q =  a[k,j] / a[i,j]
				if q != 0 :
					a[k] -= q*a[i]
                                        if withq :
                                                qs[k] -= q*qs[i]
                i += 1
                j += 1

        if withq :     # berpikir simpel saja
                return (a,qs)
        return a

# convert type into long
# input : a = array like 2d struture, s : index to start
def toLong(a, s=0) :
        i = s
        for r in a :
                a[i] = map(lambda x: long( round(x) ), r)
                i += 1

def gcd(a,b) :
        if a < b :
                (a,b) = (b,a)
        r = a % b
        while r :
                a = b
                b = r
                r = a % b
        return b

def lcm(a,b) :
        c = gcd(a,b)
        return c*(a/c)*(b/c)        
