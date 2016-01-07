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

def findNZ1(a) :
        i = 0
        while i < len(a) and not a[i] :
                i += 1
        return i

# pencarian solusi linear adalah topik yg luas
# kita tidak mungkin merangkum nya dalam satu prosedur/fungsi
# dalam prosedur ini kita asumsi matrix disusun oleh vektor kolom,
# tetapi implementasi oleh numpy utk mensupoort operasi baris
# input: A dan B dimensi sama, p = presisi numerik
#        A mtx kolom, B mtx baris
# return: solusi sistem linear
def solve_linear(A, B, p=15):     # mencari transformasi linear A.T = B
        # kita "terpaksa" memakai eliminasi Gauss karena sifat sistem
        # bawaan dri algoritma Kannan adalah tidak full rank. Menjadikan
        # adanya kolom atau baris nol dan tidak memungkinkan penggunaan
        # faktorisasi LU yg lebih efisien
        # WARN : kompleksitas algoritma ini adalah O(n^3)
        (m,n) = A.shape
        jmax = min(m,n)
        # ukuran T haruskah jmax atau n ???
        T = np.zeros((n ,len(B)))     # asumsi konsisten
        tswp = {}
        c = 0                         # indeks kolom dri T
        while c < len(B) :
                b = B[c]
                # loop
                X = np.append( A, np.resize(b,(m,1)), axis=1 )
                # cek utk row scaling
                # cari quantile
                # mulai eliminasi Gauss dg partial pivot
                i = j = 0
                while i < m and j < n :
                        # cek jika kolom dibawah sdh nol
                        if all(map(lambda x: not x, X[i:m, j])) :
                                j += 1    # geser indeks kolom
                                continue
                        while not (X[i,j] > 0 and all(map(lambda x: not x, X[i+1:m,j]) )) :
                                # s nilai maximum dalam satu kolom j
                                # lihat partial pivot, buku lin5a p.36
                                s = max(filter(lambda x: x, map(abs, X[i:m,j])))
			        rows = range(i,m)
			        tabs = np.array([rows, abs(X[rows,j])]) .transpose()
                                # index nilai minimum
			        k = min( map(lambda r: r[0], filter(lambda r: r[1]==s, tabs)) )
                                if i != k :
                                        t = np.copy(X[i])       # swap baris
                                        X[i] = X[k]
                                        X[k] = t
                                if X[i,j] < 0 :
                                        X[i]*=-1
                                for k in range(i+1,m) :
                                        q =  X[k,j] / X[i,j]
                                        if q != 0 :
                                                X[k] -= q*X[i]
                                                X[k] = map(lambda x: round(x,p),X[k])
                        # inc kolom dan baris
                        i += 1
                        j += 1
                        
                # cek karakteristik sistem, perhatikan asumsi kontek Kannan
                # 1) kasus tidak konsisten bisa muncul (tak ada solusi)
                # 2) kasus linear dependen
                # susun ulang kolom jika bentuk pivot tidak diagonal !!!
                i = m
                while i > jmax :
                        if all( map(lambda x: not x, X[i-1,:n]) ) : # baris nol
                                if X[i-1,n] :
                                        return [ ]
                        i -= 1
                
                # mulai back substitute
                for i in range(jmax-1, -1, -1) :
                        s = 0
                        for j in range(jmax-1, i, -1) :
                                s += X[i,j]*T[j, c]
                        # baris nol
                        if all( map(lambda x: not x, X[i,:n]) ) :
                                # ruas kanan harus dicek nol juga
                                if X[i,n] :
                                        return [ ]
                                T[i,c] = 0                      
                        else :
                                if X[i,i]==0 :
                                        # swap kolom
                                        j = i
                                        while X[i,j]==0 and j < n :
                                                j += 1
                                        for it in range(0,m) :
                                                t = X[it,i]
                                                X[it,i] = X[it,j]
                                                X[it,j] = t
                                        # update daftar T yg mesti diswap
                                        tswp[c] = (i,j)
                                T[i,c] = (X[i,n] - s) / X[i,i]
                c += 1  # increment kolom T
        # ubah urutan solusi karena pertukaran kolom sebelumnya (jika ada)
        # TODO : ganti konsep table swap jadi matrix elementer
        for c in tswp :
                (a,b) = tswp[c]
                t = T[a,c]
                T[a,c] = T[b,c]
                T[b,c] = t
              
        return T        # solusi msih dalam persepsi matrix kolom
