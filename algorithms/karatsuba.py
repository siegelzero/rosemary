import operator
import itertools

def convolve(M, N):
    m = len(M)
    n = len(N)
    P = [0] * (m + n - 1)
    for i in xrange(m):
        for j in xrange(n):
            P[i + j] += M[i] * N[j]
    return P

def karatsuba(M, N, kthresh=50):
    m = len(M)
    n = len(N)

    if m < kthresh or n < kthresh:
        return convolve(M, N)

    #final degree
    deg = m + n - 1
    k = max(m, n) // 2

    #
    x0 = M[0:k]
    x1 = M[k:]

    y0 = N[0:k]
    y1 = N[k:]

    z0 = karatsuba(x0, y0, kthresh)
    z2 = karatsuba(x1, y1, kthresh)

    L = itertools.izip_longest(x1, x0, fillvalue = 0) 
    x1_x0 = [ sum(e) for e in L ]
    
    L = itertools.izip_longest(y1, y0, fillvalue = 0) 
    y1_y0 = [ sum(e) for e in L ]

    xxyy = karatsuba(x1_x0, y1_y0, kthresh)

    L = itertools.izip_longest(xxyy, z2, z0, fillvalue = 0) 
    z1 = [ e[0] - e[1] - e[2] for e in L ]


    t2 = [0] * (2*k) + z2
    t1 = [0] * (k) + z1

    L = itertools.izip_longest(t2, t1, z0, fillvalue = 0) 
    xy = [ sum(e) for e in L ]

    return xy
