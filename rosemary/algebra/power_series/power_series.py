################################################################################
# Functions associated to formal power series
################################################################################

import itertools

def delta_series(n):
    eta = eta_series(n)
    delta = ps_power(eta, 24)
    return delta

def eta_series(prec = 1000):
    """
    eta_series_qexp(z, prec = 1000):
    Returns the q-expansion for eta(z) to precision p, without the q^(1/24).

    Examples:

    """
    def exps():
        n = 1
        while True:
            yield n * (3*n - 1) // 2
            yield n * (3*n + 1) // 2
            n += 1

    coeff_list = [0] * prec
    coeff_list[0] = 1
    for term, sign in itertools.izip(exps(), itertools.cycle((-1, -1, 1, 1))):
        if term >= prec:
            break
        coeff_list[term] = sign
    return coeff_list

def ps_invert(A):
    """
    power_series_invert(A):

    Given a power series p(x) = a0 + a1*x + a2*x^2 + ..., whose coefficients
    are given as the list A, this computes the coefficients of p(x)^(-1).
    """
    n = len(A)
    B = [0] * n
    a0 = A[0]
    B[0] = 1 / a0
    for k in xrange(1, n):
        ss = 0
        for i in xrange(k + 1):
            ss += A[i] * B[k - i]
        B[k] = -ss / a0
    return B

def ps_power(A, k):
    """
    power_series_power(A, k):
    Given a power series p(x) = a0 + a1*x + a2*x^2 + ..., whose coefficients
    are given as the list A, this computes the coefficients of p(x)^k using
    the J.C.P pure power recurrence.
    """
    n = len(A)
    P = [0] * n
    a0 = A[0]
    P[0] = a0**k
    for i in xrange(1, n):
        ss = 0
        for j in xrange(1, i + 1):
            ss += A[j] * ((k + 1) * j - i) * P[i - j]
        P[i] = ss // (i * a0)
    return P

def ps_mult_classical(M, N):
    m = len(M)
    n = len(N)
    P = [0] * (m + n - 1)
    for i in xrange(m):
        for j in xrange(n):
            P[i + j] += M[i] * N[j]
    return P

def ps_mult_karatsuba(M, N, kthresh=40):
    m = len(M)
    n = len(N)

    # Classical Multiplication is faster until a certain point.
    if m < kthresh or n < kthresh:
        return ps_mult_classical(M, N)

    # The algorithm works best when we split into halves
    k = max(m, n) // 2

    # write M = f1*x^k + f0, N = g1*x^k + g0
    f0 = M[0:k]
    f1 = M[k:]
    g0 = N[0:k]
    g1 = N[k:]

    # We want to compute MN = (f1*x^k + f0)(g1*x^k + g0).
    # We let z0 = f0*g0, z2 = f1*g1
    z0 = ps_mult_karatsuba(f0, g0, kthresh)
    z2 = ps_mult_karatsuba(f1, g1, kthresh)

    f1_f0 = [sum(e) for e in itertools.izip_longest(f1, f0, fillvalue = 0)]
    g1_g0 = [sum(e) for e in itertools.izip_longest(g1, g0, fillvalue = 0)]

    xxyy = ps_mult_karatsuba(f1_f0, g1_g0, kthresh)

    z1 = [e[0] - e[1] - e[2] for e in itertools.izip_longest(xxyy, z2, z0,
        fillvalue = 0)]

    t2 = [0] * (2*k) + z2
    t1 = [0] * (k) + z1
    xy = [sum(e) for e in itertools.izip_longest(t2, t1, z0, fillvalue = 0)]
    return xy

