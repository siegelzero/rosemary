# Arithmetic Functions

from rosemary.number_theory.prime_list import PRIME_LIST

import rosemary.number_theory.core
import rosemary.number_theory.factorization
import rosemary.number_theory.primality
import rosemary.number_theory.sieves

########################################################################################################################
# classical arithmetic functions
########################################################################################################################

def euler_phi(n):
    """
    euler_phi(n):
    Returns euler's totient function of n; i.e., returns the number of positive
    integers <= n that are coprime to n.

    Examples:
    >>> euler_phi(100)
    40
    >>> euler_phi(11213)
    11212
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n == 0:
            raise ValueError("Zero argument in an arithmetic function")
        n_fac = rosemary.number_theory.factorization.factor(n)
    elif isinstance(n, list):
        n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    pp = 1
    for (p, e) in n_fac:
        pp *= p**(e - 1)*(p - 1)
    return pp

def factorial(n):
    """
    factorial(n):
    Given a non-negative integer n, this returns n! = n*(n - 1)*(n - 2)...2*1
    The algorithm used is a binary splitting method.

    Examples:
    >>> factorial(10)
    3628800
    >>> factorial(40)
    815915283247897734345611269596115894272000000000
    """
    D = {0:1, 1:1, 2:2, 3:6, 4:24}
    if n < 0:
        raise ValueError
    elif n <= 4:
        return D[n]

    # bn is the number of binary digits of n
    pp = 1
    k = 1

    for k in xrange(1, n.bit_length()):
        lower = n >> k
        upper = n >> (k - 1)

        j = lower + 1
        if j % 2 == 0:
            j += 1

        t_pp = 1
        while j <= upper:
            t_pp *= j
            j += 2

        pp *= (t_pp**k)

    bc = rosemary.number_theory.core.bit_count(n)

    return 2**(n - bc) * pp

def moebius(n):
    """
    moebius(n):
    Given a positive integer n, this returns the value of the Moebius function
    mu of n, the multiplicative function satisfying:
        mu(1) = 1
        mu(p) = -1
        mu(p**e) = 0 for e > 1
        mu(p1*p2*...*pk) = (-1)**k

    Examples:
    >>> moebius(35)
    1
    >>> moebius(100)
    0
    >>> moebius(2)
    -1
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n == 0:
            raise ValueError("Zero argument in an arithmetic function")
        n_fac = rosemary.number_theory.factorization.factor(abs(n))

    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    if rosemary.number_theory.factorization.is_squarefree(n_fac):
        v = len(n_fac)
        return (-1)**v
    else:
        return 0

def primorial(n):
    """
    primorial(n):
    This returns the product of the first n primes p1, p2, ..., pn.
    """
    pp = 1
    N = len(PRIME_LIST)
    if n >= N:
        raise ValueError('n too large')
    for i in xrange(n):
        pp *= PRIME_LIST[i]
    return pp

def sigma(n, k = 1):
    """
    sigma(n, k = 1):
    Given integers n and k, this returns the sum of th k-th powers of the
    divisors of n.

    Examples:
    >>> sigma(9)
    13
    >>> sigma(10, 2)
    130
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n == 0:
            raise ValueError("Zero argument in an arithmetic function")
        n_fac = rosemary.number_theory.factorization.factor(abs(n))

    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    pp = 1
    for (p, e) in n_fac:
        pk = p**k
        pp *= (pk**(e + 1) - 1) // (pk - 1)

    return pp

def tau(n):
    """
    tau(n):
    Given an integer n, this returns the number of divisors of n.
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n == 0:
            raise ValueError("Zero argument in an arithmetic function")
        n_fac = rosemary.number_theory.factorization.factor(abs(n))

    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    pp = 1
    for (p, e) in n_fac:
        pp *= (e + 1)

    return pp

def carmichael_lambda(n):
    """
    carmichael_lambda(n):
    Returns the exponent of the multiplicative group of residues modulo n; i.e.
    the smallest positive integer m such that a^m = 1 (mod n) for all integers 
    a coprime to n.

    Examples:
    >>> carmichael_lambda(100)
    20
    >>> carmichael_lambda(113)
    112
    """
    if isinstance(n, (int, long)):
        n_fac = rosemary.number_theory.factorization.factor(abs(n))

    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    L = []
    for (p, e) in n_fac:
        if p != 2:
            L += [ p**(e - 1) * (p - 1) ]
        else:
            if e == 2:
                L += [ 2 ]
            elif e >= 3:
                L += [ 2**(e - 2) ]

    return rosemary.number_theory.core.lcm(L)

########################################################################################################################
# lists of values of arithmetic functions
########################################################################################################################

def euler_phi_list(n):
    """
    Returns a list of values of euler_phi(k), for 1 <= k <= n.

    Input:
        * n: int

    Output:
        * block: list
            This is a list of values of euler_phi(k) for 1 <= k <= n. The list begins with 0, so that block[k] holds the
            value of euler_phi(k).
    """
    block = [1]*(n + 1)
    block[0] = 0
    p_list = rosemary.number_theory.sieves.primes(n)

    for p in p_list:
        for mul in xrange(p, n + 1, p):
            block[mul] *= (p - 1)
        pk = p*p
        while pk <= n:
            for mul in xrange(pk, n + 1, pk):
                block[mul] *= p
            pk *= p
    return block

totient_list = euler_phi_list

########################################################################################################################
# sums of values of arithmetic functions
########################################################################################################################

def euler_phi_sum(n):
    """
    Returns the value of the sum of euler_phi(k) for k = 1..n.

    Input:
        * n: int

    Output:
        * sum: int
            This is the sum of euler_phi(k) for k = 1..n.

    Details:
        Let S(n) = sum_{k = 1}^{n} \phi(k). We use the identity S(n) = n*(n + 1)/2 - \sum_{d = 2}^{n} S(n/d) to compute
        the value in sublinear time by caching the recursion.
    """
    k = int(n**(0.5))
    totients = euler_phi_list(k)

    # for i <= sqrt(n), compute the sum directly
    cache = {1:1}
    for i in xrange(2, k):
        cache[i] = cache[i - 1] + totients[i]

    def S(n):
        if n in cache:
            return cache[n]
        srn = int(n**(0.5))
        S1 = sum([S(n//j) for j in xrange(srn, 1, -1)])
        S2 = sum([totients[j]*(n//j) for j in xrange(1, srn + 1)])
        S3 = srn*S(srn)

        cache[n] = n*(n + 1)//2 - S1 - S2 + S3
        return cache[n]

    return S(n)

totient_sum = euler_phi_sum

########################################################################################################################
# Miscellaneous
########################################################################################################################

def euler_phi_inverse(m):
    """
    Returns a sorted list of all integers n such that euler_phi(n) = m.
    """
    # We know euler_phi(n) is even for m >= 3
    if m % 2 == 1 and m > 1:
        return []

    mfact = rosemary.number_theory.factorization.factor(m)
    
    if m % 2 == 0:
        twopows = [ 2**i for i in range(mfact[0][1] + 1) ]
    else:
        twopows = [ 1 ]

    D = rosemary.number_theory.factorization.divisors(mfact)
    P = []

    for d in D:
        if d == 1:
            continue
        if rosemary.number_theory.primality.is_prime(d + 1):
            P.append(d + 1)

    S = [(1, m)]
    for p in reversed(P):
        T = []
        for (s0, s1) in S:
            if s1 == 1:
                continue
            pk = p
            d = s1 // (p - 1)
            mmod = s1 % (p - 1)
            while mmod == 0:
                if d % 2 == 0 or d == 1:
                    T.append((pk * s0, d))
                pk *= p
                mmod = d % p
                d = d // p
        S.extend(T)

    R = set([])
    for s in S:
        if s[1] in twopows:
            j = twopows.index(s[1])
            R.add(2**(j + 1) * s[0])
            if j == 0:
                R.add(s[0])

    return sorted(R)

