# Arithmetic Functions

from rosemary.number_theory.prime_list import PRIME_LIST

from rosemary.number_theory.core import (
    bit_count,
    lcm,
)

import rosemary.number_theory.factorization
import rosemary.number_theory.primality
import rosemary.number_theory.sieves

########################################################################################################################
# classical arithmetic functions
########################################################################################################################

def euler_phi(n):
    """
    Returns the number of positive integers <= n that are coprime to n.

    Input:
        * n: int (n > 0)

    Output:
        * r: int

    Details:
        This routine works by factoring n, and using the standard product
        definition of phi(n).

    Examples:
        >>> euler_phi(100)
        40
        >>> euler_phi(11213)
        11212
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n <= 0:
            raise ValueError("euler_phi: Input must be positive.")
        nFactors = rosemary.number_theory.factorization.factor(n)
    elif isinstance(n, list) and n[0][0] > 0:
        nFactors = n
    else:
        raise ValueError("euler_phi: Input must be a positive integer or a factorization.")

    prod = 1
    for (p, e) in nFactors:
        prod *= p**(e - 1)*(p - 1)
    return prod

def factorial(n):
    """
    Returns the factorial of n.

    Input:
        * n: int (n >= 0)

    Output:
        * prod: int

    Details:
        The algorithm used is a binary splitting method.

    Examples:
        >>> factorial(10)
        3628800
        >>> factorial(40)
        815915283247897734345611269596115894272000000000
        """
    if n < 0:
        raise ValueError("factorial: Input must be nonnegative.")
    elif n <= 1:
        return 1

    prod = 1
    for k in xrange(1, n.bit_length()):
        lower = (n >> k) + 1
        upper = (n >> (k - 1)) + 1

        if lower % 2 == 0:
            lower += 1

        partial = 1
        for j in xrange(lower, upper, 2):
            partial *= j
        prod *= (partial**k)
    return prod << (n - bit_count(n))

def moebius(n):
    """
    Returns the value of the Moebius function mu(n).

    Input:
        * n: int (n > 0)

    Output:
        * mu: int

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
        elif n <= 0:
            raise ValueError("moebius: Input must be positive.")
        nFactors = rosemary.number_theory.factorization.factor(n)
    elif isinstance(n, list) and n[0][0] > 0:
        nFactors = n
    else:
        raise ValueError("moebius: Input must be a positive integer or a factorization.")

    if rosemary.number_theory.factorization.is_squarefree(nFactors):
        k = len(nFactors)
        return (-1)**k
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

def sigma(n, k=1):
    """
    Returns the sum of th k-th powers of the divisors of n.

    Input:
        * n: int (n > 0)
        * k: int (k > 0)

    Output:
        * s: int

    Examples:
        >>> sigma(9)
        13
        >>> 1 + 3 + 9
        13
        >>> sigma(10, 2)
        130
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n <= 0:
            raise ValueError("sigma: Input must be positive.")
        nFactors = rosemary.number_theory.factorization.factor(n)
    elif isinstance(n, list) and n[0][0] > 0:
        nFactors = n
    else:
        raise ValueError("sigma: Input must be a positive integer or a factorization")

    prod = 1
    for (p, e) in nFactors:
        pk = p**k
        prod *= (pk**(e + 1) - 1)//(pk - 1)
    return prod

def tau(n):
    """
    Given an integer n, this returns the number of divisors of n.

    Input:
        * n: int (n > 0)

    Output:
        * num: int
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n <= 0:
            raise ValueError("tau: Input must be positive.")
        nFactors = rosemary.number_theory.factorization.factor(n)
    elif isinstance(n, list) and n[0][0] > 0:
        nFactors = n
    else:
        raise ValueError("tau: Input must be a positive integer or a factorization")

    prod = 1
    for (p, e) in nFactors:
        prod *= (e + 1)
    return prod

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
        n_fac = rosemary.number_theory.factorization.factor(n)

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

    return lcm(L)

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

def moebius_list(n):
    """
    Return an iterator over values of moebius(k) for 1 <= k <= n.
    """
    sr = int(n**(0.5))
    p_list = rosemary.number_theory.sieves.primes(sr)
    block = [1]*(n + 1)
    vals = [1]*(n + 1)
    block[0] = 0

    for p in p_list:
        for i in xrange(p, n + 1, p):
            if i % (p*p) == 0:
                block[i] = 0
            else:
                block[i] *= -1
                vals[i] *= p

    for i in xrange(n + 1):
        if block[i] and vals[i] < i:
            block[i] *= -1

    return block


def sigma_list(n, k=1):
    """
    Returns a list of the values of sigma(i, k), for 1 <= i <= n.
    """
    assert k >= 1

    block = [1]*(n + 1)
    block[0] = 0
    p_list = rosemary.number_theory.sieves.primes(n)

    for p in p_list:
        pk = p
        mul = p**k
        term = mul**2
        last = mul
        while pk <= n:
            for idx in xrange(pk, n + 1, pk):
                block[idx] *= (term - 1)
                block[idx] /= (last - 1)
            pk *= p
            last = term
            term *= mul

    return block

def tau_list(n, k=1):
    """
    Returns a list of the values of tau(i), for 1 <= i <= n.
    """
    block = [0]*(n + 1)
    bound = int(n**(0.5)) + 1

    for j in xrange(1, bound):
        block[j*j] += 1
        for k in xrange(j + 1, n//j + 1):
            block[k*j] += 2
    return block


########################################################################################################################
# Summatory functions
########################################################################################################################

def euler_phi_sum(n):
    """
    Returns the value of the sum of euler_phi(k) for k = 1, 2, ..., n.

    Input:
        * n: int

    Output:
        * sum: int
            This is the sum of euler_phi(k) for k = 1..n.

    Details:
        Let S(n) = \sum_{k = 1}^{n} \phi(k). We use the identity S(n) = n*(n + 1)/2 - \sum_{d = 2}^{n} S(n/d) to compute
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
        s1 = sum([S(n//j) for j in xrange(srn, 1, -1)])
        s2 = sum([totients[j]*(n//j) for j in xrange(1, srn + 1)])
        s3 = srn*S(srn)

        cache[n] = n*(n + 1)//2 - (s1 + s2 - s3)
        return cache[n]

    return S(n)

totient_sum = euler_phi_sum

def euler_phi_weighted_sum(n):
    sr = int(n**(0.5))
    totients = euler_phi_list(sr)

    # for i <= sqrt(n), compute the sum directly
    cache = {1:1}
    for i in xrange(2, sr + 1):
        cache[i] = cache[i - 1] + i*totients[i]

    def T(n):
        if n in cache:
            return cache[n]

        sn = int(n**(0.5)) + 1
        s1 = sum(d*totients[d]*(n//d - d + 1)*(n//d + d)//2 for d in xrange(1, sn))
        s2 = sum(d*(T(n//d) - T(d - 1)) for d in xrange(2, sn))
        s3 = sum(d*d*totients[d] for d in xrange(1, sn))

        cache[n] = n*(n + 1)*(2*n + 1)//6 - (s1 + s2 - s3)
        return cache[n]

    return T(n)

def moebius_sum(n):
    """
    Returns the value of the sum of moebius(k) for k = 1, 2, ..., n.

    Input:
        * n: int

    Output:
        * sum: int
    """
    sn = int(n**(0.5))
    mu = moebius_list(sn)

    cache = {1:1}
    for i in xrange(2, sn + 1):
        cache[i] = cache[i - 1] + mu[i]

    def M(n):
        if n in cache:
            return cache[n]

        sr = int(n**(0.5)) + 1
        s1 = sum(M(n//k) - M(k - 1) for k in xrange(2, sr))
        s2 = sum(mu[l]*(n//l - l + 1) for l in xrange(1, sr))
        s3 = sum(mu[k] for k in xrange(1, sr))

        cache[n] = 1 - (s1 + s2 - s3)
        return cache[n]

    return M(n)

mertens = moebius_sum

def sigma_sum(n):
    """
    Returns the value of the sum of sigma(k) for k = 1, 2, ..., n.

    Input:
        * n: int

    Output:
        * sum: int
    """
    m = int(n**(0.5))
    ss = -m*(m + 1)

    for k in xrange(1, m + 1):
        nk = n//k
        tt = nk - k + 1
        ss += 2*k*tt
        ss += tt*(nk + k)

    return ss//2

########################################################################################################################
# Miscellaneous
########################################################################################################################

def euler_phi_inverse(m):
    """
    Returns a list of all positive integers n such that euler_phi(n) = m.

    Input:
        * m: int

    Output:
        * values: list
            A sorted list of all positive integers n with euler_phi(n) = m.

    Details:
        This uses the algorithm outlined in Discovering Mathematics with Magma
        by Bosma and Cannon.

    Examples:
        >>> euler_phi_inverse(40)
        [41, 55, 75, 82, 88, 100, 110, 132, 150]
        >>> euler_phi_inverse(103)
        []
    """
    # We know euler_phi(n) is even for m >= 3
    if m % 2 == 1 and m > 1:
        return []

    mFactors = rosemary.number_theory.factorization.factor(m)
    
    powersOfTwo = set([1])
    if m % 2 == 0:
        for i in xrange(1, mFactors[0][1] + 1):
            powersOfTwo.add(2**i)

    # Odd primes p that divide n must have the property that p - 1 divides m.
    primeList = []
    for d in rosemary.number_theory.factorization.divisors(mFactors)[1:]:
        if rosemary.number_theory.primality.is_prime(d + 1):
            primeList.append(d + 1)

    # Here, we store pairs (a, b) such that a is odd and phi(a) = m/b, with b
    # even or 1. Every pair contributes at least one solution. When b = 1, the
    # pair contributes two solutions.
    pairs = [(1, m)]
    for p in reversed(primeList):
        newPairs = []
        for (a, b) in pairs:
            if b == 1:
                continue
            pk = p
            d = b//(p - 1)
            mmod = b % (p - 1)
            while mmod == 0:
                if d % 2 == 0 or d == 1:
                    newPairs.append((pk*a, d))
                pk *= p
                mmod = d % p
                d = d//p
        pairs.extend(newPairs)

    # When b = 2^k, we have the solution 2*b*a.
    values = []
    for (a, b) in pairs:
        if b in powersOfTwo:
            values.append(2*b*a)
            if b == 1:
                values.append(a)

    values.sort()
    return values

