# Sieves

from rosemary.number_theory.core import integer_sqrt
from bisect import bisect_left

def primes(n):
    """
    Returns a list of all primes <= n.

    Input:
        * n: int (n > 0)

    Output:
        * P: list

    Examples:
        >>> primes(20)
        [2, 3, 5, 7, 11, 13, 17, 19]
        >>> len(primes(10**7))
        664579

    Details:
        This program uses an optimized sieve of Eratosthenes to generate a list
        of all primes <= n. The sieve only looks at numbers <= n that are
        coprime to 6. The code is based on some found on stackoverflow.
    """
    n += 1
    offset = (n%6 > 1)
    n = {0:n, 1:n - 1, 2:n + 4, 3:n + 3, 4:n + 2, 5:n + 1}[n % 6]

    sieve = [True] * (n//3)
    sieve[0] = False
    sr = integer_sqrt(n)

    for i in xrange(sr//3 + 1):
        if sieve[i]:
            k = (3*i + 1)|1
            kk = k*k
            sieve[kk//3::2*k] = [False]*((n//6 - kk//6 - 1)//k + 1)
            sieve[(kk + 4*k - 2*k*(i&1))//3::2*k] = [False]*((n//6 - (kk + 4*k - 2*k*(i&1))//6 - 1)//k + 1)

    return [2, 3] + [(3*i + 1)|1 for i in xrange(1, n//3 - offset) if sieve[i]]

def eratosthenes(n):
    """
    Returns a list of all primes <= n.

    This program uses the sieve of Eratosthenes to generate a list of all primes <= n.

    Input:
        * n - A positive integer.

    Output:
        * L - a list of primes.

    Details:
        This slightly optimized sieve only looks at the odd numbers <= n. The
        implementation is included mainly for reference.

    Examples:
        >>> eratosthenes(100)
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
            67, 71, 73, 79, 83, 89, 97]
        >>> len(eratosthenes(10**7))
        664579
    """
    m = n//2 + 1
    b_list = [True]*m
    N = integer_sqrt(m) + 1

    for i in xrange(1, N):
        if b_list[i]:
            k = 2*i + 1
            start = (k*k - 1)/ 2
            how_many = (m - start)//k + ((m - start) % k > 0)
            b_list[start::k] = [False] * how_many

    return [2] + [2*i + 1 for i in xrange(1, n // 2 + n%2) if b_list[i]]


def sieve_interval(a, b):
    """
    Returns an iterator over all primes in the interval [a, b].

    Given positive integers a and b with a*a > b, this returns an interator
    over all primes in the interval [a, b].

    Input:
        * a - A positive integer.
        * b - A positive integer.

    Output:
        * X - An iterator over all primes in [a, b].

    Details:
        This algorithm is from Section 3.2.2 of "Prime Numbers - A
        Computational Perspective" by Crandall and Pomerance. They call it the
        "Practical Eratosthenes sieve".

    Examples:
        >>> list(sieve_interval(100, 200))
        [101, 103, 107, 109, 113, 127, 131, 137, 139, 149]
        >>> len(list(sieve_interval(10**8, 2*10**8)))
        5317482
    """
    # We need both endpoints to be even
    if a % 2 == 1:
        a -= 1
    if b % 2 == 1:
        b += 1

    bd = integer_sqrt(b)
    assert a > bd

    p_list = primes(bd)[1:]
    block_size= min((b - a) // 2, 8*bd)

    offsets = {p: -(a + 1 + p)//2 % p for p in p_list}

    for start in xrange(a, b, 2*block_size):
        block = [1]*block_size
        for p in p_list:
            offset = offsets[p]
            diff = block_size - offset
            how_many = diff//p + (diff % p > 0)
            block[offset::p] = [0]*how_many
            offsets[p] = (-diff) % p
        for j in xrange(block_size):
            if block[j]:
                m = start + 2*j + 1
                if a <= m <= b:
                    yield m

def prime_xrange(a, b=None):
    if b is None:
        b = a
        a = 2

    block_size = integer_sqrt(b)
    p_list = primes(block_size)

    if a <= p_list[-1]:
        idx = bisect_left(p_list, a)
        for i in xrange(idx, len(p_list)):
            p = p_list[i]
            if a <= p <= block_size:
                yield p

    for start in xrange(block_size + 1, b, block_size):
        block = [1]*block_size
        for p in p_list:
            offset = ((p*(start//p + 1) - 1) % block_size) % p
            for idx in xrange(offset, block_size, p):
                block[idx] = 0

        for i in xrange(block_size):
            if block[i] and a <= start + i <= b:
                yield start + i

################################################################################
# Sieves for generating values of arithmetical functions
################################################################################

def moebius_xrange(a, b=None):
    """
    Return an iterator over values of moebius(k) for 1 <= k <= n.
    """
    if b is None:
        b, a = a, 1

    blockSize = integer_sqrt(b)
    primeList = primes(blockSize)

    for start in xrange(a, b, blockSize):
        block = [1]*blockSize
        values = [1]*blockSize

        for p in primeList:
            offset = ((p*(start//p + 1) - a) % blockSize) % p
            for i in xrange(offset, blockSize, p):
                if (start + i) % (p*p) == 0:
                    block[i] = 0
                else:
                    block[i] *= -1
                    values[i] *= p

        for i in xrange(blockSize):
            if start + i >= b:
                return
            if block[i] and values[i] < start + i:
                block[i] *= -1
            yield (start + i, block[i])

def factored_xrange(a, b=None):
    if b is None:
        b, a = a, 1

    blockSize = integer_sqrt(b)
    primeList = primes(blockSize)

    if a == 1:
        yield (1, [(1, 1)])
        a += 1

    for start in xrange(a, b, blockSize):
        block = range(start, start + blockSize)
        factorizations = [[] for _ in xrange(blockSize)]

        for p in primeList:
            offset = ((p*(start//p + 1) - a) % blockSize) % p
            for i in xrange(offset, blockSize, p):
                k = 0
                while block[i] % p == 0:
                    block[i] /= p
                    k += 1
                factorizations[i].append((p, k))

        for i in xrange(blockSize):
            if start + i >= b:
                return
            if block[i] != 1:
                factorizations[i].append((block[i], 1))
            yield (start + i, factorizations[i])

