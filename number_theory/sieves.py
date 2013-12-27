# Sieves

from rosemary.number_theory.core import integer_sqrt
from bisect import bisect_left

################################################################################
# Sieves for generating primes
################################################################################

def primes(n):
    """
    Returns a list of all primes p <= n.

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

    block = [True]*(n//3)
    block[0] = False
    sqrt = int(n**(0.5))

    for i in xrange(sqrt//3 + 1):
        if block[i]:
            k = (3*i + 1)|1
            kk = k*k
            block[kk//3::2*k] = [False]*((n//6 - kk//6 - 1)//k + 1)
            block[(kk + 4*k - 2*k*(i&1))//3::2*k] = [False]*((n//6 - (kk + 4*k - 2*k*(i&1))//6 - 1)//k + 1)

    return [2, 3] + [(3*i + 1)|1 for i in xrange(1, n//3 - offset) if block[i]]


def eratosthenes(n):
    """
    Returns a list of all primes p <= n.

    Input:
        * n: int (n > 0)

    Output:
        * P: list

    Examples:
        >>> eratosthenes(20)
        [2, 3, 5, 7, 11, 13, 17, 19]
        >>> len(eratosthenes(10**7))
        664579

    Details:
        This slightly optimized sieve only looks at the odd numbers <= n. The
        implementation is included mainly for reference.
    """
    m = n//2 + 1
    block = [True]*m
    sqrt = int(m**(0.5))

    for i in xrange(1, sqrt + 1):
        if block[i]:
            k = 2*i + 1
            start = (k*k - 1)/ 2
            count = (m - start)//k + ((m - start) % k > 0)
            block[start::k] = [False]*count

    return [2] + [2*i + 1 for i in xrange(1, n//2 + n % 2) if block[i]]


def prime_xrange(a, b=None):
    """
    Returns iterator over primes in interval [a, b).
    """
    def sieve_interval(a, b, sqrt, prime_list):
        """
        This algorithm is from Section 3.2.2 of "Prime Numbers - A Computational
        Perspective" by Crandall and Pomerance. They call it the "Practical
        Eratosthenes sieve". We assume that a > sqrt(b) here.
        """
        # aa is the largest even number <= a.
        # bb is the smallest even number >= b.
        aa = a - (a % 2)
        bb = b + (b % 2)
        block_size= min((bb - aa)//2, 8*sqrt)
        offsets = {p: -(aa + 1 + p)//2 % p for p in prime_list}

        for start in xrange(aa, bb, 2*block_size):
            block = [1]*block_size
            for p in prime_list:
                offset = offsets[p]
                diff = block_size - offset
                how_many = diff//p + (diff % p > 0)
                block[offset::p] = [0]*how_many
                offsets[p] = (-diff) % p
            for j in xrange(block_size):
                if block[j]:
                    m = start + 2*j + 1
                    if a <= m < b:
                        yield m

    cutoff = 10**7
    if b is None:
        b = a
        a = 2

    # If the upper bound is less than some appropriate cutoff, we generate all
    # primes up to this limit yield the ones in the interval.
    if b <= cutoff:
        prime_list = primes(b)
        idx = bisect_left(prime_list, a)
        for i in xrange(idx, len(prime_list)):
            yield prime_list[i]
    else:
        # Otherwise use a segmented sieve on the interval.
        sqrt = int(b**(0.5))
        prime_list = primes(sqrt)
        # Our segmented sieve code assumes that a > sqrt(b). If this is not the
        # case, then we yield the primes <= sqrt(b) first.
        if a <= sqrt:
            # If a*a <= b, iterate over the primes a <= p <= sqrt(b) first.
            idx = bisect_left(prime_list, a)
            for i in xrange(idx, len(prime_list)):
                yield prime_list[i]
            # Then iterate over the primes sqrt(b) < p <= b.
            for p in sieve_interval(sqrt + 1, b, sqrt, prime_list[1:]):
                yield p
        else:
            # Otherwise, we use our segmented sieve function.
            for p in sieve_interval(a, b, sqrt, prime_list[1:]):
                yield p

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
    """
    Returns an iterator over the factorizations of the numbers in [a, b).

    Given positive integers a and b with a < b, this returns an iterator over
    all pairs (n, n_factorization) with a <= n < b, and n_factorization is the
    factorization of n into prime powers. If the optional parameter b is None,
    then a is taken to be 1, and b = a.

    Input:
        * a: int (a > 0)
        * b: int (b > a) (default=None)

    Output:
        * P: generator
            The values output by this generator are tuples (n, n_factorization),
            where n is an integer in [a, b), and n_factorization is the prime
            factorization of n.

    Examples:
        >>> list(factored_xrange(10, 20))
        [(10, [(2, 1), (5, 1)]),
         (11, [(11, 1)]),
         (12, [(2, 2), (3, 1)]),
         (13, [(13, 1)]),
         (14, [(2, 1), (7, 1)]),
         (15, [(3, 1), (5, 1)]),
         (16, [(2, 4)]),
         (17, [(17, 1)]),
         (18, [(2, 1), (3, 2)]),
         (19, [(19, 1)])]
    """
    if b is None:
        b, a = a, 1

    block_size = int(b**(0.5))
    prime_list = primes(block_size)

    if a == 1:
        yield (1, [(1, 1)])
        a += 1

    for start in xrange(a, b, block_size):
        block = range(start, start + block_size)
        factorizations = [[] for _ in xrange(block_size)]

        for p in prime_list:
            offset = ((p*(start//p + 1) - a) % block_size) % p
            for i in xrange(offset, block_size, p):
                k = 0
                while block[i] % p == 0:
                    block[i] /= p
                    k += 1
                factorizations[i].append((p, k))

        for i in xrange(block_size):
            if start + i >= b:
                return
            if block[i] != 1:
                factorizations[i].append((block[i], 1))
            yield (start + i, factorizations[i])

