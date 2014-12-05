# Sieves

from rosemary.number_theory.core import integer_sqrt
from rosemary.number_theory.prime_list import _PRIME_LIST

from bisect import bisect, bisect_left
from heapq import heappush, heappop
from math import log


################################################################################
# Sieves for generating primes
################################################################################


def primes(n):
    """Returns a list of all primes <= n.

    Input:
        * n: int

    Returns:
        * primes: list
            A list of the primes <= n.

    Examples:
        >>> primes(20)
        [2, 3, 5, 7, 11, 13, 17, 19]
        >>> len(primes(10**6))
        78498

    Details:
        For values of n <= 999983, the precomputed list of primes is used.
        Otherwise, an optimized version of the sieve of Eratosthenes is used.
        See the method `eratosthenes2` for more details.
    """
    if n < _PRIME_LIST[-1]:
        k = bisect(_PRIME_LIST, n)
        return _PRIME_LIST[:k]
    else:
        return eratosthenes2(n)


def primes_first_n(n):
    """Returns a list of the first n primes.

    Input:
        * n: int

    Returns:
        * primes: list
            A list of the first n primes.

    Examples:
        >>> primes_first_n(10)
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        >>> sum(primes_first_n(100))
        24133

    Details:
        This method uses an upper bound for the nth prime, generates all primes
        to that point, and returns the first n of them. The estimate used is p_n
        < n*(log(n) + log(log(n))) for n >= 6, where p_n denotes the nth prime.
        See Theorem 8.8.4 of "Algorithmic Number Theory - Efficient Algorithms"
        by Bach and Shallit for details.
    """
    if n <= 0:
        return []
    elif n < len(_PRIME_LIST):
        return _PRIME_LIST[:n]
    else:
        m = int(n*(log(n) + log(log(n))))
        return eratosthenes2(m)[:n]


def eratosthenes1(n):
    """Returns a list of all primes p <= n.

    Input:
        * n: int

    Returns:
        * primes: list
            A list of the primes <= n.

    Examples:
        >>> eratosthenes1(20)
        [2, 3, 5, 7, 11, 13, 17, 19]
        >>> len(eratosthenes(10**6))
        78498

    Details:
        This algorithm is an implementation of the first extension of the sieve
        of Eratosthenes by removing the multiples of 2 from the set of
        candidates. This has the effect of halving the required space and time
        of the algorithm. See the paper "A Practical Sieve Algorithm for Finding
        Prime Numbers" by Luo for details. See the commentary "Additional Notes
        on a Practical Sieve Algorithm" by Quesada for more details and
        extensions.
    """
    if n <= 1:
        return []

    m = n//2 + 1
    block = [True]*m
    sqrt = int(m**(0.5))

    for i in xrange(1, sqrt + 1):
        if block[i]:
            k = 2*i + 1
            start = k*k//2
            count = (m - start)//k + ((m - start) % k > 0)
            block[start::k] = [False]*count

    return [2] + [2*i + 1 for i in xrange(1, n//2 + n % 2) if block[i]]


def eratosthenes2(n):
    """Returns a list of all primes p <= n.

    Input:
        * n: int

    Returns:
        * primes: list
            A list of the primes <= n.

    Examples:
        >>> eratosthenes2(20)
        [2, 3, 5, 7, 11, 13, 17, 19]
        >>> len(eratosthenes2(10**6))
        78498

    Details:
        This algorithm is an implementation of the second extension of the sieve
        of Eratosthenes by removing the multiples of 2 and 3 from the set of
        candidates. This has the effect of decreasing the required space and
        time to 1/3 of the naive Eratosthenes algorithm. See the paper "A
        Practical Sieve Algorithm for Finding Prime Numbers" by Luo for details.
        See the commentary "Additional Notes on a Practical Sieve Algorithm" by
        Quesada for more details and extensions.
    """
    if n <= 1:
        return []
    elif n <= 3:
        return range(2, n + 1)

    n += 1
    offset = (n % 6 > 1)
    n = {0: n, 1: n - 1, 2: n + 4, 3: n + 3, 4: n + 2, 5: n + 1}[n % 6]

    m = n//3
    block = [True]*(m)
    block[0] = False
    sqrt = int(n**(0.5))

    for i in xrange(sqrt//3 + 1):
        if block[i]:
            k = (3*i + 1) | 1
            kk = k*k
            a = kk//3
            b = (kk + 4*k - 2*k*(i & 1))//3
            block[a::2*k] = [False]*((m - a)//(2*k) + ((m - a) % (2*k) > 0))
            block[b::2*k] = [False]*((m - b)//(2*k) + ((m - b) % (2*k) > 0))

    return [2, 3] + [(3*i + 1) | 1 for i in xrange(1, n//3 - offset) if block[i]]


def eratosthenes3(n):
    """Returns a list of all primes p <= n.

    Input:
        * n: int

    Returns:
        * primes: list
            A list of the primes <= n.

    Examples:
        >>> eratosthenes3(20)
        [2, 3, 5, 7, 11, 13, 17, 19]
        >>> len(eratosthenes3(10**6))
        78498

    Details:
        This algorithm is an implementation of the third extension of the sieve
        of Eratosthenes by removing the multiples of 2, 3, and 5 from the set of
        candidates. This has the effect of decreasing the required space and
        time to 8/30 of the naive Eratosthenes algorithm. Despite this lowered
        complexity, the second extension seems to be faster in practice due to
        the simpler implementation. See the paper "A Practical Sieve Algorithm
        for Finding Prime Numbers" by Luo for details.  See the commentary
        "Additional Notes on a Practical Sieve Algorithm" by Quesada for more
        details and extensions.
    """
    def offsets(posp):
        q = posp // 8
        r = posp % 8
        D = [0]*8
        D[0] = 2*posp + 1

        if r == 0:
            D[3] = D[0]
            D[2] = 2*D[3] - 1
            D[4] = D[2]
            D[1] = D[2] + D[3] - 1
        elif r == 1:
            D[3] = D[0] + 1
            D[4] = D[3] + D[0]
            D[2] = D[3] + D[0]
            D[1] = D[2] + D[3] + 1
        elif r == 2:
            D[3] = D[0] + 1
            D[2] = 2*D[3]
            D[4] = D[2] - 1
            D[1] = D[2] + D[3]
        elif r == 3:
            D[3] = D[0]
            D[2] = 2*D[3]
            D[4] = D[2] - 1
            D[1] = D[2] + D[3]
        elif r == 4:
            D[3] = D[0]
            D[2] = 2*D[3]
            D[4] = D[2] + 1
            D[1] = D[2] + D[3]
        elif r == 5:
            D[3] = D[0] - 1
            D[2] = 2*D[3]
            D[4] = D[2] + 1
            D[1] = D[2] + D[3]
        elif r == 6:
            D[3] = D[0] - 1
            D[2] = 2*D[3] + 1
            D[4] = D[2]
            D[1] = D[2] + D[3] - 1
        elif r == 7:
            D[3] = D[0]
            D[2] = 2*D[3] + 1
            D[4] = D[2]
            D[1] = D[2] + D[3] + 1

        D[7] = D[1]
        D[6] = D[2]
        D[5] = D[3]

        return q, r, D

    if n <= 1:
        return []
    elif n <= 5:
        return {2: [2], 3: [2, 3], 4: [2, 3], 5: [2, 3, 5]}[n]

    nn = n
    while n % 30 != 29:
        n += 1

    k = (n - 29)//30
    posp = 1

    S = [False]
    i = 0
    for i in xrange(k + 1):
        for r in [1, 7, 11, 13, 17, 19, 23, 29]:
            t = 30*i + r
            if 1 < t <= n:
                S.append(True)
            elif t > n:
                break

    f = 8*int(n**(0.5))//30
    B = [1, 7, 11, 13, 17, 19, 23, 29]
    BS = [0, 13, 32, 45, 77, 96, 141, 224]
    m = len(S)

    for posp in xrange(1, f + 1):
        if not S[posp]:
            continue

        q, r, D = offsets(posp)
        p = 30*q + B[r]

        pos = 8*q*(p + B[r]) + BS[r]
        i = (r + 1) % 8
        diff = 8*p

        for j in xrange(8):
            S[pos::diff] = [False]*((m - pos)//diff + ((m - pos) % diff > 0))
            pos += D[i]
            i = (i + 1) % 8

    return [2, 3, 5] + [30*(j//8) + B[j % 8] for j in xrange(m) if S[j] and 30*(j//8) + B[j % 8] <= nn]


def luo(n):
    """Returns a list of all primes p <= n.

    Input:
        * n: int

    Returns:
        * primes: list
            A list of the primes <= n.

    Examples:
        >>> luo(20)
        [2, 3, 5, 7, 11, 13, 17, 19]
        >>> len(luo(10**6))
        78498

    Details:
        This algorithm is an implementation of the second extension of the sieve
        of Eratosthenes by removing the multiples of 2 and 3 from the set of
        candidates. This has the effect of decreasing the required space and
        time to 1/3 of the naive Eratosthenes algorithm. This implementation is
        based on the the exposition in the paper "A Practical Sieve Algorithm
        for Finding Prime Numbers" by Luo, as well as the optimizations made by
        Pritchard in the article "Additional Notes on a Practical Sieve
        Algorithm".
    """
    if n <= 1:
        return []
    elif n <= 6:
        return {2: [2], 3: [2, 3], 4: [2, 3], 5: [2, 3, 5], 6: [2, 3, 5]}[n]

    n += 1
    offset = (n % 6 > 1)
    n = {0: n, 1: n - 1, 2: n + 4, 3: n + 3, 4: n + 2, 5: n + 1}[n % 6]

    m = n//3
    sqrt = int(n**(0.5))
    block = [True]*(n//3)
    a, k, t = 0, 1, 2

    for i in xrange(1, sqrt + 1):
        k = 3 - k
        a += 4*k*i
        t += 4*k

        if block[i]:
            b = a + 2*i*(3 - k) + 1
            block[a::t] = [False]*((m - a)//t + ((m - a) % t > 0))
            block[b::t] = [False]*((m - b)//t + ((m - b) % t > 0))

    return [2, 3] + [(3*i + 1) | 1 for i in xrange(1, m - offset) if block[i]]


def chartres(n):
    """Returns a list of all primes <= n.

    Input:
        * n: int

    Returns:
        * primes: list
            A list of primes <= n.

    Examples:
        >>> chartres(30)
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        >>> len(chartres(1000))
        168

    Details:
        This method uses the algorithm in Section 5.2.3, Exercise 15 of "The Art
        of Computer Programming Volume 3" by Knuth. The algorithm is originally
        due to Chartres. See "Algorithm 310 - Prime Number Generator" of the
        Communications of the ACM, Volume 10, Number 9, September 1967 for more
        information.

        Unlike the standard implementation of the sieve of Eratosthenes, this
        algorithm uses a heap and requires O(n*log(n)) steps, while the sieve of
        Eratosthenes runs in time O(n*log(log(n))).
    """
    primes = [0, 2, 3]
    num_primes = 2
    candidate = 5
    diff = 2
    size = 1
    q = t = 25

    # Entries in the heap have the form (u, v, 6*p), where p is the smallest
    # prime divisor of u, v = 2*p or 4*p, and u + v is not a multiple of 3.
    # The heap size is at most the number of primes <= sqrt(N).
    heap = [(25, 10, 30)]

    add_prime = primes.append

    while heap:
        while q < candidate:
            (q, qp, qpp) = heappop(heap)
            heappush(heap, (q + qp, qpp - qp, qpp))

        if candidate > n:
            return [e for e in primes if e and e <= n]

        while candidate < q:
            num_primes += 1
            add_prime(candidate)
            candidate += diff
            diff = 6 - diff

        if candidate == t:
            size += 1
            u = primes[size + 2]
            t = u*u

            if u % 3 == 2:
                heappush(heap, (t, 2*u, 6*u))
            elif u % 3 == 1:
                heappush(heap, (t, 4*u, 6*u))

        candidate += diff
        diff = 6 - diff


def prime_xrange(a, b=None):
    """Returns iterator over primes in interval [a, b).

    Input:
        a: int (a > 0)
        b: int (b > a) (default=None)

    Returns:
        * P: generator
            The values output by this generator are the primes in the interval
            [a, b).

    Examples:
        >>> list(prime_xrange(10, 20))
        [11, 13, 17, 19]

    Details:
        This method uses a variety of techniques to efficiently find the primes
        in the given interval. For b <= 10**7, the primes <= b are computed, and
        the primes in [a, b) are yielded. Note that for b <= 10**6, the
        precomptued list of primes is used here. For b > 10**7, a segmented
        sieve of Eratosthenes is used. See Section 3.2.2 of "Prime Numbers - A
        Computational Perspective" by Crandall and Pomerance for details.
    """
    def sieve_interval(a, b, sqrt, prime_list):
        """This algorithm is from Section 3.2.2 of "Prime Numbers - A
        Computational Perspective" by Crandall and Pomerance. They call it the
        "Practical Eratosthenes sieve". We assume that a > sqrt(b) here.
        """
        # aa is the largest even number <= a.
        # bb is the smallest even number >= b.
        aa = a - (a % 2)
        bb = b + (b % 2)
        block_size = min((bb - aa)//2, 8*sqrt)
        offsets = {p: -(aa + 1 + p)//2 % p for p in prime_list}

        for start in xrange(aa, bb, 2*block_size):
            block = [1]*block_size

            for p in prime_list:
                offset = offsets[p]
                diff = block_size - offset
                block[offset::p] = [0]*(diff//p + (diff % p > 0))
                offsets[p] = (-diff) % p

            for j in xrange(block_size):
                if block[j]:
                    m = start + 2*j + 1
                    if a <= m < b:
                        yield m

    cutoff = 10**7
    if a < 2:
        a = 2

    if b is None:
        b = a
        a = 1

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
            idx = bisect(prime_list, a)
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
    """Return an iterator over values of moebius(k) for a <= k < b.

    Input:
        * a: int (a > 0)
        * b: int (b > a) (default=None)

    Returns:
        * P: generator
            The values output by this generator are tuples (n, moebius(n)) where
            n is an integer in [a, b), and moebius(n) is the value of the
            moebius function.

    Examples:
        >>> list(moebius_xrange(10, 20))
        [(10, 1), (11, -1), (12, 0), (13, -1), (14, 1), (15, 1), (16, 0), (17, -1), (18, 0), (19, -1)]

    Details:
        All primes <= sqrt(b) are computed, and a segmented sieve is used to
        build the values using the multiplicative property of the moebius
        function. See the paper "Computing the Summation of the Moebius
        Function" by Deleglise and Rivat for more information.
    """
    if a < 2:
        a = 1

    if b is None:
        b, a = a, 1

    if b <= a:
        return

    block_size = integer_sqrt(b)
    prime_list = primes(block_size)

    for start in xrange(a, b, block_size):
        block = [1]*block_size
        values = [1]*block_size

        for p in prime_list:
            offset = ((p*(start//p + 1) - a) % block_size) % p
            for i in xrange(offset, block_size, p):
                if (start + i) % (p*p) == 0:
                    block[i] = 0
                else:
                    block[i] *= -1
                    values[i] *= p

        for i in xrange(block_size):
            if start + i >= b:
                return
            if block[i] and values[i] < start + i:
                block[i] *= -1
            yield (start + i, block[i])


def factored_xrange(a, b=None):
    """Returns an iterator over the factorizations of the numbers in [a, b).

    Given positive integers a and b with a < b, this returns an iterator over
    all pairs (n, n_factorization) with a <= n < b, and n_factorization is the
    factorization of n into prime powers. If the optional parameter b is None,
    then a is taken to be 1, and b = a.

    Input:
        * a: int (a > 0)
        * b: int (b > a) (default=None)

    Returns:
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

    Details:
        All primes <= sqrt(b) are computed, and a segmented sieve is used to
        construct the factorizations of the integers in the interval [a, b).
    """
    if a < 2:
        a = 1

    if b is None:
        b, a = a, 1

    if b <= a:
        return

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
