# Sieves

from bisect import bisect, bisect_left
from heapq import heappush, heappop
from math import log

from rosemary.number_theory.core import integer_sqrt
from rosemary.number_theory.prime_list import _PRIME_LIST


################################################################################
# Lists of primes
################################################################################


def primes(n):
    r"""Returns a list of all primes <= n.

    Parameters
    ----------
    n : int

    Returns
    -------
    primes : list
        A list of the primes <= n.

    See Also
    --------
    chartres, eratosthenes1, eratosthenes2, luo

    Notes
    -----
    For values of n <= 10**6, the precomputed list of primes is used.
    Otherwise, an optimized version of the sieve of Eratosthenes is used.

    Examples
    --------
    >>> primes(20)
    [2, 3, 5, 7, 11, 13, 17, 19, 23]
    >>> len(primes(10**6))
    78498
    """
    if n < _PRIME_LIST[-1]:
        k = bisect(_PRIME_LIST, n)
        return _PRIME_LIST[:k]
    else:
        return eratosthenes2(n)


def primes_first_n(n):
    r"""Returns a list of the first n primes.

    Parameters
    ----------
    n : int

    Returns
    -------
    primes : list
        A list of the first n primes.

    Notes
    -----
    This method uses an upper bound for the nth prime, generates all primes
    to that point, and returns the first n of them. The estimate used is
    p_n < n*(log(n) + log(log(n))) for n >= 6, where p_n denotes the nth
    prime. See [1] for information.

    References
    ----------
    .. [1] Bach, Shallitt, "Algorithmic Number Theory", MIT Press,
    Cambridge, MA, 1996

    Examples
    --------
    >>> primes_first_n(10)
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    >>> sum(primes_first_n(100))
    24133
    """
    if n <= 0:
        return []
    elif n < len(_PRIME_LIST):
        return _PRIME_LIST[:n]
    else:
        m = int(n*(log(n) + log(log(n))))
        return eratosthenes2(m)[:n]


def eratosthenes1(n):
    r"""Returns a list of all primes p <= n.

    Parameters
    ----------
    n : int

    Returns
    -------
    primes : list
        A list of the primes <= n.

    See Also
    --------
    chartres, eratosthenes2, luo

    Notes
    -----
    This algorithm is an implementation of the first extension of the sieve
    of Eratosthenes by removing the multiples of 2 from the set of
    candidates. This has the effect of halving the required space and time
    of the algorithm. See the paper [1] for details. See the commentary [2]
    for more details and extensions.

    All versions of the Sieve of Eratosthenes require O(n*log(log(n)))
    operations to find the primes up to n. Using a simple wheel like this
    cuts down on the required time, but the overall complexity is not
    changed.

    References
    ----------
    .. [1] X. Luo, "A Practical Sieve Algorithm for Finding Prime Numbers",
    Communications of the ACM, Volume 32, Number 3, March 1989.

    .. [2] Quesada, Pritchard, James, "Additional Notes on a Practical Sieve
    Algorithm", Communications of the ACM, Volume 35, Issue 3, March 1992.

    Examples
    --------
    >>> eratosthenes1(20)
    [2, 3, 5, 7, 11, 13, 17, 19]
    >>> len(eratosthenes(10**6))
    78498
    """
    if n <= 1:
        return []

    m = n//2 + 1
    block = bytearray([1])*m
    sqrt = integer_sqrt(m)

    for i in xrange(1, sqrt + 1):
        if block[i]:
            k = 2*i + 1
            start = k*k//2
            count = (m - start)//k + ((m - start) % k > 0)
            block[start::k] = bytearray([0])*count

    return [2] + [2*i + 1 for i in xrange(1, n//2 + n % 2) if block[i]]


def eratosthenes2(n):
    r"""Returns a list of all primes p <= n.

    Parameters
    ----------
    n : int

    Returns
    -------
    primes : list
        A list of the primes <= n.

    See Also
    --------
    chartres, eratosthenes1, luo

    Notes
    -----
    This algorithm is an implementation of the second extension of the
    sieve of Eratosthenes by removing the multiples of 2 and 3 from the set
    of candidates. This has the effect of decreasing the required space and
    time to 1/3 of the naive Eratosthenes algorithm. See the paper [1] for
    details. See the commentary in [2] for more details and extensions.
    Note that [1] and [3] both claim that the second extension of the sieve
    of Eratosthenes is the most efficient in practice.

    References
    ----------
    .. [1] X. Luo, "A Practical Sieve Algorithm for Finding Prime Numbers",
    Communications of the ACM, Volume 32, Number 3, March 1989.

    .. [2] Quesada, Pritchard, James, "Additional Notes on a Practical Sieve
    Algorithm", Communications of the ACM, Volume 35, Issue 3, March 1992.

    .. [3] J. Sorenson, "An Analysis of Two Prime Number Sieves", Technical
    Report 1028, University of Wisonsin, Computer Sciences Department, June
    1991.

    Examples
    --------
    >>> eratosthenes2(20)
    [2, 3, 5, 7, 11, 13, 17, 19]
    >>> len(eratosthenes2(10**6))
    78498
    """
    if n <= 1:
        return []
    elif n <= 3:
        return range(2, n + 1)

    n += 1
    offset = (n % 6 > 1)
    n = {0: n, 1: n - 1, 2: n + 4, 3: n + 3, 4: n + 2, 5: n + 1}[n % 6]

    m = n//3
    block = bytearray([1])*m
    block[0] = 0
    sqrt = integer_sqrt(n)

    for i in xrange(sqrt//3 + 1):
        if block[i]:
            k = (3*i + 1) | 1
            kk = k*k
            a = kk//3
            b = (kk + 4*k - 2*k*(i & 1))//3
            block[a::2*k] = bytearray([0])*((m - a)//(2*k) + ((m - a) % (2*k) > 0))
            block[b::2*k] = bytearray([0])*((m - b)//(2*k) + ((m - b) % (2*k) > 0))

    return [2, 3] + [(3*i + 1) | 1 for i in xrange(1, n//3 - offset) if block[i]]


def luo(n):
    r"""Returns a list of all primes p <= n.

    Parameters
    ----------
    n : int

    Returns
    -------
    primes : list
        A list of the primes <= n.

    See Also
    --------
    chartres, eratosthenes1, eratosthenes2

    Notes
    -----
    This algorithm is an implementation of the second extension of the
    sieve of Eratosthenes by removing the multiples of 2 and 3 from the set
    of candidates. This has the effect of decreasing the required space and
    time to 1/3 of the naive Eratosthenes algorithm. This implementation is
    based on the the exposition in the paper [1], as well as the
    optimizations made by Pritchard in [2]. Note that [1] and [3] both
    claim that the second extension of the sieve of Eratosthenes is the
    most efficient in practice.

    References
    ----------
    .. [1] X. Luo, "A Practical Sieve Algorithm for Finding Prime Numbers",
    Communications of the ACM, Volume 32, Number 3, March 1989.

    .. [2] Quesada, Pritchard, James, "Additional Notes on a Practical Sieve
    Algorithm", Communications of the ACM, Volume 35, Issue 3, March 1992.

    .. [3] J. Sorenson, "An Analysis of Two Prime Number Sieves", Technical
    Report 1028, University of Wisonsin, Computer Sciences Department, June
    1991.

    Examples
    --------
    >>> luo(20)
    [2, 3, 5, 7, 11, 13, 17, 19]
    >>> len(luo(10**6))
    78498
    """
    if n <= 1:
        return []
    elif n <= 6:
        return {2: [2], 3: [2, 3], 4: [2, 3], 5: [2, 3, 5], 6: [2, 3, 5]}[n]

    n += 1
    offset = (n % 6 > 1)
    n = {0: n, 1: n - 1, 2: n + 4, 3: n + 3, 4: n + 2, 5: n + 1}[n % 6]

    m = n//3
    sqrt = integer_sqrt(n)
    block = bytearray([1])*(n//3)
    a, k, t = 0, 1, 2

    for i in xrange(1, sqrt + 1):
        k = 3 - k
        a += 4*k*i
        t += 4*k

        if block[i]:
            b = a + 2*i*(3 - k) + 1
            block[a::t] = bytearray([0])*((m - a)//t + ((m - a) % t > 0))
            block[b::t] = bytearray([0])*((m - b)//t + ((m - b) % t > 0))

    return [2, 3] + [(3*i + 1) | 1 for i in xrange(1, m - offset) if block[i]]


def chartres(n):
    r"""Returns a list of all primes <= n.

    Parameters
    ----------
    n : int

    Returns
    -------
    primes : list
        A list of primes <= n.

    See Also
    --------
    eratosthenes1, eratosthenes2, luo

    Notes
    -----
    This method uses the algorithm in Section 5.2.3, Exercise 15 of [1].
    The algorithm is originally due to Chartres. See [2] for more
    information.

    Unlike the standard implementation of the sieve of Eratosthenes, this
    algorithm uses a heap and requires O(n*log(n)) steps, while the sieve
    of Eratosthenes runs in time O(n*log(log(n))).

    References
    ----------
    .. [1] D.E. Knuth, "The Art of Computer Programming, Volume 3: Sorting
    and Searching", Addison-Wesley Longman Publishing Co., Inc, Redwood
    City, CA, 1998.

    .. [2] B.A. Chartres, "Algorithm 310: Prime Number Generator 1",
    Communications of the ACM, Volume 10, Number 9, Sept. 1967.

    Examples
    --------
    >>> chartres(30)
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    >>> len(chartres(1000))
    168
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


def pritchard(n):
    r"""Returns a list of all primes <= n.

    Parameters
    ----------
    n : int

    Returns
    -------
    primes : list
        A list of the primes <= n.

    See Also
    --------
    chartres, eratosthenes1, eratosthenes2, luo

    Notes
    -----
    This method uses the linear sieve by Pritchard to find the primes up to
    n. Our implementation is based on exposition in the paper [1].  For
    additional details, see the original paper [2].

    Despite asymptotically faster runtime than the Sieve of Eratosthenes,
    this algorithm is slower in practice. See the paper [3] for details.

    References
    ----------
    .. [1] Dunten, Jones, Sorenson, "A Space-Efficient Fast Prime Number
    Sieve", Information Processing Letters Volume 59, Issue 2, 1996.

    .. [2] P. Pritchard, "Linear Prime-Number Sieves: A Family Tree",
    Science of Computer Programming, Volume 9, Issue 1, 1987.

    .. [3] J. Sorenson, "An Analysis of Two Prime Number Sieves", Technical
    Report 1028, University of Wisonsin, Computer Sciences Department, June
    1991.

    Examples
    --------
    >>> pritchard(20)
    [2, 3, 5, 7, 11, 13, 17, 19]
    >>> len(pritchard(10**6))
    78498
    """
    if n <= 1:
        return []

    prime_list = primes(integer_sqrt(n))

    block = [1]*(n + 1)
    block[0] = 0
    block[1] = 0

    for f in xrange(2, n//2 + 1):
        for p in prime_list:
            if p*f > n:
                break

            block[p*f] = 0

            if f % p == 0:
                break

    return [i for i in xrange(2, n + 1) if block[i]]


################################################################################
# Generators of primes
################################################################################


def prime_xrange(a, b=None):
    r"""Returns generator of primes in interval [a, b).

    Parameters
    ----------
    a : int (a > 0)
    b : int (b > a) (default=None)

    Returns
    -------
    primes : generator
        Yields the primes in the interval [a, b).

    Notes
    -----
    This method uses a variety of techniques to efficiently find the primes
    in the given interval. For b <= 10**7, the primes <= b are computed,
    and the primes in [a, b) are yielded. Note that for b <= 10**6, the
    precomptued list of primes is used here. For b > 10**7, a segmented
    sieve of Eratosthenes is used. See Section 3.2.2 of [1] for details.

    References
    ----------
    .. [1] R. Crandall, C. Pomerance, "Prime Numbers: A Computational
    Perspective", Springer-Verlag, New York, 2001.

    Examples
    --------
    >>> list(prime_xrange(10, 20))
    [11, 13, 17, 19]
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
                # offset = -((start - start % 2) + 1 + p)//2 % p
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
