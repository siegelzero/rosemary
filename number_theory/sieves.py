# Sieves

from rosemary.number_theory.core import integer_sqrt
import rosemary.number_theory.misc

def primes(n):
    """
    Returns a list of all primes <= n.

    This program uses the sieve of Eratosthenes to generate a list of all primes <= n.

    Input:
        * n: int (n > 0)

    Output:
        * p_list: list
            A list of the primes <= n.

    Details:
        This highly optimized sieve only looks at the numbers <= n that are coprime to 6. The code is based on some
        found on stackoverflow.

    Examples:
        >>> primes(100)
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
        >>> len(primes(10**7))
        664579
    """
    n += 1
    offset = (n%6 > 1)
    n = {0:n, 1:n - 1, 2:n + 4, 3:n + 3, 4:n + 2, 5:n + 1}[n % 6]

    sieve = [True] * (n // 3)
    sieve[0] = False
    sr = integer_sqrt(n)

    for i in xrange(sr // 3 + 1):
        if sieve[i]:
            k = (3*i + 1)|1
            kk = k*k
            sieve[kk // 3::2*k] = [False]*((n // 6 - kk // 6 - 1) // k + 1)
            sieve[(kk + 4*k - 2*k*(i&1)) // 3::2*k] = [False]*((n // 6 - (kk + 4*k - 2*k*(i&1)) // 6 - 1) // k + 1)

    return [2, 3] + [(3*i + 1)|1 for i in xrange(1, n // 3 - offset) if sieve[i]]


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


def segmented_sieve(n):
    """
    Returns an iterator over all primes <= n.

    Given a positive integer n, this finds all primes <= n by breaking up the
    interval from 1 to n into intervals of size sqrt(d) and sieving them
    separately.

    Input:
        * n - A positive integer.

    Output:
        * X - an iterator over all primes <= n.

    Details:
        The algorithm used is a standard segmented sieve. The main improvement
        is that this algorithm uses O(sqrt(n)) space, a significant improvement
        over the O(n) space used in th standard sieve of Eratosthenes.  See
        chapter 9 of "Algorithmic Number Theory I - Efficient Algorithms" by
        Bach and Shallit for details.

    Examples:
        >>> len(segmented_sieve(10**7))
        664579
        >>> X = segmented_sieve(1000)
        >>> len([p for p in X if p % 4 == 1])
        80
    """
    # compute the primes <= sqrt(n)
    delta = integer_sqrt(n)
    initial_primes = primes(delta)
    psquares = [p*p for p in initial_primes]
    r = len(initial_primes)

    for m in xrange(0, n, delta):
        d = [1] * delta
        for i in xrange(r):
            while psquares[i] < delta:
                d[psquares[i]] = 0
                psquares[i] += initial_primes[i]
            psquares[i] -= delta
        for i in xrange(delta):
            if d[i] and 1 < m + i <= n:
                yield m + i

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

    p_list = primes(bd)
    num_primes = len(p_list)
    max_size = min((b - a) // 2, 8*bd)
    block_size = rosemary.number_theory.misc.largest_divisor(b - a, max_size)

    q_list = [ -(a + 1 + p)//2 % p for p in p_list ]

    for t in xrange(a, b, 2*block_size):
        b_list = [1] * block_size
        for k in xrange(1, num_primes):
            pk = p_list[k]
            qk = q_list[k]
            diff = block_size - qk
            how_many = diff // pk + (diff % pk > 0)
            b_list[qk:block_size:pk] = [0] * how_many
            q_list[k] = (-diff) % pk
        for j in xrange(block_size):
            if b_list[j]:
                yield t + 2*j + 1

def number_of_primes_in_interval(a, b):
    """
    Returns the number of primes p, a <= p <= b.

    Given positive integers a and b with a*a > b, this returns the number of
    primes in the interval [a, b].

    Input:
        * a - A positive integer.
        * b - A positive integer.

    Output:
        * c - The number of primes in the interval [a, b].

    Details:
        This algorithm is from Section 3.2.2 of "Prime Numbers - A
        Computational Perspective" by Crandall and Pomerance. They call it the
        "Practical Eratosthenes sieve".

    Examples:
        >>> number_of_primes_in_interval(10**8, 2*10**8)
        5317482
        >>> number_of_primes_in_interval(10**9, 2*10**9)
        47389752
    """
    # We need both endpoints to be even
    if a % 2 == 1:
        a -= 1
    if b % 2 == 1:
        b += 1

    bd = integer_sqrt(b)
    assert a > bd

    p_list = primes(bd)
    num_primes = len(p_list)
    max_size = min((b - a) // 2, 8*bd)
    block_size = rosemary.number_theory.misc.largest_divisor(b - a, max_size)

    q_list = [ -(a + 1 + p)//2 % p for p in p_list ]

    count = 0
    for _ in xrange(a, b, 2*block_size):
        b_list = [1] * block_size
        for k in xrange(1, num_primes):
            pk = p_list[k]
            qk = q_list[k]
            diff = block_size - qk
            how_many = diff // pk + (diff % pk > 0)
            b_list[qk:block_size:pk] = [0] * how_many
            q_list[k] = (-diff) % pk
        count += sum(b_list)
    return count


def prime_xrange(a, b=None):
    if b is None:
        b = a
        a = 2

    block_size = integer_sqrt(b)
    p_list = primes(block_size)

    if a <= p_list[-1]:
        for p in p_list:
            if a <= p <= b:
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



def sieve_interval2(a, b):
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


def prime_xrange2(a, b=None):
    if b is None:
        b = a
        a = 2

    block_size = integer_sqrt(b)
    p_list = primes(block_size)

    if a <= p_list[-1]:
        for p in p_list:
            if a <= p <= b:
                yield p
    
    block_size *= 20

    for start in xrange(block_size + 1, b, block_size):
        block = [1]*block_size
        for p in p_list:
            offset = ((p*(start//p + 1) - 1) % block_size) % p
            count = (block_size - offset)//p + ((block_size - offset) % p > 0)
            block[offset::p] = [0]*count

        for i in xrange(block_size):
            if block[i] and a <= start + i <= b:
                yield start + i


