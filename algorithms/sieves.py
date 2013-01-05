# Sieves

from itertools import islice, count, izip

def erat2():
    D = {}
    yield 2
    for q in islice(count(3), 0, None, 2):
        p = D.pop(q, None)
        if p is None:
            D[q * q] = q
            yield q
        else:
            x = p + q
            while x in D or not (x & 1):
                x += p
            D[x] = p

def primes(n):
    """
    primes(n):
    This returns a list of primes <= n.
    """
    n += 1
    correction = (n%6 > 1)

    n = { 0:n,
          1:n - 1,
          2:n + 4,
          3:n + 3,
          4:n + 2,
          5:n + 1}[n % 6]

    sieve = [True] * (n // 3)
    sieve[0] = False
    sr = integer_sqrt(n)

    for i in xrange(sr // 3 + 1):
        if sieve[i]:
            k = (3*i + 1)|1
            sieve[(k*k) // 3::2*k] \
                    = [False]*((n // 6 - (k*k) // 6 - 1) // k + 1)
            sieve[(k*k + 4*k - 2*k*(i&1)) // 3::2*k] \
                    = [False]*((n // 6 - (k*k + 4*k - 2*k*(i&1)) // 6 - 1) // k + 1)

    return [2, 3] + [ (3*i + 1)|1 for i in xrange(1, n // 3 - correction) if sieve[i] ]

def eratosthenes2(n):
    """
    eratosthenes(n):
    This uses the sieve of Eratosthenes to generate a list of all prime <= n.
    """
    m = n // 2 + 1
    b_list = [True] * m
    N = integer_sqrt(m) + 1

    for i in xrange(1, N):
        if b_list[i]:
            k = 2*i + 1
            start = (k*k - 1) // 2
            how_many = (m - start) // k + ((m - start) % k > 0)
            b_list[start::k] = [False] * how_many

    p_list = [2] + [ 2*i + 1 for i in xrange(1, n // 2 + n%2) if b_list[i] ]
    return p_list

def eratosthenes(n):
    """
    eratosthenes(n):
    This uses the sieve of Eratosthenes to generate a list of all prime <= n.
    """
    b_list = [True] * (n + 1)
    b_list[0] = False
    b_list[1] = False
    b_list[4::2] = [False] * (n // 2 - 1)
    N = integer_sqrt(n) + 1

    for i in xrange(3, N, 2):
        if b_list[i]:
            b_list[i*i::i] = [False] * (n // i + 1 - i)

    p_list = [ i for i in xrange(n + 1) if b_list[i] ]
    return p_list

def segmented_sieve(n):
    """
    segmented_sieve(n, d):
    Given positive integers n and d, this finds all primes < n by breaking up
    the interval from 1 to n into ceil(n / d) intervals of size d and sieving
    them separately.
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
    Returns an iterator over all primes p, a <= p <= b.

    This algorithm finds all primes in an interval [a, b].

    Input:
        * "a" - A positive integer
        * "b" - A positive integer

    Output:
        * "X" - An iterator over all primes in [a, b].

    Details:
        This algorithm is from Section 3.2.2 of "Prime Numbers - A Computational
        Perspective" by Crandall and Pomerance. They call it the "Practical
        Eratosthenes sieve".

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
    block_size = largest_divisor(b - a, max_size)

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

    This algorithm counts the number of primes in the interval [a, b].

    Input:
        * "a" - A positive integer
        * "b" - A positive integer

    Output:
        * "c" - The number of primes in the interval [a, b].

    Details:
        This algorithm is from Section 3.2.2 of "Prime Numbers - A Computational
        Perspective" by Crandall and Pomerance. They call it the "Practical
        Eratosthenes sieve".

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
    block_size = largest_divisor(b - a, max_size)

    q_list = [ -(a + 1 + p)//2 % p for p in p_list ]

    count = 0
    for t in xrange(a, b, 2*block_size):
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

