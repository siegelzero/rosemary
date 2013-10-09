# Miscellaneous combinatorial functions
# Kenneth Brown, 3/06/11

from rosemary.number_theory.arithmetic_functions import factorial
from fractions import Fraction

def bell_number(n):
    """
    bell_number(n):
    Returns the Bell number B_n; i.e. the number of partitions of an n-set.
    """
    B = [1] + [ 0 for k in range(n) ]
    for i in range(1, n + 1):
        B[i] = sum([ binomial(i - 1, k - 1) * B[i - k] for k in range(1, i + 1) ])
    return B[-1]

def bernoulli_number(n):
    """
    This returns the nth bernoulli number, using the algorithm due to
    akiyama-tanigawa.
    """
    if n == 0:
        return Fraction(1, 1)
    elif n == 1:
        return Fraction(-1, 2)
    elif n == 2:
        return Fraction(1,6)
    elif n % 2 == 1:
        return Fraction(0, 1)

    A = [0] * (n + 1)
    for m in xrange(0, n + 1):
        A[m] = Fraction(1, m + 1)
        for j in xrange(m, 0, -1):
            A[j - 1] = j * (A[j - 1] - A[j])
    return A[0]

def bernoulli_list(n):
    """
    This returns a list of of bernoulli numbers B_0, B_1, B_2, ..., B_n
    """
    b_list = [0] * (n + 1)
    b_list[0] = Fraction(1, 1)
    b_list[1] = Fraction(-1, 2)

    for m in xrange(2, n + 1, 2):
        ss = sum([ binomial(m + 1, k) * b_list[k] for k in range(m) ])
        b_list[m] = - ss / (m + 1)
    return b_list

def binomial(n, k):
    """
    Returns the binomial coefficient.

    Input:
        n: a nonnegative integer
        k: a nonnegative integer

    Returns:
        The binomial coefficient n choose k; i.e. the number of ways to choose k
        elements from a set of n elements.

    Examples:
        >>> binomial(5, 2)
        10
        >>> binomial(2, 0)
        1
    """
    pp = 1
    k = min(k, n - k)
    for e in xrange(n - k + 1, n + 1):
        pp = pp * e
    for e in xrange(2, k + 1):
        pp = pp / e
    return pp

def multinomial(n, k):
    """
    Returns the multinomial coefficient.

    Input:
        n: A positive integer
        k: A list or tuple of positive integers
    Returns:
        This returns the multinomial coefficient (n, k1, k2, ..., km) defined by
        n! / (k1! * k2! * ... * km!), where k1 + k2 + ... + km == n. If sum(k)
        != n, then the other elements of k are all assumed to be 1.

    Examples:
        >>> multinomial(6, [3, 3])
        20
        >>> multinomial(6, [2, 2, 2])
        90
    """
    if sum(k) > n:
        raise ValueError('Too many ks')
    elif len(k) == 0:
        raise ValueError('Not enough ks')

    num = factorial(n)
    den = 1
    for e in k:
        den *= factorial(e)
    return num / den

def combinations(A, t):
    """
    lex_combinations_set(A, t):
    Given a sequence A of 
    """
    n = len(A)
    c = [ 0 ] + [ j - 1 for j in range(1, t+1) ] + [ n, 0 ]

    while True:
        yield [ A[i] for i in c[1:t+1] ]

        j = 1
        while c[j] + 1 == c[j + 1]:
            c[j] = j - 1
            j += 1
        if j > t:
            return
        c[j] += 1

def permutations(A):
    """
    permutations(A):
    Given a sequence A of n elements a_1, a_2, ..., a_n, initially sorted so
    a_1 <= a_2 <= ... <= a_n, this algorithm generates all permutations of
    {a_1, a_2, ..., a_n} in lexicographic order.

    This algorithm comes from knuth 4 (Algorithm L)

    Examples:
    >>> L = [1, 2, 2, 3]
    >>> P = list(permutations(L))
    >>> P
    [[1, 2, 2, 3],
     [1, 2, 3, 2],
     [1, 3, 2, 2],
     [2, 1, 2, 3],
     [2, 1, 3, 2],
     [2, 2, 1, 3],
     [2, 2, 3, 1],
     [2, 3, 1, 2],
     [2, 3, 2, 1],
     [3, 1, 2, 2],
     [3, 2, 1, 2],
     [3, 2, 2, 1]]
    """
    B = sorted(A)
    n = len(B) - 1

    while True:
        yield B[:]

        # Find the largest index j such that A[j] < A[j + 1]
        j = n - 1
        while j >= 0 and B[j] >= B[j + 1]:
            j -= 1

        if j < 0:
            return

        # Find the largest index l such that A[j] < A[l]. Then j < l.
        l = n
        while B[j] >= B[l]:
            l -= 1

        # Swap A[j] with A[l]
        B[j], B[l] = B[l], B[j]

        # Reverse the sequence from A[j + 1] to the end A[n]
        l = n
        k = j + 1
        while k < l:
            B[k], B[l] = B[l], B[k]
            k += 1
            l -= 1

def n_tuples(n, M = None):
    """
    n_tuples(n M):
    This algorithm generates all n-tuples a_0, a_1, ..., a_{n - 1} that satisfy
    0 <= a_j < M[j] for 0 <= j <= n - 1. If M is None, this returns all 0,1-tuples
    of length n.

    Examples:
    >>> L = n_tuples(2)
    >>> list(L)
    [[0, 0], [0, 1], [1, 0], [1, 1]]
    >>> L = n_tuples(3, [3, 2, 2])
    >>> list(L)
    [[0, 0, 0],
     [0, 0, 1],
     [0, 1, 0],
     [0, 1, 1],
     [1, 0, 0],
     [1, 0, 1],
     [1, 1, 0],
     [1, 1, 1],
     [2, 0, 0],
     [2, 0, 1],
     [2, 1, 0],
     [2, 1, 1]]
    >>> L = n_tuples(2, 3)
    >>> list(L)
    [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
    """
    A = [ 0 for j in range(n + 1) ]

    if isinstance(M, list):
        if len(M) == n:
            M = [ 2 ] + M
    else:
        if isinstance(M, (int, long)):
            bd = M
        else:
            bd = 2
        M = [ 2 ] + [ bd for j in range(n) ]

    while True:
        yield A[1:]
        j = n
        while A[j] == M[j] - 1:
            A[j] = 0
            j -= 1
        if j == 0:
            return
        else:
            A[j] = A[j] + 1

