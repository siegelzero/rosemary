# Algorithms related to integer factorization

import collections
import itertools
from random import randint
from math import exp, log, sqrt

import rosemary.number_theory.primes.sieves as sieves
import rosemary.number_theory.primes.primality as primality

from rosemary.number_theory.prime_list import _PRIME_LIST
from rosemary.number_theory.core import (
    gcd,
    integer_nth_root,
    integer_sqrt,
    is_power,
    jacobi_symbol,
)


###########################################################################
# Basic factorization algorithms
###########################################################################


def trial_division(n, bound=None):
    r"""Returns the smallest prime divisor of `n`, or proves primality.

    Given an integer `n` > 1, this algorithm attempts to find a factor of
    `n` by using trial division.

    Parameters
    ----------
    n : int (n > 1)

    bound : int, optional (default=None)
        If a non-None value is passed, then only trial divisors <= `bound`
        are considered. Otherwise, all trial divisors up the the square
        root of `n` are considered.

    Returns
    -------
    d : int
        The smallest prime factor of `n`, or `n` itself if `n` is prime.
        Note that if `bound` less than the square root of `n` is passed,
        then there are no guarantees of the primality of `n`.

    Notes
    -----
    This method uses trial division to find the smallest prime factor of n,
    requiring O(sqrt(n)) arithmetic operations. The algorithm runs through
    the precomputed list of primes first, and then uses a simple modulo 30
    wheel after that. This algorithm is very useful to find small prime
    factors of a number, but serves as a poor primality test for numbers
    more than 12 digits.

    References
    ----------
    .. [1] H. Cohen, "A Course in Computational Algebraic Number Theory",
    Springer-Verlag, New York, 2000.

    .. [2] R. Crandall, C. Pomerance, "Prime Numbers: A Computational
    Perspective", Springer-Verlag, New York, 2001.

    .. [3] V. Shoup, "A Computational Introduction to Number Theory and
    Algebra", Cambridge University Press, New York, 2009.

    Examples
    --------
    >>> trial_division(100)
    2
    >>> trial_division(10000004400000259)
    100000007
    >>> trial_division(10000004400000259, 100)
    10000004400000259
    """
    if bound is None:
        bound = integer_sqrt(n) + 1

    # First trial divide by the primes in our stored list.
    for d in _PRIME_LIST:
        if d > bound:
            return n
        if n % d == 0:
            return d

    # Next we use a wheel for the rest. All primes >= 7 fall into one of
    # eight residue classes modulo 30. This uses that fact to avoid trial
    # dividing by numbers divisible by 2, 3, or 5.
    offset = {1: 6, 7: 4, 11: 2, 13: 4, 17: 2, 19: 4, 23: 6, 29: 2}
    d += offset[d % 30]

    while d <= bound:
        if n % d == 0:
            return d
        d += offset[d % 30]

    return n


def fermat(n):
    r"""Returns a nontrivial divisor of `n`, or proves primality.

    Given an integer `n` > 1, this algorithm attempts to find a factor of
    `n` using Fermat's method.

    Parameters
    ----------
    n : int (n > 1)

    Returns
    -------
    d : int
        If n == d, then n is proven prime. Otherwise, d is a nontrivial
        divisor of n.

    Notes
    -----
    This algorithm runs in time O((n + 1) / 2 - sqrt(n)), where n is the
    number to be factored. In the worst case, n = 2*p, this algorithm is
    much worse than trial division. Note that the worse case for Fermat's
    method is the best case for trial division.

    References
    ----------
    .. [1] R. Crandall, C. Pomerance, "Prime Numbers: A Computational
    Perspective", Springer-Verlag, New York, 2001.

    .. [2] D.E. Knuth, "The Art of Computer Programming, Volume 2:
    Seminumerical Algorithms", Addison-Wesley Longman Publishing Co., Inc,
    Boston, MA, 1997.

    Examples
    --------
    >>> fermat(100)
    2
    >>> fermat(1112470797641561909)
    1052788969
    """
    if n % 2 == 0:
        return 2

    for a in xrange(integer_sqrt(n) + 1, (n + 9)//6 + 1):
        t = a*a - n
        b = integer_sqrt(t)
        if b*b == t:
            return a - b

    return n


def lehman(n):
    r"""Returns a nontrivial divisor of `n`, or proves primality.

    Given an integer `n` > 1, this algorithm attempts to find a factor of
    `n` using Lehman's method.

    Parameters
    ----------
    n : int (n > 1)

    Returns
    -------
    d : int
        If d == n, then n is proven prime. Otherwise, d is a nontrivial
        divisor of n.

    Notes
    -----
    This algorithm runs in time O(n^(1/3)). This is substantially better
    than O(n^(1/2)) trial division for reasonably small value of n, but
    this algorithm is not suited for large values of n. See section 8.4 of
    [1] or section 5.1.2 of [2] for more details.

    References
    ----------
    .. [1] H. Cohen, "A Course in Computational Algebraic Number Theory",
    Springer-Verlag, New York, 2000.

    .. [2] R. Crandall, C. Pomerance, "Prime Numbers: A Computational
    Perspective", Springer-Verlag, New York, 2001.

    Examples
    --------
    >>> lehman(100)
    2
    >>> lehman(1112470797641561909)
    1056689261L
    """
    # first, we trial divide up to floor(n^(1/3))
    bound = integer_nth_root(3, n)
    d = trial_division(n, bound)
    if d < n:
        return d

    for k in xrange(1, bound + 1):
        if k % 2 == 0:
            r = 1
            m = 2
        else:
            r = k + n
            m = 4
        # we want to iterate over a, where 4*k*n <= a^2 <= 4*k*n + bound^2
        # and a = r (mod m)
        fkn = 4*k*n
        a = integer_sqrt(fkn)
        # now, increase a until a = r (mod m)
        rm = r % m
        while a % m != rm:
            a += 1
        ub = fkn + bound**2
        while fkn <= a*a <= ub:
            c = a*a - fkn
            b = integer_sqrt(c)
            if b*b == c:
                return gcd(a + b, n)
            a += m
    return n


def pollard_p_minus_1(n, limit=100000):
    r"""Attempts to find a nontrivial factor of `n`.

    Given a composite number `n` and search bound `limit`, this algorithm
    attempts to find a nontrivial factor of `n` using Pollard's p - 1
    method.

    Parameters
    ----------
    n : int (n > 1)
        Integer to be factored.

    limit : int (limit > 0) (default=100000)
        Search limit for possible prime divisors of p - 1.

    Returns
    -------
    d : int
        Divisor of n.

    Notes
    -----
    The algorithm used here is a one-stage form of Pollard's p - 1 method.
    The idea used is that if p - 1 divides Q, then p divides a^Q - 1 if (a,
    p) == 1. So if p is a prime factor of an integer n, then p divides
    gcd(a^Q - 1, n). The idea here is to choose values Q with many divisors
    of the form p - 1, and so search for many primes p as possible divisors
    at once. This is not a general purpose factoring algorithm, but this
    technique works well when n has a prime factor P such that P - 1 is
    divisible by only small primes. See chapter 6 of [2] for details.

    References
    ----------
    .. [1] R. Crandall, C. Pomerance, "Prime Numbers: A Computational
    Perspective", Springer-Verlag, New York, 2001.

    .. [2] H. Riesel, "Prime Numbers and Computer Methods for
    Factorization", 2nd edition, Birkhauser Verlag, Basel, Switzerland,
    1994.

    Examples
    --------
    >>> pollard_p_minus_1(1112470797641561909)
    1056689261L
    """
    p_list = sieves.primes(limit)
    pow_list = []

    for p in p_list:
        pk = p
        while pk <= limit:
            pow_list.append((pk, p))
            pk *= p

    pow_list.sort()
    c = 13

    for (count, (pk, p)) in enumerate(pow_list):
        c = pow(c, p, n)
        if count % 100 == 0:
            g = gcd(c - 1, n)
            if 1 < g < n:
                return g

    return gcd(c - 1, n)


def pollard_rho(n, iterations=None):
    r"""Attempts to find a nontrivial factor of `n`.

    Given a composite number `n`, this algorithm attempts to find a
    nontrivial factor of `n` using Pollard's rho algorithm.

    Parameters
    ----------
    n : int (n > 1)
        Integer to be factored.

    iterations : int, optional (default=None)
        Maximum number of iterations. If None, then the algorithm runs
        until a factor is found.

    Returns
    -------
    d : int
        Factor of `n`.

    Notes
    -----
    This method should return a nontrivial factor of n in O(sqrt(p)) steps,
    where p is the least prime factor of n. Because of this dependence on
    the smallest prime dividing n and not n itself, this method is
    especially useful for large composites with small prime factors out of
    range of trial division. For more details, see section 5.2.1 of [2] and
    section 8.5 of [1].

    References
    ----------
    .. [1] H. Cohen, "A Course in Computational Algebraic Number Theory",
    Springer-Verlag, New York, 2000.

    .. [2] R. Crandall, C. Pomerance, "Prime Numbers: A Computational
    Perspective", Springer-Verlag, New York, 2001.

    Examples
    --------
    >>> pollard_rho(1112470797641561909)
    1052788969L
    >>> pollard_rho(2175282241519502424792841)
    513741730823L
    """
    # Instead of computing a gcd in each iteration, we accumulate the
    # products and take the gcd only when the number of terms is a multiple
    # of `step`
    step = 100
    prod = 1
    a = randint(1, n - 3)
    u = randint(0, n - 1)
    v = u
    d = 1

    if iterations is None:
        X = itertools.count()
    else:
        X = xrange(iterations)

    for k in X:
        u = (u*u + a) % n
        v = (v*v + a) % n
        v = (v*v + a) % n
        prod = prod*(u - v) % n

        if k % step == 0:
            d = gcd(prod, n)
            prod = 1
            if d > 1:
                return d

    return n


def pollard_rho_brent(n):
    r"""Attempts to find a nontrivial factor of `n`.

    Given a composite number `n`, this algorithm attempts to find a
    nontrivial factor of `n` using Brent's improvement to the Pollard-rho
    algorithm.

    Parameters
    ----------
    n : int (n > 1)
        Integer to factor.

    Returns
    -------
    d : int
        Factor of `n`.

    Notes
    -----
    This method should return a nontrivial factor of n in O(sqrt(p)) steps,
    where p is the least prime factor of n. Because of this dependence on
    the smallest prime dividing n and not n itself, this method is
    especially useful for large composites with small prime factors out of
    range of trial division. For more details, see [1] and section 8.5 of
    [2].

    References
    ----------
    .. [1] R.P. Brent, "An Improved Monte Carlo Factorization Algorithm",
    BIT, Vol. 20, 1980.

    .. [2] H. Cohen, "A Course in Computational Algebraic Number Theory",
    Springer-Verlag, New York, 2000.

    .. [3] J.M. Pollard, "A Monte Carlo Method for Factorization", BIT,
    Vol. 15, 1975.

    Examples
    --------
    >>> pollard_rho_brent(1112470797641561909)
    1052788969L
    >>> pollard_rho_brent(2175282241519502424792841)
    513741730823L
    """
    c = randint(1, n - 3)
    y = randint(0, n - 1)
    step = 100
    prod = 1
    r = 1

    while True:
        x = y
        for i in xrange(r):
            y = (y*y + c) % n

        for k in itertools.count(0, step):
            ys = y
            for i in xrange(min(step, r - k)):
                y = (y*y + c) % n
                prod = prod*(x - y) % n
            g = gcd(prod, n)
            if k >= r or g > 1:
                break

        r *= 2
        if g > 1:
            break

    if g == n:
        while True:
            ys = (ys*ys + c) % n
            g = gcd(x - ys, n)
            if g > 1:
                break
    return g


def cfrac(n, k=None):
    r"""Returns a nontrivial divisor of `n`.

    Given a composite integer `n` > 1, this algorithm attempts to find a
    factor of `n` using Morrison and Brillhart's continued fraction method
    CFRAC.

    Parameters
    ----------
    n : int (n > 1)
        Number to be factored.

    k : int, optional (default=None)
        Multiplier to use in case the period of sqrt(n) is too short.

    Returns
    -------
    d : int
        Divisor of n.

    Notes
    -----
    Morrison and Brillhart's continued fraction method was the first
    factorization algorithm of subexponential running time. The idea is to
    find a nontrivial solution to the congrunce x^2 = y^2 (mod n), and
    extracting a factor from gcd(x + y, n). See [1] for details, and see
    [2] for implementation details.

    References
    ----------
    .. [1] H. Cohen, "A Course in Computational Algebraic Number Theory",
    Springer-Verlag, New York, 2000.

    .. [2] M.A. Morrison, J. Brillhart, "A Method of Factoring and the
    Factorization of F7", Mathematics of Computation, Vol. 29, Num. 129,
    Jan. 1975.

    .. [3] H. Riesel, "Prime Numbers and Computer Methods for
    Factorization", 2nd edition, Birkhauser Verlag, Basel, Switzerland,
    1994.

    Examples
    --------
    >>> cfrac(12007001)
    4001
    >>> cfrac(1112470797641561909)
    1052788969L
    >>> cfrac(2175282241519502424792841)
    513741730823L
    """
    # B is our smoothness bound.
    B = int(exp(0.5*sqrt(log(n)*log(log(n))))) + 1
    prime_list = sieves.primes(B)

    # Choose a multiplier if none is provided.
    if k is None:
        k = _cfrac_multiplier(n)
    kn = k*n

    # Perform simple trial division by the primes we computed to find any
    # small prime divisors.
    for p in prime_list:
        if n % p == 0:
            return p

    # Our factor base needs to include -1 and 2, and the odd primes p for
    # which (kN|p) = 0 or 1.
    factor_base = [-1]
    for p in prime_list:
        if p == 2 or jacobi_symbol(kn, p) >= 0:
            factor_base.append(p)

    num_primes = len(factor_base)

    # Compute the product of the elements in our factor base for smoothness
    # checking computations later.
    prod = 1
    for p in factor_base:
        prod *= p

    # Instead of using trial division to check each value for smoothness
    # individually, we use a batch smoothness test, processing batches of
    # size `batch_size` at once.

    # Set e as the least positive integer with n <= 2**e.
    e = 1
    while 2**e < n:
        e *= 2

    aq_pairs = _cfrac_aq_pairs(kn)

    num_smooths_found = 0
    exponent_matrix = []
    a_list = []

    while num_smooths_found <= num_primes:
        (i, a, q) = aq_pairs.next()

        # This is from the batch smoothness test given as Algorithm 3.3.1
        # in [1].
        if gcd(q, pow(prod, e, q)) != q:
            continue

        if i % 2 == 1:
            q *= -1

        # At this point, we know q is smooth, and we can factor it
        # completely using our factor base.
        exponent_vector = smooth_factor(q, factor_base)
        exponent_matrix.append(exponent_vector)
        num_smooths_found += 1
        a_list.append(a)

    kernel = _z2_gaussian_elimination(exponent_matrix)

    for i in xrange(len(kernel)):
        y = 1
        x2_exponents = [0]*num_primes
        for j in xrange(len(kernel[i])):
            if kernel[i][j]:
                y = (a_list[j]*y) % n
                for f in xrange(num_primes):
                    x2_exponents[f] += exponent_matrix[j][f]

        x = 1
        for j in xrange(num_primes):
            x *= factor_base[j]**(x2_exponents[j]//2)

        for val in [x - y, x + y]:
            d = gcd(val, n)
            if 1 < d < n:
                return d


def _z2_gaussian_elimination(exponents):
    r"""Returns the zero linear combinations of the given vectors.

    Given a list of vectors, this returns the linear combinations of the
    vectors that sum to zero (modulo 2).

    Parameters
    ----------
    exponents : list
        List of vectors, each vector given as a list.

    Returns
    -------
    relations : list
        List of lists, each sublist has the coefficients of a linear
        combination of the input vectors that sums to zero modulo 2.

    Notes
    -----
    This is a subroutine useful in several factorization algorithms. Our
    implementation is basically straightforward Gaussian elimination with a
    few modifications, as outlined in [1].

    References
    ----------
    .. [1] M.A. Morrison, J. Brillhart, "A Method of Factoring and the
    Factorization of F7", Mathematics of Computation, Vol. 29, Num. 129,
    Jan. 1975.

    Examples
    --------
    >>> vectors = [
        [1, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0],
        [0, 0, 1, 0, 0, 1, 0],
        [1, 1, 0, 0, 0, 1, 0],
        [0, 1, 0, 1, 0, 0, 1],
        [1, 1, 0, 0, 0, 1, 0],
        [0, 1, 0, 0, 1, 0, 1]
    ]
    >>> _z2_gaussian_elimination(vectors)
    [[1, 0, 1, 1, 0, 0, 0], [1, 0, 1, 0, 0, 1, 0], [0, 1, 0, 0, 1, 0, 1]]
    """
    # Arithmetic done over Z2, so we reduce the exponent vectors modulo 2.
    num_rows = len(exponents)
    num_cols = len(exponents[0])

    reduced = []
    for vector in exponents:
        value = 0
        for i in xrange(num_cols):
            if vector[i] % 2 != 0:
                value += 1 << (num_cols - i - 1)
        reduced.append(value)

    # history is the identity matrix with 'num_rows' rows and columns.
    history = []
    value = 1 << (num_rows - 1)
    for i in xrange(num_rows):
        history.append(value)
        value >>= 1

    # Now we have the exponent vectors encoded in bits, so we can use
    # efficient bit operations.

    # Starting in the rightmost column, find the first vector whose
    # rightmost 1 is in column j.
    bit = 1
    max_size = 1 << (num_cols - 1)
    while bit <= max_size:
        pivot_found = False
        for i in xrange(num_rows):
            entry = reduced[i]
            if entry & bit and entry % bit == 0:
                pivot_found = True
                his = history[i]
                break
        # If we haven't found a pivot yet, move to the next column left and
        # look again.
        if pivot_found:
            for m in xrange(i + 1, num_rows):
                if reduced[m] & bit:
                    reduced[m] ^= entry
                    history[m] ^= his

        bit <<= 1

    # Now convert the binary encoded vectors back to exponent lists.
    vectors = []
    for i in xrange(num_rows):
        value = history[i]
        if reduced[i] != 0 or value == 0:
            continue
        vector = [0]*num_rows
        j = num_rows - 1
        while j >= 0:
            if value % 2:
                vector[j] = 1
            value >>= 1
            j -= 1

        vectors.append(vector)

    return vectors


def smooth_factor(n, factor_base):
    r"""Factors `n` over the factor base.

    Returns the exponents in the prime factorization of `n` if `n` factors
    completely using the primes in 'factor_base'.

    Parameters
    ----------
    n : int
        Integer to factor.

    factor_base : list
        A sorted list of primes, possibly including -1.

    Returns
    -------
    exponents: list
        If `n` factors completely over `factor_base`, then this is a list
        of the exponents appearing in the prime factorization; exponents[i]
        holds the power to which factor_base[i] divides `n`.

    Raises
    ------
    ValueError
        If `n` does not factor over `factor_base`.

    Notes
    -----
    This is a straightforward implementation of factorization using trial
    division, using only the primes in the given factor base.

    References
    ----------
    .. [1] R. Crandall, C. Pomerance, "Prime Numbers: A Computational
    Perspective", Springer-Verlag, New York, 2001.

    Examples
    --------
    >>> smooth_factor(100, [2, 5])
    [2, 2]
    >>> smooth_factor(100, [2, 3, 5, 7])
    [2, 0, 2, 0]
    >>> smooth_factor(100, [3, 5])
    ValueError: smooth_factor: n does not factor over the given factor base
    """
    num_primes = len(factor_base)
    exponents = [0]*num_primes
    start = 0

    if factor_base[0] == -1:
        start = 1
        if n < 0:
            exponents[0] = 1
            n *= -1

    for i in xrange(start, num_primes):
        p = factor_base[i]
        if n % p == 0:
            e = 1
            n = n//p
            while n % p == 0:
                n = n//p
                e += 1
            exponents[i] = e
        if n == 1:
            return exponents

    raise ValueError("smooth_factor: n does not factor over the given factor base")


def _cfrac_aq_pairs(n):
    r"""Yields triples (i, A_{i - 1}, Q_i) for i > 0, where A_{i - 1}^2 =
    (-1)^i Q_i (mod n)

    Parameters
    ----------
    n : int

    Yields
    ------
    (i, A_{i - 1}, Q_i) : tuple
        Yields integer tuples consisting of the index `i`, the numerator of
        the (i - 1)th convergent A({i - 1} (mod n), and the quantity Q_{i}
        coming from the continued fraction expansion of sqrt(n).

    Notes
    -----
    The continued fraction expansion of sqrt(n) is computed, using the
    formulas given in section 2 of [1]. See also appendix 8 of [2].

    References
    ----------
    .. [1] M.A. Morrison, J. Brillhart, "A Method of Factoring and the
    Factorization of F7", Mathematics of Computation, Vol. 29, Num. 129,
    Jan. 1975.

    .. [2] H. Riesel, "Prime Numbers and Computer Methods for
    Factorization", 2nd edition, Birkhauser Verlag, Basel, Switzerland,
    1994.

    Examples
    --------
    >>> X = _cfrac_aq_pairs(1000009)
    >>> X.next()
    (1, 1000, 9)
    >>> X.next()
    (2, 222001, 445)
    >>> X.next()
    (3, 889004, 873)
    >>> X.next()
    (4, 1000000, 81)
    """
    g = integer_sqrt(n)
    A0, A1 = 0, 1
    Q0, Q1 = n, 1
    P0 = 0
    r0 = g

    for i in itertools.count():
        q = (g + P0)//Q1
        r1 = g + P0 - q*Q1
        A2 = (q*A1 + A0) % n
        P1 = g - r1
        Q2 = Q0 + q*(r1 - r0)

        if i > 0:
            yield (i, A1, Q1)

        A0, A1 = A1, A2
        Q0, Q1 = Q1, Q2
        P0 = P1
        r0 = r1


def _cfrac_multiplier(n):
    r"""Computes a multiplier for input into the cfrac algorithm.

    Parameters
    ----------
    n : int
        Integer to be factored.

    Returns
    -------
    k : int
        Multiplier for the cfrac algorithm.

    Notes
    -----
    The continued fraction of sqrt(n) is always periodic. In the cases
    where the period of sqrt(n) is too short, it is necessary to expand
    sqrt(k*n) for some k > 1. We choose the multiplier `k` which allows th
    elargest number of primes <= 31 to be in the factor base. The method
    follows Remark 5.3 of [1].

    References
    ----------
    .. [1] M.A. Morrison, J. Brillhart, "A Method of Factoring and the
    Factorization of F7", Mathematics of Computation, Vol. 29, Num. 129,
    Jan. 1975.

    Examples
    --------
    >>> _cfrac_multiplier(5**77 - 1)
    781
    """
    prime_list = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    choices = {}
    # Look for multpliers in the range [1, 1000).
    for k in xrange(1, 1000):
        if jacobi_symbol(k*n, 3) >= 0 and jacobi_symbol(k*n, 5) >= 0:
            # Find the multiplier k that allows the largest number of
            # primes <= 31 into the factor base.
            count = sum(1 for p in prime_list if jacobi_symbol(k*n, p) >= 0)
            if count not in choices:
                choices[count] = [k]
            else:
                choices[count].append(k)
    # If several values of k allow this maximal number, we simply choose the
    # smallest of them.
    max_count = max(choices)
    return min(choices[max_count])



