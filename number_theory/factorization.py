# Algorithms related to integer factorization

import collections
import itertools
from random import randint
from math import exp, log, sqrt

import rosemary.number_theory.sieves
import rosemary.number_theory.primality

from rosemary.number_theory.prime_list import PRIME_LIST
from rosemary.number_theory.core import (
    gcd,
    integer_log,
    integer_nth_root,
    integer_sqrt,
    jacobi_symbol,
)

###############################################################################
# Classical algorithms
###############################################################################


def trial_division(n, bound=None):
    """
    Returns the smallest prime divisor of n, or proves primality.

    Given an integer n > 1, this algorithm attempts to find a factor of n by
    using trial division.

    Input:
        * n: int (n > 1)
        * bound: int (default=None)
            The bound for trial division.

    Output:
        * d: int
            The smallest prime factor of n, or n itself if n is prime.

    Examples:
        >>> trial_division(100)
        2
        >>> trial_division(10000004400000259)
        100000007

    Details:
        This method uses trial division to find the smallest prime factor of n.
        The algorithm runs through the precomputed list of primes first, and
        then uses a simple wheel modulo 30 after that. This algorithm is very
        useful to find small prime factors of a number, but serves a poor
        primality test for numbers more than 12 digits.
    """
    if bound is None:
        bound = integer_sqrt(n) + 1

    # First trial divide by the primes in our stored list.
    for d in PRIME_LIST:
        if d > bound:
            return n
        if n % d == 0:
            return d

    # Next we use a wheel for the rest. All primes >= 7 fall into one of eight
    # residue classes modulo 30. This uses that fact to avoid trial dividing by
    # numbers divisible by 2, 3, or 5.
    offset = {1: 6, 7: 4, 11: 2, 13: 4, 17: 2, 19: 4, 23: 6, 29: 2}
    d += offset[d % 30]

    while d <= bound:
        if n % d == 0:
            return d
        d += offset[d % 30]
    return n


def fermat(n):
    """
    Returns a nontrivial divisor of n, or proves primality.

    Given an integer n > 1, this algorithm attempts to find a factor of n using
    Fermat's method.

    Input:
        * n: int (n > 1)

    Output:
        * d: int
            If n == d, then n is proven prime. Otherwise, d is a nontrivial
            divisor of n.

    Examples:
        >>> m = 1112470797641561909
        >>> fermat(m)
        1052788969

    Details:
        The algorithm used here is Fermat's method. This algorithm runs in time
        O((n + 1) / 2 - sqrt(n)), where n is the number to be factored.  In the
        worst case, this algorithm is much worse than trial division.  See
        section 5.1.1 of "Prime Numbers - A Computation Perspective" by Crandall
        and Pomerance for more details.
    """
    if n % 2 == 0:
        return 2

    a = integer_sqrt(n) + 1
    while a <= (n + 9)//6:
        t = a**2 - n
        b = integer_sqrt(t)
        if b*b == t:
            return a - b
        a += 1
    return n


def lehman(n):
    """
    Returns a nontrivial divisor of n, or proves primality.

    Given an integer n > 1, this algorithm finds a nontrivial factor of n if
    n is not prime, or returns n if n is prime.

    Input:
        * n: int (n >= 3)

    Output:
        * d: int
            If d == n, then n is proven prime. Otherwise, d is a nontrivial
            divisor of n.

    Examples:
        >>> l = 1112470797641561909
        >>> lehman(l)
        1056689261L

    Details:
        The algorithm used is Lehman's Method. This algorithm runs in time
        O(n^(1/3)), where n is the number to be factored / proven prime. This is
        substantially better than O(n^(1/2)) trial division for reasonably small
        value of n, but this algorithm is not suited for large values of n. See
        section 8.4 of "A Course in Computational Algebraic Number Theory" by
        Cohen or section 5.1.2 of "Prime Numbers - A Computational Perspective"
        by Crandall and Pomerance for more details.


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


def pollard_p_minus_1(n, B=20000):
    """
    Attempts to find a nontrivial factor of n.

    Given a composite number n and search bound B, this algorithm attempts to
    find a nontrivial factor of n.

    Input:
        * n: int (n > 1)
            A composite integer.

        * B: int (B > 0) (default=20000)
            Search bound.

    Output:
        * d: int
            Nontrivial divisor of n.

    Examples:
        >>> m = 1112470797641561909
        >>> pollard_p_minus_1(m)
        1056689261L

    Details:
        The algorithm used is Pollard's p - 1 method. We know that if p is an
        odd prime, then 2^(p - 1) = 1 (mod p), and 2^M = 1 (mod p) if (p - 1) |
        M. So if p is a prime factor of an integer n, then p divides gcd(2^M -
        1, n). The idea behind this algorithm is to choose M with many divisors
        of the form p - 1, and then search for many primes p as possible
        divisors of n at once. For more information, see section 5.4 of "Prime
        Numbers - A Computational Perspective" by Crandall and Pomerance.
    """
    c = randint(2, 20)
    p_list = rosemary.number_theory.sieves.primes(B)
    for p in p_list:
        a = integer_log(B, p)
        for _ in xrange(a):
            c = pow(c, p, n)
    g = gcd(c - 1, n)
    return g


def pollard_rho(n):
    """
    Attempts to find a nontrivial factor of n.

    Given a composite number n, this algorithm attempts to find a nontrivial
    factor of n. The algorithm used is Pollard's rho algorithm.

    Input:
        * n: int (n > 1)

    Output:
        * d: int

    Examples:
        >>> m = 1112470797641561909
        >>> pollard_rho(m)
        1052788969L
        >>> m = 2175282241519502424792841
        >>> pollard_rho(m)
        513741730823L

    Details:
        The algorithm used is Pollard's rho factorization method. This method
        should return a nontrivial factor of n in O(sqrt(p)) steps, where p is
        the least prime factor of n. Because of this dependence on the smallest
        prime dividing n and not n itself, this method is especially useful for
        large composites with small prime factors out of range of trial
        division. For more details, see section 5.2.1 of "Prime Numbers - A
        Computational Perspective" by Crandall and Pomerance and section 8.5 of
        "A Course in Computational Algebraic Number Theory" by Cohen.
    """
    # Instead of computing a gcd in each iteration, we accumulate the products
    # and take the gcd only when the number of terms is a multiple of 'step'
    step = 20
    prod = 1
    terms = 0

    # a is the constant term of the polynomial x^2 + a
    a = randint(1, n - 3)

    # u is the random seed
    u = randint(0, n - 1)
    v = u

    while True:
        u = (u*u + a) % n
        v = (v*v + a) % n
        v = (v*v + a) % n
        prod = prod*(u - v) % n
        terms += 1

        if terms == step:
            g = gcd(prod, n)
            if 1 < g < n:
                return g
            elif g == n:
                return n
            prod = 1
            terms = 0


def pollard_rho_brent(n):
    """
    Attempts to find a nontrivial factor of n.

    Given a composite number n, this algorithm attempts to find a nontrivial
    factor of n. The method used is Brent's improvement to the Pollard-Rho
    algorithm.

    Input:
        * n: int (n > 1)

    Output:
        * d: int

    Examples:
        >>> m = 1112470797641561909
        >>> pollard_rho_brent(m)
        1052788969L
        >>> m = 2175282241519502424792841
        >>> pollard_rho_brent(m)
        513741730823L

    Details:
        The algorithm used is Brent's improvement to Pollard's rho factorization
        method. This method should return a nontrivial factor of n in O(sqrt(p))
        steps, where p is the least prime factor of n. Because of this
        dependence on the smallest prime dividing n and not n itself, this
        method is especially useful for large composites with small prime
        factors out of range of trial division. For more details, see the paper
        "An Improved Monte Carlo Factorization Algorithm" by R.P. Brent and
        section 8.5 of "A Course in Computational Algebraic Number Theory" by
        Cohen.
    """
    # c is the constant term of the polynomial
    c = randint(1, n - 3)

    # y is the random seed
    y = randint(0, n - 1)

    # m is the number of terms to multiply together before taking a gcd
    m = 100
    prod = 1
    r = 1

    while True:
        x = y
        for _ in xrange(r):
            y = (y*y + c) % n

        k = 0
        while True:
            ys = y
            for _ in xrange(m):
                y = (y*y + c) % n
                prod = prod*(x - y) % n
            g = gcd(prod, n)
            k += m
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


def z2_gaussian_elimination(exponents):
    """

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

    # Starting in the rightmost column, find the first vector whose rightmost 1
    # is in column j.
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
        # If we haven't found a pivot yet, move to the next column left and look
        # again.
        if pivot_found:
            for m in xrange(i + 1, num_rows):
                if reduced[m] & bit:
                    reduced[m] ^= entry
                    history[m] ^= his

        bit <<= 1

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


def smooth_factor(n, factor_base, num_primes=None):
    """
    Returns the exponents in the prime factorization of n if n factors
    completely using the primes in 'factor_base', and returns False otherwise.

    Input:
        * n: int

        * factor_base: list
            A sorted list of primes, possibly including -1.

        * num_primes: int (default=None)
            Number of elements in 'factor_base'.

    Output cases:
        * exponents: list
            If n factors completely over the factor_base, then this is a list of
            the exponents appearing in the prime factorization; exponents[i]
            holds the power to which factor_base[i] divides n.

        * Returns False if n does not factor over factor_base

    Examples:
        >>> smooth_factor(100, [2, 5])
        [2, 2]
        >>> smooth_factor(100, [2, 3, 5, 7])
        [2, 0, 2, 0]
        >>> smooth_factor(100, [3, 5])
        False
    """
    if num_primes is None:
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

    return False


def dixon(n, num_trials=10000):
    B = int(exp((log(n)*log(log(n)))**(0.5))**(0.5)) + 1
    factor_base = rosemary.number_theory.sieves.primes(B)
    num_primes = len(factor_base)

    for p in factor_base:
        if n % p == 0:
            return p

    num_smooths_found = 0
    max_element = n
    power = 1
    while 2**power < max_element:
        power *= 2

    prod = 1
    for p in factor_base:
        prod *= p

    for trial in xrange(num_trials):
        exponent_matrix = []
        a_list = []
        while num_smooths_found <= num_primes:
            a = randint(1, n)
            s = a*a % n
            pp = pow(prod, power, s)
            g = gcd(s, pp)
            if g != s:
                continue

            exponent_vector = smooth_factor(s, factor_base, num_primes)
            exponent_matrix.append(exponent_vector)
            num_smooths_found += 1
            print num_smooths_found, num_primes
            a_list.append(a)

        kernel = z2_gaussian_elimination(exponent_matrix)

        for i in xrange(len(kernel)):
            y = 1
            x2_exponents = [0]*num_primes
            for j in xrange(len(kernel[i])):
                if kernel[i][j]:
                    y = (a_list[j]*y) % n
                    for f in xrange(num_primes):
                        x2_exponents[f] += exponent_matrix[j][f]

            x = 1
            for f in xrange(num_primes):
                x *= factor_base[f]**(x2_exponents[f]//2)

            for val in [x - y, x + y]:
                d = gcd(val, n)
                if 1 < d < n:
                    return d


def cfrac(n, k=1):
    """
    """
    def find_multiplier():
        """
        Computes a multiplier for input into the cfrac factorization algorithm.

        The method is based on Remark 5.3 from the paper "A Method of Factoring
        and the Factorization of F7" by Morrison and Brillhart.
        """
        choices = {}
        # Look for multpliers in the range [1, 1000).
        for k in xrange(1, 1000):
            if jacobi_symbol(k*n, 3) >= 0 and jacobi_symbol(k*n, 5) >= 0:
                # Find the multiplier k that allows the largest number of primes
                # <= 31 into the factor base.
                count = 0
                for p in prime_list:
                    # We've already looked at p = 3 and p = 5.
                    if p <= 5:
                        continue
                    if p > 31:
                        break
                    if jacobi_symbol(k*n, p) >= 0:
                        count += 1
                if count not in choices:
                    choices[count] = [k]
                else:
                    choices[count].append(k)
        # If several values of k allow this maximal number, we simply choose the
        # smallest of them.
        max_count = max(choices)
        return min(choices[max_count])

    def cfrac_aq_pairs(n):
        """
        Yields tripes (i, A_{i - 1}, Q_i) for i > 0, where
            A_{i - 1}^2 = (-1)^i Q_i (mod n)

        Input:
            * n: int

        Output:
            * X: generator

        Details
            This algorithm expands sqrt(n) into a simple continued fraction. The
            values (i, A_{i - 1}, Q_i) output by this algorithm correspond to
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

    # B is our smoothness bound.
    B = int(exp(0.5*sqrt(log(n)*log(log(n))))) + 1
    prime_list = rosemary.number_theory.sieves.primes(B)

    # Choose a multiplier if none is provided.
    if k == 0:
        k = find_multiplier()
    kn = k*n

    # Perform simple trial division by the primes we computed to find any small
    # prime divisors.
    for p in prime_list:
        if n % p == 0:
            return p

    # Our factor base needs to include -1 and 2, and the odd primes p for which
    # (kN|p) = 0 or 1.
    factor_base = [-1]
    for p in prime_list:
        if p == 2 or jacobi_symbol(kn, p) >= 0:
            factor_base.append(p)

    # We compute the product of the elements in our factor base for smoothness
    # checking computations later.
    prod = 1
    for p in factor_base:
        prod *= p

    num_primes = len(factor_base)

    # Instead of using trial division to check each value for smoothness
    # individually, we use a batch smoothness test, processing batches of size
    # 'batch_size' at once.

    exponent_matrix = []
    num_smooths_found = 0
    a_list = []

    max_element = n
    power = 1
    while 2**power < max_element:
        power *= 2

    aq_pairs = cfrac_aq_pairs(kn)

    while num_smooths_found <= num_primes:
        (i, a, q) = aq_pairs.next()

        pp = pow(prod, power, q)
        g = gcd(q, pp)
        if g != q:
            continue

        if i % 2 == 1:
            q *= -1

        exponent_vector = smooth_factor(q, factor_base, num_primes)
        exponent_matrix.append(exponent_vector)
        num_smooths_found += 1
        a_list.append(a)

    kernel = z2_gaussian_elimination(exponent_matrix)

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


################################################################################
# Methods related to factorization and divisors
################################################################################

def factor(n):
    """
    factorization(n):
    Returns the factorization of n. Currently, this routine calls a number of
    different algorithms.

    Examples:
    >>> factor(100)
    [(2, 2), (5, 2)]
    >>> factor(537869)
    [(37, 1), (14537, 1)]
    """
    D = collections.defaultdict(int)

    if n == 0:
        raise ValueError("Prime factorization of 0 not defined")

    # Take care of the sign
    if n < 0:
        D[-1] = 1
        n = -1 * n

    # First, strip off all small factors on n
    for p in PRIME_LIST:
        if n == 1:
            break
        elif p*p > n:
            D[n] += 1
            n = 1
            break

        while n % p == 0:
            D[p] += 1
            n = n // p

    # Next, use pollard rho
    while n > 1:
        if rosemary.number_theory.primality.is_probable_prime(n):
            D[n] += 1
            n = 1
        else:
            d = pollard_rho_brent(n)
            while not rosemary.number_theory.primality.is_probable_prime(d):
                d = pollard_rho_brent(d)

            while n % d == 0:
                D[d] += 1
                n = n // d

    p_list = sorted([(p, D[p]) for p in D])
    return p_list


def factor_back(F):
    """
    factor_back(F):
    Given a factorization F, this return the factored integer.

    Examples:
    >>> factor_back([(2, 2), (5, 2)])
    100
    """
    if not isinstance(F, list):
        raise ValueError("Not a factorization in factor_back")

    pp = 1
    for (p, e) in F:
        pp *= p**e

    return pp


def divisors(n):
    """
    divisors(n):
    Returns a sorted list of the positive integer divisors of n. The argument
    n can be an integer, or the factorization of an integer.

    Examples:
    >>> divisors(100)
    [1, 2, 4, 5, 10, 20, 25, 50, 100]
    >>> divisors([(2, 2), (5, 2)])
    [1, 2, 4, 5, 10, 20, 25, 50, 100]
    >>> divisors(426497)
    [1, 71, 6007, 426497]
    """
    if isinstance(n, (int, long)):
        n_fac = factor(abs(n))
    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    p_divs = [p for (p, e) in n_fac]
    div_list = []
    iter_list = (xrange(e + 1) for (p, e) in n_fac)

    for tup in itertools.product(*iter_list):
        pp = 1
        for i in xrange(len(tup)):
            pp *= p_divs[i]**tup[i]
        div_list += [pp]

    div_list.sort()
    return div_list


def prime_divisors(n):
    """
    prime_divisors(n):
    Returns a list of the primes dividing n

    Examples:
    >>> prime_divisors(120)
    [2, 3, 5]
    >>> prime_divisors(5272)
    [2, 659]
    """
    if isinstance(n, (int, long)):
        n_fac = factor(abs(n))
    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    p_divs = [p for (p, e) in n_fac]
    return p_divs


def xdivisors(n):
    """
    xdivisors(n):
    Returns an iterator over the positive integer divisors of n.
    The divisors are not yielded in increasing order.
    """
    if isinstance(n, (int, long)):
        n_fac = factor(n)
    elif isinstance(n, list):
        n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    p_divs = [p for (p, e) in n_fac]
    iter_list = (xrange(e + 1) for (p, e) in n_fac)
    for tup in itertools.product(*iter_list):
        pp = 1
        for i in xrange(len(tup)):
            pp *= p_divs[i]**tup[i]
        yield pp
    return


def is_squarefree(n):
    """
    Determines if n is squarefree.

    Given a nonnegative integer n, this return True iff n is not divisible by
    the square of an integer > 1.

    Input:
        * n: int (n >= 0)

    Output:
        * b: bool

    Details:
        If n is a nonnegative integer, this factors n and checks if n is
        divisible by the square of a prime.  If n is in factored form, this
        directly checks the prime factorization.

    Examples:
        >>> is_squarefree(35)
        True
        >>> is_squarefree(100)
        False
    """
    if isinstance(n, (int, long)):
        n_fac = factor(abs(n))

    elif isinstance(n, list):
        if n[0][0] == -1:
            n_fac = n[1:]
        else:
            n_fac = n
    else:
        raise ValueError("Input must be an integer or a factorization")

    for (p, e) in n_fac:
        if e > 1:
            return False

    return True
