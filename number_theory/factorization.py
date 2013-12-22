import collections
import itertools
from random import randint

import rosemary.number_theory.sieves
import rosemary.number_theory.primality

from rosemary.number_theory.prime_list import PRIME_LIST
from rosemary.number_theory.core import (
    gcd,
    integer_log,
    integer_nth_root,
    integer_sqrt,
)

################################################################################
# Factorization algorithms
################################################################################

def fermat(n):
    """
    Returns a nontrivial divisor of n, or proves primality.

    Given an odd integer n > 1, this algorithm attempts to find a factor of n
    using Fermat's method.

    Input:
        * n: int (n > 1)

    Output:
        * d: int
            If n == d, then n is proven prime. Otherwise, d is a nontrivial
            divisor of n.

    Details:
        The algorithm used here is Fermat's method. This algorithm runs in time
        O((n + 1) / 2 - sqrt(n)), where n is the number to be factored.  In the
        worst case (n = 3*p), this algorithm is much worse than trial division.
        See section 5.1.1 of "Prime Numbers - A Computation Perspective" by
        Crandall and Pomerance for more details.

    Examples:
        >>> m = 1112470797641561909
        >>> fermat(m)
        1052788969
    """
    if n % 2 == 0:
        return 2
    a = integer_sqrt(n) + 1
    while a <= (n + 9) // 6:
        t = a**2 - n
        b = integer_sqrt(t)
        if b*b == t:
            return a - b
        a += 1
    return n


def lehman(n):
    """
    Returns a nontrivial divisor of n, or proves primality.

    Given an integer n >= 3, this algorithm finds a nontrivial factor of n if n
    is not prime, or returns n if n is prime.

    Input:
        * n: int (n >= 3)

    Output:
        * d: int
            If d == n, then n is proven prime. Otherwise, d is a nontrivial
            divisor of n.

    Details:
        The algorithm used is Lehman's Method. This algorithm runs in time
        O(n^(1/3)), where n is the number to be factored / proven prime. This is
        substantially better than O(n^(1/2)) trial division for reasonably small
        value of n, but this algorithm is not suited for large values of n. See
        section 8.4 of "A Course in Computational Algebraic Number Theory" by
        Cohen or section 5.1.2 of "Prime Numbers - A Computational Perspective"
        by Crandall and Pomerance for more details.

    Examples:
        >>> l = 1112470797641561909
        >>> lehman(l)
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


def pollard_p_1(n, B=20000):
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

    Details:
        The algorithm used is Pollard's p - 1 method. We know that if p is an
        odd prime, then 2^(p - 1) = 1 (mod p), and 2^M = 1 (mod p) if (p - 1) |
        M. So if p is a prime factor of an integer n, then p divides gcd(2^M -
        1, n). The idea behind this algorithm is to choose M with many divisors
        of the form p - 1, and then search for many primes p as possible
        divisors of n at once. For more information, see section 5.4 of "Prime
        Numbers - A Computational Perspective" by Crandall and Pomerance.

    Examples:
        >>> m = 1112470797641561909
        >>> pollard_p_1(m)
        1056689261L
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
    factor of n.

    Input:
        * n: int (n > 1)

    Outout:
        * d: int
            Nontrivial divisor of n.

    Details:
        The algorithm used is Pollard's rho factorization method. This method
        should return a nontrivial factor of n in O(sqrt(p)) steps, where p is
        the least prime factor of n. For more details, see section 5.2.1 of
        "Prime Numbers - A Computational Perspective" by Crandall and Pomerance
        and section 8.5 of "A Course in Computational Algebraic Number Theory"
        by Cohen.

    Examples:
        >>> m = 1112470797641561909
        >>> pollard_rho(m)
        1052788969L
        >>> m = 2175282241519502424792841
        >>> pollard_rho(m)
        513741730823L
    """
    a = randint(1, n - 3)
    s = randint(0, n - 1)
    u = s
    v = s
    while True:
        u = (u*u + a) % n
        v = (v*v + a) % n
        v = (v*v + a) % n
        g = gcd(u - v, n)
        if 1 < g < n:
            return g
        elif g == n:
            return n


def one_line_factor(k, M=10000000):
    #bound = integer_nth_root(3, k)
    #d = trial_division(k, bound)
    #if d < k:
    #    return d

    n = 480*k
    for i in xrange(1, M + 1):
        s = int((n*i)**(0.5)) + 1
        m = s*s % n
        t = int(m**(0.5))
        if t*t == m:
            g = gcd(k, s - t)
            return g
    return k


def trial_division(n, b=None):
    """
    trial_division(n, b):
    This algoritm performs trial division on n with divisors <= b. The smallest
    prime dividing n is returned if found, otherwise n is returned. We have a
    list of primes to check through first.
    """
    if b is None:
        b = integer_sqrt(n) + 1

    for d in PRIME_LIST:
        if d > b:
            # no prime divisors found <= b, so return n
            return n
        if n % d == 0:
            # return any prime divisors found
            return d

    # Next we use a segmented sieve for the rest
    for p in rosemary.number_theory.sieves.sieve_interval(PRIME_LIST[-1], b):
        if n % p == 0:
            return p
    return n


def factor_pollard_rho_brent(n):
    """
    Attempts to find a nontrivial factor of n.

    Given a composite number n, this algorithm attempts to find a nontrivial
    factor of n. The method used is Brent's improvement to the Pollard-Rho
    algorithm.
    """
    x0 = randint(0, n - 1)
    (r, q, y, g, m, c) = (1, 1, x0, 1, 3, 1)

    while True:
        x = y
        for _ in xrange(r):
            y = (y*y + c) % n
        k = 0
        while True:
            ys = y
            for _ in xrange(min(m, r - k)):
                y = (y*y + c) % n
                q = q * abs(x - y) % n
            g = gcd(q, n)
            k += m
            if k >= r or g > 1:
                break
        
        r *= 2
        if g > 1:
            break

    if g == n:
        while True:
            ys = (ys*ys + c) % n
            g = gcd(abs(x - ys), n)
            if g > 1:
                break
    return g

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
            d = factor_pollard_rho_brent(n)
            while not rosemary.number_theory.primality.is_probable_prime(d):
                d = factor_pollard_rho_brent(d)

            while n % d == 0:
                D[d] += 1
                n = n // d

    p_list = sorted([ (p, D[p]) for p in D ])
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

    p_divs = [ p for (p, e) in n_fac ]
    div_list = []
    iter_list = ( xrange(e + 1) for (p, e) in n_fac )

    for tup in itertools.product(*iter_list):
        pp = 1
        for i in xrange(len(tup)):
            pp *= p_divs[i]**tup[i]
        div_list += [ pp ]

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

    p_divs = [ p for (p, e) in n_fac ]
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

