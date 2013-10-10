# Fp arithmetic

from rosemary.number_theory.core import (
    gcd,
    integer_sqrt,
    inverse_mod,
    jacobi_symbol,
    valuation,
)

from rosemary.number_theory.factorization import prime_divisors
from random import randint

########################################################################################################################
# primitive roots
########################################################################################################################

def is_primitive_root(a, p, p_divs=None):
    """
    Returns True if a is a primitive root mod p and False otherwise.

    Input:
        * a: int
            Potential primitive root.

        * p: int
            A prime.

        * p_divs: list (default=None)
            List of prime divisors of p - 1.

    Output:
        * b: bool
            b is True if a is a primitive root mod p and False otherwise.

    Details:
        The algorithm is based on the observation that a is a primitive root modulo p if and only if a^((p - 1) / q) !=
        1 (mod p) for all primes q dividing p - 1. The primality of p is not verified.
    """
    p_minus_one = p - 1

    if p_divs is None:
        p_divs = prime_divisors(p_minus_one)

    for pk in p_divs:
        if pow(a, p_minus_one // pk, p) == 1:
            return False
    return True


def primitive_root_mod_p(p):
    """
    Returns a primitive root modulo p.

    Input:
        * p: int
            A prime number.

    Output:
        * a: int
            This is a primitive root modulo p; i.e. a is a generator for the multiplicative group of nonzero residues
            modulo p. The primality of p is not verified.
    """
    p_minus_one = p - 1
    p_divs = prime_divisors(p_minus_one)

    for a in xrange(2, p):
        if is_primitive_root(a, p, p_divs):
            return a


def fibonacci_primitive_roots(p):
    """
    Returns a list of the Fibonacci primitive roots mod p.

    A Fibonacci primitive root mod p is a primitive root satisfying g^2 = g + 1. These can only exist for primes p = 1,
    9 (mod 10).

    Input:
        * p: int
            A prime number.
    
    Output:
        * roots: list
            This is a list of the fibonacci primitive roots mod p. This list will contain no more than 2 elements.
    """
    if p == 5:
        return [3]
    if p < 5 or jacobi_symbol(5, p) != 1:
        return []

    sqrt5 = sqrt_mod_p(5, p)
    inv = inverse_mod(2, p)
    r1 = (1 + sqrt5) * inv % p
    r2 = (1 - sqrt5) * inv % p
    roots = []
    
    p_divs = prime_divisors(p - 1)

    for r in (r1, r2):
        if is_primitive_root(r, p, p_divs):
            roots.append(r)

    return roots

########################################################################################################################
# mod p root extraction
########################################################################################################################

def sqrt_mod_p(a, p):
    """
    Returns a solution x to x^2 = a (mod p).

    Given an odd prime p and an integer a with (a|p) = 1, this algorithm returns a solution x to x^2 = a (mod p). The
    algorithm used is due to Tonelli and Shanks.

    Input:
        * a: int
            A square modulo p.

        * p: int
            An odd prime.

    Output:
        * x: int
            This is an integer in [0, p) such that x^2 = a (mod p).

    Examples:
        >>> sqrt_mod_p(10, 13)
        7
        >>> pow(7, 2, 13)
        10
    """
    assert jacobi_symbol(a, p) == 1, "a must be a square modulo p"

    a = a % p
    if p % 8 in (3, 7):
        x = pow(a, (p + 1)//4, p)
    elif p % 8 == 5:
        x = pow(a, (p + 3)// 8, p)
        if x*x % p != a:
            x = x*pow(2, (p - 1)//4, p) % p
    else:
        # Find a quadratic nonresidue
        d = randint(2, p - 1)
        while jacobi_symbol(d, p) != -1:
            d = randint(2, p - 1)

        # Write p - 1 = 2^s * t with t odd
        s = valuation(p - 1, 2)
        t = (p - 1) >> s

        A = pow(a, t, p)
        D = pow(d, t, p)
        m = 0
        for i in xrange(s):
            if pow(A*D**m, 2**(s - 1 - i), p) == p - 1:
                m += 2**i
        x = pow(a, (t + 1)//2, p) * pow(D, m//2, p) % p
    return x


def discrete_log(a, b, p):
    """
    Finds a positiveinteger k such that a^k = b (mod p).

    Input:
        * a: int
            Base of the exponential. Typically, this will be a primitive root modulo p.

        * b: int
            The target. Must be
        * p: int

    Output:
        * k: int
            We have 0 <= k < p is an integer satisfying a^k = b (mod p).

    Details:
        The algorithm is based on the Baby-Step Giant-Step method due to Shanks. Using dictionary data type, this
        algoritm runs in O(sqrt(n)) time and space. We are following the method as described in Algorithm 2.4.1 in
        Number Theory for Computing by Yan.

    Examples:
        >>> discrete_log(2, 6, 19)
        14
        >>> pow(2, 14, 19)
        6
        >>> discrete_log(59, 67, 113)
        11
        >>> pow(59, 11, 113)
        67
    """
    m = integer_sqrt(p) + 1
    am = pow(a, m, p)

    cache = {}
    val = am
    for k in xrange(1, m + 1):
        cache[val] = k
        val = (val*am) % p

    val = b
    for k in xrange(m + 1):
        if val in cache:
            return cache[val]*m - k
        val = (val*a) % p


def nth_roots_of_unity_mod_p(n, p, g=None):
    """
    Returns the nth roots of unity modulo p.

    Given a positive integer n and a prime p, this returns a list of the nth roots of unity modulo p. The primality of p
    is assumed and not verified.

    Input:
        * n: int
            The index of the root.

        * p: int
            Prime modulus.

        * g: int (default=None)
            Primitive root modulo p.

    Output:
        * roots: list
            List of all nth roots of unity.
    """
    if p == 2:
        return [1]
    if g is None:
        g = primitive_root_mod_p(p)

    d = gcd(n, p - 1)
    # This is one solution. We find all others by looking at powers of this one.
    h = pow(g, (p - 1)//d, p)

    roots = []
    val = 1
    for _ in xrange(d):
        roots.append(val)
        val = (val*h) % p

    roots.sort()
    return roots


def nth_roots_of_minus1_mod_p(n, p, g=None):
    """
    Returns the nth roots of -1 modulo p.

    Given a positive integer n and a prime p, this returns a list of all nth roots of -1 modulo p. The primality of p is
    not verified.

    Input:
        * n: int
            The index of the root.

        * p: int
            Prime modulus.

        * g: int (default=None)
            Primitive root modulo p.

    Output:
        * roots: list
            List of all nth roots of -1 modulo p.
    """
    if p == 2:
        return [1]

    d = gcd(n, p - 1)

    if (p - 1)//d % 2 != 0:
        return []

    if g is None:
        g = primitive_root_mod_p(p)

    # This is one solution to x^n = -1 (mod p).
    root = pow(g, (p - 1)//(2*d), p)
    roots_of_unity = nth_roots_of_unity_mod_p(n, p, g)

    return sorted([root*h % p for h in roots_of_unity])

