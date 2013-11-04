# Fp arithmetic

import itertools
from random import randint
from rosemary.number_theory.core import (
    chinese,
    chinese_preconditioned,
    crt_preconditioning_data,
    gcd,
    integer_sqrt,
    inverse_mod,
    jacobi_symbol,
    valuation,
)
from rosemary.number_theory.factorization import factor, prime_divisors

########################################################################################################################
# primitive roots
########################################################################################################################

def is_primitive_root(a, p, primeDivisors=None):
    """
    Returns True if a is a primitive root mod p and False otherwise.

    Input:
        * a: int
            Potential primitive root.

        * p: int
            A prime.

        * primeDivisors: list (default=None)
            List of prime divisors of p - 1.

    Output:
        * b: bool
            b is True if a is a primitive root mod p and False otherwise.

    Details:
        The algorithm is based on the observation that a is a primitive root
        modulo p if and only if a^((p - 1) / q) != 1 (mod p) for all primes q
        dividing p - 1. The primality of p is not verified.
    """
    pMinusOne = p - 1
    if primeDivisors is None:
        primeDivisors = prime_divisors(pMinusOne)

    for pk in primeDivisors:
        if pow(a, pMinusOne//pk, p) == 1:
            return False
    return True

def primitive_root_mod_p(p, primeDivisors=None):
    """
    Returns a primitive root modulo p.

    Input:
        * p: int
            A prime number. The primality of p is not verified.

        * primeDivisors: list (default=None)
            List of prime divisors of p - 1.

    Output:
        * a: int
            This is a primitive root modulo p; i.e. a is a generator for the
            multiplicative group of nonzero residues modulo p.
    """
    pMinusOne = p - 1
    if primeDivisors is None:
        primeDivisors = prime_divisors(pMinusOne)

    for a in xrange(2, p):
        if is_primitive_root(a, p, primeDivisors):
            return a

def fibonacci_primitive_roots(p):
    """
    Returns a list of the Fibonacci primitive roots mod p.

    A Fibonacci primitive root mod p is a primitive root satisfying g^2 = g + 1.
    These can only exist for primes p = 1, 9 (mod 10).

    Input:
        * p: int
            A prime number. The primality of p is not verified.
    
    Output:
        * roots: list
            This is a list of the fibonacci primitive roots mod p. This list
            will contain no more than 2 elements.
    """
    if p == 5:
        return [3]
    if p < 5 or jacobi_symbol(5, p) != 1:
        return []

    sqrt5 = sqrt_mod_p(5, p)
    inverse = inverse_mod(2, p)
    r1 = (1 + sqrt5) * inverse % p
    r2 = (1 - sqrt5) * inverse % p
    primeDivisors = prime_divisors(p - 1)
    roots = []

    for r in (r1, r2):
        if is_primitive_root(r, p, primeDivisors):
            roots.append(r)
    return roots

########################################################################################################################
# modular root extraction
########################################################################################################################

def sqrt_mod_p(a, p):
    """
    Returns a solution x to x^2 = a (mod p).

    Given an odd prime p and an integer a with (a|p) = 1, this algorithm returns
    a solution x to x^2 = a (mod p). The algorithm used is due to Tonelli and
    Shanks.

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
        x = pow(a, (t + 1)//2, p)*pow(D, m//2, p) % p
    return x

def sqrt_mod_pk(a, p, k):
    """
    Return all integers b, 1 <= b <= p^k, such that b^2 = a (mod p^k).

    Input:
        * a: int
            A square modulo p.

        * p: int
            The prime modulus.

        * k: int
            Exponent on the prime modulus.

    Output:
        * solutions: list
            A list of all solutions to b^2 = a (mod p^k).

    Details:
        For odd primes p, the procedure finds a square root, and then uses
        Hensel lifting to find a root modulo p^k. For p = 2, the method proceeds
        by cases for lifts past modulus 16.

    """
    def hensel_lift(a, b, p, e):
        """
        Lifts a solutions b^2 = a (mod p) to a solution modulo p^e.

        Input:
            * a: int
                A square modulo p.

            * b: int
                An integer satisfying b^2 = a (mod p^e).
            
            * p: int
                The odd prime modulus.

            * e: int
                The exponent of the prime power modulus.

        Output:
            * c: int
                An integer satifsying c^2 = a (mod p^e).
        """
        pk = p
        for k in xrange(1, e):
            tt = (a - b*b)//pk
            h = tt*inverse_mod(2*b, p)
            b = (b + h*pk) % (pk*p)
            pk *= p
        return b

    if p % 2 == 1:
        b = sqrt_mod_p(a, p)
        c = hensel_lift(a, b, p, k)
        solutions = [c]
    else:
        if a % 2 == 0:
            raise ValueError("a cannot be even")

        if k == 1:
            solutions = [1]
            c = 1
        elif k == 2:
            if a % 4 != 1:
                raise ValueError("Must have a = 1 (mod 4)")
            solutions = [1]
        else:
            if a % 8 != 1:
                raise ValueError("Must have a = 1 (mod 8)")
            solutions = [1, 3]
            if k > 3:
                for e in xrange(3, k + 1):
                    liftedSolutions = solutions[:]
                    for (i, xi) in enumerate(solutions):
                        h = (xi**2 - a) % 2**k
                        diff = h//(2**e)
                        if diff % 2 == 1:
                            liftedSolutions[i] += 2**(e - 1)
                    solutions = liftedSolutions[:]

    allSolutions = solutions + [p**k - e for e in solutions]
    allSolutions.sort()
    return allSolutions

def sqrt_mod_n(a, n, nFactorization=None):
    """
    Returns all solutions x, 1 <= x <= n, to the congruence x^2 = a (mod n).

    Input:
        * a: int
            A square modulo n.

        * n: int
            The modulus.

        * nFactorization: list (default=None)
            The factorization of n.

    Output:
        * roots: list
            A list of all square roots of a modulo n.
    """
    if nFactorization is None:
        nFactorization = factor(n)

    congruences = []
    for (p, k) in nFactorization:
        roots = sqrt_mod_pk(a, p, k)
        pairs = [(r, p**k) for r in roots]
        congruences.append(pairs)

    numCongruences = len(congruences)

    if numCongruences == 1:
        values = [r for (r, pk) in congruences[0]]
    else:
        values = []
        for system in itertools.product(*congruences):
            values.append(chinese(system))

    values.sort()
    return values

def discrete_log(a, b, p):
    """
    Finds a positiveinteger k such that a^k = b (mod p).

    Input:
        * a: int
            Base of the exponential. Typically, this will be a primitive root
            modulo p.

        * b: int
            The target. Must be
        * p: int

    Output:
        * k: int
            We have 0 <= k < p is an integer satisfying a^k = b (mod p).

    Details:
        The algorithm is based on the Baby-Step Giant-Step method due to Shanks.
        Using a dictionary data type, this algoritm runs in O(sqrt(n)) time and
        space. We are following the method as described in Algorithm 2.4.1 in
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

    Given a positive integer n and a prime p, this returns a list of the nth
    roots of unity modulo p. The primality of p is assumed and not verified.

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

    Given a positive integer n and a prime p, this returns a list of all nth
    roots of -1 modulo p. The primality of p is not verified.

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
    rootsOfUnity = nth_roots_of_unity_mod_p(n, p, g)
    allRoots = [root*h % p for h in rootsOfUnity]
    allRoots.sort()
    return allRoots

def quadratic_roots_mod_n(coeffList, n, nFactorization=None):
    """
    Returns the roots of the quadratic equation modulo n.

    Input:
        * coeffList: list
            A list of the coefficients of the quadratic.

        * n: int
            The modulus.

        * nFactorization: list (default=None)
            The factorization of the modulus n.

    Output:
        * roots: list
            A list of the roots of the quadratic modulo n.

    Details:
        The algorithm proceeds by finding roots modulo p^k for each prime power
        p^k in the factorization of the modulus n. These roots are then combined
        using the Chinese Remainder Theorem.

        The quadratic formula is used to find the roots, with the square roots
        computed modulo p^k, for each p^k dividing n. For odd primes p, the
        Tonelli-Shanks algorithm is used. For p = 2, a simple search is
        performed.

    Examples:
        >>> quadratic_roots_mod_n([1, 3, -18], 1000)
        [3, 378, 619, 994]
    """
    if nFactorization is None:
        nFactorization = factor(n)

    (a, b, c) = coeffList
    discriminant = b*b - 4*a*c
    allRoots = []

    for (p, k) in nFactorization:
        pk = p**k
        # Use dumb search for powers of 2
        modpRoots = []
        if p == 2:
            for x in xrange(pk):
                if (a*x*x + b*x + c) % pk == 0:
                    modpRoots.append((x, pk))
        else:
            squareRoots = sqrt_mod_pk(discriminant, p, k)
            inverse = inverse_mod(2*a, pk)

            for root in squareRoots:
                r1 = (-b + root)*inverse % n
                modpRoots.append((r1, pk))
        allRoots.append(modpRoots)

    combinedRoots = []
    for system in itertools.product(*allRoots):
        combinedRoots.append(chinese(system))

    combinedRoots.sort()
    return combinedRoots

def idempotents_mod_n(n, nFactorization=None):
    """
    Returns a list of idempotents modulo n; i.e. elements such that a^2 = a.
    """
    if n == 1:
        return [0]

    if nFactorization is None:
        nFactorization = factor(n)

    allRoots = []
    moduli = []
    for (p, k) in nFactorization:
        pk = p**k
        moduli.append(pk)
        modpRoots = [(0, pk), (1, pk)]
        allRoots.append(modpRoots)

    preconditioningData = crt_preconditioning_data(moduli)
    combinedRoots = []
    for combination in itertools.product(*allRoots):
        combinedRoots.append(chinese_preconditioned(combination, preconditioningData))

    combinedRoots.sort()
    return combinedRoots

