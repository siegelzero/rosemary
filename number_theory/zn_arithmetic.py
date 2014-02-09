# Modular arithmetic

import itertools
import rosemary.number_theory.factorization

from random import randint
from rosemary.number_theory.core import (
    chinese_preconditioned,
    crt_preconditioning_data,
    ext_gcd,
    gcd,
    integer_sqrt,
    inverse_mod,
    jacobi_symbol,
    valuation,
)

########################################################################################################################
# primitive roots
########################################################################################################################

def is_primitive_root(a, p, prime_divisors=None):
    """
    Returns True if a is a primitive root mod p and False otherwise.

    Input:
        * a: int
        * p: int
        * prime_divisors: list (default=None)
            List of prime divisors of p - 1.

    Output:
        * b: bool
            b is True if a is a primitive root mod p and False otherwise.

    Examples:
        >>> is_primitive_root(2, 11)
        True
        >>> is_primitive_root(3, 11)
        False

    Details:
        The algorithm is based on the observation that a is a primitive root
        modulo p if and only if a^((p - 1)/q) != 1 (mod p) for all primes q
        dividing p - 1. The primality of p is not verified.
    """
    p_minus_one = p - 1
    if prime_divisors is None:
        prime_divisors = rosemary.number_theory.factorization.prime_divisors(p_minus_one)

    for pk in prime_divisors:
        if pow(a, p_minus_one//pk, p) == 1:
            return False
    return True


def primitive_root_mod_p(p, prime_divisors=None):
    """
    Returns a primitive root modulo p.

    Input:
        * p: int
            A prime number. The primality of p is not verified.

        * prime_divisors: list (default=None)
            List of prime divisors of p - 1.

    Output:
        * a: int
            This is a primitive root modulo p; i.e. a is a generator for the
            multiplicative group of nonzero residues modulo p.

    Examples:
        >>> primitive_root_mod_p(7)
        3
        >>> primitive_root_mod_p(11)
        2
    """
    p_minus_one = p - 1
    if prime_divisors is None:
        prime_divisors = rosemary.number_theory.factorization.prime_divisors(p_minus_one)

    for a in xrange(2, p):
        if is_primitive_root(a, p, prime_divisors):
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

    Examples:
        >>> fibonacci_primitive_roots(11)
        [8]
        >>> fibonacci_primitive_roots(19)
        [15]
        >>> fibonacci_primitive_roots(41)
        [7, 35]
    """
    if p % 2 == 0:
        return []
    if p == 5:
        return [3]
    if p < 5 or jacobi_symbol(5, p) != 1:
        return []

    sqrts = sqrts_mod_p(5, p)
    inverse = inverse_mod(2, p)
    roots = [(1 + root)*inverse % p for root in sqrts]
    prime_divisors = rosemary.number_theory.factorization.prime_divisors(p - 1)
    primitive_roots = []

    for r in roots:
        if is_primitive_root(r, p, prime_divisors):
            primitive_roots.append(r)
    return primitive_roots

########################################################################################################################
# root extraction
########################################################################################################################

def sqrts_mod_p(a, p):
    """
    Returns a solution x to x^2 = a (mod p).

    Given a prime p and an integer a, this function returns all solutions to the
    congruence x^2 = a (mod p).

    Input:
        * a: int
        * p: int

    Output:
        * solutions: list

    Examples:
        >>> sqrts_mod_p(10, 13)
        [6, 7]

    Details:
        For odd primes p, this function uses the algorithm due to Tonelli and
        Shanks. See Algorithm 2.3.8 in "Prime Numbers - A Computational
        Perspective" by Crandall and Pomerance for details of the algorithm. For
        p = 2, simple checking is done to find solutions.
    """
    if p == 2:
        return [a % p]

    a = a % p
    if a == 0:
        return [0]

    if jacobi_symbol(a, p) != 1:
        return []

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
        s = valuation(2, p - 1)
        t = (p - 1) >> s

        A = pow(a, t, p)
        D = pow(d, t, p)
        m = 0
        for i in xrange(s):
            if pow(A*D**m, 2**(s - 1 - i), p) == p - 1:
                m += 2**i
        x = pow(a, (t + 1)//2, p)*pow(D, m//2, p) % p

    solutions = [x, p - x]
    solutions.sort()
    return solutions


def sqrts_mod_pk(a, p, k):
    def hensel(a, r, p, k):
        """
        Lifts solution r^2 = a (mod p^k) to solutions modulo p^{k + 1}.
        """
        solutions = []
        if 2*r % p != 0:
            # In this case, r lifts to a unique solution.
            inverse = inverse_mod(2*r, p)
            # Note that r^2 = a (mod p^k) so the division below is exact.
            t = (-inverse*((r*r - a)//p**k)) % p
            root = r + t*p**k
            solutions.append(root)
        else:
            # Here, r lifts to p incongruent solutions modulo p^k.
            if (r*r - a) % (p**(k + 1)) == 0:
                for t in xrange(p):
                    root = r + t*p**k
                    solutions.append(root)
        solutions.sort()
        return solutions

    if k == 1:
        return sqrts_mod_p(a, p)

    if p == 2:
        if k == 2:
            solutions = [e for e in xrange(4) if e*e % 4 == a % 4]
        else:
            solutions = [e for e in xrange(8) if e*e % 8 == a % 8]
            if k > 3:
                for e in xrange(3, k):
                    lifted_solutions = []
                    for root in solutions:
                        lifts = hensel(a, root, p, e)
                        lifted_solutions.extend(lifts)
                    solutions = list(lifted_solutions)
    else:
        solutions = sqrts_mod_p(a, p)
        for e in xrange(1, k):
            lifted_solutions = []
            for root in solutions:
                lifts = hensel(a, root, p, e)
                lifted_solutions.extend(lifts)
            solutions = list(lifted_solutions)

    solutions.sort()
    return solutions


def sqrts_mod_n(a, n, n_factorization=None):
    """
    Returns all solutions x, 1 <= x <= n, to the congruence x^2 = a (mod n).

    Input:
        * a: int
            A square modulo n.

        * n: int
            The modulus.

        * n_factorization: list (default=None)
            The factorization of n.

    Output:
        * roots: list
            A list of all square roots of a modulo n.
    """
    if n_factorization is None:
        n_factorization = rosemary.number_theory.factorization.factor(n)

    congruences = []
    moduli = []
    for (p, k) in n_factorization:
        roots = sqrts_mod_pk(a, p, k)
        pk = p**k
        pairs = [(r, pk) for r in roots]
        congruences.append(pairs)
        moduli.append(pk)

    if len(congruences) == 1:
        values = [r for (r, pk) in congruences[0]]
    else:
        preconditioning_data = crt_preconditioning_data(moduli)
        values = []

        for system in itertools.product(*congruences):
            values.append(chinese_preconditioned(system, preconditioning_data))

    values.sort()
    return values


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
    roots_of_unity = nth_roots_of_unity_mod_p(n, p, g)
    all_roots = [root*h % p for h in roots_of_unity]
    all_roots.sort()
    return all_roots


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

################################################################################
# solutions to congruences
################################################################################

def linear_congruence(a, b, n):
    """
    Returns a list of the solutions x to the congruence a*x = b (mod n).

    Input:
        * a: int
        * b: int
        * n: int (n > 1)

    Output:
        * solutions: list

    Examples:
        >>> linear_congruence(10, 6, 12)
        [3, 9]
        >>> linear_congruence(12, 9, 15)
        [2, 7, 12]
        >>> linear_congruence(10, 3, 12)
        []

    Details:
        The linear congruence a*x = b (mod n) has a solution if and only if d =
        gcd(a, n) divides b. This function uses a straightforward application of
        this theorem. See Theorem 3.7 from "Elementary Number Theory" By Jones
        and Jones for details.
    """
    if n < 2:
        raise ValueError("linear_congruence: Must have n >= 2.")

    (u, v, d) = ext_gcd(a, n)
    # The congruence has a solution if and only if gcd(a, n) | b.
    if b % d != 0:
        return []

    # x0 is our particular solution.
    # There will be exactly d incongruent solutions modulo n.
    x0 = b*u//d
    solutions = [(x0 + k*n//d) % n for k in xrange(d)]

    solutions.sort()
    return solutions


def quadratic_congruence(coeff_list, n, n_factorization=None):
    """
    Returns the roots of the quadratic equation modulo n.

    Input:
        * coeff_list: list
            A list of the coefficients of the quadratic.

        * n: int
            The modulus.

        * n_factorization: list (default=None)
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
        >>> quadratic_congruence([1, 3, -18], 1000)
        [3, 378, 619, 994]
    """

    if n_factorization is None:
        n_factorization = rosemary.number_theory.factorization.factor(n)

    (a, b, c) = coeff_list
    discriminant = b*b - 4*a*c
    all_roots = []

    if gcd(2*a, n) == 1:
        # This is the easy case where we can complete the square.
        # First, we solve the congruence y^2 = b^2 - 4*a*c (mod n) for y.
        discriminant_roots = sqrts_mod_n(discriminant, n, n_factorization)

        # Next, we solve the congruence 2*a*x = y - b (mod n) for x to obtain
        # the solutions to the quadratic. Since gcd(2*a, n) == 1, this is
        # simple, and each each value of y leads to one value of x.
        inverse = inverse_mod(2*a, n)
        for y in discriminant_roots:
            x = (y - b)*inverse % n
            all_roots.append(x)
    else:
        # Here, gcd(4*a, n) != 1, so we can't complete the square as usual.
        # Write 4*a = a1*a2, with a2 coprime to n.
        a1 = 1
        a2 = 4*a
        d = gcd(n, a2)
        while d > 1:
            a1 *= d
            a2 /= d
            d = gcd(n, a2)

        # We solve the congruence y^2 = b^2 - 4*a*c (mod a1*n) for y.
        discriminant_roots = sqrts_mod_n(discriminant, a1*n)

        # For each solution y, we solve 2*a*x = y - b (mod a1*n) for x. Since
        # gcd(2*a, n) > 1, each solution y leads to multiple values of x.
        for y in discriminant_roots:
            roots = linear_congruence(2*a, y - b, a1*n)
            all_roots.extend(roots)

        # Eliminate repeated solutions, and reduce modulo the original modulus.
        distinct_roots = {x % n for x in all_roots}
        all_roots = list(distinct_roots)

    all_roots.sort()
    return all_roots
