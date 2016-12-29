# Modular arithmetic

import itertools
import rosemary.number_theory.factorization.factorization as factor
import rosemary.number_theory.arithmetic_functions.functions as functions

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


################################################################################
# primitive roots
################################################################################


def is_primitive_root(a, n, phi=None, phi_divisors=None):
    """Returns True if a is a primitive root mod n and False otherwise.

    Given an integer a and modulus n, this method determines whether or not a is
    a primitive root modulo n.

    Input:
        * a: int

        * n: int (n > 0)
            The modulus.

        * phi: int (default=None)
            Value of euler_phi(n).

        * phi_divisors: list (default=None)
            List of prime divisors of euler_phi(n).

    Returns:
        * b: bool
            b is True if a is a primitive root mod n and False otherwise. Note
            that b will be False when no primitive roots exist modulo n.

    Raises:
        * ValueError: if n <= 1.

    Examples:
        >>> is_primitive_root(2, 11)
        True
        >>> is_primitive_root(3, 11)
        False
        >>> is_primitive_root(3, -2)
        Traceback (most recent call last):
        ...
        ValueError: is_primitive_root: n must be >= 2.

    Details:
        An integer a is called a primitive root mod n if a generates the
        multiplicative group of units modulo n. Equivalently, a is a primitive
        root modulo n if the order of a modulo m is euler_phi(n).

        The algorithm is based on the observation that an integer a coprime to n
        is a primitive root modulo n if and only if a**(euler_phi(n)/d) != 1
        (mod n) for each divisor d of euler_phi(n). See Lemma 6.4 of "Elementary
        Number Theory" by Jones and Jones for a proof.

        Primitive roots exist only for the moduli n = 1, 2, 4, p**a, and 2*p**a,
        where p is an odd prime and a >= 1. These cases are captured in the
        above Lemma, however some simple divisibility tests allow us to exit
        early for certain moduli. See Chapter 6 of "Elementary Number Theory" by
        Jones and Jones for more details. See also Chapter 10 of "Introduction
        to Analytic Number Theory" by Apostol.
    """
    if n <= 1:
        raise ValueError("is_primitive_root: n must be >= 2.")

    if n > 4 and n % 4 == 0:
        return False

    if gcd(a, n) != 1:
        return False

    if phi is None:
        phi = functions.euler_phi(n)

    if phi_divisors is None:
        phi_divisors = factor.prime_divisors(phi)

    for d in phi_divisors:
        if pow(a, phi//d, n) == 1:
            return False

    return True


def primitive_root(n, n_factorization=None, phi=None, phi_divisors=None):
    """Returns a primitive root modulo n.

    Given a positive integer n of the form n = 2, 4, p**a, or 2*p**a for p an
    odd prime and a >= 1, this function returns the least primitive root of n.

    Input:
        * n: int (n > 1)

        * n_factorization: list (default=None)
            The prime factorization of n.

        * phi: int (default=None)
            The value of euler_phi(n).

        * phi_divisors: list (default=None)
            A list of the prime divisors of euler_phi(n).

    Returns:
        * g: int
            The least primitive root of n.

    Raises:
        * ValueError: if primitive roots don't exist modulo n.

    Examples:
        >>> primitive_root(7)
        3
        >>> primitive_root(11)
        2
        >>> primitive_root(14)
        3
        >>> primitive_root(12)
        Traceback (most recent call last):
        ...
        ValueError: primitive_root: No primitive root for n.

    Details:
        An integer a is called a primitive root mod n if a generates the
        multiplicative group of units modulo n. Equivalently, a is a primitive
        root modulo n if the order of a modulo m is euler_phi(n).

        As noted above, primitive roots only exist for moduli of the form n = 2,
        4, p**a, 2*p**a, where p > 2 is prime, and a >= 1. See Theorem 6.11 of
        "Elementary Number Theory" by Jones and Jones for a proof of this fact.
        See also Chapter 10 of "Introduction to Analytic Number Theory" by
        Apostol.

        This method uses the fact that an integer a coprime to n is a primitive
        root if and only if a**(euler_phi(n)/p) != 1 (mod n) for each prime p
        dividing euler_phi(n). See Lemma 6.4 in Jones and Jones for a proof of
        this.

        Note also that there is another way to construct primitive roots for
        composite moduli. To find a primitive root modulu p**a for p an odd
        prime and a >= 2, first compute g, a primitive root modulo p using the
        above criteria. Next, compute g1 = g**(p - 1) (mod p**2). If g1 != 1,
        then g is a primitive root modulo p**a for every a >= 1. Otherwise, g +
        p is. Finally, note that when p is an odd prime, if g is a primitive
        root modulo p**a, then either g or g + p**a (whichever is odd) is a
        primitive root modulo 2*p**a. See Lemma 1.4.5 of "A Course in
        Computational Algebraic Number Theory" by Cohen for more details.
    """
    if n_factorization is None:
        n_factorization = factor.factor(n)

    if n % 4 == 0 and n != 4:
        raise ValueError("primitive_root: No primitive root for n.")

    if len(n_factorization) > 2:
        raise ValueError("primitive_root: No primitive root for n.")

    if n % 2 == 1 and len(n_factorization) > 1:
        raise ValueError("primitive_root: No primitive root for n.")

    if phi is None:
        phi = functions.euler_phi(n_factorization)

    if phi_divisors is None:
        phi_divisors = factor.prime_divisors(phi)

    for g in xrange(1, n):
        if is_primitive_root(g, n, phi=phi, phi_divisors=phi_divisors):
            return g


def fibonacci_primitive_roots(p):
    """Returns a sorted list of the Fibonacci primitive roots mod p.

    Input:
        * p: int
            A prime number. The primality of p is not verified.

    Returns:
        * roots: list
            This is a list of the Fibonacci primitive roots mod p. This list
            will contain no more than 2 elements.

    Examples:
        >>> fibonacci_primitive_roots(11)
        [8]
        >>> fibonacci_primitive_roots(19)
        [15]
        >>> fibonacci_primitive_roots(41)
        [7, 35]

    Details:
        A Fibonacci primitive root mod p is a primitive root satisfying g**2 = g
        + 1. These can only exist for primes p = 1, 9 (mod 10). See the paper
        "Fibonacci Primitive Roots" by Shanks for details.
    """
    if p % 2 == 0:
        return []
    if p == 5:
        return [3]
    if p < 5 or jacobi_symbol(5, p) != 1:
        return []

    sqrts = sqrts_mod_p(5, p)
    inverse = inverse_mod(2, p)
    quad_roots = [(1 + root)*inverse % p for root in sqrts]
    phi_divisors = factor.prime_divisors(p - 1)
    roots = []

    for r in quad_roots:
        if is_primitive_root(r, p, phi=p - 1, phi_divisors=phi_divisors):
            roots.append(r)

    return sorted(roots)


################################################################################
# root extraction
################################################################################


def sqrts_mod_p(a, p):
    """Returns the solutions x to x**2 = a (mod p).

    Given a prime p and an integer a, this function returns a sorted list of all
    solutions to the congruence x**2 = a (mod p).

    Input:
        * a: int

        * p: int (p >= 2)
            The prime modulus. The primality of p is not verified.

    Returns:
        * solutions: list
            A sorted list of solutions to x**2 = a (mod p).

    Raises:
        * ValueError: if p < 2.

    Examples:
        >>> sqrts_mod_p(10, 13)
        [6, 7]
        >>> sqrts_mod_p(3615, 2**16 + 1)
        [367, 65170]
        >>> a = 552512556430486016984082237
        >>> sqrts_mod_p(a, 2**89 - 1)
        [1000000000000000000L, 618970018642690137449562111L]
        >>> sqrts_mod_p(10, 1)
        Traceback (most recent call last):
        ...
        ValueError: sqrts_mod_p: Must have p >= 2 be prime.

    Details:
        This function uses the algorithm due to Tonelli and Shanks. One problem
        that must be solved in this algorithm is that of finding a quadratic
        nonresidue d (mod p). There is no known efficient deterministic
        algorithm to do this, so we use a randomized algorithm to find such a d.
        This gives the algorithm of Tonelli and Shanks an expected runtime of
        O(log(p)**4).

        Our implementation follows the description given as Algorithm 2.3.8 in
        "Prime Numbers - A Computational Perspective" by Crandall and Pomerance.
        See also Section 7.1 of "Algorithmic Number Theory - Efficient
        Algorithms" by Bach and Shallit. For information about some of the
        optimizations made, see Section 1.5 of "A Course in Computational
        Algebraic Number Theory" by Cohen.
    """
    # We don't test the primality of p, but we do some simple error detection.
    if p < 2:
        raise ValueError("sqrts_mod_p: Must have p >= 2 be prime.")

    # The easy case is when p == 2.
    if p == 2:
        return [a % p]

    # Take care of the case when p | a.
    a = a % p
    if a == 0:
        return [0]

    # No solutions if a is not a square modulo p.
    if jacobi_symbol(a, p) != 1:
        return []

    # Below, p is odd, so we use Tonelli-Shanks.
    # Half of the time, there is an easy solution. Namely, for primes
    # p = 3 (mod 4) a solution is given by x = a**((p + 1)/4) (mod p).
    if p % 8 in (3, 7):
        x = pow(a, (p + 1)//4, p)

    # A slightly less trivial solution works for half of the remaining primes.
    # For p = 5 (mod 8), one can verify that a**((p - 1)/4) = 1, -1 (mod p). In
    # the positive case, then x = a**((p + 3)//8) (mod p) is a solution.
    elif p % 8 == 5:
        x = pow(a, (p + 3)//8, p)
        if x*x % p != a:
            x = x*pow(2, (p - 1)//4, p) % p

    # The remaining case is p = 1 (mod 8).
    else:
        # Find a quadratic nonresidue. Each choice of d will be a quadratic
        # nonresidue with probability close to 1/2, so we expect to find one
        # very quickly.
        d = randint(2, p - 1)
        while jacobi_symbol(d, p) != -1:
            d = randint(2, p - 1)

        # Write p - 1 = 2**s * t with t odd
        s = valuation(2, p - 1)
        t = (p - 1) >> s

        A = pow(a, t, p)
        D = pow(d, t, p)
        m = 0

        for i in xrange(s):
            if pow(A*D**m, 2**(s - 1 - i), p) == p - 1:
                m += 2**i

        x = pow(a, (t + 1)//2, p)*pow(D, m//2, p) % p

    # We want both solutions.
    solutions = sorted([x, p - x])
    return solutions


def _quadratic_lift(a, r, p, k):
    """Lifts a solution v to the congruence r**2 = a (mod p**k) to solutions
    modulo p**(k + 1).

    Suppose that for some k >= 1 we have an integer r satisfying
    r**2 = a (mod p**k). This function returns all solutions q satisfying
    q**2 = a (mod p**(k + 1)).

    Input:
        * a: int

        * r: int
            The root the the congruence r**2 = a (mod p**k).

        * p: int (p >= 2)
            A prime p, the primality of which is not verified.

        * k: int (k >= 1)
            The exponent such that r**2 = a (mod p**k).

    Returns:
        * solutions: list
            A list of all solutions c to the congruence
            c**2 = a (mod p**(k + 1))

    Raises:
        * ValueError: if k < 1 or p < 2 or r**2 != a (mod p**k).

    Examples:
        >>> 5154**2 % 11213
        119
        >>> _quadratic_lift(119, 5154, 11213, 1)
        [5869553]
        >>> 5869553**2 % 11213**2
        119
        >>> _quadratic_lift(0, 0, 7, 1)
        [0, 7, 14, 21, 28, 35, 42]

    Details:
        This function uses a simple form of Hensel lifting to obtain solutions
        modulo higher powers. See Section 12.5.2 of "A Computational
        Introduction to Number Theory and Algebra" by Shoup for details. See
        also Section 15.1 of "Introduction to Number Theory" by Hua and Theorem
        2.24 in "Fundamental Number Theory with Applications" by Mollin.
    """
    if k < 1:
        raise ValueError("_quadratic_lift: Must have k >= 1.")

    if p < 2:
        raise ValueError("_quadratic_lift: Must have p >= 2.")

    if (r*r - a) % p**k != 0:
        raise ValueError("_quadratic_lift: Must have r**2 = a (mod p**k).")

    solutions = []
    if 2*r % p != 0:
        # In this case, r lifts to a unique solution.
        inverse = inverse_mod(2*r, p)
        # Note that r**2 = a (mod p**k) so the division below is exact.
        t = (-inverse*((r*r - a)//p**k)) % p
        root = r + t*p**k
        solutions.append(root)
    else:
        # Here, r lifts to p incongruent solutions modulo p**k.
        if (r*r - a) % (p**(k + 1)) == 0:
            for t in xrange(p):
                root = r + t*p**k
                solutions.append(root)

    solutions.sort()
    return solutions


def _sqrts_mod_pk(a, p, k):
    """Returns a sorted list of the solutions x to the congruence
    x**2 = a (mod p**k).
    """
    if k == 1:
        return sqrts_mod_p(a, p)

    if p == 2:
        if k == 2:
            # For k = 2, there are two or no roots depending on whether
            # a = 1 or 3 (mod 4).
            solutions = [e for e in xrange(4) if e*e % 4 == a % 4]
        else:
            # For k > 2, there are four or no roots depending on whether
            # a = 1 (mod 8) or a != 1 (mod 8). We find all solutions (mod 8),
            # and lift to higher powers.
            solutions = [e for e in xrange(8) if e*e % 8 == a % 8]
            if k > 3 and solutions:
                for e in xrange(3, k):
                    lifted_solutions = []
                    for root in solutions:
                        lifts = _quadratic_lift(a, root, p, e)
                        lifted_solutions.extend(lifts)
                    solutions = list(lifted_solutions)
    else:
        solutions = sqrts_mod_p(a, p)
        for e in xrange(1, k):
            lifted_solutions = []
            for root in solutions:
                lifts = _quadratic_lift(a, root, p, e)
                lifted_solutions.extend(lifts)
            solutions = list(lifted_solutions)

    solutions.sort()
    return solutions


def sqrts_mod_n(a, n, n_factorization=None):
    """Returns the solutions x to the congruence x**2 = a (mod n).

    Given integers a and n, this function returns all solutions to the
    congruence x**2 = a (mod n).

    Input:
        * a: int
            A square modulo n.

        * n: int
            The modulus.

        * n_factorization: list (default=None)
            The factorization of n.

    Returns:
        * roots: list
            A sorted list of all square roots of a modulo n.

    Examples:
        >>> sqrts_mod_n(3, 11**3*23**3)
        [16152263, 28026114, 150110933, 161984784]
        >>> sqrts_mod_n(49, 53**3*61**4)
        [7, 721465236980, 1339862033577, 2061327270550]
        >>> sqrts_mod_n(-1, 5**3*7**2)
        []
    """
    if n_factorization is None:
        n_factorization = factor.factor(n)

    congruences = []
    moduli = []
    for (p, k) in n_factorization:
        roots = _sqrts_mod_pk(a, p, k)
        pk = p**k
        pairs = [(r, pk) for r in roots]
        congruences.append(pairs)
        moduli.append(pk)

    if len(congruences) == 1:
        values = [r for (r, _) in congruences[0]]
    else:
        preconditioning_data = crt_preconditioning_data(moduli)
        values = []

        for system in itertools.product(*congruences):
            values.append(chinese_preconditioned(system, preconditioning_data))

    values.sort()
    return values


def nth_roots_of_unity_mod_p(n, p, g=None):
    """Returns the nth roots of unity modulo p.

    Given a positive integer n and a prime p, this returns a list of the
    solutions to the congruence x**n = 1 (mod p).

    Input:
        * n: int
            The index of the root.

        * p: int
            Prime modulus.

        * g: int (default=None)
            Primitive root modulo p.

    Returns:
        * roots: list
            List of all nth roots of unity.

    Examples:
        >>> nth_roots_of_unity_mod_p(4, 101)
        [1, 10, 91, 100]
        >>> nth_roots_of_unity_mod_p(8, 17)
        [1, 2, 4, 8, 9, 13, 15, 16]

    Details:
        The congruence x**n == 1 (mod p) has gcd(n, p - 1) roots. We find one
        solution given by g**((p - 1)/d), where d = gcd(n, p - 1), and g is a
        primitive root modulo p. Given this, we find all other solutions by
        looking at the powers of this one solution. See Theorem 7.1 of
        "Introduction to Number Theory" by Hua for more information.
    """
    if p == 2:
        return [1]

    if g is None:
        g = primitive_root(p, phi=p-1)

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
    """Returns the nth roots of -1 modulo p.

    Given a positive integer n and a prime p, this returns a list of the
    solutions to the congruence x**n == -1 (mod p).

    Input:
        * n: int
            The index of the root.

        * p: int
            Prime modulus.

        * g: int (default=None)
            Primitive root modulo p.

    Returns:
        * roots: list
            List of all nth roots of -1 modulo p.

    Examples:
        >>> nth_roots_of_minus1_mod_p(5, 11)
        [2, 6, 7, 8, 10]
        >>> nth_roots_of_minus1_mod_p(2, 101)
        [10, 91]

    Details:
        The congruence x**n == -1 (mod p) has gcd(n, p - 1) roots. We find one
        solution to the congruence by considering g**((p - 1)/(2*d)), where d =
        gcd(n, p - 1), and g is a primitive root modulo p. We then find all
        other solutions by multiplying this solution by the nth roots of unity
        modulo p. See Section 3.7 of "Introduction to Number Theory" by Hua for
        more details.
    """
    if p == 2:
        return [1]

    d = gcd(n, p - 1)

    if (p - 1)//d % 2 != 0:
        return []

    if g is None:
        g = primitive_root(p, phi=p - 1)

    # This is one solution to x^n = -1 (mod p).
    root = pow(g, (p - 1)//(2*d), p)
    roots_of_unity = nth_roots_of_unity_mod_p(n, p, g)
    all_roots = [root*h % p for h in roots_of_unity]
    all_roots.sort()

    return all_roots


def discrete_log(a, b, p):
    """Finds a positive integer k such that a**k = b (mod p).

    Input:
        * a: int
            Base of the exponential. Typically, this will be a primitive root
            modulo p.

        * b: int
            The target.

        * p: int
            The modulus. This must be a prime. The primality of p is not
            verified.

    Returns:
        * k: int
            A positive integer such that a**k = b (mod p). If no solution
            exists, then None is returned.

    Raises:
        * ValueError: If no solution is found.

    Examples:
        >>> discrete_log(2, 6, 19)
        14
        >>> pow(2, 14, 19)
        6
        >>> discrete_log(59, 67, 113)
        11
        >>> pow(59, 11, 113)
        67
        >>> discrete_log(5, 3, 2017)
        1030
        >>> pow(5, 1030, 2017)
        3
        >>> discrete_log(2, 76, 100)
        120
        >>> discrete_log(2, 77, 100)
        Traceback (most recent call last):
        ...
        ValueError: discrete_log: No solution found.

    Details:
        The algorithm is based on the Baby-Step Giant-Step method due to Shanks.
        Using a dictionary data type, this algoritm runs in O(sqrt(n)) time and
        space. We follow the method as described in Algorithm 2.4.1 in Number
        Theory for Computing by Yan. See Chapter 5 of "Prime Numbers - A
        Computational Perspective" by Crandall and Pomerance for variants. See
        "A Computational Introduction to Number Theory and Algebra" by Shoup for
        more in depth analysis.
    """
    m = integer_sqrt(p) + 1
    am = pow(a, m, p)

    cache = {}
    val = am
    for k in xrange(1, m + 1):
        cache[val] = k
        val = (val*am) % p

    val = b
    x = 1
    for k in xrange(m + 1):
        if val in cache:
            x = cache[val]*m - k
            break
        val = (val*a) % p

    if pow(a, x, p) == b % p:
        return x
    else:
        raise ValueError("discrete_log: No solution found.")


################################################################################
# solutions to congruences
################################################################################


def linear_congruence(a, b, n):
    """Returns a list of the solutions x to the congruence a*x = b (mod n).

    Input:
        * a: int
        * b: int
        * n: int (n > 1)

    Returns:
        * solutions: list

    Raises:
        * ValueError: if n <= 1.

    Examples:
        >>> linear_congruence(10, 6, 12)
        [3, 9]
        >>> linear_congruence(12, 9, 15)
        [2, 7, 12]
        >>> linear_congruence(10, 3, 12)
        []
        >>> linear_congruence(10, 3, 0)
        Traceback (most recent call last):
        ...
        ValueError: linear_congruence: Must have n >= 2.

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
    """Returns a list of the roots of the quadratic equation modulo n.

    Input:
        * coeff_list: list
            A list of the coefficients of the quadratic.

        * n: int
            The modulus.

        * n_factorization: list (default=None)
            The factorization of the modulus n.

    Returns:
        * roots: list
            A sorted list of the roots of the quadratic modulo n.

    Examples:
        >>> quadratic_congruence([1, 3, -18], 1000)
        [3, 378, 619, 994]
        >>> quadratic_congruence([1, -31, -12], 36)
        [7, 15, 16, 24]
        >>> quadratic_congruence([11, 5, 18], 29)
        [22, 25]

    Details:
        The algorithm proceeds by finding roots modulo p**k for each prime power
        p**k in the factorization of the modulus n. These roots are then
        combined using the Chinese Remainder Theorem. See Chapter 5 of "The
        Theory of Numbers - A Text and Source Book of Problems" by Adler and
        Coury for detailed information.
    """
    if n_factorization is None:
        n_factorization = factor.factor(n)

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
