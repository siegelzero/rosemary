from rosemary.number_theory.prime_list import PRIME_LIST
from rosemary.number_theory.elementary import (gcd, integer_sqrt,
        integer_nth_root, integer_log, trial_division, power_mod, randint,
        primes, jacobi_symbol, is_square)

def fermat(n):
    """
    Returns a nontrivial factor, or proves primality.

    Given an odd integer n > 1, this algorithm attempts to find a factor of n
    using Fermat's method.

    Input:
        * "n" - An odd integer > 1

    Output:
        * "d" - If n == d, then n is proven prime. Otherwise, d is a nontrivial
          factor of n.

    Details:
        The algorithm used here is Fermat's method. This algorithm runs in time
        O((n + 1) / 2 - sqrt(n)), where n is the number to be factored.  In the
        worst case (n = 3*p), this algorithm is much worse than trial division.
        See section 5.1.1 of "Prime Numbers - A Computation Perspective" by
        Crandall and Pomerance for more details.

    Examples:
        >>> m = 1112470797641561909
        >>> factor_fermat(m)
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
    Returns a nontrivial factor, or proves primality.

    Given an integer n >= 3, this algorithm finds a nontrivial factor of n if n
    is not prime, or returns n if n is prime.

    Input:
        * "n" - An integer >= 3

    Output:
        * "d" - If d == n, then n is proven prime. Otherwise, d is a nontrivial
          factor of n.

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
        >>> lehman_factor(l)
        1056689261L
    """
    # first, we trial divide up to floor(n^(1/3))
    bound = integer_nth_root(3, n)
    d = trial_division(n, bound)
    # if a non-trivlal divisor is found, return it
    if d < n:
        return d

    # loop until k > floor(n^(1/3))
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
        * "n" - A composite integer.
        * "B" - A positive integer (default = 20000).

    Output:
        * "d" - A (hopefully) nontrivial factor of n.

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
        >>> factor_pollard_p_1(m)
        1056689261L
    """
    c = randint(2, 20)
    p_list = primes(B)
    for p in p_list:
        a = integer_log(B, p)
        for i in xrange(a):
            c = power_mod(c, p, n)
    g = gcd(c - 1, n)
    return g

def pollard_rho(n):
    """
    Attempts to find a nontrivial factor of n.

    Given a composite number n, this algorithm attempts to find a nontrivial
    factor of n.

    Input:
        * "n" - A composite integer.

    Outout:
        * "d" - A (hopefully) nontrivial factor of n.

    Details:
        The algorithm used is Pollard's rho factorization method. This method
        should return a nontrivial factor of n in O(sqrt(p)) steps, where p is
        the least prime factor of n. For more details, see section 5.2.1 of
        "Prime Numbers - A Computational Perspective" by Crandall and Pomerance
        and section 8.5 of "A Course in Computational Algebraic Number Theory"
        by Cohen.

    Examples:
        >>> m = 1112470797641561909
        >>> factor_pollard_rho(m)
        1052788969L
        >>> m = 2175282241519502424792841
        >>> factor_pollard_rho(m)
        513741730823L
    """
    # Choose seeds.
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
            # We've found a nontrivial factor.
            return g
        elif g == n:
            # In this case, we choose new seeds.
            return n

def one_line_factor(k, M):
    # first, we trial divide up to floor(n^(1/3))
    bound = integer_nth_root(3, k)
    print "trial_division to {0}".format(bound)
    d = trial_division(k, bound)
    # if a non-trivlal divisor is found, return it
    if d < k:
        print "found"
        return d

    n = 480*k
    for i in xrange(1, M + 1):
        s = int((n*i)**(0.5)) + 1
        m = s*s % n
        if is_square(m):
            t = integer_sqrt(m)
            g = gcd(k, s - t)
            return g
    return k

