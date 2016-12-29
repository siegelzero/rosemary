################################################################################
# Arithmetic Functions
################################################################################

from rosemary.number_theory.prime_list import _PRIME_LIST

from rosemary.number_theory.core import (
    bit_count,
    lcm_list,
)

import rosemary.number_theory.factorization.factorization as factor
import rosemary.number_theory.primes.primality as primality


################################################################################
# Multiplicative functions
################################################################################


def euler_phi(n=None, factorization=None):
    r"""Returns the Euler totient of `n`.

    For positive integers `n`, this method returns the number of positive
    integers less than or equal to `n` which are relatively prime to `n`.

    Parameters
    ----------
    n : int, optional

    factorization : list, optional
        Factorization of `n` given as a list of (prime, exponent) pairs.

    Returns
    -------
    phi : int
        The value euler_phi(n)

    Raises
    ------
    ValueError
        If ``n <= 0``.

    See Also
    --------
    euler_phi_sum, euler_phi_list

    Notes
    -----
    The Euler phi function is a multiplicative function satisfying
    .. math::
        \phi(p^k) = (p - 1) p^{k - 1}

    for prime powers ``p^k``. For positive integers `n`, this method
    computes the factorization of `n`, and then uses this product
    definition. If instead a factorization of `n` is given, then the
    product formula is used directly. See section 2.3 of [1] for details,
    and see section 3.7 for information related to the growth rate.

    References
    ----------
    .. [1] T.M. Apostol, "Introduction to Analytic Number Theory",
    Springer-Verlag, New York, 1976.

    Examples
    --------
    >>> euler_phi(1)
    1
    >>> euler_phi(7)
    6
    >>> euler_phi(100)
    40
    >>> euler_phi(factorization=[(2, 2), (5, 2)])
    40
    >>> euler_phi(0)
    Traceback (most recent call last):
    ...
    ValueError: euler_phi: Must have n > 0.
    """
    if n is None and factorization is None:
        raise TypeError("euler_phi: Expected at least one argument.")
    elif factorization is None:
        if n == 1:
            return 1
        elif n <= 0:
            raise ValueError("euler_phi: Must have n > 0.")
        factorization = factor.factor(n)
    else:
        if factorization[0][0] == -1:
            raise ValueError("euler_phi: Must have n > 0.")

    prod = 1
    for (p, e) in factorization:
        prod *= p**(e - 1)*(p - 1)
    return prod


def moebius(n=None, factorization=None):
    r"""Returns the value of the Moebius function ``mu(n)``.

    Parameters
    ----------
    n : int, optional

    factorization : list, optional
        Factorization of `n` given as a list of (prime, exponent) pairs.

    Returns
    -------
    mu : int
        The value moebius(n).

    Raises
    ------
    ValueError
        If ``n <= 0``.

    See Also
    --------
    moebius_sum, moebius_list

    Notes
    -----
    The Moebius mu function is a multiplicative function satisfying
    .. math::
        mu(p_1 p_2 \cdots p_k) = (-1)^k

    and ``mu(p**e) = 0`` for ``e >= 2``. For positive integers `n`, this
    method computes the factorization of `n`, and then uses this product
    definition. If instead a factorization of `n` is given, then the
    product formula is used directly. See section 2.2 of [1] for details,
    and see section 3.9 for information related to the growth rate.

    References
    ----------
    .. [1] T.M. Apostol, "Introduction to Analytic Number Theory",
    Springer-Verlag, New York, 1976.

    Examples
    --------
    >>> moebius(1)
    1
    >>> moebius(35)
    1
    >>> moebius(factorization=[(5, 1), (7, 1)])
    1
    >>> moebius(100)
    0
    >>> moebius(2)
    -1
    >>> moebius(-1)
    Traceback (most recent call last):
    ...
    ValueError: moebius: Must have n > 0.
    """
    if n is None and factorization is None:
        raise TypeError("euler_phi: Expected at least one argument.")
    elif factorization is None:
        if n == 1:
            return 1
        elif n <= 0:
            raise ValueError("moebius: Must have n > 0.")
        factorization = factor.factor(n)
    else:
        if factorization[0][0] == -1:
            raise ValueError("moebius: Must have n > 0.")

    if factor.is_squarefree(factorization=factorization):
        k = len(factorization)
        return (-1)**k
    else:
        return 0


def sigma(n, k=1):
    """Returns the sum of the kth powers of the divisors of n.

    Parameters
        * n: int or list (n > 0)
            The value of n can be an int or a factorization.

        * k: int (k >= 0) (default=1)

    Returns:
        * s: int

    Raises:
        * ValueError: If n <= 0 or k < 0.

    Examples:
        >>> sigma(9)
        13
        >>> 1 + 3 + 9
        13
        >>> sigma([(3, 2)])
        13
        >>> sigma(10, 2)
        130
        >>> 1**2 + 2**2 + 5**2 + 10**2
        130
        >>> sigma(10, 0)
        4
        >>> sigma(-1)
        Traceback (most recent call last):
        ...
        ValueError: sigma: Must have n > 0.
        >>> sigma(10, -1)
        Traceback (most recent call last):
        ...
        ValueError: sigma: Must have k >= 0.

    Details:
        The divisor sigma function is a multiplicative function satisfying
        sigma(p**e, k) = 1 + p**k + p**(2*k) + ... + p**(e*k) for prime powers
        p**e. For positive integers n, this method computes the factorization of
        n, and then uses this product definition. If instead a factorization
        of n is given, then the product formula is used directly.
    """
    if k < 0:
        raise ValueError("sigma: Must have k >= 0.")
    elif k == 0:
        return tau(n)

    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n <= 0:
            raise ValueError("sigma: Must have n > 0.")
        n_factorization = factor.factor(n)
    elif isinstance(n, list) and n[0][0] > 0:
        n_factorization = n
    else:
        raise ValueError("sigma: Input must be a positive integer or a factorization.")

    prod = 1

    for (p, e) in n_factorization:
        pk = p**k
        prod *= (pk**(e + 1) - 1)//(pk - 1)

    return prod


def tau(n):
    """Returns the number of divisors of n.

    Parameters
        * n: int (n > 0)
            The value of n can be an int or a factorization.

    Returns:
        * s: int

    Raises:
        * ValueError: If n <= 0 or if n is not an int or factorization.

    Examples:
        >>> tau(9)
        3
        >>> tau([(2, 2), (5, 2)])
        9
        >>> tau(100)
        9
        >>> tau('cat')
        Traceback (most recent call last):
        ...
        ValueError: tau: Input must be a positive integer or a factorization.

        >>> tau(-1)
        Traceback (most recent call last):
        ...
        ValueError: tau: Must have n > 0.

    Details:
        The divisor tau function is a multiplicative function satisfying
        tau(p**k) = 1 + k for prime powers p**k. For positive integers n, this
        method computes the factorization of n, and then uses this product
        definition. If instead a factorization of n is given, then the product
        formula is used directly.
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n <= 0:
            raise ValueError("tau: Must have n > 0.")
        n_factorization = factor.factor(n)
    elif isinstance(n, list) and n[0][0] > 0:
        n_factorization = n
    else:
        raise ValueError("tau: Input must be a positive integer or a factorization.")

    prod = 1
    for (p, e) in n_factorization:
        prod *= (e + 1)
    return prod


################################################################################
# Other arithmetical functions
################################################################################


def carmichael_lambda(n):
    """Returns the smallest positive integer m such that a**m = 1 (mod n) for
    all integers a coprime to n.

    Parameters
        * n: int or list (n > 1)
            The modulus. This value can be an int or a factorization.

    Returns:
        * m: int

    Raises:
        * ValueError: If n <= 0 or n is not an integer or factorization.

    Examples:
        >>> carmichael_lambda(100)
        20
        >>> carmichael_lambda([(2, 2), (5, 2)])
        20
        >>> carmichael_lambda(113)
        112
        >>> carmichael_lambda(0)
        Traceback (most recent call last):
        ...
        ValueError: carmichael_lambda: Must have n > 0.

    Details:
        The Carmichael lambda function gives the exponent of the multiplicative
        group of integers modulo n. For each n, lambda(n) divides euler_phi(n).
        This method uses the standard definition of the Carmichael lambda
        function. In particular, lambda(n) is the least common multiple of the
        values of lambda(p**e) for each prime power in the factorization of n.
        See "Fundamental Number Theory with Applications" by Mollin for details.
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n <= 0:
            raise ValueError("carmichael_lambda: Must have n > 0.")
        n_factorization = factor.factor(n)
    elif isinstance(n, list) and n[0][0] > 0:
        n_factorization = n
    else:
        raise ValueError("carmichael_lambda: Input must be a positive integer or a factorization.")

    terms = []
    for (p, e) in n_factorization:
        if p != 2:
            term = p**(e - 1)*(p - 1)
        else:
            if e <= 2:
                term = 1 << (e - 1)
            else:
                term = 1 << (e - 2)
        terms.append(term)

    value = lcm_list(terms)
    return value


def factorial(n):
    """Returns the factorial of n.

    Parameters
        * n: int (n >= 0)

    Returns:
        * prod: int

    Raises:
        * ValueError: If n < 0.

    Examples:
        >>> factorial(10)
        3628800
        >>> factorial(40)
        815915283247897734345611269596115894272000000000L
        >>> factorial(-1)
        Traceback (most recent call last):
        ...
        ValueError: factorial: Must have n >= 0.

    Details:
        The algorithm used is a binary splitting method. See Section 10.2.3 of
        "Pi and the AGM" by Borwein and Borwein for information about this
        method.
    """
    if n < 0:
        raise ValueError("factorial: Must have n >= 0.")
    elif n in (0, 1):
        return 1

    prod = 1
    for k in xrange(1, n.bit_length()):
        lower = (n >> k) + 1
        upper = (n >> (k - 1)) + 1

        if lower % 2 == 0:
            lower += 1

        partial = 1
        for j in xrange(lower, upper, 2):
            partial *= j
        prod *= (partial**k)

    return prod << (n - bit_count(n))


def primorial(n):
    """Returns the product of the first n primes p1, p2, ..., pn.

    Parameters
        * n: int (n > 0)

    Returns:
        * prod: int

    Raises:
        * ValueError: If n < 1 or n >= 78498.

    Examples:
        >>> primorial(3)
        30
        >>> primorial(10)
        6469693230
        >>> primorial(0)
        Traceback (most recent call last):
        ...
        ValueError: primorial: Must have n >= 1.
        >>> primorial(10000000)
        Traceback (most recent call last):
        ...
        ValueError: primorial: Must have n < 78498.

    Details:
        The current implementation uses the precomputed list of primes in the
        computation. Because of this, the input n is limited to be no larger
        than the number of precomputed primes. This shouldn't be a problem in
        most cases, since this number is quite large.
    """
    if n <= 0:
        raise ValueError("primorial: Must have n >= 1.")

    num_primes = len(_PRIME_LIST)
    if n >= num_primes:
        raise ValueError("primorial: Must have n < {}.".format(num_primes))

    prod = 1
    for (i, prime) in enumerate(_PRIME_LIST):
        if i == n:
            return prod
        prod *= prime


def euler_phi_inverse(n):
    """Returns a sorted list of all positive integers k such that euler_phi(k) = n.

    Parameters
        * n: int or list (n > 0)
            The value of n can be an int or a factorization.

    Returns:
        * values: list

    Raises:
        * ValueError: If n <= 0 or n is not an int or factorization.

    Examples:
        >>> euler_phi_inverse(40)
        [41, 55, 75, 82, 88, 100, 110, 132, 150]
        >>> euler_phi_inverse(100)
        [101, 125, 202, 250]
        >>> euler_phi_inverse([(2, 2), (5, 2)])
        [101, 125, 202, 250]
        >>> euler_phi_inverse(103)
        []
        >>> euler_phi_inverse(-1)
        Traceback (most recent call last):
        ...
        ValueError: euler_phi_inverse: Must have n > 0.

    Details:
        This uses the algorithm outlined in "Discovering Mathematics with Magma"
        by Bosma and Cannon. For a different method, see the paper "Euler's
        Totient Function and its Inverse" by Gupta.
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return [1, 2]
        elif n > 1 and n % 2 == 1:
            # We know euler_phi(n) is even for n >= 3
            return []
        elif n <= 0:
            raise ValueError("euler_phi_inverse: Must have n > 0.")
        n_factorization = factor.factor(n)
    elif isinstance(n, list) and n[0][0] > 0:
        n_factorization = n
        n = factor.factor_back(n_factorization)
    else:
        raise ValueError("euler_phi_inverse: Input must be a positive integer or a factorization.")

    powers_of_two = set([1])
    if n % 2 == 0:
        for i in xrange(1, n_factorization[0][1] + 1):
            powers_of_two.add(2**i)

    # Odd primes p that divide n must have the property that p - 1 divides m.
    prime_list = []
    for d in factor.divisors(factorization=n_factorization)[1:]:
        if primality.is_prime(d + 1):
            prime_list.append(d + 1)

    # Here, we store pairs (a, b) such that a is odd and phi(a) = m/b, with b
    # even or 1. Every pair contributes at least one solution. When b = 1, the
    # pair contributes two solutions.
    pairs = [(1, n)]
    for p in reversed(prime_list):
        new_pairs = []
        for (a, b) in pairs:
            if b == 1:
                continue
            pk = p
            d = b//(p - 1)
            mmod = b % (p - 1)
            while mmod == 0:
                if d % 2 == 0 or d == 1:
                    new_pairs.append((pk*a, d))
                pk *= p
                mmod = d % p
                d = d//p
        pairs.extend(new_pairs)

    # When b = 2^k, we have the solution 2*b*a.
    values = []
    for (a, b) in pairs:
        if b in powers_of_two:
            values.append(2*b*a)
            if b == 1:
                values.append(a)

    values.sort()
    return values
