# Arithmetic Functions

import inspect

import rosemary.number_theory.factorization.factorization
import rosemary.number_theory.primes.sieves as sieves
import rosemary.number_theory.classification

from rosemary.number_theory.core import (
    bit_count,
    lcm_list,
)


################################################################################
# Utilities
################################################################################


def _validate_input(n):
    r"""Validates input of a multiplicative arithmetical function.

    This method ensures that the input `n` of a multiplicative arithmetical
    function is either an integer type or the factorization of an integer,
    given as a list of (prime, exponent) pairs.

    Parameters
    ----------
    n : integer or factorization (n > 0)
        The value of `n` can be an integer or the factorization of an
        integer given as a list of (prime, exponent) pairs.

    Returns
    -------
    factorization : list
        The validated factorization of `n`.

    Raises
    ------
    ValueError
        If ``n <= 0`` or `n` is not a valid factorization.

    OneFactorization
        If ``n == 1``, since the factoring doesn't make sense in this case.
    """
    # Name of calling function
    caller_name = inspect.stack()[1][3]

    if isinstance(n, (int, long)):
        if n == 1:
            return [(1, 1)]
        elif n <= 0:
            message = "{}: Must have n > 0.".format(caller_name)
            raise ValueError(message)
        factorization = rosemary.number_theory.factorization.factorization.factor(n)
    elif isinstance(n, list) and n[0][0] > 0:
        factorization = n
    else:
        message = "{}: Input must be a positive integer or a factorization.".format(caller_name)
        raise ValueError(message)

    return factorization


################################################################################
# Multiplicative functions
################################################################################


def euler_phi(n):
    r"""Returns the Euler totient of `n`.

    For positive integers `n`, this method returns the number of positive
    integers less than or equal to `n` which are relatively prime to `n`.

    Parameters
    ----------
    n : integer or factorization (n > 0)
        The value of `n` can be an integer or the factorization of an
        integer given as a list of (prime, exponent) pairs.

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

    for prime powers ``p**k``. For positive integers `n`, this method
    computes the factorization of `n`, and then uses this product
    definition. If instead a factorization of `n` is given, then the
    product formula is used directly. See section 2.3 of [1] and section
    16.1 of [2] for details, and see section 3.7 for information related to
    the growth rate.

    References
    ----------
    .. [1] T.M. Apostol, "Introduction to Analytic Number Theory",
    Springer-Verlag, New York, 1976.

    .. [2] G.H. Hardy, E.M. Wright, "An Introduction to the Theory of
    Numbers", 6th edition, Oxford University Press, 2008.

    Examples
    --------
    >>> euler_phi(1)
    1
    >>> euler_phi(7)
    6
    >>> euler_phi(100)
    40
    >>> euler_phi([(2, 2), (5, 2)])
    40
    >>> euler_phi(0)
    Traceback (most recent call last):
    ...
    ValueError: euler_phi: Must have n > 0.
    """
    factorization = _validate_input(n)

    if factorization == [(1, 1)]:
        return 1

    prod = 1
    for (p, e) in factorization:
        prod *= p**(e - 1)*(p - 1)
    return prod


def moebius(n):
    r"""Returns the value of the Moebius function ``mu(n)``.

    Parameters
    ----------
    n : integer or factorization (n > 0)
        The value of `n` can be an integer or the factorization of an
        integer given as a list of (prime, exponent) pairs.

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
    product formula is used directly. See section 2.2 of [1] and section
    16.3 of [2] for details, and see section 3.9 of [1] for information
    related to the growth rate.

    References
    ----------
    .. [1] T.M. Apostol, "Introduction to Analytic Number Theory",
    Springer-Verlag, New York, 1976.

    .. [2] G.H. Hardy, E.M. Wright, "An Introduction to the Theory of
    Numbers", 6th edition, Oxford University Press, 2008.

    Examples
    --------
    >>> moebius(1)
    1
    >>> moebius(35)
    1
    >>> moebius([(5, 1), (7, 1)])
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
    factorization = _validate_input(n)

    if factorization == [(1, 1)]:
        return 1

    if rosemary.number_theory.classification.is_squarefree(factorization):
        k = len(factorization)
        return (-1)**k
    else:
        return 0


def sigma(n, k=1):
    r"""Returns the sum of the `k`th powers of the divisors of `n`.

    Parameters
    ----------
    n : integer or factorization (n > 0)
        The value of `n` can be an integer or the factorization of an
        integer given as a list of (prime, exponent) pairs.

    k : integer (k >= 0) (default=1)
        The exponent in the divisor sum.

    Returns
    -------
    sum : int

    Raises
    ------
    ValueError
        If ``n <= 0`` or ``k < 0``.

    See Also
    --------
    sigma_sum, sigma_list, tau

    Notes
    -----
    The divisor sigma function is a multiplicative function satisfying
    .. math::
        \sigma(p^e, k) &= 1 + p^k + p^{2k} + \cdots + p^{ek} \\
                       &= \frac{p^{k(e + 1)} - 1}{p^{k} - 1}

    for prime powers ``p**e``. For positive integers `n`, this method
    computes the factorization of `n`, and then uses this product
    definition. If instead a factorization of `n` is given, then the
    product formula is used directly. See section 2.13 of [1] and section
    16.7 of [2] for details, and see section 3.6 of [1] for information
    related to the growth rate.

    References
    ----------
    .. [1] T.M. Apostol, "Introduction to Analytic Number Theory",
    Springer-Verlag, New York, 1976.

    .. [2] G.H. Hardy, E.M. Wright, "An Introduction to the Theory of
    Numbers", 6th edition, Oxford University Press, 2008.

    Examples
    --------
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
    """
    if k < 0:
        raise ValueError("sigma: Must have k >= 0.")
    elif k == 0:
        return tau(n)

    factorization = _validate_input(n)

    prod = 1
    for (p, e) in factorization:
        pk = p**k
        prod *= (pk**(e + 1) - 1)//(pk - 1)
    return prod


sum_of_divisors = sigma


def tau(n):
    r"""Returns the number of divisors of n.

    Parameters
    ----------
    n : integer or factorization (n > 0)
        The value of `n` can be an integer or the factorization of an
        integer given as a list of (prime, exponent) pairs.

    Returns
    -------
    tau : int
        The number of divisors of `n`.

    Raises
    ------
    ValueError
        If ``n <= 0`` or `n` is not a valid factorization.

    See Also
    --------
    sigma, tau_list

    Notes
    -----
    The divisor tau function is a multiplicative function satisfying
    .. math::
        tau(p**k) = 1 + k

    for prime powers ``p**k``. For positive integers `n`, this method
    computes the factorization of `n`, and then uses this product
    definition. If instead a factorization of `n` is given, then the
    product formula is used directly. See section 2.13 of [1] and section
    16.7 of [2] for details, and see section 3.6 of [1] for information
    related to the growth rate.

    References
    ----------
    .. [1] T.M. Apostol, "Introduction to Analytic Number Theory",
    Springer-Verlag, New York, 1976.

    .. [2] G.H. Hardy, E.M. Wright, "An Introduction to the Theory of
    Numbers", 6th edition, Oxford University Press, 2008.

    Examples
    --------
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
    """
    factorization = _validate_input(n)

    prod = 1
    for (p, e) in factorization:
        prod *= (e + 1)
    return prod


################################################################################
# Other arithmetical functions
################################################################################


def carmichael_lambda(n):
    r"""Returns the value of the Carmichael lambda function at `n`.

    This function returns the smallest positive integer `m` such that
    ``a**m = 1 (mod n)`` for all integers `a` coprime to `n`. This integer
    `m` is also called the least universal exponent for `n`.

    Parameters
    ----------
    n : integer or factorization (n > 1)
        The value of `n` can be an integer or the factorization of an
        integer given as a list of (prime, exponent) pairs.


    Returns
    -------
    m : int
        The least universal exponent for `n`.

    Raises
    ------
    ValueError
        If ``n <= 0`` or if `n` is a valid factorization.

    Notes
    -----
    The Carmichael lambda function gives the exponent of the multiplicative
    group of integers modulo n. For each n, lambda(n) divides euler_phi(n).
    This method uses the standard definition of the Carmichael lambda
    function. In particular, lambda(n) is the least common multiple of the
    values of lambda(p**e) for each prime power in the factorization of n.
    See "Fundamental Number Theory with Applications" by Mollin for
    details.

    Examples
    --------
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
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n <= 0:
            raise ValueError("carmichael_lambda: Must have n > 0.")
        n_factorization = rosemary.number_theory.factorization.factorization.factor(n)
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

    prod = 1
    for (i, prime) in sieves.primes_first_n(n):
        prod *= prime
    return prod
