# Primality tests

from rosemary.number_theory.core import jacobi_symbol
from rosemary.number_theory.factorization.algorithms import trial_division
from rosemary.number_theory.prime_list import _PRIME_LIST

from bisect import bisect_left
from math import log
from random import randint


################################################################################
# Algorithms
################################################################################


def solovay_strassen(n, a=None):
    r"""Determines if `n` is a probable prime.

    Parameters
    ----------
    n : positive integer

    Returns
    -------
    b : boolean
        Returns the value True if `n` is a probable prime, and False
        otherwise.

    See Also
    --------
    miller_rabin, solovay_strassen_deterministic, is_probable_prime

    Notes
    -----
    This method uses the Solovay Strassen test to determine the
    compositeness of `n`. The test is based on the relationship
    ``a^((p - 1)/2) = (a | p) (mod p)`` proven by Euler, where ``(a | p)``
    is the Jacobi symbol, `p` is prime, and ``gcd(a, p) == 1``.

    The Solovay Strassen test computes this for a random base `a` Since
    this is an identity for primes, this method always returns True when
    `n` is prime. On the other hand, for odd composite `n`, this test
    returns False for at least half of all `a` coprime to `n`. Thus, the
    likelihood of returning True for a composite number is less than 1/2.

    This idea may be used to create a more accurate compositeness test, and
    may even be generalized to a deterministic primality test, pending the
    Extended Riemann Hypothesis.

    References
    ----------
    .. [1] Bach, Shallitt, "Algorithmic Number Theory", MIT Press,
    Cambridge, MA, 1996

    Examples
    --------
    >>> solovay_strassen(341, 2)
    False
    >>> solovay_strassen(341, 29)
    True
    """
    if a is None:
        a = randint(2, n - 2)

    jacobi = jacobi_symbol(a, n)
    if jacobi == 0 or jacobi % n != pow(a, (n - 1)//2, n):
        return False
    return True


def miller_rabin(n, a=None):
    r"""Determines if `n` is a probable prime.

    Parameters
    ----------

    Returns
    -------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------
    """
    if a is None:
        a = randint(2, n - 2)

    t = n - 1
    s = 0
    while t % 2 == 0:
        s += 1
        t = t // 2

    b = pow(a, t, n)
    if b == 1 or b == n - 1:
        return True

    for _ in xrange(1, s):
        b = (b*b) % n
        if b == n - 1:
            return True
    return False


def miller_rabin_deterministic(n):
    """
    Determines if n is prime (pending ERH).

    Given a positive integer n, this function returns True if n is a prime;
    otherwise returns False.

    Input:
        * "n" - A positive integer.

    Output:
        * "b" - A boolean value. True if n is prime, else False.


    """
    if n < 2:
        return False

    W = min(int(2*log(n)**2), n - 1)

    for a in xrange(2, W + 1):
        if not miller_rabin(n, a):
            return False
    return True


################################################################################
# Primality tests for specific types of numbers
################################################################################


def lucas_lehmer(p):
    """
    lucas_lehmer(p):
    Given an odd prime p, this algorithm determines whether 2^p - 1 is prime.
    """
    v = 4
    mp = 2**p - 1
    for _ in xrange(1, p - 1):
        v = (v*v - 2) % mp
    return (v == 0)


################################################################################
# Generic algorithms for testing primality
################################################################################


def is_probable_prime(n):
    """Determines if n is a probable prime.

    Given a positive integer `n`, this function returns True if `n` is a
    probable prime, otherwise returns False.

    Input:
        * n: odd positive integer (n > 1)

    Output:
        * primality: bool

    Examples:
        >>> is_probable_prime(17)
        True

    Details:
        The algorithm used is based on a method popularized by J. Selfridge. See
        algorithm 3.5.2 of "Prime Numbers: A Computational Perspective" by
        Crandall and Pomerance for details.
    """
    # return all(solovay_strassen(n, a) for a in xrange(2, min(n - 1, 20)))
    return all(miller_rabin(n, a) for a in xrange(2, min(n - 1, 20)))


def is_prime(n):
    r"""Determines if `n` is a probable prime.

    Parameters
    ----------

    Returns
    -------

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------
    """
    if n < 2:
        return False
    elif n <= _PRIME_LIST[-1]:
        # check to see if the number is in our list of precomputed primes
        index = bisect_left(_PRIME_LIST, n)
        return _PRIME_LIST[index] == n
    else:
        if not is_probable_prime(n):
            return False
        # otherwise, use trial division to ensure primality
        d = trial_division(n)
        return d == n


################################################################################
# Algorithms for finding specific prime values
################################################################################


def next_prime(n):
    """
    next_prime(n):
    This returns the smallest prime > n.

    Examples:
    >>> next_prime(100)
    101
    >>> next_prime(10**10)
    10000000019
    """
    if n == 1:
        return 2

    n += 1
    if n % 2 == 0:
        n += 1

    while not is_prime(n):
        n += 2
    return n


def next_probable_prime(n):
    """
    next_probable_prime(n):
    This returns the smallest probable prime > n.

    Examples:
    >>> next_probable_prime(100)
    101
    >>> next_probable_prime(10**10)
    10000000019
    """
    if n == 1:
        return 2

    n += 1
    if n % 2 == 0:
        n += 1

    while not is_probable_prime(n):
        n += 2
    return n


def random_prime(a, b):
    """
    Returns a random prime in [a, b].
    """
    while True:
        n = randint(a, b)
        if is_prime(n):
            return n
