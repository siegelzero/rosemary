# General algorithms related to integer factorization

import collections
import itertools

import rosemary.number_theory.factorization.algorithms as algorithms
import rosemary.number_theory.classification

from rosemary.number_theory.prime_list import _PRIME_LIST
from rosemary.number_theory.primes.primality import is_probable_prime


###########################################################################
# General methods related to factorization and divisors
###########################################################################


def factor(n, skip_trial_division=False, skip_power_detection=True, use_cfrac=False):
    r"""Returns the factorization of the integer `n`.

    Parameters
    ----------
    n : int
        Integer to be factored.

    skip_trial_division : bool, optional (default=False)
        If False, performs a preliminary trial division using the
        precomputed list of primes. If True, then this step is skipped.

    skip_power_detection : bool, optional (default=True)

    use_cfrac : bool, optional (default=False)
        If True, applies the cfrac method to search for factors instead of
        the Pollard rho method. If False, then the Pollard rho method is
        used to search for factors.

    Returns
    -------
    factorization : list
        A list of pairs (p, e) consisting of the primes p and exponents e
        such that ``p**e`` divides `n`.

    Notes
    -----
    This method uses a variety of techniques to obtain a factorization of
    `n`. First, we perform trial division using a precomputed list of
    primes to find all small prime factors. After this, we apply Pollard's
    rho method to find any remaining small prime divisors.

    References
    ----------
    .. [1] H. Riesel, "Prime Numbers and Computer Methods for
    Factorization", 2nd edition, Birkhauser Verlag, Basel, Switzerland,
    1994.

    Examples
    --------
    >>> factor(100)
    [(2, 2), (5, 2)]
    >>> factor(537869)
    [(37, 1), (14537, 1)]
    >>> factor(2**91 - 1)
    [(127, 1), (911, 1), (8191, 1), (112901153L, 1), (23140471537L, 1)]
    >>> factor(10**22 + 1)
    [(89, 1), (101, 1), (1052788969L, 1), (1056689261L, 1)]
    """
    fac = collections.defaultdict(int)

    _is_prime = is_probable_prime
    if use_cfrac:
        _factor = algorithms.cfrac
    else:
        _factor = algorithms.pollard_rho_brent

    if n == 0:
        raise ValueError("Prime factorization of 0 not defined")

    # Detect perfect powers.
    if not skip_power_detection:
        ab = rosemary.number_theory.classification.is_power(n)
        if ab:
            (n, multiplier) = ab
        else:
            (n, multiplier) = (n, 1)
    else:
        (n, multiplier) = (n, 1)

    # Take care of the sign
    if n < 0:
        fac[-1] = 1
        n *= -1

    # Strip off all small factors of n
    for p in _PRIME_LIST:
        if n == 1:
            break
        elif p*p > n:
            fac[n] += 1
            n = 1
            break

        while n % p == 0:
            fac[p] += 1
            n //= p

    # Next, strip off remaining factors.
    while n > 1:
        if _is_prime(n):
            fac[n] += 1
            n = 1
        else:
            d = _factor(n)
            while not _is_prime(d):
                d = _factor(d)
            while n % d == 0:
                fac[d] += 1
                n = n // d

    factorization = sorted([(p, fac[p]*multiplier) for p in fac])
    return factorization


def factor_back(fac):
    r"""Given the factorization of an integer, this returns the integer.

    Parameters
    ----------
    fac : list
        A list of (p, e) prime, exponent pairs corresponding to the
        factorization of some integer.

    Returns
    -------
    n : int
        The product p**e for (p, e) in `fac`. This gives back the
        unfactored integer.

    Examples
    --------
    >>> factor(537869)
    [(37, 1), (14537, 1)]
    >>> factor_back([(37, 1), (14537, 1)])
    537869
    >>> factor(100)
    [[(2, 2), (5, 2)]
    >>> factor_back([(2, 2), (5, 2)])
    100
    """
    if not isinstance(fac, list):
        raise ValueError("Not a factorization in factor_back")

    pp = 1
    for (p, e) in fac:
        pp *= p**e
    return pp


def divisors(n=None, factorization=None):
    r"""Returns a sorted list of the positive integer divisors of `n`.

    Parameters
    ----------
    n : int, optional

    factorization : list, optional
        Factorization of `n` given as a list of (prime, exponent) pairs.

    Returns
    -------
    div_list : list
        Sorted list of the positive divisors of `n`.

    See Also
    --------
    prime_divisors, unitary_divisors, xdivisors

    Examples
    --------
    >>> divisors(100)
    [1, 2, 4, 5, 10, 20, 25, 50, 100]
    >>> divisors(factorization=[(2, 2), (5, 2)])
    [1, 2, 4, 5, 10, 20, 25, 50, 100]
    """
    if n is None and factorization is None:
        raise TypeError("divisors: Expected at least one argument.")
    elif factorization is None:
        factorization = factor(abs(n))
    else:
        if factorization[0][0] == -1:
            factorization = factorization[1:]

    p_divisors = [p for (p, e) in factorization]
    div_list = []
    iter_list = (xrange(e + 1) for (p, e) in factorization)

    for exp_tuple in itertools.product(*iter_list):
        prod = 1
        for i in xrange(len(exp_tuple)):
            prod *= p_divisors[i]**exp_tuple[i]
        div_list += [prod]

    return sorted(div_list)


def divisor_pairs(n=None, factorization=None):
    r"""Returns a sorted list of the positive integer divisors of `n`.

    Parameters
    ----------
    n : int, optional

    factorization : list, optional
        Factorization of `n` given as a list of (prime, exponent) pairs.

    Returns
    -------
    div_list : list
        Sorted list of the positive divisors of `n`.

    See Also
    --------
    prime_divisors, unitary_divisors, xdivisors

    Examples
    --------
    >>> divisors(100)
    [1, 2, 4, 5, 10, 20, 25, 50, 100]
    >>> divisors(factorization=[(2, 2), (5, 2)])
    [1, 2, 4, 5, 10, 20, 25, 50, 100]
    """
    if n is None and factorization is None:
        raise TypeError("divisors: Expected at least one argument.")
    elif factorization is None:
        factorization = factor(abs(n))
    else:
        if factorization[0][0] == -1:
            factorization = factorization[1:]

    div_list = divisors(n, factorization=factorization)

    pairs = []
    for d in div_list:
        if d > n//d:
            break
        pairs.append((d, n//d))

    return pairs


def prime_divisors(n=None, factorization=None):
    r"""Returns a sorted list of the primes dividing `n`.

    Parameters
    ----------
    n : int, optional

    factorization : list, optional
        Factorization of `n` given as a list of (prime, exponent) pairs.

    Returns
    -------
    p_list : list
        Sorted list of the prime divisors of `n`.

    See Also
    --------
    divisors, unitary_divisors, xdivisors

    Examples
    --------
    >>> prime_divisors(120)
    [2, 3, 5]
    >>> prime_divisors(factorization=[(2, 2), (5, 2)])
    [2, 5]
    """
    if n is None and factorization is None:
        raise TypeError("prime_divisors: Expected at least one argument.")
    elif factorization is None:
        factorization = factor(abs(n))
    else:
        if factorization[0][0] == -1:
            factorization = factorization[1:]

    return [p for (p, e) in factorization]


def xdivisors(n=None, factorization=None):
    r"""Yields the positive integer divisors of `n`.

    Parameters
    ----------
    n : int, optional

    factorization : list, optional
        Factorization of `n` given as a list of (prime, exponent) pairs.

    Yields
    -------
    divisors : iterator
        Iterator over the positive integer divisors of `n`.

    Examples
    --------
    >>> list(xdivisors(10))
    [1, 5, 2, 10]
    >>> list(xdivisors(factorization=[(2, 1), (5, 1)]))
    [1, 5, 2, 10]
    """
    if n is None and factorization is None:
        raise TypeError("xdivisors: Expected at least one argument.")
    elif factorization is None:
        factorization = factor(abs(n))
    else:
        if factorization[0][0] == -1:
            factorization = factorization[1:]

    p_divisors = [p for (p, e) in factorization]
    iter_list = (xrange(e + 1) for (p, e) in factorization)

    for exp_tuple in itertools.product(*iter_list):
        prod = 1
        for (p, e) in itertools.izip(p_divisors, exp_tuple):
            prod *= p**e
        yield prod


def unitary_divisors(n=None, factorization=None):
    if n is None and factorization is None:
        raise TypeError("unitary_divisors: Expected at least one argument.")
    elif factorization is None:
        factorization = factor(abs(n))
    else:
        if factorization[0][0] == -1:
            factorization = factorization[1:]

    pk_divisors = [p**k for (p, k) in factorization]
    div_list = [1, n]

    for i in xrange(1, len(pk_divisors)):
        for subset in itertools.combinations(pk_divisors, i):
            prod = 1
            for pk in subset:
                prod *= pk
            div_list.append(prod)

    div_list.sort()
    return div_list
