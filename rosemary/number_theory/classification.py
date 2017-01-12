# Functions for classification of different types of numbers

import rosemary.number_theory.factorization.factorization

from rosemary.number_theory.core import (
    gcd_list,
    integer_nth_root,
    integer_sqrt
)

import rosemary.number_theory.arithmetic_functions.functions


###########################################################################
# Power detection
###########################################################################


class _IsSquare(object):
    """Class for square-detection. The algorithm here is based on the Section
    1.7.2 in "A Course in Computational Algebraic Number Theory" by Cohen.
    """
    # The following tables are built when the module is imported.
    # We keep track of the squares modulo 64 and 45045 = 11*63*65.
    q64 = [0]*64
    for k in xrange(32):
        q64[k*k % 64] = 1

    q45045 = [0]*45045
    for k in xrange(22522):
        q45045[k*k % 45045] = 1

    @classmethod
    def is_square(cls, n):
        """Returns False if n is not a perfect square, and returns the square
        root of n otherwise.
        """
        # The number of squares modulo 64, 63, 65, and 11 is 12, 16, 21, and 6,
        # respectively. So if n is not a square, then the probability that this
        # will no be detected in the table lookups is 12/64*16/63*21/65*6/11 =
        # 6/715, which is less than 1 percent.
        if cls.q64[n % 64] == 0:
            return False

        if cls.q45045[n % 45045] == 0:
            return False
        else:
            q = integer_sqrt(n)
            if q*q == n:
                return q
            return False


def is_square(n):
    """Determines if n is a perfect square.

    Given a nonnegative integer n, this returns False if n is not a perfect
    square, and returns the square root of n otherwise.

    Input:
        * n: int or list (n >= 0)
            The value of n can be an int or a factorization.

    Returns:
        * sqrt: int
            If n is a perfect square, the square root of n is returned.

        * False if n is not a perfect square.

    Examples:
        >>> is_square(16)
        4
        >>> is_square(15)
        False
        >>> is_square([(2, 2), (5, 2)])
        10
        >>> is_square([(2, 3), (5, 2)])
        False

    Details:
        If n is a nonnegative integer, the function computes the integer square
        root of n, and checks if this root squared is n. If n is in factored
        form, this checks if the exponents in the prime factorization are all
        even.
    """
    if isinstance(n, list):
        if n[0][0] == -1:
            return False
        n_factorization = n[:]
        sqrt = 1
        for (p, e) in n_factorization:
            if e % 2 == 1:
                return False
            else:
                sqrt *= p**(e//2)
        return sqrt
    else:
        if n < 0:
            return False
        elif n == 0:
            return True
        return _IsSquare.is_square(n)


def is_power(n, k=None):
    """Determines if n is a perfect power.

    If n is a perfect power, this returns (b, k) where n = b^k with k maximal.
    If the optional parameter k is given, this returns (b, k) if n = b^k for
    some b. Otherwise, returns False.

    Input:
        * n: int or list (n >= 1)
            The value of n can be an int or a factorization.

        * k: int (k >= 1) (default=None)

    Returns:
        * (b, k): tuple of ints
            Values such that n = b^k if such values exist.

        * Returns False if no such values exist.

    Raises:
        * ValueError: If n <= 0 or k <= 0.

    Examples:
        >>> is_power(36)
        (6, 2)
        >>> is_power(81)
        (3, 4)
        >>> is_power(81, 2)
        (9, 2)
        >>> is_power(1330)
        False
        >>> is_power([(2, 4), (5, 4)])
        (10, 4)
        >>> is_power([(2, 4), (5, 4)], 2)
        (100, 2)
        >>> is_power([(2, 3), (5, 2)])
        False
        >>> is_power(17, 0)
        Traceback (most recent call last):
        ...
        ValueError: is_power: Must have k >= 1.
        >>> is_power(0)
        Traceback (most recent call last):
        ...
        ValueError: is_power: Must have n >= 1.

    Details:
        If n is given as an int and a value of k is given, then this function
        computes the integer kth root r of n and checks if r^k = n. In the case
        where no such k is given, this looks at each k >= 2, extracts the
        integer kth root of n, and uses this to determine if n is a kth power.

        If n is given as a factorization and a value of k is given, this
        function checks if each exponent in the prime factorization of n is
        divisible by k. If no such k is given, this computes the gcd d of the
        exponents in the prime factorization and checks if d > 1.
    """
    if k is not None and k < 1:
        raise ValueError("is_power: Must have k >= 1.")

    if isinstance(n, list):
        if n[0][0] == -1:
            raise ValueError("is_power: Must have n >= 1.")
        n_factorization = n[:]

        kth_root = 1
        if k is not None:
            for (p, e) in n_factorization:
                if e % k == 0:
                    kth_root *= p**(e//k)
                else:
                    return False
            return (kth_root, k)
        else:
            exponents = [e for (p, e) in n_factorization]
            d = gcd_list(exponents)
            if d == 1:
                return False
            else:
                kth_root = 1
                for (p, e) in n_factorization:
                    kth_root *= p**(e//d)
                return (kth_root, d)
    else:
        if n < 1:
            raise ValueError("is_power: Must have n >= 1.")

        if k is not None:
            if k == 1:
                return (n, 1)
            root = integer_nth_root(k, n)
            if root**k == n:
                return (root, k)
            return False

        for k in xrange(2, n.bit_length() + 1):
            base = integer_nth_root(k, n)
            if base**k == n:
                exponent = k

                while True:
                    next_step = is_power(base)
                    if next_step:
                        base = next_step[0]
                        exponent *= next_step[1]
                    else:
                        return (base, exponent)

        return False


def is_squarefree(n):
    r"""Determines if n is squarefree.

    Given a nonnegative integer n, this return True iff n is not divisible
    by the square of an integer > 1.

    Parameters
    ----------
    n : int, optional

    factorization : list, optional
        Factorization of `n` given as a list of (prime, exponent) pairs.

    Returns
    -------
    b : bool
        True if `n` is squarefree; False otherwise.

    Notes
    -----
    If `n` is a nonnegative integer, this factors `n` and checks if `n` is
    divisible by the square of a prime.  If `n` is in factored form, this
    directly checks the prime factorization.

    Examples
    --------
    >>> is_squarefree(35)
    True
    >>> is_squarefree(100)
    False
    >>> is_squarefree(factorization=[(2, 2), (5, 2)])
    False
    """
    n = abs(n)
    factorization = rosemary.number_theory.arithmetic_functions.functions._validate_input(n)

    return all(e == 1 for (p, e) in factorization)
