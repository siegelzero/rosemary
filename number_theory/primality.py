# Primality tests

from rosemary.number_theory.prime_list import PRIME_LIST

import rosemary.number_theory.factorization

from bisect import bisect_left
from math import log
from random import randint


def is_probable_prime(n, a=None):
    """
    Determines if n is a probable prime.

    Given positive integers n and a, 1 < a < n - 1, this function returns True if n is a probable prime to the base a,
    otherwise returns False. If a is omitted, then a random value is chosen.

    Input:
        * n: odd positive integer (n > 1)
        * a: positive integer (default=None)

    Output:
        * b: bool

    Details:
        The algorithm used is based on a method popularized by J. Selfridge. See algorithm 3.5.2 of "Prime Numbers: A
        Computational Perspective" by Crandall and Pomerance for details.

    Examples:
        >>> is_probable_prime(17)
        True
    """
    if a is None:
        a = randint(2, n - 2)

    t = n - 1
    s = 0
    while t % 2 == 0:
        s += 1
        t = t // 2

    if t % 2 == 0:
        return False
    if t == 3:
        return True

    b = pow(a, t, n)
    if b == 1 or b == n - 1:
        return True

    for _ in xrange(1, s):
        b = (b*b) % n
        if b == n - 1:
            return True
    return False


def is_prime(n):
    """
    is_prime(n):
    Returns True if n is prime, otherwise returns False
    """
    if n < 2:
        return False
    elif n <= PRIME_LIST[-1]:
        # check to see if the number is in our list of precomputed primes
        index = bisect_left(PRIME_LIST, n)
        return PRIME_LIST[index] == n
    else:
        if not is_probable_prime(n):
            return False
        # otherwise, use trial division to ensure primality
        d = rosemary.number_theory.factorization.trial_division(n)
        return d == n


def is_prime_miller(n):
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
        if not is_probable_prime(n, a):
            return False
    return True


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
