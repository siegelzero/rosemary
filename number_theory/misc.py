# Miscellaneous number theoretic functions

import rosemary.number_theory.factorization

from rosemary.number_theory.prime_list import _PRIME_LIST
#from fractions import Fraction


def largest_divisor(n, b):
    """
    Returns the largest divisor d of n with d <= b.

    Given a positive integer n, this returns the largest divisor d of n such
    that d <= b.

    Input:
        * "n" - A positive integer.
        * "b" - A positive integer <= n.

    Output:
        * "d" - The largest positive integer d | n with d <= b.

    Details:
        This uses a sort of "meet in the middle" approach to find the largest
        such divisor. This method gives a nice tradeoff between time and space.

    Examples:
    """
    # If b divides n, we're done.
    target = b
    if n % target == 0:
        return target

    n_factorization = rosemary.number_theory.factorization.factor(n)
    L = []
    for (p, e) in n_factorization:
        pp = p
        for _ in xrange(e):
            L.append(pp)
            pp *= p

    k = len(L) // 2
    P1 = L[:k]
    P2 = L[k:]

    # Form all of the products in each half
    prods1 = [1]
    for p in P1:
        prods1 += [x*p for x in prods1 if x*p <= b and n % (x*p) == 0]
    prods2 = [1]
    for p in P2:
        prods2 += [x*p for x in prods2 if x*p <= b and n % (x*p) == 0]

    prods1.sort()
    prods2.sort(reverse=True)

    l1 = len(prods1)
    l2 = len(prods2)
    i = 0
    j = 0
    best = 0

    # this is a meet-in-the-middle approach
    while i < l1 and j < l2:
        current = prods1[i] * prods2[j]
        if current > target:
            j += 1
        elif current < target:
            i += 1
            if current > best:
                best = current
    return best


def maximize_divisors(B):
    """
    Computes a positive integer n <= B with the maximal number of divisors,
    along with this number of divisors.

    Input:
        * B: int (B > 0)

    Output:
        * (num_divs, n): tuple
            A tuple of positive integers. The value n is the positive integer n
            with the maximal number of divisors, and the value num_divs is this
            number of divisors.

    Examples:
        >>> maximize_divisors(1000)
        (32, 840)
        >>> maximize_divisors(10**30)
        (13271040, 950542574818669103079134726400L)

    Details:
        The algorithm uses the following fact: If n <= B has the maximal number
        of divisors and both p_i**e_i | n and p_j**e_j | n with p_i < p_j, then
        e_i >= e_j. It follows that if p_j | n, then p_i | n for all primes p_i
        < p_j. We perform an exhaustive search of this solution space.
    """
    p_list = []
    pp = 1
    for prime in _PRIME_LIST:
        pp *= prime
        if pp <= B:
            p_list.append(prime)
        else:
            break

    num_primes = len(p_list)
    max_tau = [1, 1]

    def backtrack(last_i, last_e, total, tau):
        if tau > max_tau[0]:
            max_tau[0] = tau
            max_tau[1] = total

        if last_i < num_primes - 1:
            next_p = p_list[last_i + 1]
            k = 1
            t_pk = next_p
            while k <= last_e:
                t_total = total * t_pk
                if t_total <= B:
                    t_tau = tau * (k + 1)
                    backtrack(last_i + 1, k, t_total, t_tau)
                else:
                    break
                k += 1
                t_pk *= next_p

    k = 1
    while 2**k <= B:
        backtrack(0, k, 2**k, k + 1)
        k += 1
    return tuple(max_tau)


#def bernoulli_list(n):
#    """
#    Returns a list of of Bernoulli numbers B_0, B_1, B_2, ..., B_n, with the
#    usual convention B_1 = -1/2.

#    Input:
#        * n: int (n >= 0)

#    Output:
#        * L: list
#            A list of the Bernoulli numbers.

#    Examples:
#        >>> bernoulli_list(10)
#        [Fraction(1, 1), Fraction(-1, 2), Fraction(1, 6), 0, Fraction(-1, 30),
#         0, Fraction(1, 42), 0, Fraction(-1, 30), 0, Fraction(5, 66)]

#    Details:
#        The Bernoulli numbers are defined by the recursion
#            (m + 1)*B_m = -\sum_{k = 0}^{m - 1} \binom{m + 1}{k} B_k.
#        We use this recursion directly.
#    """
#    values = [0]*(n + 1)
#    values[0] = 1
#    values[1] = Fraction(-1, 2)

#    for m in xrange(2, n + 1, 2):
#        total = 1
#        # binom holds the value binomial(m + 1, k) below.
#        binom = 1
#        for k in xrange(1, m):
#            binom *= (m + 2 - k)
#            binom /= k
#            if values[k]:
#                total += binom*values[k]
#        values[m] = -total/(m + 1)
#    return values
