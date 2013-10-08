from rosemary.number_theory.prime_list import PRIME_LIST

import rosemary.number_theory.factorization

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
        >>> pp = primorial(15) * primorial(10)
        >>> pp
        3978148263608934751190154300L
        >>> tau(pp)
        1889568
        >>> b = 5413374899
        >>> time max([d for d in xdivisors(pp) if d <= b ])
        CPU times: user 7.82 s, sys: 0.01 s, total: 7.83 s
        Wall time: 7.83 s
        5413363329L
        >>> time largest_divisor(pp, b)
        CPU times: user 0.01 s, sys: 0.00 s, total: 0.01 s
        Wall time: 0.01 s
        5413363329
    """
    # If b divides n, we're done.
    target = b
    if n % target == 0:
        return target

    n_fac = rosemary.number_theory.factorization.factor(n)
    L = []
    for (p, e) in n_fac:
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
    maximize_divisors(B):
    This returns the integer k, 1 <= k <= B, with the most divisors.
    """
    p_list = []
    pp = 1
    for prime in PRIME_LIST:
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
    return max_tau



