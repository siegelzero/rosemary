# Preimages of Arithmetic Functions

from collections import defaultdict

import rosemary.number_theory.factorization.factorization as factor
import rosemary.number_theory.primes.primality as primality


###########################################################################
# Euler Phi Inversion
###########################################################################

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


def inverse_phi(n):
    prime_list = [a + 1 for a in factor.divisors(n) if primality.is_prime(a + 1)]
    phi_cache = defaultdict(list)

    for p in prime_list[::-1]:
        a = p - 1
        pk = p

        while a <= n:
            if n % a == 0:
                phi_cache[a] += [pk]

            for b in phi_cache.keys():
                if a*b > n:
                    continue
                for x in phi_cache[b]:
                    if x % p != 0 and n % (a*b) == 0:
                        phi_cache[a*b] += [x*pk]

            a *= p
            pk *= p

    return sorted(phi_cache[n])
