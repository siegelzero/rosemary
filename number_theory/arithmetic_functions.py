# Arithmetic Functions

from rosemary.number_theory.prime_list import PRIME_LIST

from rosemary.number_theory.core import (
    bit_count,
    lcm_list,
)

import rosemary.number_theory.factorization
import rosemary.number_theory.primality
import rosemary.number_theory.sieves


################################################################################
# multiplicative functions
################################################################################


def euler_phi(n):
    """Returns the number of positive integers <= n that are coprime to n.

    Input:
        * n: int or list (n > 0)
            The value of n can be an int or a factorization.

    Returns:
        * r: int

    Raises:
        * ValueError: if n <= 0 or n is not an int or factorization.

    Examples:
        >>> euler_phi(100)
        40
        >>> euler_phi([(2, 2), (5, 2)])
        40
        >>> euler_phi(11213)
        11212
        >>> euler_phi(0)
        Traceback (most recent call last):
        ...
        ValueError: euler_phi: Must have n > 0.

    Details:
        The Euler phi function is a multiplicative function satisfying phi(p**k)
        = (p - 1)*p**(k - 1) for prime powers p**k. For positive integers n,
        this method computes the factorization of n, and then uses this product
        definition. If instead a factorization of n is given, then the product
        formula is used directly. See "Fundamental Number Theory with
        Applications" by Mollin for details.
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n <= 0:
            raise ValueError("euler_phi: Must have n > 0.")
        n_factorization = rosemary.number_theory.factorization.factor(n)
    elif isinstance(n, list) and n[0][0] > 0:
        n_factorization = n
    else:
        raise ValueError("euler_phi: Input must be a positive integer or a factorization.")

    prod = 1
    for (p, e) in n_factorization:
        prod *= p**(e - 1)*(p - 1)
    return prod


def moebius(n):
    """Returns the value of the Moebius function mu(n).

    Input:
        * n: int or list (n > 0)
            The value of n can be an int or a factorization.

    Returns:
        * mu: int

    Raises:
        * ValueError: If n <= 0 or n is not an int or factorization.

    Examples:
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

    Details:
        The Moebius mu function is a multiplicative function satisfying
        mu(p1*p2*...pk) = (-1)**k, and mu(p**e) = 0 for e >= 2. For positive
        integers n, this method computes the factorization of n, and then uses
        this product definition. If instead a factorization of n is given, then
        the product formula is used directly.
    """
    if isinstance(n, (int, long)):
        if n == 1:
            return 1
        elif n <= 0:
            raise ValueError("moebius: Must have n > 0.")
        n_factorization = rosemary.number_theory.factorization.factor(n)
    elif isinstance(n, list) and n[0][0] > 0:
        n_factorization = n
    else:
        raise ValueError("moebius: Input must be a positive integer or a factorization.")

    if rosemary.number_theory.factorization.is_squarefree(n_factorization):
        k = len(n_factorization)
        return (-1)**k
    else:
        return 0


def sigma(n, k=1):
    """Returns the sum of the kth powers of the divisors of n.

    Input:
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
        n_factorization = rosemary.number_theory.factorization.factor(n)
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

    Input:
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
        n_factorization = rosemary.number_theory.factorization.factor(n)
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

    Input:
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
        n_factorization = rosemary.number_theory.factorization.factor(n)
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

    Input:
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

    Input:
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

    num_primes = len(PRIME_LIST)
    if n >= num_primes:
        raise ValueError("primorial: Must have n < {}.".format(num_primes))

    prod = 1
    for (i, prime) in enumerate(PRIME_LIST):
        if i == n:
            return prod
        prod *= prime


################################################################################
# Lists of values of multiplicative functions
################################################################################


def euler_phi_list(n):
    """Returns a list of values of euler_phi(k) for 1 <= k <= n.

    Input:
        * n: int (n > 0)

    Returns:
        * values: list
            This is a list of values of euler_phi(k) for 1 <= k <= n. The list
            begins with 0, so that L[k] holds the value of euler_phi(k).

    Raises:
        * ValueError: If n <= 0.

    Examples:
        >>> euler_phi_list(10)
        [0, 1, 1, 2, 2, 4, 2, 6, 4, 6, 4]
        >>> [0] + [euler_phi(k) for k in xrange(1, 11)]
        [0, 1, 1, 2, 2, 4, 2, 6, 4, 6, 4]
        >>> euler_phi_list(0)
        Traceback (most recent call last):
        ...
        ValueError: euler_phi_list: Must have n > 0.

    Details:
        This function initially creates a list of (n + 1) zeros, and fills the
        list by sieving, using the product definition of euler_phi(n).
    """
    if n <= 0:
        raise ValueError("euler_phi_list: Must have n > 0.")

    block = [1]*(n + 1)
    block[0] = 0
    prime_list = rosemary.number_theory.sieves.primes(n)

    for p in prime_list:
        for mul in xrange(p, n + 1, p):
            block[mul] *= (p - 1)
        pk = p*p
        while pk <= n:
            for mul in xrange(pk, n + 1, pk):
                block[mul] *= p
            pk *= p

    return block


totient_list = euler_phi_list


def moebius_list(n):
    """Return a list of the values of moebius(k) for 1 <= k <= n.

    Input:
        * n: int (n > 0)

    Returns:
        * L: list
            This is a list of values of moebius(k) for 1 <= k <= n. The list
            begins with 0, so that L[k] holds the value of moebius(k).

    Raises:
        * ValueError: If n <= 0.

    Examples:
        >>> moebius_list(10)
        [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1]
        >>> [0] + [moebius(k) for k in xrange(1, 11)]
        [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1]
        >>> moebius_list(0)
        Traceback (most recent call last):
        ...
        ValueError: moebius_list: Must have n > 0.

    Details:
        This function creates a list of (n + 1) elements, and fills the list by
        sieving and using the product definition of moebius(n).
    """
    if n <= 0:
        raise ValueError("moebius_list: Must have n > 0.")

    sqrt = int(n**(0.5))
    prime_list = rosemary.number_theory.sieves.primes(sqrt)

    values = [1]*(n + 1)
    block = [1]*(n + 1)
    block[0] = 0

    for p in prime_list:
        for i in xrange(p, n + 1, p):
            if i % (p*p) == 0:
                block[i] = 0
            else:
                block[i] *= -1
                values[i] *= p

    for i in xrange(n + 1):
        if block[i] and values[i] < i:
            block[i] *= -1

    return block


def sigma_list(n, k=1):
    """Return a list of the values of sigma(i, k) for 1 <= i <= n.

    Input:
        * n: int (n > 0)
        * k: int (k > 0) (default=1)

    Returns:
        * L: list
            This is a list of values of sigma(i, k) for 1 <= i <= n. The list
            begins with 0, so that L[i] holds the value of sigma(i, k).

    Raises:
        * ValueError: If n <= 1 or k < 0.

    Examples:
        >>> sigma_list(10)
        [0, 1, 3, 4, 7, 6, 12, 8, 15, 13, 18]
        >>> [0] + [sigma(k) for k in xrange(1, 11)]
        [0, 1, 3, 4, 7, 6, 12, 8, 15, 13, 18]
        >>> sigma_list(10, 2)
        [0, 1, 5, 10, 21, 26, 50, 50, 85, 91, 130]
        >>> [0] + [sigma(k, 2) for k in xrange(1, 11)]
        [0, 1, 5, 10, 21, 26, 50, 50, 85, 91, 130]
        >>> sigma_list(-1)
        Traceback (most recent call last):
        ...
        ValueError: sigma_list: Must have n > 0.
        >>> sigma_list(10, -1)
        Traceback (most recent call last):
        ...
        ValueError: sigma_list: Must have k >= 0.

    Details:
        This function creates a list of (n + 1) elements, and fills the list by
        sieving and using the product definition of sigma(n, k). For k == 1 and
        small values of n, a simpler but faster algorithm is used. See Section
        9.8 in "Algorithmic Number Theory - Efficient Algorithms" by Bach and
        Shallit for details.
    """
    if n <= 0:
        raise ValueError("sigma_list: Must have n > 0.")

    if k == 0:
        return tau_list(n)
    elif k < 0:
        raise ValueError("sigma_list: Must have k >= 0.")

    # Use a special algorithm when k == 1 for small values of n.
    if k == 1 and n <= 10**4:
        block = [0]*(n + 1)
        sqrt = int(n**(0.5))

        for j in xrange(1, sqrt + 1):
            block[j*j] += j
            for k in xrange(j + 1, n//j + 1):
                block[k*j] += j + k

    else:
        block = [1]*(n + 1)
        block[0] = 0
        prime_list = rosemary.number_theory.sieves.primes(n)

        for p in prime_list:
            pk = p
            mul = p**k
            term = mul**2
            last = mul
            while pk <= n:
                for idx in xrange(pk, n + 1, pk):
                    block[idx] *= (term - 1)
                    block[idx] /= (last - 1)
                pk *= p
                last = term
                term *= mul

    return block


def tau_list(n):
    """Returns a list of the values of tau(i), for 1 <= i <= n.

    Input:
        * n: int (n > 0)

    Returns:
        * L: list
            This is a list of values of tau(k) for 1 <= k <= n. The list begins
            with 0, so that L[k] holds the value of tau(k).

    Raises:
        * ValueError: If n <= 0.

    Examples:
        >>> tau_list(10)
        [0, 1, 2, 2, 3, 2, 4, 2, 4, 3, 4]
        >>> [0] + [tau(k) for k in xrange(1, 11)]
        [0, 1, 2, 2, 3, 2, 4, 2, 4, 3, 4]
        >>> tau_list(-1)
        Traceback (most recent call last):
        ...
        ValueError: tau_list: Must have n > 0.

    Details:
        The algorithm used here comes from Section 9.8 in "Algorithmic Number
        Theory - Efficient Algorithms" by Bach and Shallit.
    """
    if n <= 0:
        raise ValueError("tau_list: Must have n > 0.")

    block = [0]*(n + 1)
    sqrt = int(n**(0.5))

    for j in xrange(1, sqrt + 1):
        block[j*j] += 1
        for k in xrange(j + 1, n//j + 1):
            block[k*j] += 2
    return block


################################################################################
# Summatory functions of multiplicative functions
################################################################################


def euler_phi_sum(n):
    """Returns the value of the sum of euler_phi(k) for k = 1, 2, ..., n.

    Input:
        * n: int (n > 0)

    Returns:
        * s: int

    Raises:
        * ValueError: If n <= 0.

    Examples:
        >>> euler_phi_sum(10)
        32
        >>> sum(euler_phi_list(10))
        32
        >>> euler_phi_sum(10**5)
        3039650754
        >>> euler_phi_sum(0)
        Traceback (most recent call last):
        ...
        ValueError: euler_phi_sum: Must have n > 0.

    Details:
        Let S(n) = \sum_{k = 1}^{n} \phi(k). We use the identity
        S(n) = n*(n + 1)/2 - \sum_{d = 2}^{n} S(n/d) to compute the value. By
        applying a variant of Dirichlet's Hyperbola Method, we are able to cut
        the limits of summation, and by memoizing the recursion, we are able to
        achieve a sublinear runtime. See Section 4.2 of "The Prime Numbers and
        Their Distribution" by Mendes, et al for the Hyperbola Method. See
        section 4.9 equation 4.60 in "Concrete Mathematics" by Graham, et al for
        a proof of the identity for S(n).
    """
    if n <= 0:
        raise ValueError("euler_phi_sum: Must have n > 0.")

    sqrt = int(n**(0.5))
    phi_list = euler_phi_list(sqrt)

    # for i <= sqrt(n), compute the sum directly
    cache = {1: 1}
    for i in xrange(2, sqrt + 1):
        cache[i] = cache[i - 1] + phi_list[i]

    def S(n):
        if n in cache:
            return cache[n]

        sqrt_n = int(n**(0.5))
        value = n - sqrt_n*cache[sqrt_n]

        for d in xrange(2, sqrt_n + 1):
            value += (S(n//d) + phi_list[d]*(n//d))

        cache[n] = n*(n + 1)//2 - value
        return cache[n]

    return S(n)


totient_sum = euler_phi_sum


def euler_phi_weighted_sum(n):
    """Returns the value of the sum of k*euler_phi(k) for k = 1, 2, ..., n.

    Input:
        * n: int (n > 0)

    Returns:
        * s: int

    Raises:
        * ValueError: If n <= 0.

    Examples:
        >>> euler_phi_weighted_sum(100)
        203085
        >>> sum(k*euler_phi(k) for k in xrange(1, 101))
        203085
        >>> euler_phi_weighted_sum(0)
        Traceback (most recent call last):
        ...
        ValueError: euler_phi_weighted_sum: Must have n > 0.

    Details:
        Let T(n) = \sum_{k = 1}^{n} k*euler_phi(k). We use the relationship
        \sum_{d = 1}^{n} d*T(n/d) = n*(n + 1)*(2*n + 1) / 6 to form a recurrence
        for T(n). Using a form of Dirichlet's hyperbola method and memoizing the
        recursion gives us a sublinear runtime.
    """
    if n <= 0:
        raise ValueError("euler_phi_weighted_sum: Must have n > 0.")

    sqrt = int(n**(0.5))
    totients = euler_phi_list(sqrt)

    # for i <= sqrt(n), compute the sum directly
    cache = {1: 1}
    for i in xrange(2, sqrt + 1):
        cache[i] = cache[i - 1] + i*totients[i]

    def T(n):
        if n in cache:
            return cache[n]

        sqrt_n = int(n**(0.5)) + 1
        s1 = sum(d*totients[d]*(n//d - d + 1)*(n//d + d)//2 for d in xrange(1, sqrt_n))
        s2 = sum(d*(T(n//d) - T(d - 1)) for d in xrange(2, sqrt_n))
        s3 = sum(d*d*totients[d] for d in xrange(1, sqrt_n))

        cache[n] = n*(n + 1)*(2*n + 1)//6 - (s1 + s2 - s3)
        return cache[n]

    return T(n)


def moebius_sum(n):
    """Returns the value of the sum of moebius(k) for k = 1, 2, ..., n.

    Input:
        * n: int (n > 0)

    Returns:
        * s: int

    Raises:
        * ValueError: If n <= 0.

    Examples:
        >>> moebius_sum(10)
        -1
        >>> sum(moebius_list(10))
        -1
        >>> moebius_sum(0)
        Traceback (most recent call last):
        ...
        ValueError: moebius_sum: Must have n > 0.

    Details:
        The formula used here is based on the relationship
        \sum_{d = 1}^{n} M(n/d) = 1, where M(n) = \sum_{k = 1}^{n} moebius(k).
        We convert this to a recursive formula, whereby using some form of
        Dirichlet's hyperbola method and memoizing the recursion gives us a
        sublinear runtime.
    """
    if n <= 0:
        raise ValueError("moebius_sum: Must have n > 0.")

    sqrt = int(n**(0.5))
    mu_list = moebius_list(sqrt)

    # for i <= sqrt(n), compute the sum directly
    cache = {1: 1}
    for i in xrange(2, sqrt + 1):
        cache[i] = cache[i - 1] + mu_list[i]

    def M(n):
        if n in cache:
            return cache[n]
        sqrt_n = int(n**(0.5)) + 1
        value = n - 1
        for d in xrange(2, sqrt_n):
            value += M(n//d) - M(d - 1)
            value += mu_list[d]*(n//d - d)
        cache[n] = 1 - value
        return cache[n]

    return M(n)


mertens = moebius_sum


def sigma_sum(n):
    """Returns the value of the sum of sigma(k) for k = 1, 2, ..., n.

    Input:
        * n: int (n > 0)

    Returns:
        * sum: int

    Raises:
        * ValueError: If n <= 0.

    Examples:
        >>> sigma_sum(100)
        8299
        >>> sum(sigma_list(100))
        8299
        >>> sigma_sum(-1)
        Traceback (most recent call last):
        ...
        ValueError: sigma_sum: Must have n > 0.

    Details:
        The formula used is obtained by interchanging the order of summation and
        using some form of Dirichlet's hyperbola method to achieve a sublinear
        runtime.
    """
    if n <= 0:
        raise ValueError("sigma_sum: Must have n > 0.")

    sqrt = int(n**(0.5))
    value = -sqrt*(sqrt + 1)

    for k in xrange(1, sqrt + 1):
        nk = n//k
        tt = nk - k + 1
        value += 2*k*tt
        value += tt*(nk + k)

    return value//2


################################################################################
# Miscellaneous
################################################################################


def euler_phi_inverse(n):
    """Returns a sorted list of all positive integers k such that euler_phi(k) = n.

    Input:
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
        n_factorization = rosemary.number_theory.factorization.factor(n)
    elif isinstance(n, list) and n[0][0] > 0:
        n_factorization = n
        n = rosemary.number_theory.factorization.factor_back(n_factorization)
    else:
        raise ValueError("euler_phi_inverse: Input must be a positive integer or a factorization.")

    powers_of_two = set([1])
    if n % 2 == 0:
        for i in xrange(1, n_factorization[0][1] + 1):
            powers_of_two.add(2**i)

    # Odd primes p that divide n must have the property that p - 1 divides m.
    prime_list = []
    for d in rosemary.number_theory.factorization.divisors(n_factorization)[1:]:
        if rosemary.number_theory.primality.is_prime(d + 1):
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
