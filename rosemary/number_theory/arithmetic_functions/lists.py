import rosemary.number_theory.primes.sieves as prime_sieves

################################################################################
# Lists of values of multiplicative functions
################################################################################


def euler_phi_list(n):
    """Returns a list of values of euler_phi(k) for 1 <= k <= n.

    Parameters
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
    prime_list = prime_sieves.primes(n)

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

    Parameters
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
    prime_list = prime_sieves.primes(sqrt)

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

    Parameters
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
        prime_list = prime_sieves.primes(n)

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

    Parameters
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
