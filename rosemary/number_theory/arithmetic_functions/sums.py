################################################################################
# Summatory functions of multiplicative functions
################################################################################

import rosemary.number_theory.arithmetic_functions.sieves as sieves


def euler_phi_sum(n):
    r"""Returns the value of the sum of euler_phi(k) for k = 1, 2, ..., n.

    Parameters
    ----------
    n : int (n > 0)

    Returns
    -------
    s : int

    Raises
    ------
    ValueError
        If n <= 0.

    Notes
    -----
    Let S(n) = \sum_{k = 1}^{n} \phi(k). We use the identity S(n) = n*(n +
    1)/2 - \sum_{d = 2}^{n} S(n/d) to compute the value. By applying a
    variant of Dirichlet's Hyperbola Method, we are able to cut the limits
    of summation, and by memoizing the recursion, we are able to achieve a
    sublinear runtime. See Section 4.2 of "The Prime Numbers and Their
    Distribution" by Mendes, et al for the Hyperbola Method. See section
    4.9 equation 4.60 in "Concrete Mathematics" by Graham, et al for a
    proof of the identity for S(n).

    Examples
    --------
    >>> euler_phi_sum(10)
    32
    >>> sum(lists.euler_phi_list(10))
    32
    >>> euler_phi_sum(10**5)
    3039650754
    >>> euler_phi_sum(0)
    Traceback (most recent call last):
    ...
    ValueError: euler_phi_sum: Must have n > 0.
    """
    if n <= 0:
        raise ValueError("euler_phi_sum: Must have n > 0.")

    sqrt = int(n**(0.5))
    phi_list = sieves.euler_phi_list(sqrt)

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

    Parameters
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
    totients = sieves.euler_phi_list(sqrt)

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
    r"""Returns the value of the sum of moebius(k) for k = 1, 2, ..., n.

    Parameters
    ----------
    n : int (n > 0)

    Returns
    -------
    s : int

    Raises
    ------
    ValueError
        If n <= 0.

    Notes
    -----
    The formula used here is based on the relationship \sum_{d = 1}^{n}
    M(n/d) = 1, where M(n) = \sum_{k = 1}^{n} moebius(k).  We convert this
    to a recursive formula, whereby using some form of Dirichlet's
    hyperbola method and memoizing the recursion gives us a sublinear
    runtime.

    Examples
    --------
    >>> moebius_sum(10)
    -1
    >>> sum(lists.moebius_list(10))
    -1
    >>> moebius_sum(0)
    Traceback (most recent call last):
    ...
    ValueError: moebius_sum: Must have n > 0.
    """
    if n <= 0:
        raise ValueError("moebius_sum: Must have n > 0.")

    sqrt = int(n**(0.5))
    mu_list = sieves.moebius_list(sqrt)

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

    Parameters
        * n: int (n > 0)

    Returns:
        * sum: int

    Raises:
        * ValueError: If n <= 0.

    Examples:
        >>> sigma_sum(100)
        8299
        >>> sum(lists.sigma_list(100))
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
