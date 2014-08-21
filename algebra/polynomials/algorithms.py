# Algorithms on polynomials


def mul_classical(M, N):
    """
    Returns the product of M and N.

    Given lists M and N of numbers, this returns the coefficients of the product
    f(x) * g(x), where f(x) = M[0] + M[1]*x + ... + M[k - 1]*x^{k - 1} and g(x)
    = N[0] + N[1]*x + ... + N[l - 1]*x^{l - 1}.

    Input:
        * M - A list or tuple of numbers.
        * N - A list or tuple of numbers.

    Output:
        * L - A list of numbers. These are the coefficients of the product of
        f(x) * g(x), where M and N are the lists of their coefficients,
        respectively.

    Details:
        This function uses the classical "gradeschool" algorithm for
        multiplication. This algorithm requires O(k*l) operations to compute the
        product of two polynomials with degrees k >= 1 and l >= 1.

    Examples:
        >>> M = [1, 1]
        >>> mul_classical(M, M)
        [1, 2, 1]
    """
    m = len(M)
    n = len(N)
    if m == 0 or n == 0:
        return []
    P = [0] * (m + n - 1)
    for i in xrange(m):
        for j in xrange(n):
            P[i + j] += M[i] * N[j]
    return P


def mul_karatsuba(M, N):
    """
    Returns the product of M and N.

    Given lists M and N of numbers, this returns the coefficients of the product
    f(x) * g(x), where f(x) = M[0] + M[1]*x + ... + M[n]*x^n and g(x) = N[0] +
    N[1]*x + ... + N[k]*x^k.

    Input:
        * M - A list or tuple of numbers.
        * N - A list or tuple of numbers.

    Output:
        * L - A list of numbers. These are the coefficients of the product of
        f(x) * g(x), where M and N are the lists of their coefficients,
        respectively.

    Details:
        This function uses the Karatsuba algorithm for multiplication. For
        multiplying two degree n polynomials, this method has time complexity
        O(n^c), where c = log_2(3) = 1.58496..., so it is much faster than the
        classical O(n^2) method. For more information, see section 4.3.3 of "The
        Art of Computer Programming, Volume II" by D. Knuth.

    Examples:
        >>> M = [1]*2000
        >>> N = [1]*3000
        >>> A = mul_classical(M, N)
        >>> B = mul_karatsuba(M, N)
        >>> A == B
        True
    """
    m = len(M)
    n = len(N)
    kthresh = 30

    # Classical Multiplication is faster until a certain point.
    if m < kthresh or n < kthresh:
        return mul_classical(M, N)

    # The algorithm works best when we split into halves
    k = max(m, n) // 2

    # write M = f1*x^k + f0, N = g1*x^k + g0
    f0 = M[0:k]
    f1 = M[k:]
    g0 = N[0:k]
    g1 = N[k:]

    # We want to compute MN = (f1*x^k + f0)(g1*x^k + g0).
    # We let z0 = f0*g0, z2 = f1*g1
    z0 = mul_karatsuba(f0, g0)
    z2 = mul_karatsuba(f1, g1)
    f1_f0 = list_add(f1, f0)
    g1_g0 = list_add(g1, g0)
    xxyy = mul_karatsuba(f1_f0, g1_g0)
    z1 = list_subtract(list_subtract(xxyy, z2), z0)
    t2 = [0] * (2*k) + z2
    t1 = [0] * (k) + z1
    xy = list_add(list_add(t2, t1), z0)
    return xy


def power_jcp(A, k):
    """
    Returns the kth power

    Given a nonnegative integer k and a list A, this returns the coefficients of
    f(x)^k, where f(x) is the polynomial f(x) = A[0] + A[1]*x + ... + A[n]*x^n.

    Input:
        * A - A list or tuple of numbers
        * k - a nonnegative integer

    Output:
        * P - A list of numbers. These are the coefficients of f(x)^k, where
        f(x) is the polynomial with coefficients given by the list A.

    Details:
        The algorithm used is the J.C.P. Miller Pure-Power Recurrence. For a
        fixed polynomial of degree L, this algorithm has time complexity
        O(L^2*m) = O(m) and space complexity O(L) = O(1). For more details, look
        at the paper titled "The J.C.P. Miller Recurrence for Exponentiating a
        Polynomial, and its q-Analog" by D. Zeilberger.

    Examples:
        >>> L = [1, 3, 2]
        >>> power_jcp(L, 2)
        [1, 6, 13, 12, 4]
    """
    shift = 0
    while A[shift] == 0:
        shift += 1

    B = A[shift:]
    n = len(B)
    new_deg = (n - 1) * k + 1
    P = [0] * new_deg
    a0 = B[0]
    P[0] = a0**k
    for i in xrange(1, new_deg):
        ss = 0
        for j in xrange(1, min(n, i + 1)):
            ss += B[j] * ((k + 1) * j - i) * P[i - j]
        P[i] = ss // (i * a0)
    return [0]*shift*k + P


def list_add(A, B):
    """
    Adds two lists.

    Let A be a list of length m, and let B be a list of length n. This returns
    the list [A[i] + B[i] | 0 <= i <= max(m, n)], where A[i], B[i] = 0 if i >=
    m, n, respectively. In effect, this works the same as polynomial addition,
    where we assume that coefficients not appearing in a polynomial are zero.

    Input:
        * A - A list of numbers.
        * B - A list of numbers.

    Output:
        * C - A list of numbers.

    Examples:
        >>> list_add([0, 1, 2], [1, -1, 1])
        [1, 0, 3]
        >>> list_add([1, 1, 1, 1], [2, 3, 5])
        [3, 4, 6, 1]
    """
    if len(A) >= len(B):
        C = list(A)
        for (i, c) in enumerate(B):
            C[i] += c
    else:
        C = list(B)
        for (i, c) in enumerate(A):
            C[i] += c
    return C


def list_subtract(A, B):
    """
    Subtract two lists.

    Let A be a list of length m, and let B be a list of length n. This returns
    the list [A[i] - B[i] | 0 <= i <= max(m, n)], where A[i], B[i] = 0 if i >=
    m, n, respectively. In effect, this works the same as polynomial addition,
    where we assume that coefficients not appearing in a polynomial are zero.

    Input:
        * A (list)
        * B (list)

    Output:
        * C (list)

    Examples:
        >>> list_subtract([0, 1, 2], [1, -1, 1])
        [-1, 2, 1]
        >>> list_subtract([1, 1, 1, 1], [2, 3, 5])
        [-1, -2, -4, 1]
    """
    if len(A) >= len(B):
        C = list(A)
        for (i, c) in enumerate(B):
            C[i] -= c
    else:
        C = [-e for e in B]
        for (i, c) in enumerate(A):
            C[i] += c
    return C
