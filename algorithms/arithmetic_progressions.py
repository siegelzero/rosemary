def has_three_term_ap(A, **kwargs):
    """
    Determines if A contains a three term arithmetic progression.

    Given an array of positive integers A in increasing order, this determines
    if A has a three term arithmetic progression. If one is found, this returns
    either the elements themselves (the default) or the indices of the elements
    in A.

    Input:
        * A - A sorted list of positive integers.

    Output:
        * (i, j, k) - A tuple of nonnegative integers.

    Details:
        The algorithm used is based on an algorithm published by J. Erickson.

    Examples:
        >>> has_three_term_ap([2, 3, 5, 7, 11, 13, 17, 19])
        (1, 2, 3)
        >>> has_three_term_ap([2, 3, 5, 7, 11, 13, 17, 19], indices=True)
        >>> has_three_term_ap([2, 3, 7, 13])
        False
    """
    return_indices = kwargs.get('indices', False)

    n = len(A)
    L = []
    for j in range(1, n):
        i = j - 1
        k = j + 1
        while i >= 0 and k < n:
            if A[i] + A[k] < 2 * A[j]:
                k += 1
            elif A[i] + A[k] > 2 * A[j]:
                i -= 1
            else:
                if return_indices:
                    return (i, j, k)
                else:
                    return (A[i], A[j], A[k])
    return False

def longest_ap(A):
    """
    Returns the longest arithmetic progressions in A.

    Given an array of positive integers A in increasing order, this returns a
    tuple of the elements of A that form the longest arithmetic progression. If
    there are multiple progressions of this same length, all are returned in a
    list.

    Input:
        * A - A sorted list of positive integers.

    Output:
        * L - A sorted list of tuples.

    Details:
        This algorithm is based on one published by J. Erickson.

    Examples:
        >>> longest_ap([2, 3, 5, 7, 11])
        [(3, 5, 7), (5, 7, 11)]
        >>> longest_ap([2, 3, 5, 7, 11, 17, 23])
        [(5, 11, 17, 23])
        >>> p_list = primes(3000)
        >>> longest_ap(p_list)
        [(199, 409, 619, 829, 1039, 1249, 1459, 1669, 1879, 2089)]
    """
    L = {}
    n = len(A)
    longest = 2

    for j in xrange(n - 2, 0, -1):
        i = j - 1
        k = j + 1
        while i >= 0 and k <= n - 1:
            if A[i] + A[k] < 2 * A[j]:
                k += 1
            elif A[i] + A[k] > 2 * A[j]:
                L[i, j] = 2
                i -= 1
            else:
                L[i, j] = L.get((j, k), 2) + 1
                longest = max(longest, L[i, j])
                i -= 1
                k += 1
        while i >= 0:
            L[i, j] = 2
            i -= 1

    progs = []
    for (i, j) in L:
        if L[i, j] == longest:
            diff = A[j] - A[i]
            progs.append(tuple(A[i] + diff*k for k in xrange(longest)))

    progs.sort()
    return progs
