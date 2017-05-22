###########################################################################
# Algorithms for Generating Permutations
###########################################################################


def phillips(elements):
    """Yields all arrangements of the given multiset.

    Given a sequence of elements, this generates all permutations of the
    elements in lexicographic order. The elements do not have to be
    distinct.

    Parameters
    ----------
    elements : list
        List of elements

    Yields
    ------
    arrangement : list
        Arrangement of the sequence.

    See Also
    --------
    knuth_algorithm_l

    Notes
    -----
    This is an implementation of the sped-up version of Knuth's Algorithm
    L outlined in Exercise 1 of 7.2.1.2 [1]. The algorithm is based on the
    description given in [2]

    The idea is that in step L2 of Algorithm L, we have `j = n - 1` half of
    the time when the elements are distinct, because exactly one half of
    the permutations have `a[n - 1] < a[n]`. This special case is
    recognized, as well as a few other cases.

    References
    ----------
    .. [1] D.E. Knuth, "The Art of Computer Programming, Volume 4, Fascicle
    2: Generating all Tuples and Permutations", Addison-Wesley
    Professional, 2005.

    .. [2] J.P.N. Phillips, "Permutations of the Elements of a Vector in
    Lexicographic Order (Algorithm 28)", The Computer Journal, Volume 10,
    Number 3, 1967.

    Examples
    --------
    >>> L = [1, 2, 2, 3]
    >>> list(permutations(L))
    [[1, 2, 2, 3],
     [1, 2, 3, 2],
     [1, 3, 2, 2],
     [2, 1, 2, 3],
     [2, 1, 3, 2],
     [2, 2, 1, 3],
     [2, 2, 3, 1],
     [2, 3, 1, 2],
     [2, 3, 2, 1],
     [3, 1, 2, 2],
     [3, 2, 1, 2],
     [3, 2, 2, 1]]
    """
    # We need the elements to be in sorted order. We don't want to modify
    # the input elements, so we make a copy here.
    entries = sorted(elements)
    n = len(entries) - 1

    # The Phillips algorithm assumes that there at least three elements in
    # the list. We handle the other cases separately.
    if n == 0:
        yield entries[:]
        return

    if n == 1:
        yield [entries[0], entries[1]]
        yield [entries[1], entries[0]]
        return

    while True:
        # L1 [Visit]
        # Yield a copy of the arrangement
        yield entries[:]

        # L2' [Easiest case?]
        y = entries[n - 1]
        z = entries[n]
        if y < z:
            entries[n - 1] = z
            entries[n] = y
            continue

        # L2.1' [Next easiest case?]
        x = entries[n - 2]
        if x < y:
            if x < z:
                entries[n - 2] = z
                entries[n - 1] = x
                entries[n] = y
            else:
                entries[n - 2] = y
                entries[n - 1] = z
                entries[n] = x
            continue
        else:
            # L2.2' [Find j]
            j = n - 3
            y = entries[j]
            while y >= x:
                j -= 1
                x = y
                y = entries[j]

            if j < 0:
                return

            while True:
                # L3' [Easy increase?]
                if y < z:
                    entries[j] = z
                    entries[j + 1] = y
                    entries[n] = x
                    break

                # L3.1' [Increase a[j]]
                l = n - 1
                while y >= entries[l]:
                    l -= 1

                entries[j] = entries[l]
                entries[l] = y

                # L4' [Begin to reverse]
                entries[n] = entries[j + 1]
                entries[j + 1] = z
                break

            # L4.1' [Reverse a[j + 1], ..., a[n - 1]]
            k = j + 2
            l = n - 1
            while k < l:
                entries[k], entries[l] = entries[l], entries[k]
                k += 1
                l -= 1


def knuth_algorithm_l(elements):
    """Yields all permutations of the given multiset.

    Given a sequence of elements, this generates all permutations of the
    elements in lexicographic order. The elements do not have to be
    distinct.

    Parameters
    ----------
    elements : list
        List of elements

    Yields
    ------
    permutation : list
        Permutation of the sequence.

    See Also
    --------
    phillips

    Notes
    -----
    This method uses Algorithm L from 7.2.1.2 in [1] to generate all
    permutations of the given sequence. This is a fairly standard algorithm
    to generate permutations of the given sequence of elements. Knuth's
    version in [1] generalizes the version of 5.1.1 in [2] to handle the
    case where the elements may not be distinct.

    References
    ----------
    .. [1] D.E. Knuth, "The Art of Computer Programming, Volume 4, Fascicle
    2: Generating all Tuples and Permutations", Addison-Wesley
    Professional, 2005.

    .. [2] E.M. Reingold, J. Nievergelt, N. Deo, "Combinatorial Algorithms:
    Theory and Practice", Prentice Hall, 1977

    Examples
    --------
    >>> L = [1, 2, 2, 3]
    >>> list(knuth_algorithm_l(L))
    [[1, 2, 2, 3],
     [1, 2, 3, 2],
     [1, 3, 2, 2],
     [2, 1, 2, 3],
     [2, 1, 3, 2],
     [2, 2, 1, 3],
     [2, 2, 3, 1],
     [2, 3, 1, 2],
     [2, 3, 2, 1],
     [3, 1, 2, 2],
     [3, 2, 1, 2],
     [3, 2, 2, 1]]
    """
    # We need the elements to be in sorted order. We don't want to modify
    # the input elements, so we make a copy here.
    entries = sorted(elements)
    n = len(entries) - 1

    while True:
        # L1 [Visit]
        # Yield a copy of the permutation
        yield entries[:]

        # L2 [Find j]
        # Knuth uses `a` for his permutation, so we'll follow.
        # Find the largest index j such that a[j] < a[j + 1]
        j = n - 1
        while j >= 0 and entries[j] >= entries[j + 1]:
            j -= 1

        if j < 0:
            # If we reach this point, then all permutations have been
            # visited.
            return

        # L3 [Increase a[j]]
        # Find the smallest element to the right of a[j] and greater than a[j]
        l = n
        while entries[j] >= entries[l]:
            l -= 1

        # Swap A[j] with A[l]
        entries[j], entries[l] = entries[l], entries[j]

        # L4 [Reverse a[j + 1], ..., a[n]]
        # Reverse the sequence from A[j + 1] to the end A[n]
        l = n
        k = j + 1
        while k < l:
            entries[k], entries[l] = entries[l], entries[k]
            k += 1
            l -= 1


###########################################################################
# Algorithms for Generating Variations
###########################################################################


def knuth_algorithm_l_variations(elements):
    """Yields all variations of the given multiset.

    Given a sequence of elements, this generates all variations of the
    elements in lexicographic order. The variations are the permutations of
    all its submultisets.

    Parameters
    ----------
    elements : list
        List of elements

    Yields
    ------
    variation : list
        Variation of the sequence.

    See Also
    --------
    knuth_algorithm_l_rvariations

    Notes
    -----
    This method is based on Algorithm L from 7.2.1.2 in [1] to generate all
    permutations of the given sequence. In Exercise 8 of the section, an
    algorithm is described to generate all variations of the given
    multiset.

    References
    ----------
    .. [1] D.E. Knuth, "The Art of Computer Programming, Volume 4, Fascicle
    2: Generating all Tuples and Permutations", Addison-Wesley
    Professional, 2005.

    Examples
    --------
    >>> L = [1, 2, 2, 3]
    >>> list(knuth_algorithm_l_variations(L))
    [[], [1], [1, 2], [1, 2, 2], [1, 2, 2, 3], [1, 2, 3], [1, 2, 3, 2],
     [1, 3], [1, 3, 2], [1, 3, 2, 2], [2], [2, 1], [2, 1, 2], [2, 1, 2, 3],
     [2, 1, 3], [2, 1, 3, 2], [2, 2], [2, 2, 1], [2, 2, 1, 3], [2, 2, 3],
     [2, 2, 3, 1], [2, 3], [2, 3, 1], [2, 3, 1, 2], [2, 3, 2], [2, 3, 2, 1],
     [3], [3, 1], [3, 1, 2], [3, 1, 2, 2], [3, 2], [3, 2, 1], [3, 2, 1, 2],
     [3, 2, 2], [3, 2, 2, 1]]
    """
    # We need the elements to be in sorted order. We don't want to modify
    # the input elements, so we make a copy here.
    entries = sorted(elements)
    n = len(entries) - 1
    j = 0
    yield []

    while True:
        # L1 [Visit]
        # Yield a copy of the variation
        while j <= n:
            yield entries[:j + 1]
            j += 1

        # L2 [Find j]
        # Knuth uses `a` for his permutation, so we'll follow.
        # Find the largest index j such that a[j] < a[j + 1]
        j = n - 1
        while j >= 0 and entries[j] >= entries[j + 1]:
            j -= 1

        if j < 0:
            # If we reach this point, then all permutations have been
            # visited.
            return

        # L3 [Increase a[j]]
        # Find the smallest element to the right of a[j] and greater than a[j]
        l = n
        while entries[j] >= entries[l]:
            l -= 1

        # Swap A[j] with A[l]
        entries[j], entries[l] = entries[l], entries[j]

        # L4 [Reverse a[j + 1], ..., a[n]]
        # Reverse the sequence from A[j + 1] to the end A[n]
        l = n
        k = j + 1
        while k < l:
            entries[k], entries[l] = entries[l], entries[k]
            k += 1
            l -= 1


def knuth_algorithm_l_rvariations(elements, r):
    """Yields all r-variations of the given multiset.

    Given a sequence of elements, this generates all r-variations of the
    elements in lexicographic order. The r-variations are the permutations of
    its r-element submultisets.

    Parameters
    ----------
    elements : list
        List of elements

    Yields
    ------
    variation : list
        Variation of the sequence.

    See Also
    --------
    knuth_algorithm_l_variations

    Notes
    -----
    This method is based on Algorithm L from 7.2.1.2 in [1] to generate all
    permutations of the given sequence. In Exercise 9 of the section, an
    algorithm is described to generate all r-variations of the given
    multiset.

    References
    ----------
    .. [1] D.E. Knuth, "The Art of Computer Programming, Volume 4, Fascicle
    2: Generating all Tuples and Permutations", Addison-Wesley
    Professional, 2005.

    Examples
    --------
    >>> L = [1, 2, 2, 3]
    >>> list(knuth_algorithm_l_rvariations(L, 2))
    [[1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2]]
    """
    entries = sorted(elements)
    n = len(entries) - 1
    r = r - 1

    while True:
        # R1 [Visit]
        # Yield a copy of the arrangement
        yield entries[:r + 1]

        # R2 [Find j]
        # Knuth uses `a` for his permutation, so we'll follow.
        # Find the largest index j such that a[j] < a[j + 1]
        if entries[r] < entries[n]:
            j = r + 1
            while entries[j] <= entries[r]:
                j += 1
            entries[r], entries[j] = entries[j], entries[r]
            continue

        # R3
        l = n
        k = r + 1
        while k < l:
            entries[k], entries[l] = entries[l], entries[k]
            k += 1
            l -= 1

        # R4
        j = r - 1
        while j >= 0 and entries[j] >= entries[j + 1]:
            j -= 1

        if j < 0:
            # If we reach this point, then all permutations have been
            # visited.
            return

        # R5 [Increase a[j]]
        # Find the smallest element to the right of a[j] and greater than a[j]
        l = n
        while entries[j] >= entries[l]:
            l -= 1

        # Swap A[j] with A[l]
        entries[j], entries[l] = entries[l], entries[j]

        # R6 [Reverse a[j + 1], ..., a[n]]
        # Reverse the sequence from A[j + 1] to the end A[n]
        l = n
        k = j + 1
        while k < l:
            entries[k], entries[l] = entries[l], entries[k]
            k += 1
            l -= 1
