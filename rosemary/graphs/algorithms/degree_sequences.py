def is_graphic(S):
    r"""
    Returns True if S is a graphic sequence and False otherwise.

    A sequence of non-negative integers is called graphic if there exists a
    graph whose degree sequence is precisely that sequence.

    Parameters
    ----------
    S : iterable
        A iterable of the degrees in the sequence.

    Returns
    -------
    bool : boolean
        True if the sequence is graphic and False otherwise.

    Notes
    -----
    This method uses the theorem of Havel and Hakimi to determine if the
    sequence is graphic. The theorem states that a sequence `d` is graphic
    if and only if `d'` is graphic, where `d'` is obtained from `d` by
    deleting its largest element `r`, and subtracting 1 from the next `r`
    largest elements.

    See Theorem 1.1.2 of [1] for proof of the theorem. See Algorithm 2.5.7
    of [2] for a description and details of the implementation. See also
    Theorem 1.3.31 of [3]

    References
    ----------
    .. [1] N. Hartsfield, G. Ringel, "Pearls in Graph Theory - A
    Comprehensive Introduction", Revised and Augmented Edition, Academic
    Press, Inc, 1994.

    .. [2] K.M. Koh, et al, "Graph Theory - Undergraduate Mathematics",
    World Scientific Publishing, 2015.

    .. [3] D.B. West, "Introduction to Graph Theory", 2nd Edition, Prentice
    Hall, 2000.

    Examples
    --------
    >>> is_graphic([5, 4, 4, 3, 1, 1])
    False
    >>> is_graphic([4, 4, 3, 2, 2, 1])
    True
    """
    # Can't have any degrees >= than the number of vertices.
    n = len(S)
    for d in S:
        if d >= n:
            return False

    T = list(S)
    while True:
        # Detect if the sequence is all-zero, or return False if there are
        # any negative entries.
        all_zero = True
        for e in T:
            if e < 0:
                return False
            elif e > 0:
                all_zero = False

        if all_zero:
            return True

        T = sorted(T, reverse=True)
        print T
        d = T.pop(0)

        # Subtract 1 from the next d largest elements.
        for i in xrange(d):
            T[i] -= 1


def is_graphic_erdos_gallai(S):
    T = sorted(S, reverse=True)
    d_sum = 0

    for k in xrange(len(T)):
        d_sum += T[k]
        rhs = k*(k + 1) + sum(min(k + 1, T[i]) for i in xrange(k + 1, len(T)))

        if d_sum > rhs:
            return False

    return True
