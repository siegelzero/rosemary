def n_tuples(n, M=None):
    """
    n_tuples(n M):
    This algorithm generates all n-tuples a_0, a_1, ..., a_{n - 1} that satisfy
    0 <= a_j < M[j] for 0 <= j <= n - 1. If M is None, this returns all 0,1-tuples
    of length n.

    Examples:
    >>> L = n_tuples(2)
    >>> list(L)
    [[0, 0], [0, 1], [1, 0], [1, 1]]
    >>> L = n_tuples(3, [3, 2, 2])
    >>> list(L)
    [[0, 0, 0],
     [0, 0, 1],
     [0, 1, 0],
     [0, 1, 1],
     [1, 0, 0],
     [1, 0, 1],
     [1, 1, 0],
     [1, 1, 1],
     [2, 0, 0],
     [2, 0, 1],
     [2, 1, 0],
     [2, 1, 1]]
    >>> L = n_tuples(2, 3)
    >>> list(L)
    [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
    """
    A = [0 for j in range(n + 1)]

    if isinstance(M, list):
        if len(M) == n:
            M = [2] + M
    else:
        if isinstance(M, (int, long)):
            bd = M
        else:
            bd = 2
        M = [2] + [bd for j in range(n)]

    while True:
        yield A[1:]
        j = n
        while A[j] == M[j] - 1:
            A[j] = 0
            j -= 1
        if j == 0:
            return
        else:
            A[j] = A[j] + 1
