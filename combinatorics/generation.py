def combinations(n, t):
    """
    generates all combinations of 0, 1, ..., n - 1 of size t.
    """
    c = [0] + [j - 1 for j in range(1, t + 1)] + [n, 0]

    while True:
        yield c[t:0:-1]
        j = 1
        while c[j] + 1 == c[j + 1]:
            c[j] = j - 1
            j += 1
        if j > t:
            return
        c[j] += 1


def permutations(A):
    """
    permutations(A):
    Given a sequence A of n elements a_1, a_2, ..., a_n, initially sorted so
    a_1 <= a_2 <= ... <= a_n, this algorithm generates all permutations of
    {a_1, a_2, ..., a_n} in lexicographic order.

    This algorithm comes from knuth 4 (Algorithm L)

    Examples:
    >>> L = [1, 2, 2, 3]
    >>> P = list(permutations(L))
    >>> P
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
    B = sorted(A)
    n = len(B) - 1

    while True:
        yield B[:]

        # Find the largest index j such that A[j] < A[j + 1]
        j = n - 1
        while j >= 0 and B[j] >= B[j + 1]:
            j -= 1

        if j < 0:
            return

        # Find the largest index l such that A[j] < A[l]. Then j < l.
        l = n
        while B[j] >= B[l]:
            l -= 1

        # Swap A[j] with A[l]
        B[j], B[l] = B[l], B[j]

        # Reverse the sequence from A[j + 1] to the end A[n]
        l = n
        k = j + 1
        while k < l:
            B[k], B[l] = B[l], B[k]
            k += 1
            l -= 1


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


def integer_partitions(n, **kwargs):
    """
    Returns an generator over the integer partitions of n, subject to the
    constraints given by the keywords.

    Given a nonnegative integer n and optional constraints, this returns a
    generator over the partitions of n satisfying these constraints. Each
    partition n = a_1 + a_2 + ... + a_m satisfyies a_1 >= a_2 >= ... a_m.

    Input:
        * n: int (n >= 0)

        * constraints:
            * distinct: bool (default=False)
                If True, only the partitions made up of distinct parts are
                yielded. Otherwise, parts may be used multiple times in each
                partition.

            * num_parts: int (default=0)
                If positive, only the partitions made up of exactly `num_parts`
                parts are yielded. If this value is zero, then no restrictions
                are placed on the number of parts.

            * parts: list (default=range(1, n + 1))
                Only the partitions that are comprised of parts from this list
                are yielded.

    Output:
        * partitions (generator)
            This is an generator of the partitions of n that satisfy the given
            constraints (if any).

    Examples:
        >>> list(integer_partitions(4))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
        >>> list(integer_partitions(4, distinct=True))
        [[4], [3, 1]]
        >>> list(integer_partitions(4, parts=[1, 3]))
        [[3, 1], [1, 1, 1, 1]]
        >>> list(integer_partitions(4, num_parts=2))
        [[3, 1], [2, 2]]
        >>> list(integer_partitions(4, num_parts=2, distinct=True))
        [[3, 1]]

    Details:
        The Art of Computer Programming Volume 4, section 7.2.1.4 is a great
        reference for the theory and algorithms here. Specifically, see Algorithm
        P for generating partitions with no constraints and Algorithm H for
        partitions into m parts. Another good reference here is Combinatorial
        Algorithms by Kreher and Stinson, specifically chapter 3.

        The function here uses Knuth's algorithms for the cases listed above,
        and uses a backtracking algorithm to generate partitions satisfying
        multiple constraints.
    """
    if 'parts' in kwargs:
        parts = sorted(kwargs['parts'], reverse=True)
        custom_parts = True
    else:
        parts = range(n, 0, -1)
        custom_parts = False
    total_number = len(parts)

    if 'distinct' in kwargs and kwargs['distinct']:
        distinct = 1
    else:
        distinct = 0

    if 'num_parts' in kwargs:
        num_parts = kwargs['num_parts']
        if num_parts > n:
            yield []
            return
    else:
        num_parts = 0

    def algorithm_p(n):
        """
        Generates all partitions of n. This is Algorithm P from 7.2.1.4 of
        Knuth, Vol. 4.
        """
        partition = [0]*n
        last_replaced = 0
        partition[last_replaced] = n
        idx = last_replaced - (n == 1)

        while True:
            yield partition[0:last_replaced + 1]
            if idx < 0:
                return
            if partition[idx] == 2:
                partition[idx] = 1
                idx -= 1
                last_replaced += 1
                partition[last_replaced] = 1
            else:
                replacement = partition[idx] - 1
                partition[idx] = replacement
                n = last_replaced - idx + 1
                last_replaced = idx + 1
                while n > replacement:
                    partition[last_replaced] = replacement
                    last_replaced += 1
                    n -= replacement
                partition[last_replaced] = n
                idx = last_replaced - (n == 1)

    def algorithm_h(n, m):
        """
        Generates all partitions of n into m parts. This is Algorithm H from
        7.2.1.4 of Knuth, Vol. 4.
        """
        partition = [1]*m
        partition[0] = n - m + 1

        while True:
            yield partition[:]
            if partition[1] < partition[0] - 1:
                partition[0] -= 1
                partition[1] += 1
            else:
                j = 2
                s = partition[0] + partition[1] - 1
                while j < m and partition[j] >= partition[0] - 1:
                    s += partition[j]
                    j += 1
                if j >= m:
                    return
                replacement = partition[j] + 1
                partition[j] = replacement
                j -= 1
                while j > 0:
                    partition[j] = replacement
                    s -= replacement
                    j -= 1
                partition[0] = s

    def backtrack(partial_sum, used, num_used, last_idx):
        if partial_sum == n:
            if not num_parts or (num_parts and num_used == num_parts):
                yield used
        elif partial_sum < n:
            if num_parts and num_used >= num_parts:
                return
            idx = 0
            if last_idx != 0:
                idx = last_idx + distinct
            for i in xrange(idx, total_number):
                part = parts[i]
                for partition in backtrack(partial_sum + part,
                                           used + [part], num_used + 1, i):
                    yield partition

    if distinct or custom_parts:
        partition_gen = backtrack(0, [], 0, 0)
    elif not distinct and not custom_parts and num_parts != 0:
        partition_gen = algorithm_h(n, num_parts)
    else:
        partition_gen = algorithm_p(n)

    for partition in partition_gen:
        yield partition
