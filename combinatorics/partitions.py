################################################################################
# Functions related to counting integer partitions
################################################################################

from rosemary.algebra import power_series
from rosemary.utilities import cached_function
import math

def bipartitions(m, n):
    """
    bipartitions(b, w):
    This returns the number of partitions of a multiset containing two distinct
    elements - m of one kind, and n of the other.

    Examples:

    >>> bipartitions(3, 1)
    7
    >>> bipartitions(60, 40)
    83735848679360680
    >>> bipartitions(4, 0)  # this is the usual integer partition function
    5
    """
    p = [[1] * (n + 1) for i in xrange(m + 1)]

    for i in xrange(0, m + 1):
        for j in xrange(0, n + 1):
            if i + j > 1:
                for k in xrange(i, m + 1):
                    for l in xrange(j, n + 1):
                       p[k][l] += p[k - i][l - j]

    return p[m][n]

################################################################################

@cached_function
def partition_number(n):
    """
    partition_number(n):
    Given a non-negative integer n, this returns the number of integer partitions
    of n; i.e. the number of ways to write n as a non-increasing sum of positive
    integers

    Examples:
    >>> partition_number(4)
    5
    >>> len(integer_partitions(5))
    5
    """
    if n in (0, 1):
        return 1
    if n < 0:
        return 0

    max_k = int((1 + math.sqrt(1 + 24*n)) / 6) + 1
    s = (-1)**(max_k + 1)
    ss = 0
    for k in xrange(max_k, 0, -1):
        t1 = partition_number(n - k * (3*k - 1) // 2)
        t2 = partition_number(n - k * (3*k + 1) // 2)
        ss += s * (t1 + t2)
        s *= (-1)

    return ss

################################################################################

def partition_list(n):
    """
    partition_list(n):
    This returns a list of the values p(0), p(1), ..., p(n)
    """
    p_vals = [0] * (n + 1)
    p_vals[0] = 1
    p_vals[1] = 1

    for m in xrange(2, n + 1):
        max_k = int((1 + math.sqrt(1 + 24*m)) / 6)

        # this is the sign in the alternating sum
        s = (-1)**(max_k + 1)

        # this comes from the pentagonal number theorem
        ss = 0
        for k in xrange(max_k, 0, -1):
            v1 = m - k * (3*k - 1) // 2
            t1 = p_vals[v1]

            v2 = m - k * (3*k + 1) // 2
            if v2 >= 0:
                t2 = p_vals[v2]
            else:
                t2 = 0

            ss += s * (t1 + t2)
            s *= (-1)

        p_vals[m] = ss

    return p_vals

################################################################################

def rth_power_partition_number(r, n):
    """
    rth_power_partition_number(r, n):
    This returns p_r(n); the nth coefficient of the rth power of the generating
    function for p(n).
    """
    if n == 0:
        return 1

    if r == 1:
        return partition_number(n)

    Pr = [0] * (n + 1)
    Pr[0] = 1

    # it's faster to compute these in descending order
    P = [ partition_number(k) for k in xrange(0, n + 1) ]

    # this bit uses the j.c.p. pure power recurrence to compute the rth
    # power coefficients of the partition generating function
    for k in xrange(1, n + 1):
        ss = 0
        for i in xrange(1, k + 1):
            ss += P[i] * ((r + 1) * i - k) * Pr[k - i]
        Pr[k] = ss // k

    return Pr[n]

################################################################################

def rth_power_partition_list(r, n):
    """
    rth_power_partition_number(r, n):
    This returns p_r(n); the nth coefficient of the rth power of the generating
    function for p(n).
    """
    if n == 0:
        return 1

    if r == 1:
        return partition_number(n)

    P = partition_list(n)
    Pr = power_series.power_series_power(P, r)

    return Pr

################################################################################
# Functions associated to constructing integer partitions
################################################################################

def integer_partitions(n):
    """
    integer_partitions(n):
    (Partitions in reverse lexicographic order).
    This algorithm returns all partitions a_1 >= a_2 >= ... >= a_m >= 1 with
    a_1 + a_2 + ... + a_m = n and 1 <= m <= n, assuming that n >= 1.

    Examples:
    >>> list(integer_partitions(4))
    [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    >>> L = integer_partitions(4)
    >>> for e in L: print e
      [4]
      [3, 1]
      [2, 2]
      [2, 1, 1]
      [1, 1, 1, 1]
    """
    a = [0] * (n + 1)
    m = 1
    a[1] = n
    q = m - (n == 1)
    L = []

    while True:
        L.append(a[1:m + 1])

        if a[q] != 2:
            if q == 0:
                return L
            x = a[q] - 1
            a[q] = x
            n = m - q + 1
            m = q + 1
        else:
            a[q] = 1
            q -= 1
            m += 1
            a[m] = 1
            continue

        while n > x:
            a[m] = x
            m += 1
            n -= x
        else:
            a[m] = n
            q = m - (n == 1)

################################################################################

def integer_partitions_gen(n):
    """
    integer_partitions_gen(n):
    (Partitions in reverse lexicographic order).
    This algorithm returns a generator over all partitions
    a_1 >= a_2 >= ... >= a_m >= 1 with a_1 + a_2 + ... + a_m = n
    and 1 <= m <= n, assuming that n >= 1.

    Examples:
    >>> list(integer_partitions(4))
    [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    >>> L = integer_partitions(4)
    >>> for e in L: print e
      [4]
      [3, 1]
      [2, 2]
      [2, 1, 1]
      [1, 1, 1, 1]
    """
    a = [0] * (n + 1)
    m = 1
    a[1] = n
    q = m - (n == 1)

    while True:
        yield a[1:m + 1]

        if a[q] != 2:
            if q == 0:
                return
            x = a[q] - 1
            a[q] = x
            n = m - q + 1
            m = q + 1
        else:
            a[q] = 1
            q -= 1
            m += 1
            a[m] = 1
            continue
        while n > x:
            a[m] = x
            m += 1
            n -= x
        else:
            a[m] = n
            q = m - (n == 1)

################################################################################

def integer_partitions_m_parts(n, m):
    """
    integer_partitions_m_parts(n, m):
    This algorithm generates all integer m-tuples a_1, ..., a_m such that
    a_1 >= ... a_m >= 1 and a_1 + ... + a_m = n, assuming that n >= m >= 2.
    """
    a = [ n - m + 1 ] + [1] * (m - 1)
    L = []

    while True:
        L.append(a[:])

        if a[1] < a[0] - 1:
            a[0] -= 1
            a[1] += 1
            continue

        j = 2
        s = a[0] + a[1] - 1
        while j < m and a[j] >= a[0] - 1:
            s += a[j]
            j += 1

        if j >= m:
            return L

        x = a[j] + 1
        a[j] = x
        j -= 1

        while j > 0:
            a[j] = x
            s -= x
            j -= 1

        a[0] = s

################################################################################

def integer_partitions_m_parts_gen(n, m):
    """
    integer_partitions_m_parts_gen(n, m):
    This algorithm returns a generator over all integer m-tuples a_1, ..., a_m
    such that a_1 >= ... a_m >= 1 and a_1 + ... + a_m = n, assuming that
    n >= m >= 2.
    """
    a = [ n - m + 1 ] + [1] * (m - 1)

    while True:
        yield a[:]

        if a[1] < a[0] - 1:
            a[0] -= 1
            a[1] += 1
            continue

        j = 2
        s = a[0] + a[1] - 1
        while j < m and a[j] >= a[0] - 1:
            s += a[j]
            j += 1

        if j >= m:
            return

        x = a[j] + 1
        a[j] = x
        j -= 1

        while j > 0:
            a[j] = x
            s -= x
            j -= 1

        a[0] = s

################################################################################

def partitions_nonrecursive(n):
    # L is our list of possible parts
    L = range(n, 0, -1)
    # this will store all partitions of n
    good = []
    so_far = []
    # S is our stack, S[1] is the first level
    S = {}
    S[1] = L[:]
    k = 1

    while k > 0:
        while S[k]:
            # while the stack is nonempty, pick the next element
            x = S[k].pop()
            # append this element onto the current partial solution
            so_far.append(x)
            # check that the new config is still a partial solution
            t_sum = sum(so_far)
            if t_sum == n:
                good.append(so_far[:])

            k += 1
            S[k] = []
            if t_sum < n:
                for t in L:
                    if t < x:
                        break
                    if t_sum + t <= n:
                        S[k].append(t)

        if so_far:
            so_far.pop()
        k -= 1

    return good

################################################################################

def restricted_partitions(n, parts = None):
    """
    restricted_partitions(n, parts):
    Given a positive integer n, this function returns all partitions of n using
    parts from the list parts.

    Examples:
    >>> restricted_partitions(4)
    [[1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1], [4]]
    >>> restricted_partitions(10, [2, 3, 5, 7])
    [[2, 2, 2, 2, 2], [3, 3, 2, 2], [5, 3, 2], [5, 5], [7, 3]]
    """
    def can_add_one(n, parts, used, partial_sum, last):
        for part in parts:
            if part > last:
                break
            t_sum = partial_sum + part
            if t_sum <= n:
                t_used = used + [ part ]
                if t_sum == n:
                    L.append(t_used)
                else:
                    can_add_one(n, parts, t_used, t_sum, part)

    if parts is None:
        parts = range(1, n + 1)
    else:
        parts.sort()

    L = []
    can_add_one(n, parts, [], 0, parts[-1])
    return L

################################################################################

def restricted_partitions_gen(n, parts = None):
    """
    restricted_partitions_gen(n, parts):
    Given a positive integer n, this function returns a generator over all
    partitions of n using parts from the list parts.

    Examples:
    >>> L = restricted_partitions(4)
    >>> list(L)
    [[1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1], [4]]
    >>> L = restricted_partitions(10, [2, 3, 5, 7])
    >>> list(L)
    [[2, 2, 2, 2, 2], [3, 3, 2, 2], [5, 3, 2], [5, 5], [7, 3]]
    """
    if parts is None:
        parts = range(1, n + 1)
    else:
        parts.sort()

    so_far = []
    S = {}
    S[1] = parts[:]
    k = 1

    while k > 0:
        while S[k]:
            # pick the next element in the stack
            x = S[k].pop()

            # append this element onto the current partial solution
            so_far.append(x)

            # check that the new config is still a partial solution
            t_sum = sum(so_far)
            if t_sum == n:
                yield so_far[:]
                #good.append(so_far[:])

            k += 1
            S[k] = []
            for t in parts:
                if t > x:
                    break
                if t_sum + t <= n:
                    S[k].append(t)

        if so_far:
            so_far.pop()
        k -= 1

################################################################################

def restricted_partitions_distinct(n, parts = None):
    """
    restricted_partitions(n, parts):
    Given a positive integer n, this function returns all partitions of n using
    parts from the list parts.

    Examples:
    >>> restricted_partitions(4)
    [[1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1], [4]]
    >>> restricted_partitions(10, [2, 3, 5, 7])
    [[2, 2, 2, 2, 2], [3, 3, 2, 2], [5, 3, 2], [5, 5], [7, 3]]
    """
    def can_add_one(n, parts, used, partial_sum, last):
        for part in parts:
            if part >= last:
                break
            t_sum = partial_sum + part
            if t_sum <= n:
                t_used = used + [ part ]
                if t_sum == n:
                    L.append(t_used)
                else:
                    can_add_one(n, parts, t_used, t_sum, part)

    if parts is None:
        parts = range(1, n + 1)
    else:
        if sum(parts) < n:
            return []
        parts.sort()

    L = []
    can_add_one(n, parts, [], 0, parts[-1])
    return L


