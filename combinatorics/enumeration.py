# Algorithms for enumerative combinatorics: methods to count the number of
# certain structures and configurations. Every generation algorithm can be made
# in to an enumeration algorithm, but we prefer algorithms more advanced than
# simply listing and counting whenever possible.

from rosemary.number_theory.core import gcd
from rosemary.number_theory.factorization import xdivisors
from rosemary.number_theory.arithmetic_functions import euler_phi, factorial

def bell_number(n):
    """
    Returns the Bell number B(n).

    The Bell number B(n) is defined to be the number of partitions of an n-set
    (or equivalently, the number of equivalence relations on an n-set).

    Input:
        * n: int (n >= 0)

    Output:
        * B: int

    Examples:
        >>> bell_number(4)
        15
        >>> bell_number(10)
        115975

    Details:
        The algorithm uses the recurrence $B(m) = \sum_{i = 0}^{m - 1} \binom{m
        - 1}{i} B(i)$. See Theorem 3.9 of Combinatorial Algorithms by Kreher and
          Stinson for a proof of the recurrence.
    """
    if n < 0:
        raise ValueError("bell_number: Must have n >= 0.")

    B = [0]*(n + 1)
    B[0] = 1

    for i in xrange(1, n + 1):
        coeff = 1
        value = 0
        for k in xrange(1, i + 1):
            value += coeff*B[i - k]
            coeff *= (i - k)
            coeff /= k
        B[i] = value
    return B[n]

def stirling_number2(m, n):
    S = {}
    S[0, 0] = 1

    for i in xrange(1, m + 1):
        S[i, 0] = 0

    for i in xrange(m + 1):
        S[i, i + 1] = 0

    for i in xrange(1, m + 1):
        for j in xrange(1, min(i, n) + 1):
            S[i, j] = j*S[i - 1, j] + S[i - j, j - 1]

    return S[m, n]


def binomial(n, k):
    """
    Returns the binomial coefficient.

    Input:
        n: a nonnegative integer
        k: a nonnegative integer

    Returns:
        The binomial coefficient n choose k; i.e. the number of ways to choose k
        elements from a set of n elements.

    Examples:
        >>> binomial(5, 2)
        10
        >>> binomial(2, 0)
        1
    """
    pp = 1
    k = min(k, n - k)
    for e in xrange(n - k + 1, n + 1):
        pp = pp * e
    for e in xrange(2, k + 1):
        pp = pp / e
    return pp

def multinomial(n, k):
    """
    Returns the multinomial coefficient.

    Input:
        n: A positive integer
        k: A list or tuple of positive integers
    Returns:
        This returns the multinomial coefficient (n, k1, k2, ..., km) defined by
        n! / (k1! * k2! * ... * km!), where k1 + k2 + ... + km == n. If sum(k)
        != n, then the other elements of k are all assumed to be 1.

    Examples:
        >>> multinomial(6, [3, 3])
        20
        >>> multinomial(6, [2, 2, 2])
        90
    """
    if sum(k) > n:
        raise ValueError('Too many ks')
    elif len(k) == 0:
        raise ValueError('Not enough ks')

    num = factorial(n)
    den = 1
    for e in k:
        den *= factorial(e)
    return num / den


def bipartitions(m, n):
    """
    Returns the number of partitions of a multiset containing two distinct
    elements - m of one kind, and n of the other.

    Input:
        * m: int (m >= 0)
        * n: int (n >= 0)

    Output:
        * p: int
            The number of bipartitions

    Examples:
        >>> bipartitions(3, 1)
        7
        >>> bipartitions(60, 40)
        83735848679360680
        >>> bipartitions(10, 0)
        42

    Details:
        See 7.2.1.5 of TAOCP Vol 4 for more details.
    """
    p = [[1]*(n + 1) for i in xrange(m + 1)]

    for i in xrange(m + 1):
        for j in xrange(n + 1):
            if i + j > 1:
                for k in xrange(i, m + 1):
                    for l in xrange(j, n + 1):
                       p[k][l] += p[k - i][l - j]

    return p[m][n]

def partition_list(n):
    """
    partition_list(n):
    This returns a list of the values p(0), p(1), ..., p(n)
    """
    values = [0]*(n + 1)
    values[0] = 1
    values[1] = 1

    for m in xrange(2, n + 1):
        max_idx = int((1 + (1 + 24*m)**(0.5))/6)
        sign = (-1)**(max_idx + 1)
        total = 0

        for k in xrange(max_idx, 0, -1):
            v1 = m - k*(3*k - 1)//2
            t1 = values[v1]

            v2 = m - k*(3*k + 1)//2
            t2 = 0 if v2 < 0 else values[v2]

            total += sign*(t1 + t2)
            sign *= (-1)
        values[m] = total
    return values


def partition_number(n):
    return partition_list(n)[n]

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

def number_of_domino_tilings(m, n):
    """
    This returns the number of ways to tile an m x n grid with dominos.
    We can think of the board as follows:

        0       1       2   ... n-1
        n       n+1     n+2 ... 2*n-1
        .
        .
        .
        (m-1)*n ...     ...     m*n-1

    For example, the 4x6 can be though of as

        0   1   2   3   4   5
        6   7   8   9   10  11
        12  13  14  15  16  17
        18  19  20  21  22  23

    Observe that the element in position 15 has coordinates (3, 2).
    Also, observe that 15 // 6 == 2, while 15 % 4 == 3.
    """
    # create the board, and store the largest number appearing
    largest = m*n - 1
    D = {}
    def backtrack(covered, uncovered):
        t_c = tuple(covered)
        t_u = tuple(uncovered)

        if (t_c, t_u) in D:
            return D[t_c, t_u]

        ss = 0

        # Find the smallest number not covered by a domino, and cover it
        smallest = min(uncovered)

        # there are two possible domino configurations to cover this spot
        # this is the vertical configuration
        piece = (smallest, smallest + n)

        # ensure we're not extending beyond the board
        if piece[1] <= largest:
            # don't cover a square that is already covered
            if not covered.intersection(piece):
                # remove the piece from the set of uncovered squares
                new_uncovered = uncovered.difference(piece)
                if not new_uncovered:
                    # if nothing is left uncovered, then we have a good tiling
                    return 1
                else:
                    # otherwise, we continue to cover the board
                    new_covered = covered.union(piece)
                    ss += backtrack(new_covered, new_uncovered)

        # this is the horizontal configuration
        piece = (smallest, smallest + 1)

        # ensure that we're not extending beyong the board
        if piece[1] <= largest and piece[1] % n != 0:
            # don't cover a square that is already covered
            if not covered.intersection(piece):
                # remove the piece from the set of uncovered squares
                new_uncovered = uncovered.difference(piece)
                if not new_uncovered:
                    # if nothing is left uncovered, then we have a good tiling
                    return 1
                else:
                    # otherwise, we continue to cover the board
                    new_covered = covered.union(piece)
                    ss += backtrack(new_covered, new_uncovered)

        D[t_c, t_u] = ss
        return ss

    # call the function with an empty covering
    total = backtrack(set([]), set(range(m * n)))
    return total

def triomino_tilings(nrows, ncols):
    """
    triomino_tilings(nrows, ncols):
    This returns the number of ways to tile an m x n grid with triominos.
    A triomino is a shape consisting of three squares joined via the edges.
    There are six possible layouts:

        1       2       3       4       5       6

        00      00      0        0      000     0
        0        0      00      00              0
                                                0
    """
    if nrows < ncols:
        return triomino_tilings(ncols, nrows)

    def next_pieces(root):
        i, j = root
        L = []
        L += [ ((i, j), (i, j + 1), (i + 1, j)) ]
        L += [ ((i, j), (i, j + 1), (i + 1, j + 1)) ]
        L += [ ((i, j), (i + 1, j), (i + 1, j + 1)) ]
        L += [ ((i, j), (i + 1, j), (i + 1, j - 1)) ]
        L += [ ((i, j), (i, j + 1), (i, j + 2)) ]
        L += [ ((i, j), (i + 1, j), (i + 2, j)) ]

        return L

    def backtrack(uncovered, D = {}):
        t_u = tuple(uncovered)
        if t_u in D:
            return D[t_u]

        root = min(uncovered)

        ss = 0
        for piece in next_pieces(root):
            if uncovered.issuperset(piece):
                new_uncovered = uncovered.difference(piece)
                if not new_uncovered:
                    return 1
                else:
                    ss += backtrack(new_uncovered)

        D[t_u] = ss
        return ss

    # create the board, and store the largest number appearing
    board = set([ (i, j) for i in range(nrows) for j in range(ncols) ])
    total = backtrack(board)
    return total

def necklaces(n, k, weights=None):
    """
    Returns the number of necklace colorings.

    Given integers n >= 1 and k >= 2, this returns the number of ways to color a
    necklace of n beads with k colors. A necklace of length n colored by k
    colors is an equivalence class of n-character strings over an alphabet of
    size k, taking all rotations as equivalent.

    The optional argument weights can be a list or tuple of positive integers
    summing to n. If the colors used are k_0, k_1, ..., k_{n - 1}, then
    weights[i] is the number of times that color k_i is used.

    See Charalambides - Enumerative Combinatorics chapter 13 for more details.

    Input:
        n: Positive integer.
        k: Positive integer >= 2.
        weights: Optional list of weights.

    Returns:
        The number of colorings.

    Examples:
        >>> necklaces(3, 2)
        4
        >>> necklaces(10, 5)
        976887

        This is how we would count the colorings of a necklace with 6 beads
        using 3 colors, where each color is used exactly 2 times:
        >>> necklaces(6, 3, [2, 2, 2])
        16
    """
    if weights:
        j = gcd(weights)
        ss = 0
        for d in xdivisors(j):
            ll = [ e / d for e in weights ]
            ss += euler_phi(d)*multinomial(n//d, ll)
        return ss//n
    else:
        ss = 0
        for d in xdivisors(n):
            ss += euler_phi(d)*k**(n//d)
        return ss//n
