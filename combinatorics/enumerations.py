from rosemary.combinatorics.combinatorics import multinomial
from rosemary.number_theory.core import gcd
from rosemary.number_theory.factorization import xdivisors
from rosemary.number_theory.arithmetic_functions import euler_phi
from math import cos, pi

def count_domino_tilings(m, n):
    """
    count_domino_tilings(m, n):
    This returns the number of ways to cover an m x n grid with dominos. The
    formula used is due to Temperley & Fisher, and also by Kasteleyn.
    """
    pp = 1
    for l in xrange(1, m + 1):
        t1 = 4 * pow(cos(pi * l / (m + 1)), 2)
        for k in xrange(1, n + 1):
            t2 = 4 * pow(cos(pi * k / (n + 1)), 2)
            pp *= pow(t1 + t2, 0.25)
    return pp

def domino_tilings(nrows, ncols):
    """
    domino_tilings(nrows, ncols):
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
    largest = nrows * ncols - 1
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
        piece = (smallest, smallest + ncols)

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
        if piece[1] <= largest and piece[1] % ncols != 0:
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
    total = backtrack(set([]), set(range(nrows * ncols)))
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
            ss += euler_phi(d) * multinomial(n // d, ll)
        return ss // n
    else:
        ss = 0
        for d in xdivisors(n):
            ss += euler_phi(d) * k**(n // d)
        return ss // n

