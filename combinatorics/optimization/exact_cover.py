from collections import defaultdict


class ExactCover(object):
    """
    Class representing the Exact Cover problem.

    Given a collection S of subsets of the set X, does there exist a
    subcollection S' such that every element of X is contained in exactly one
    member of S'?
    """
    def __init__(self, sets):
        self.sets = sets.copy()
        universe = set()

        for (key, value) in sets.iteritems():
            universe.update(value)

        self.universe = universe
        self.num_cols = len(universe)

    def __repr__(self):
        return "Exact Cover Problem"


class ExactCoverBinary(ExactCover):
    """
    Class representing the Exact Cover problem.

    Given a matrix M with entries of 0 and 1, does there exist a subset of the
    rows of M which sum to the vector of all 1s?
    """
    def __init__(self, integers):
        vectors = defaultdict(list)
        num_bits = 0

        for n in integers:
            least_set_bit = n & (-n)
            num_bits = max(num_bits, n.bit_length())
            vectors[least_set_bit].append(n)

        self.vectors = vectors
        self.num_bits = num_bits

    def number_of_solutions(self):
        """
        Returns the number of exact covers. The algorithm is a memoized
        recursive backtracking method.
        """
        vectors = self.vectors
        all_bits = 2**self.num_bits - 1
        cache = {all_bits: 1}

        def backtrack(covered):
            if covered in cache:
                return cache[covered]

            total = 0
            least_unset_bit = ~covered & (covered + 1)
            for v in vectors[least_unset_bit]:
                if (covered & v) == 0:
                    total += backtrack(covered | v)

            cache[covered] = total
            return total

        total = backtrack(0)
        return total

    def solutions(self):
        all_bits = 2**self.num_bits - 1
        vectors = self.vectors
        stack = [(0, [])]

        while stack:
            (covered, used) = stack.pop()
            if covered == all_bits:
                yield used
            else:
                least_unset_bit = ~covered & (covered + 1)
                for v in vectors[least_unset_bit]:
                    if (covered & v) == 0:
                        stack.append((covered | v, used + [v]))


def langford_ints(n):
    L = []
    for p in xrange(1, n + 1):
        c = 1 << (3*n - p)
        l = 2*n - 1
        r = l - p - 1
        while r >= 0:
            b = c | (1 << l)
            b |= (1 << r)
            L.append(b)
            r -= 1
            l -= 1
    return L

def langford_pairings(n):
    ints = langford_ints(n)
    E = ExactCoverBinary(ints)
    num_sols = E.number_of_solutions()
    return num_sols // 2

def domino_ints(n):
    L = []
    for i in xrange(n*n):
        if i + n < n*n:
            L.append((1 << i) | (1 << (i + n)))

        if i + 1 < n*n and (i + 1) % n != 0:
            L.append((1 << i) | (1 << (i + 1)))

    return L

def domino_tilings(n):
    ints = domino_ints(n)
    E = ExactCoverBinary(ints)
    num_sols = E.number_of_solutions()
    return num_sols
