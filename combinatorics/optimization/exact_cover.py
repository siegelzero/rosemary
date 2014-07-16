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
        self.all_bits = 2**num_bits - 1


    def number_of_solutions(self):
        """
        Returns the number of exact covers. The algorithm is a memoized
        recursive backtracking method.
        """
        vectors = self.vectors
        cache = {self.all_bits: 1}

        def backtrack(bits):
            if bits in cache:
                return cache[bits]

            total = 0
            least_unset_bit = ~bits & (bits + 1)
            for v in vectors[least_unset_bit]:
                if (bits & v) == 0:
                    total += backtrack(bits | v)

            cache[bits] = total
            return total

        total = backtrack(0)
        return total


    def solutions(self):
        all_bits = self.all_bits
        vecs = self.vectors
        stack = [(0, [])]

        while stack:
            (bits, used) = stack.pop()

            if bits == all_bits:
                yield used
            else:
                bit = ~bits & (bits + 1)
                for v in vecs[bit]:
                    if (bits & v) == 0:
                        stack.append((bits | v, used + [v]))


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
