from collections import defaultdict


def to_string(entries):
    s = []
    for i in xrange(9):
        for j in xrange(9):
            s.append(str(entries[i][j]))
    return ''.join(s)


def parse(digits):
    M = [[0]*9 for _ in xrange(9)]
    for (i, d) in enumerate(digits):
        M[i//9][i % 9] = int(d)
    return M


def algorithmx(target, sets):
    solution = []

    Y = {}
    for i in xrange(len(sets)):
        Y[i] = list(sets[i])

    X = {j: set() for j in target}
    for i in Y:
        for j in Y[i]:
            X[j].add(i)

    def solve():
        if not X:
            yield [sets[i] for i in solution]
        else:
            c = min(X, key=lambda c: len(X[c]))

            for r in list(X[c]):
                solution.append(r)

                cols = []
                for j in Y[r]:
                    for i in X[j]:
                        for k in Y[i]:
                            if k != j:
                                X[k].remove(i)
                    cols.append(X.pop(j))

                for s in solve():
                    yield s

                for j in reversed(Y[r]):
                    X[j] = cols.pop()
                    for i in X[j]:
                        for k in Y[i]:
                            if k != j:
                                X[k].add(i)

                solution.pop()

    return solve()


def algorithm_dlx(target, sets):
    solution = []

    Y = {}
    for i in xrange(len(sets)):
        Y[i] = list(sets[i])

    X = {j: set() for j in target}
    for i in Y:
        for j in Y[i]:
            X[j].add(i)

    def solve():
        if not X:
            yield [sets[i] for i in solution]
        else:
            c = min(X, key=lambda c: len(X[c]))

            for r in list(X[c]):
                solution.append(r)

                cols = []
                for j in Y[r]:
                    for i in X[j]:
                        for k in Y[i]:
                            if k != j:
                                X[k].remove(i)
                    cols.append(X.pop(j))

                for s in solve():
                    yield s

                for j in reversed(Y[r]):
                    X[j] = cols.pop()
                    for i in X[j]:
                        for k in Y[i]:
                            if k != j:
                                X[k].add(i)

                solution.pop()

    return solve()


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


def latin_square_ints(n):
    ints = []
    all_blocks = {}
    for i in xrange(n):
        all_blocks[i] = []
        for j in xrange(n):
            blocks = [[0]*n for _ in xrange(n)]
            blocks[j][i] = 1
            all_blocks[i].append(blocks)

    for p in xrange(n):
        for (i, a) in enumerate(all_blocks[p]):
            for (j, b) in enumerate(all_blocks[p]):
                t = []
                for s in a + b:
                    t.extend(s)
                val = sum(2**(2*n*n - k - 1) for k in xrange(2*n*n) if t[k])
                ints.append((val << (n*n)) | (1 << (n*n - j*n - i - 1)))

    return ints


def sudoku_vectors(M):
    vectors = []
    block_idx = {
        (0, 0): 0,
        (0, 3): 1,
        (0, 6): 2,
        (3, 0): 3,
        (3, 3): 4,
        (3, 6): 5,
        (6, 0): 6,
        (6, 3): 7,
        (6, 6): 8
    }

    for row in range(9):
        for col in range(9):
            if M[row][col] == 0:
                digits = range(1, 10)
            else:
                digits = [M[row][col]]

            for digit in digits:
                entries = set()
                entries.add("{}O{}".format(row, col))
                entries.add("{}R{}".format(digit, row))
                entries.add("{}C{}".format(digit, col))

                rr = row - (row % 3)
                cc = col - (col % 3)
                entries.add("{}B{}".format(digit, block_idx[rr, cc]))

                vectors.append(entries)

    return vectors


def sudoku_dlx(a):
    M = parse(a)
    vecs = sudoku_vectors(M)
    target = set()
    for e in vecs:
        target.update(e)

    sol = algorithmx(target, vecs).next()

    for entry in sol:
        for e in entry:
            if e[1] == 'O':
                r = int(e[0])
                c = int(e[2])

            if e[1] == 'C':
                d = int(e[0])

        M[r][c] = d

    return to_string(M)
