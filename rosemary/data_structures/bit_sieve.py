class BITSieve(object):
    def __init__(self, n):
        """
        Initializes sieve of [1, 2, ..., n]
        """
        A = {}
        parent = {}

        for j in xrange(1, n + 1):
            A[0, j] = 1

        i = 1
        while n >> i:
            for j in xrange(1, (n >> i) + 1):
                A[i, j] = A[i - 1, 2*j - 1] + A[i - 1, 2*j]
                parent[i - 1, 2*j - 1] = (i, j)
                parent[i - 1, 2*j] = (i, j)

            i += 1

        self.A = A
        self.n = n
        self.parent = parent

    def entries(self):
        return [self.A[0, i] for i in xrange(1, self.n + 1)]

    def update(self):
        A = self.A
        n = self.n
        i = 1
        while n >> i:
            for j in xrange(1, (n >> i) + 1):
                A[i, j] = A[i - 1, 2*j - 1] + A[i - 1, 2*j]
            i += 1

    def mark_multiples(self, p):
        # A = self.A
        # for i in xrange(p, self.n + 1, p):
        #     A[0, i] = 0

        # self.update()

        A = self.A
        parent = self.parent
        for i in xrange(p, self.n + 1, p):
            if A[0, i] == 0:
                continue

            pos = (0, i)
            A[pos] = 0

            while pos in parent:
                A[parent[pos]] -= 1
                pos = parent[pos]

    def partial_sum(self, m):
        A = self.A
        total = 0
        i = 0

        while m:
            r = m % 2
            if r:
                total += A[i, m]
            m >>= 1
            i += 1

        return total


class BITSieveArray(object):
    def __init__(self, n):
        """
        Initializes sieve of [1, 2, ..., n]
        """
        A = [[1]*(n + 1)]
        A[0][0] = 0

        i = 1
        while n >> i:
            A.append([0]*((n >> i) + 1))
            for j in xrange(1, (n >> i) + 1):
                A[i][j] = A[i - 1][2*j - 1] + A[i - 1][2*j]

            i += 1

        self.A = A
        self.n = n

    def entries(self):
        return self.A[0]

    def mark_multiples(self, p):
        A = self.A

        for i in xrange(p, self.n + 1, p):
            if A[0][i] == 0:
                continue

            y = i
            for x in xrange(len(A)):
                if y >= len(A[x]):
                    break
                A[x][y] -= 1
                y = (y + 1)//2

    def partial_sum(self, m):
        A = self.A
        total = 0
        i = 0

        while m:
            r = m % 2
            if r:
                total += A[i][m]
            m >>= 1
            i += 1

        return total
