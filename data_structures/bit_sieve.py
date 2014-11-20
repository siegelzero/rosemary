class BITSieve(object):
    def __init__(self, n):
        self.n = n
        entries = [0]*(n + 1)

        for i in xrange(1, n + 1):
            k = i
            while k <= n:
                entries[k] += 1
                k |= k + 1

        self.entries = entries

    def add(self, i, d):
        entries = self.entries
        n = self.n

        while i <= n:
            entries[i] += d
            i |= i + 1

    def sum(self, i=None):
        if i is None:
            i = self.n

        total = 0
        entries = self.entries

        while i >= 0:
            total += entries[i]
            i &= i + 1
            i -= 1

        return total

    def mark(self, p):
        for i in xrange(p, self.n + 1, p):
            v = self.sum(i) - self.sum(i - 1)
            if v:
                self.add(i, -v)
