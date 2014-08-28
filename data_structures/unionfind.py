class NamedUnionFind(object):
    def __init__(self, node_list):
        """
        Initialize a new NamedUnionFind object containing nodes in node_list.
        """
        self.count = {}
        self.name = {}
        self.father = {}
        self.root = {}
        for u in node_list:
            self.count[u] = 1
            self.name[u] = u
            self.father[u] = None
            self.root[u] = u

    def __getitem__(self, u):
        return self.find(u)

    def find(self, u):
        """
        Returns the name of the set containing u.
        """
        path = []
        while self.father[u] is not None:
            path.append(u)
            u = self.father[u]
        for v in path:
            self.father[v] = u
        return self.name[u]

    def union(self, u, v, w):
        """
        Combines the sets named u and v into a new set named w.
        """
        if self.count[self.root[u]] > self.count[self.root[v]]:
            u, v = v, u
        large = self.root[v]
        small = self.root[u]
        self.father[small] = large
        self.count[large] += self.count[small]
        self.name[large] = w
        self.root[w] = large


class UnionFind(object):
    def __init__(self, node_list):
        """
        Initialize a new UnionFind object containing nodes in node_list.
        """
        rank = {}
        father = {}

        for u in node_list:
            rank[u] = 0
            father[u] = u

        self.rank = rank
        self.father = father

    def __getitem__(self, u):
        return self.find(u)

    def find(self, x):
        """
        Returns the name of the set containing x.
        """
        father = self.father
        path = []
        while x != father[x]:
            path.append(x)
            x = father[x]

        for v in path:
            father[v] = x

        return x

    def union(self, x, y):
        """
        Combines the sets named x and y.
        """
        father = self.father
        rank = self.rank

        fx = father[x]
        fy = father[y]

        if fx == fy:
            return

        if rank[fx] > rank[fy]:
            father[fy] = fx
        else:
            father[fx] = fy

            if rank[fx] == rank[fy]:
                rank[fy] += 1
