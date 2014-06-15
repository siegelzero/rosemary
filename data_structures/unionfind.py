class NamedUnionFind(object):
    def __init__(self, node_list):
        """
        Initialize a new UnionFind object containing nodes in node_list.
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


