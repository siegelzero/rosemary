class NamedUnionFind(object):
    def __init__(self, node_list):
        """
        Initialize a new NamedUnionFind object containing nodes in node_list.
        """
        count = {}
        name = {}
        father = {}
        root = {}

        for u in node_list:
            count[u] = 1
            name[u] = u
            father[u] = None
            root[u] = u

        self.count = count
        self.name = name
        self.father = father
        self.root = root

    def find(self, u):
        """
        Returns the name of the set containing u.
        """
        father = self.father
        path = []
        append = path.append

        while father[u] is not None:
            append(u)
            u = father[u]

        for v in path:
            father[v] = u

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
    """
    Implementation of a disjoint set data structure.

    This data structure solves the problem of maintaining a collection of
    disjoint sets under the operation of union. To determine if two elements
    live in the same set, we implement a find operation which returns the name
    of the set containing the given element.

    Our implementation follows Chapter 2 of "Data Structures and Network
    Algorithms" by Tarjan. We implement both the union by rank and the path
    compression heuristics for fast amortized runtime.
    """
    def __init__(self, elements):
        """
        Initialize a new UnionFind object containing the given elements.

        Given a list of elements, we initialize the data structure by putting
        each element into its own set. We represent each set by a rooted tree.
        The nodes of the tree are the elements of the set, and the
        representative of each set is the tree root. For each node x, we keep
        track of its parent node. By convention, we make the root element its
        own parent.
        """
        rank = {}
        parent = {}

        for u in elements:
            rank[u] = 0
            parent[u] = u

        self.rank = rank
        self.parent = parent

    def find(self, x):
        """
        Returns the name of the set containing the element x.

        Here, we use the path compression heuristic. This changes the structure
        of the tree during a find by moving nodes closer to the root. When
        carrying out find(x), after locating the root r of the tree containing
        x, we make every node on the path from x to r point directly to r. This
        heuristic increases the time of a single find by a constant factor, but
        saves enough time in later finds to more than pay for itself.
        """
        # We follow parent pointers from x to the root of the tree containing x.
        parent = self.parent
        path = []
        append = path.append

        while x != parent[x]:
            append(x)
            x = parent[x]

        # Make every node on the path from x to the root point directly to the
        # root.
        for v in path:
            parent[v] = x

        return x

    def union(self, x, y):
        """
        Combines the sets named x and y.

        Here, we use the union by rank heuristic to keep the trees shallow.
        With each node x, we store a nonnegative integer rank that is an upper
        bound on the height of x.
        """
        parent = self.parent
        rank = self.rank

        px = parent[x]
        py = parent[y]

        if rank[px] > rank[py]:
            x, y = y, x

        if rank[px] == rank[py]:
            rank[py] += 1

        parent[px] = py
