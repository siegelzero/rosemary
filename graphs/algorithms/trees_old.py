import heapq
import networkx
from rosemary.data_structures.unionfind import NamedUnionFind

def spanning_trees_ordered(G, k=None):
    """
    Returns an object that iterates over the spanning trees of a graph.

    Given a networkx graph G, this returns an iterator over the spanning trees
    of G, listed in order of increasing weight. The algorithm is based on the
    paper "Two Algorithms for Generating Weighted Spanning Trees in Order" by
    H.N.  Gabow.

    Input:
        * G: NetworkX graph

        * k: positive integer, optional (default=None)
            Number of spanning trees to return. If k is None, then all are
            returned.

    Output:
        * trees: iterator
            Iterator over the spanning trees. The trees are represented as
            lists of edges.

    Examples:
        >>> G = networkx.Graph()
        >>> G.add_weighted_edges_from([(1, 2, 3), (2, 3, 4), (3, 4, 5), \
            (1, 4, 5), (1, 3, 4)])
        >>> X = spanning_trees_ordered(G)
        >>> list(X)
        [[(2, 1), (3, 1), (4, 1)],
         [(2, 1), (3, 2), (4, 1)],
         [(2, 1), (3, 1), (4, 3)],
         [(2, 1), (3, 2), (4, 3)],
         [(2, 3), (3, 1), (4, 1)],
         [(2, 3), (3, 1), (4, 3)],
         [(2, 1), (3, 4), (4, 1)],
         [(2, 3), (3, 4), (4, 1)]]
        >>> H = networkx.Graph()
        >>> H.add_weighted_edges_from([(1, 2, 5), (1, 3, 5), (1, 4, 6),\
            (2, 3, 2), (2, 5, 7), (2, 7, 4), (3, 4, 1), (3, 5, 2), (4, 5, 3),\
            (4, 6, 2), (5, 6, 1), (5, 7, 2), (6, 7, 3)])
        >>> Y = spanning_trees_ordered(H)
        >>> len(list(Y))
        615
        >>> list(spanning_trees_ordered(H, 5))
        [[(2, 1), (3, 2), (4, 3), (5, 3), (6, 5), (7, 5)],
         [(2, 1), (3, 2), (4, 3), (5, 6), (6, 4), (7, 5)],
         [(2, 3), (3, 1), (4, 3), (5, 3), (6, 5), (7, 5)],
         [(2, 3), (3, 1), (4, 3), (5, 6), (6, 4), (7, 5)],
         [(2, 1), (3, 2), (4, 3), (5, 4), (6, 5), (7, 5)]]
        >>> K8 = networkx.complete_graph(8)
    """
    def first_common_ancestor(X, father, x, y):
        """
        Finds the first common ancestor of nodes x and y.

        In a tree with root node r, the first common ancestor of nodes x and y
        is the first node in both the path from x to r and the path from y to
        r.

        Input:
            * X: UnionFind structure

            * father: dict
                This is the father dict representing the tree we're looking at.

            * x: node
                First node.

            * y: node
                Second node.

        Output:
            * a: node
                The first common ancestor of x and y. If there is no common
                ancestor (i.e. the graph is not connected), then None is
                returned.
        """
        # First, build the path from x to the root node.
        x1 = X.find(x)
        x_path = [x1]
        while True:
            x2 = X.find(father[x1])
            if x1 == x2:
                break
            x_path.append(x2)
            x1 = x2

        # Next, build the path from y to the root node.
        y1 = X.find(y)
        y_path = [y1]
        while True:
            y2 = X.find(father[y1])
            if y1 == y2:
                break
            y_path.append(y2)
            y1 = y2

        # Finally, find the first common element of these lists.
        for ancestor in x_path:
            if ancestor in y_path:
                return ancestor

        # Otherwise, return no ancestor was found.
        return None

    def find_exchange(father, included, excluded):
        """
        Finds a minimum exchange.

        Given a spanning tree T represented by the dict father, this finds a
        minimum weight exchange to perform.

        Input:
            * father: dict
                This is the father dict of the tree T we're currently looking
                at.
            * included: list
                This is a list of edges in T that must remain in T.
            * excluded: list
                This is a list of edges out of T that must remain out of T.

        Output:
            * (r, e, f): tuple
                * r: integer
                    This is the weight of the tree.
                * e: edge
                    This is the edge to remove.
                * f: edge
                    This is the edge to add.
        """
        X = NamedUnionFind(nodes)
        (r, e, f) = (inf, None, None)

        # First, make all edges in the included list ineligible. This is done
        # by by placing x and y in the same set (since they will have the same
        # first eligible ancestor).
        for edge in included:
            # We make it so that y = father[x]
            x, y = edge
            if y != father[x]:
                x, y = y, x
            y = X.find(y)
            X.union(x, y, y)

        # Next, we mark all edges in the excluded list so they are not
        # considered for exchanges. These are later unmarked.
        for edge in excluded:
            marked.add(edge)

        for (_, x, y) in edges:
            if (x, y) in marked or (y, x) in marked:
                marked.discard((x, y))
                marked.discard((y, x))
            elif father[x] != y and father[y] != x:
                a = first_common_ancestor(X, father, x, y)
                if a is None:
                    continue
                for u in [x, y]:
                    v = X.find(u)
                    while v != a:
                        exchange_weight = weight[x, y] - weight[v, father[v]]
                        if exchange_weight < r:
                            r = exchange_weight
                            e = (v, father[v])
                            f = (x, y)
                        w = X.find(father[v])
                        X.union(v, w, w)
                        v = w
        return (r, e, f)

    inf = float('inf')
    if k is None:
        k = inf

    # Dictionary of weights, with default edge weight 1 for unweighted graphs.
    weight = {(u, v): G[u][v].get('weight', 1) for u in G for v in G[u]}

    # Set of marked edges.
    marked = set()

    # Compute a minimum spanning tree and its weight.
    T = networkx.minimum_spanning_tree(G).edges()
    t = sum(weight[u, v] for u in G for v in G[u])

    # Edges of G in order of increasing weight.
    edges = sorted([(weight[u, v], u, v) for u in G for v in G[u] if v < u])
    nodes = G.nodes()

    # father[u] is the next node after u in the path from u to the root node 1.
    father = networkx.bfs_predecessors(networkx.from_edgelist(T), 1)
    father[1] = 1

    (r, e, f) = find_exchange(father, [], [])
    P = [(t + r, e, f, father, [], [])]

    yield [(x, y) for (x, y) in father.iteritems() if x != y]

    j = 1
    while j < k:
        (t, e, f, father, included, excluded) = heapq.heappop(P)

        if t == inf:
            return

        G_new = networkx.from_edgelist(father.items())
        G_new.remove_edge(*e)
        G_new.add_edge(*f)

        father_new = networkx.bfs_predecessors(G_new, 1)
        father_new[1] = 1

        yield [(x, y) for (x, y) in father_new.iteritems() if x != y]

        ti = t - weight[f] + weight[e]

        included_i = included + [e]
        excluded_j = excluded + [e]

        (r, e, f) = find_exchange(father, included_i, excluded)
        heapq.heappush(P, (ti + r, e, f, father, included_i, excluded))

        (r, e, f) = find_exchange(father_new, included, excluded_j)
        heapq.heappush(P, (t + r, e, f, father_new, included, excluded_j))
        j += 1


def spanning_trees(G):
    """
    Returns an object that iterates over the spanning trees of a graph.

    Given a networkx graph G, this returns an iterator over the spanning trees
    of G. The method used is a not-so-smart backtracking routine enhanced by a
    UnionFind data structure. Since many dict copies are performed, we avoid
    the overhead of using the UnionFind data structure.

    Input:
        * G: NetworkX graph.

    Output:
        * trees: iterator
            Iterator over the spanning trees.

    Examples:
        >>> G = networkx.Graph()
        >>> G.add_weighted_edges_from([(1, 2, 3), (2, 3, 4), (3, 4, 5),\
                (1, 4, 5), (1, 3, 4)])
        >>> X = spanning_trees(G)
        print list(X)
        [[(1, 4), (1, 3), (1, 2)]
         [(2, 3), (1, 4), (1, 2)]
         [(2, 3), (1, 4), (1, 3)]
         [(3, 4), (1, 3), (1, 2)]
         [(3, 4), (1, 4), (1, 2)]
         [(3, 4), (2, 3), (1, 2)]
         [(3, 4), (2, 3), (1, 3)]
         [(3, 4), (2, 3), (1, 4)]]
        >>> H = networkx.Graph()
        >>> H.add_weighted_edges_from([(1, 2, 5), (1, 3, 5), (1, 4, 6),\
                (2, 3, 2), (2, 5, 7), (2, 7, 4), (3, 4, 1), (3, 5, 2),\
                (4, 5, 3), (4, 6, 2), (5, 6, 1), (5, 7, 2), (6, 7, 3)])
        >>> Y = spanning_trees(H)
        >>> len(list(Y))
        615
        >>> K8 = networkx.complete_graph(8)
        >>> len(list(spanning_trees(K8)))
        262144
    """
    def find(parents, u):
        """
        This returns the name of the set containing u.
        """
        # Find without path compression. This is faster is some tests.
        while u != parents[u]:
            u = parents[u]
        return u

        #Below is the path compression heuristic. This seems to be slower in
        #some tests (possibly due to the recursion)
        #if u != parents[u]:
        #    parents[u] = find(parents, parents[u])
        #return parents[u]

    def union(parents, weights, u, v):
        """
        This subroutine combines the sets containing u and v into a new set.
        """
        u = find(parents, u)
        v = find(parents, v)
        if weights[u] > weights[v]:
            parents[v] = u
        else:
            parents[u] = v
            if weights[u] == weights[v]:
                weights[v] += 1

    edges = G.edges()
    n = G.number_of_nodes() - 1
    m = len(edges)

    def backtrack(partial, last, k, parents, weights):
        """
        This routine works by extending a tree until it has the correct number
        of edges to be a spanning tree. We check for cycles by keeping a
        disjoint set data structure in the dicts "find" and "union".

        Input:
            * partial: list
                This is a list of edges representing a subtree of the graph. We
                extend each tree until it covers every node; i.e. becomes a
                spanning tree.

            * last: int
                This is the index of the last edge added to the tree.

            * k: int
                This is the number of edges in the current tree "partial".

            * parents: dict
                This is a dict of parent nodes used in checking for cycles.

            * weights: dict
                This is a dict of tree weights used in the disjoint set
                structure that we maintain.
        """
        # We extend the current tree by adding each edge that does not create a
        # cycle to it, then recursively try to extend each of these new trees
        # in the same manner. We only add edges of greater index than the
        # current.
        for i in xrange(last + 1, m):
            u, v = edges[i]
            # Check if the edge (u, v) will create a cycle
            if find(parents, u) != find(parents, v):
                if k + 1 == n:
                    # Any spanning tree will have n edges in it.
                    yield partial + [(u, v)]
                else:
                    # Extend the current tree and pass copies of the disjoint
                    # set information.
                    new_parents = parents.copy()
                    new_weights = weights.copy()
                    union(new_parents, new_weights, u, v)
                    for t in backtrack(partial + [(u, v)], i, k + 1,
                                       new_parents, new_weights):
                        yield t

    # A spanning tree must have exactly n edges in it, so we stop when this is
    # impossible.
    for i in xrange(m - n + 1):
        u, v = edges[i]
        parents = {n: n for n in G}
        weights = {n: 0 for n in G}
        union(parents, weights, u, v)
        for t in backtrack([(u, v)], i, 1, parents, weights):
            yield t

def kruskal(G):
    """
    kruskal(G):
    Given a weighted graph G, this algorithm returns a minimal spanning tree
    for G using Kruskal's algorithm. The algorithm builds a spanning tree by
    successively inserting edges from G in order of increasing cost. We
    insert each edge e as long as it does not create a cycle when added to
    the edges we've already inserted. If, on the other hand, inserting e
    would result in a cycle, then we simply discard e and continue.
    """
    def find(parents, u):
        """
        This subroutine returns a representative element from C that contains
        u. Using this, we can determine whether two vertices u and v belong to
        the same tree by testing whether find(C, u) = find(C, v).
        """
        if parents[u] != u:
            parents[u] = find(parents, parents[u])
        return parents[u]

    def union(parents, weights, u, v):
        """
        This subroutine
        """
        u = find(parents, u)
        v = find(parents, v)
        if weights[u] > weights[v]:
            parents[v] = u
        else:
            parents[u] = v
        if weights[u] == weights[v]:
            weights[v] += 1

    E = [(G[u][v], u, v) for u in G for v in G[u]]
    T = set()
    parents = {u: u for u in G}
    weights = {u: 0 for u in G}

    for (_, u, v) in sorted(E):
        if find(parents, u) != find(parents, v):
            T.add((u, v))
            union(parents, weights, u, v)

    return T
