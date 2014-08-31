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
