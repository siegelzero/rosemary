from heapq import heappush, heappop, heapify, merge
from rosemary.data_structures.unionfind import UnionFind


###############################################################################
# Algorithms for minimum spanning trees
###############################################################################


def prim(graph, root=None, edge_list=True):
    """
    Returns a minimum spanning tree of the given graph.

    Given a connected weighted graph and optional vertex w, this method returns
    a minimum spanning tree of the graph, rooted at vertex w. In the case that
    the graph is not connected, a minimum spanning tree of the component
    containing the root vertex is returned.

    Input:
        * graph: Graph

        * root: vertex of graph (default=None)
            Root vertex of spanning tree. If root is None, a (somewhat) random
            vertex is chosen.

        * edge_list: bool (default=True)
            If True, a list of the edges in the minimum spanning tree is
            returned. Otherwise, a predecessor dict is returned.

    Output:
        * (weight, edges): tuple
            * weight: number
                Weight of the minimum spanning tree.

            * edges: list
                If edge_list is True, this is a list of the edges in the minimum
                spanning tree. Otherwise, this is the predecessor dict of the
                spanning tree.

    Examples:
        >>> G = Graph()
        >>> G.add_edges([('a', 'b', 5), ('a', 'c', 6), ('a', 'd', 4),
                         ('b', 'c', 1), ('b', 'd', 2), ('c', 'd', 2),
                         ('c', 'e', 5), ('c', 'f', 3), ('d', 'f', 4),
                         ('e', 'f', 4)])
        >>> prim(G)
        (14, [('a', 'd'), ('b', 'c'), ('b', 'd'), ('c', 'f'), ('e', 'f')])
        >>> prim(G, edge_list=False)
        (14, {'a': None, 'b': 'd', 'c': 'b', 'd': 'a', 'e': 'f', 'f': 'c'})

    Details:
        This method uses Prim's algorithm to find a minimum spanning tree of the
        given graph. The notion of a spanning tree only makes sense if the graph
        is connected. In the case that the graph is not connected, then a
        minimum spanning tree of the component containing the root vertex is
        returned.

        Prim's algorithm begins with a tree containing only the root vertex, and
        then repeatedly grows the tree by adding the lightest edge between a
        vertex in the tree and a vertex outside of the tree. Using Python heaps,
        this algorithm runs in time O((|V| + |E|)*log(|V|)). Our implementation
        is based on the exposition in the book "Python Algorithms" by Lie
        Hetland. Another standard reference is the book "Algorithms" by
        Dasgupta, et al. See "Network Flows: Theory, Algorithms, and
        Applications" by Ahuja, et al for detailed overview of minimum spanning
        tree algorithms.
    """
    # Choose a root vertex if none is given
    if root is None:
        root = graph.graph_dict.keys()[0]

    previous = {}
    heap = [(0, None, root)]

    while heap:
        (_, p, u) = heappop(heap)

        if u in previous:
            continue

        previous[u] = p

        # Note that in typical implementations of Prim's algorithm, we relax
        # before we push to the heap. Here, we use the heap property to get
        # around this.
        for (v, w) in graph[u].iteritems():
            heappush(heap, (w, u, v))

    # We have the predecessor dict. Now we compute the edges in the tree.
    edges = []
    total = 0
    for u in previous:
        v = previous[u]
        if v is not None:
            edge = min(u, v), max(u, v)
            total += graph[u][v]
            edges.append(edge)

    if edge_list:
        return total, sorted(edges)

    return total, previous


def kruskal(graph):
    """
    Returns the edge list of a minimum spanning tree of the given graph.

    Given a connected weighted graph and optional vertex w, this method returns
    a minimum spanning tree of the graph. In the case that the graph is not
    connected, it returns a list of the edges of a minimum spanning tree for
    each connected component.

    Input:
        * graph: Graph

    Output:
        * (weight, edges): tuple
            * weight: number
                Weight of the minimum spanning tree.

            * edges: list
                A list of the edges in the minimum spanning tree.

    Examples:
        >>> G = Graph()
        >>> G.add_edges([('a', 'b', 5), ('a', 'c', 6), ('a', 'd', 4),
                         ('b', 'c', 1), ('b', 'd', 2), ('c', 'd', 2),
                         ('c', 'e', 5), ('c', 'f', 3), ('d', 'f', 4),
                         ('e', 'f', 4)])
        >>> kruskal(G)
        (14, [('a', 'd'), ('b', 'c'), ('b', 'd'), ('c', 'f'), ('e', 'f')])

    Details:
        This method uses Kruskal's algorithm to find a minimum spanning tree of
        the given graph. The notion of a spanning tree only makes sense if the
        graph is connected. In the case that the graph is not connected, then a
        list of the edges of a minimum spanning forest is returned.

        Kruskal's algorithm begins with an empty tree, and then repeatedly adds
        the next lightest edge of the graph that doesn't produce a cycle, using
        a disjoint-set data structure for cycle detection. Our version of
        Kruskal's algorithm runs in time O(|E| log(|V|)). See "Algorithms" by
        Dasgupta, et al for details. See "Network Flows: Theory, Algorithms, and
        Applications" by Ahuja, et al for detailed overview of minimum spanning
        tree algorithms.
    """
    vertices = graph.vertex_set()
    num_vertices = len(vertices)

    edges = [(graph[u][v], u, v) for u in graph for v in graph[u]]
    heapify(edges)

    tree_edges = []
    total = 0
    num_edges = 0

    X = UnionFind(vertices)
    find = X.find
    union = X.union

    # A tree on n vertices has n - 1 edges, so break when finished.
    while num_edges < num_vertices - 1:
        (weight, u, v) = heappop(edges)
        if find(u) != find(v):
            edge = min(u, v), max(u, v)
            total += weight
            tree_edges.append(edge)
            num_edges += 1
            union(u, v)

    return total, sorted(tree_edges)


def cheriton_tarjan(graph):
    import random
    vertices = graph.vertices()
    X = UnionFind(vertices)
    find = X.find
    union = X.union
    PQ = {}

    edges = []

    for u in vertices:
        Q = []
        for v in graph[u]:
            w = graph[u][v]
            Q.append((w, u, v))
        heapify(Q)
        PQ[u] = Q

    t = 0
    while t < len(vertices) - 1:
        u = random.choice(vertices)
        Tu = find(u)

        while True:
            (_, x, y) = heappop(PQ[Tu])
            Tx = find(x)
            Ty = find(y)

            if Tx == Ty:
                continue

            if Tx == Tu:
                Tv = Ty
            else:
                Tv = Tx

            break

        union(Tu, Tv)
        u = find(Tu)

        Q = list(merge(PQ[Tu], PQ[Tv]))
        heapify(Q)

        del PQ[Tv]
        del PQ[Tu]

        PQ[u] = Q

        edges.append((min(x, y), max(x, y)))

        t += 1

    total = 0
    for (u, v) in edges:
        total += graph[u][v]

    return total, sorted(edges)
