from heapq import heappush, heappop


###############################################################################
# Algorithms for minimum spanning trees
###############################################################################


def prim(graph, w=None, edge_list=True):
    """
    Returns a minimum spanning tree of graph of the given graph.

    Given a connected weighted graph and optional vertex w, this method returns
    a minimum spanning tree of the graph, rooted at vertex w. In the case that
    the graph is not connected, a minimum spanning tree of the component
    containing the root vertex is returned.

    Input:
        * graph: Graph

        * w: vertex of graph (default=None)
            Root vertex of spanning tree. If w is None, a (somewhat) random
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

        Using Python heaps, this algorithm runs in time O((m + n)*log(n)), where
        n and m are the number of vertices and edges of the graph, respectively.
        Our implementation is based on the exposition in the book "Python
        Algorithms" by Lie Hetland. Another solid reference is the book
        "Algorithms" by Dasgupta, et al.
    """
    # Choose a root vertex if none is given
    if w is None:
        w = graph.graph_dict.keys()[0]

    previous = {}
    heap = [(0, None, w)]

    while heap:
        (_, p, u) = heappop(heap)

        if u in previous:
            continue

        previous[u] = p

        for (v, w) in graph[u].iteritems():
            heappush(heap, (w, u, v))

    # We have the predecessor dict. Now we compute the edges in the tree.
    edges = []
    total = 0
    for u in previous:
        v = previous[u]
        if v is not None:
            edge = min(u, v), max(u, v)
            weight = graph[u][v]
            total += weight
            edges.append(edge)

    if edge_list:
        return total, sorted(edges)

    return total, previous
