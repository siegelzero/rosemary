from heapq import heappush, heappop, heapify
from random import choice

from rosemary.data_structures.unionfind import UnionFind, NamedUnionFind
from rosemary.data_structures.heaps import LeftistHeap

from rosemary.graphs.algorithms.traversal import breadth_first_search_tree
from rosemary.graphs.graphs import Graph


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
        is based on the exposition in the book "Introduction to Algorithms" by
        Cormen, et al. Another standard reference is the book "Algorithms" by
        Dasgupta, et al. See "Network Flows: Theory, Algorithms, and
        Applications" by Ahuja, et al. for a detailed overview of minimum
        spanning tree algorithms.
    """
    # Choose a root vertex if none is given
    if root is None:
        root = graph.graph_dict.keys()[0]

    inf = float('inf')

    previous = {u: None for u in graph}
    cost = {u: inf for u in graph}
    cost[root] = 0

    heap = [(0, root)]
    visited = set()
    add = visited.add

    while heap:
        (_, u) = heappop(heap)

        if u in visited:
            continue
        add(u)

        for (v, w) in graph[u].iteritems():
            if v not in visited and w < cost[v]:
                cost[v] = w
                previous[v] = u
                heappush(heap, (w, v))

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
        Applications" by Ahuja, et al for a detailed overview of minimum
        spanning tree algorithms.
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
    """
    Returns the edge list of a minimum spanning tree of the given graph.

    Given a connected weighted graph, this method returns a minimum spanning
    tree of the graph. In the case that the graph is not connected, it returns a
    list of the edges of a minimum spanning tree for each connected component.

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
        >>> cheriton_tarjan(G)
        (14, [('a', 'd'), ('b', 'c'), ('b', 'd'), ('c', 'f'), ('e', 'f')])

    Details:
        This method uses the algorithm due to Cheriton and Tarjan, based on the
        exposition in "Graphs, Algorithms, and Optimization" by Kocay and
        Kreher. The Cheriton-Tarjan algorithm finds a minimum spanning tree in
        O(|E| log log |V|) time, although in practice, it is beaten by simpler
        algorithms such as Prim's. See the book "Data Structures and Network
        Algorithms" by Tarjan for detailed information. For another exposition,
        see the book "Algorithms from P to NP" by Moret and Shapiro.
    """
    vertices = graph.vertices()
    num_vertices = len(vertices)

    X = UnionFind(vertices)
    find = X.find
    union = X.union
    PQ = {}

    tree_edges = []
    num_edges = 0

    for u in vertices:
        triples = []
        for v in graph[u]:
            weight = graph[u][v]
            triples.append((weight, u, v))
        PQ[u] = LeftistHeap(triples)

    while num_edges < num_vertices - 1:
        u = choice(PQ.keys())
        Tu = find(u)

        while True:
            (_, x, y) = PQ[Tu].delete_min()
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

        if u != Tu:
            Tu, Tv = u, Tu

        PQ[Tu].merge(PQ[Tv])

        del PQ[Tv]
        tree_edges.append((min(x, y), max(x, y)))

        num_edges += 1

    total = 0
    for (u, v) in tree_edges:
        total += graph[u][v]

    return total, sorted(tree_edges)


################################################################################
# Algorithms for generating all spanning trees of a graph
################################################################################


def spanning_trees_ordered(graph, k=None):
    def find_exchange(father, included, excluded):
        X = NamedUnionFind(vertices)
        find = X.find
        union = X.union

        (min_weight, e, f) = (inf, None, None)

        for (x, y) in included:
            # Make it so that y = father[x]
            if y != father[x]:
                x, y = y, x
            y = find(y)
            union(x, y, y)

        for edge in excluded:
            mark(edge)

        for (x, y) in edges:
            if (x, y) in marked or (y, x) in marked:
                unmark((x, y))
                unmark((y, x))

            elif father[x] != y and father[y] != x:
                a = find(x)
                ancestors = set()
                while a not in ancestors:
                    ancestors.add(a)
                    a = find(father[a])

                a = find(y)
                while a not in ancestors:
                    a = find(father[a])

                for u in [x, y]:
                    v = find(u)
                    while v != a:
                        fv = father[v]
                        exchange_weight = weight[x, y] - weight[v, fv]
                        if exchange_weight < min_weight:
                            min_weight = exchange_weight
                            e = (v, fv)
                            f = (x, y)
                        w = find(fv)
                        union(v, w, w)
                        v = w

        return (min_weight, e, f)

    inf = float('inf')
    if k is None:
        k = inf

    weight = {(u, v): graph[u][v] for u in graph for v in graph[u]}

    marked = set()
    mark = marked.add
    unmark = marked.discard

    edges = sorted(graph.edge_set(), key=lambda (u, v): weight[(u, v)])
    vertices = graph.vertices()

    # arbitrary root vertex
    root = vertices[0]
    tree_weight, father = prim(graph, root, edge_list=False)
    father[root] = root

    (exchange_weight, e, f) = find_exchange(father, [], [])
    heap = [(tree_weight + exchange_weight, e, f, father, [], [])]

    tree_edges = sorted([(min(x, y), max(x, y)) for (x, y) in father.items() if x != y])
    yield tree_weight, tree_edges

    j = 1
    while j < k:
        (tree_weight, e, f, father, included, excluded) = heappop(heap)

        if tree_weight == inf:
            return

        new_graph = Graph()
        new_graph.add_edges(father.items())
        new_graph.delete_edge(e)
        new_graph.add_edge(f)

        new_father = breadth_first_search_tree(new_graph, root)
        new_father[root] = root

        tree_edges = sorted([(min(x, y), max(x, y)) for (x, y) in new_father.items() if x != y])
        yield tree_weight, tree_edges

        new_tree_weight = tree_weight - weight[f] + weight[e]

        included_i = included + [e]
        excluded_j = excluded + [e]

        (exchange_weight, e, f) = find_exchange(father, included_i, excluded)
        heappush(heap, (new_tree_weight + exchange_weight, e, f, father, included_i, excluded))

        (exchange_weight, e, f) = find_exchange(new_father, included, excluded_j)
        heappush(heap, (tree_weight + exchange_weight, e, f, new_father, included, excluded_j))
        j += 1
