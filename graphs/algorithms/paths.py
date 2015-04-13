import itertools

from collections import deque, defaultdict
from heapq import heappush, heappop

from rosemary.data_structures.heaps import PairingHeap


################################################################################
# Algorithms for the single-source shortest path problem.
################################################################################


def dijkstra(graph, s):
    r"""Returns a shortest path tree of `graph`.

    Given a weighted graph with nonnegative edge weights and vertex `s`,
    this method returns a shortest path tree of graph rooted at `s`.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Weighted graph with nonnegative edge weights.

    s : vertex of graph
        Root vertex of shortest path tree.

    Returns
    -------
    (distance, previous) : tuple
        The dict `distance` contains the length of the shortest path from
        `s` to each vertex of `graph`; i.e. `distance[u]` gives the
        distance from `s` to `u` in `graph`. The dict `previous` gives the
        previous node in the tree; i.e. `previous[u]` gives the parent of
        node `u`.

    Raises
    ------
    ValueError
        If `s` is not a vertex of `graph`.

    See Also
    --------
    dijkstra_bidirectional, dijkstra_buckets, dijkstra_iterator,
    dijkstra_pairing_heap

    Notes
    -----
    This method uses Dijkstra's algorithm to compute the shortest-path
    spanning tree of the graph. Using the Python heapq module, we are able
    to obtain a O(|E|*log(|V|)) runtime.

    Dijkstra's algorithm is very similar to the minimum spanning tree
    algorithm of Prim. The major difference between the two is in the
    objective function: in the shortest path problem we minimize path
    length from the root, while in the minimum spanning tree problem we
    minimize total sum of edge lengths.

    References
    ----------
    .. [1] R. Tarjan, "Data Structures and Network Algorithms", Society for
    Industrial and Applied Mathematics, 1983.

    .. [2] S. Dasgupta, C.H. Papadimitrious, and U. Vazirani, "Algorithms",
    McGraw-Hill, 2008.

    .. [3] K. Mehlhorn, "Data Structures and Algorithms 2: Graph Algorithms
    and NP-Completeness", Springer-Verlag, 1984.

    Examples
    --------
    >>> G = Graph()
    >>> G.add_edges([('a', 'b', 2), ('a', 'c', 1), ('b', 'c', 1),
                     ('b', 'd', 2), ('b', 'e', 3), ('c', 'e', 4),
                     ('d', 'e', 2)])
    >>> dijkstra(G, 'a')
    ({'a': 0, 'b': 2, 'c': 1, 'd': 4, 'e': 5},
     {'a': None, 'b': 'a', 'c': 'a', 'd': 'b', 'e': 'c'})
    """
    graph_dict = graph.graph_dict

    if s not in graph_dict:
        raise ValueError("dijkstra: {} is not a vertex of graph.".format(s))

    inf = float('inf')
    distance = {u: inf for u in graph_dict}
    previous = {u: None for u in graph_dict}

    to_visit = [(0, s)]
    distance[s] = 0
    visited = set()
    add_to_visited = visited.add

    while to_visit:
        (w, u) = heappop(to_visit)

        if u in visited:
            continue
        add_to_visited(u)

        for v in graph_dict[u]:
            new_estimate = w + graph_dict[u][v]
            if new_estimate < distance[v]:
                distance[v] = new_estimate
                previous[v] = u
                # Since Python heaps don't support the decrease_key operation,
                # we insert the new vertex with the new priority into the heap.
                # This doesn't change the asymptotic runtime, but could lead to
                # a larger heap.
                heappush(to_visit, (new_estimate, v))

    return distance, previous


def dijkstra_iterator(graph, s):
    r"""Yields the vertices of `graph` in non-decreasing order of distance
    from the vertex `s`.

    Given a weighted graph with nonnegative edge weights and vertex `s`,
    this method returns a shortest path tree of graph rooted at `s`.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Weighted graph with nonnegative edge weights.

    s : vertex of graph
        Root vertex of shortest path tree.

    Yields
    ------
    (node, distance, path) : tuple
        Each tuple consists of a node name `node`, distance from `s` given
        by `distance`, and the path from `s` to `node`.

    Raises
    ------
    ValueError
        If `s` is not a vertex of `graph`.

    See Also
    --------
    dijkstra, dijkstra_bidirectional, dijkstra_buckets,
    dijkstra_pairing_heap

    Notes
    -----
    This method uses Dijkstra's algorithm to compute the shortest-path
    spanning tree of the graph. Using the Python heapq module, we are able
    to obtain a O(|E|*log(|V|)) runtime.

    References
    ----------
    .. [1] R. Tarjan, "Data Structures and Network Algorithms", Society for
    Industrial and Applied Mathematics, 1983.

    .. [2] R.K. Ahuja, T.L. Magnanti, J.B. Orlin, "Network Flows: Theory,
    Algorithms, and Applications", Prentice-Hall, 1993.

    Examples
    --------
    >>> G = Graph()
    >>> G.add_edges([('a', 'b', 2), ('a', 'c', 1), ('b', 'c', 1),
                     ('b', 'd', 2), ('b', 'e', 3), ('c', 'e', 4),
                     ('d', 'e', 2)])
    >>> X = dijkstra_iterator(G, 'a')
    >>> X.next()
    ('a', 0, ['a'])
    >>> X.next()
    ('c', 1, ['a', 'c'])
    >>> X.next()
    ('b', 2, ['a', 'b'])
    """
    graph_dict = graph.graph_dict

    if s not in graph_dict:
        raise ValueError("dijkstra_iterator: {} is not a vertex of graph.".format(s))

    inf = float('inf')
    estimate = {u: inf for u in graph_dict}
    estimate[s] = 0

    to_visit = [(0, s, [s])]
    visited = set()
    add_to_visited = visited.add

    while to_visit:
        (w, u, path) = heappop(to_visit)

        if u in visited:
            continue
        add_to_visited(u)

        yield u, w, path

        for v in graph_dict[u]:
            new_estimate = w + graph_dict[u][v]
            if new_estimate < estimate[v]:
                estimate[v] = new_estimate
                heappush(to_visit, (new_estimate, v, path + [v]))


def dijkstra_buckets(graph, s):
    r"""Returns a shortest path tree of `graph`.

    Given a weighted graph with nonnegative integral edge weights and
    vertex `s`, this method returns a shortest path tree of graph rooted at
    `s`.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Weighted graph with nonnegative integral edge weights.

    s : vertex of graph
        Root vertex of shortest path tree.

    Returns
    -------
    (distance, previous) : tuple
        The dict `distance` contains the length of the shortest path from
        `s` to each vertex of `graph`; i.e. `distance[u]` gives the
        distance from `s` to `u` in `graph`. The dict `previous` gives the
        previous node in the tree; i.e. `previous[u]` gives the parent of
        node `u`.

    Raises
    ------
    ValueError
        If any edge of the graph has non-integral weight.
        If `s` is not a vertex of `graph`.

    See Also
    --------
    dijkstra

    Notes
    -----
    This method uses Dijkstra's algorithm to compute the shortest-path
    spanning tree of the graph. The implementation is based on Dial's
    implementation, as given in Section 4.6 of [2]. Dial's algorithm
    maintains sets (called buckets), where bucket k stores all nodes with
    distance estimate equal to k. This algorithm runs in O(m + n*C) time,
    where m is the number of edges in the graph, n is the number of nodes,
    and C is the largest edge-length in the network.

    References
    ----------
    .. [1] R. Dial, F. Glover, D. Karney, D. Klingman, "A Computational
    Analysis of Alternative Algorithms and Labeling Techniques for Finding
    Shortest Path Trees", Networks, Vol. 9, 1979.

    .. [2] R.K. Ahuja, T.L. Magnanti, J.B. Orlin, "Network Flows: Theory,
    Algorithms, and Applications", Prentice-Hall, 1993.

    .. [3] R. Tarjan, "Data Structures and Network Algorithms", Society for
    Industrial and Applied Mathematics, 1983.

    Examples
    --------
    >>> G = Graph()
    >>> G.add_edges([('a', 'b', 2), ('a', 'c', 1), ('b', 'c', 1),
                     ('b', 'd', 2), ('b', 'e', 3), ('c', 'e', 4),
                     ('d', 'e', 2)])
    >>> dijkstra_buckets(G, 'a')
    ({'a': 0, 'b': 2, 'c': 1, 'd': 4, 'e': 5},
     {'a': None, 'b': 'a', 'c': 'a', 'd': 'b', 'e': 'c'})
    """
    graph_dict = graph.graph_dict

    if s not in graph_dict:
        raise ValueError("dijkstra_buckets: {} is not a vertex of graph.".format(s))

    inf = float('inf')
    estimate = {v: inf for v in graph_dict}
    estimate[s] = 0

    previous = {s: None}
    buckets = defaultdict(list)
    buckets[0] = [s]

    min_weight = 0
    visited = set()
    add_to_visited = visited.add

    while buckets:
        for w in itertools.count(min_weight):
            if w in buckets:
                min_weight = w
                break

        vertices = buckets[min_weight]

        for u in vertices:
            if u in visited:
                continue
            add_to_visited(u)

            for v in graph_dict[u]:
                # Ensure that all weights are integral.
                if not isinstance(graph_dict[u][v], (int, long)):
                    raise ValueError("dijkstra_buckets: Weights must be integral.")
                d = min_weight + graph_dict[u][v]
                if d < estimate[v]:
                    estimate[v] = d
                    previous[v] = u
                    buckets[d] += [v]

        del buckets[min_weight]

    return estimate, previous


def dijkstra_pairing_heap(graph, s):
    r"""Returns a shortest path tree of `graph`.

    Given a weighted graph with nonnegative edge weights and vertex `s`,
    this method returns a shortest path tree of graph rooted at `s`.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Weighted graph with nonnegative edge weights.

    s : vertex of graph
        Root vertex of shortest path tree.

    Returns
    -------
    (distance, previous) : tuple
        The dict `distance` contains the length of the shortest path from
        `s` to each vertex of `graph`; i.e. `distance[u]` gives the
        distance from `s` to `u` in `graph`. The dict `previous` gives the
        previous node in the tree; i.e. `previous[u]` gives the parent of
        node `u`.

    Raises
    ------
    ValueError
        If `s` is not a vertex of `graph`.

    See Also
    --------
    dijkstra

    Notes
    -----
    This method uses Dijkstra's algorithm to compute the shortest-path
    spanning tree of the graph. This implementation uses pairing heaps
    instead of Python's heapq module. While harder to analyze, pairing
    heaps offer good performance in practice, making them a simpler
    alternative to Fibonacci heaps.

    References
    ----------
    .. [1] R. Tarjan, "Data Structures and Network Algorithms", Society for
    Industrial and Applied Mathematics, 1983.

    .. [2] M.L. Fredman, R. Sedgewick, D.D. Sleator, and R.E. Tarjan, "The
    Pairing Heap: A New Form of Self-Adjusting Heap", Algorithmica, Volume
    1 Issue 1, Jan. 1986.

    .. [3] M.A. Weiss, "Data Structures and Algorithms Analysis in C++",
    Pearson, 2014.

    .. [4] K. Mehlhorn, "Data Structures and Algorithms 2: Graph Algorithms
    and NP-Completeness", Springer-Verlag, 1984.

    Examples
    --------
    >>> G = Graph()
    >>> G.add_edges([('a', 'b', 2), ('a', 'c', 1), ('b', 'c', 1),
                     ('b', 'd', 2), ('b', 'e', 3), ('c', 'e', 4),
                     ('d', 'e', 2)])
    >>> dijkstra_pairing_heap(G, 'a')
    ({'a': 0, 'b': 2, 'c': 1, 'd': 4, 'e': 5},
     {'a': None, 'b': 'a', 'c': 'a', 'd': 'b', 'e': 'c'})
    """
    graph_dict = graph.graph_dict

    if s not in graph_dict:
        raise ValueError("dijkstra_pairing_heap: {} is not a vertex of graph.".format(s))

    inf = float('inf')
    distance = {v: inf for v in graph_dict}
    previous = {s: None}

    distance[s] = 0
    heap = PairingHeap()
    heap.insert(0, s)

    delete_min = heap.delete_min
    decrease_key = heap.decrease_key
    insert = heap.insert
    lookup = heap.lookup

    while True:
        node = delete_min()
        if node is None:
            break

        u = node.value
        w = node.key

        for v in graph_dict[u]:
            d = w + graph_dict[u][v]
            if d < distance[v]:
                distance[v] = d
                previous[v] = u
                if v in lookup:
                    decrease_key(v, d)
                else:
                    insert(d, v)

    return distance, previous


def bellman_ford(graph, s):
    r"""Returns a shortest path tree of `graph`.

    Given a weighted graph (negative weight edges allowed) and vertex `s`,
    this method returns a shortest path tree of `graph` rooted at `s`.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Weighted graph. Negative edges are allowed.

    s : vertex of graph
        Root vertex of shortest path tree.

    Returns
    -------
    (distance, previous) : tuple
        The dict `distance` contains the length of the shortest path from
        `s` to each vertex of `graph`; i.e. `distance[u]` gives the
        distance from `s` to `u` in `graph`. The dict `previous` gives the
        previous node in the tree; i.e. `previous[u]` gives the parent of
        node `u`.

    Raises
    ------
    ValueError
        If there is a negative cycle detected in the graph.
        If `s` is not a vertex of `graph`.

    Notes
    -----
    This method uses the Bellman-Ford algorithm to solve the general
    single-source shortest path problem in time O(|V|*|E|). As opposed to
    Dijkstra's algorithm which uses a shortest-first approach to scanning
    the vertices, this method uses a breadth-first approach.

    While negative edge weights are allowed, the notion of a shortest path
    only makes sense if there are no negative cycles in the graph. By
    counting the number of times that each vertex is removed from the queue,
    we are able to detect the presence of negative cycles, raising a
    ValueError in this case.

    Our implementation follows that given in section IV.7.3 of [2]. See
    also [1] for more details.

    References
    ----------
    .. [1] R. Tarjan, "Data Structures and Network Algorithms", Society for
    Industrial and Applied Mathematics, 1983.

    .. [2] K. Mehlhorn, "Data Structures and Algorithms 2: Graph Algorithms
    and NP-Completeness", Springer-Verlag, 1984.

    Examples
    --------
    >>> G = Graph()
    >>> G.add_edges([('a', 'b', 2), ('a', 'c', 1), ('b', 'c', 1),
                     ('b', 'd', 2), ('b', 'e', 3), ('c', 'e', 4),
                     ('d', 'e', 2)])
    >>> bellman_ford(G, 'a')
    ({'a': 0, 'b': 2, 'c': 1, 'd': 4, 'e': 5},
     {'a': None, 'b': 'a', 'c': 'a', 'd': 'b', 'e': 'c'})
    """
    graph_dict = graph.graph_dict

    if s not in graph_dict:
        raise ValueError("bellman_ford: {} is not a vertex of graph.".format(s))

    n = len(graph_dict)
    inf = float('inf')
    count = {v: 0 for v in graph_dict}
    previous = {v: None for v in graph_dict}

    cost = {v: inf for v in graph_dict}
    cost[s] = 0

    stack = deque([s])
    pop = stack.popleft
    append = stack.append

    stack_elements = set([s])
    discard = stack_elements.discard
    add = stack_elements.add

    while stack:
        u = stack[0]
        count[u] += 1

        if count[u] >= n + 1:
            raise ValueError('bellman_ford: Negative cycle detected.')

        pop()
        discard(u)

        for v in graph_dict[u]:
            if cost[u] + graph_dict[u][v] < cost[v]:
                cost[v] = cost[u] + graph_dict[u][v]
                previous[v] = u
                if v not in stack_elements:
                    append(v)
                    add(v)

    return cost, previous


################################################################################
# Algorithms for the single-pair shortest path problem.
################################################################################


def dijkstra_bidirectional(graph, s, t):
    r"""Returns the shortest path between `s` and `t` in `graph`.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Weighted graph with nonnegative edge weights.

    s : vertex of graph

    t : vertex of graph

    Returns
    -------
    (distance, path) : tuple
        `distance` is the length of a shortest path from `s` to `t` in
        `graph`, and `path` is the sequence of vertices in this shortest
        path.

    Raises
    ------
    ValueError
        If `s` or `t` is not a vertex of `graph`.

    Notes
    -----
    In the bidirectional version of Dijkstra's algorithm, we simultaneously
    apply Dijkstra's algorithm forward from vertex `s` and backward from
    vertex `t`, stopping when both searches have labeled the same vertex.
    At this point, we find the optimal way to combine both paths in to one
    shortest path.

    The runtime of this algorithm is the same as unidirectional
    implementation, but the practical running time of this version is
    typically better. A simple heuristic allows us to estimate this
    speedup: a circle of a certain diameter has twice the area of two
    circles of half the diameter.

    References
    ----------
    .. [1] R.K. Ahuja, T.L. Magnanti, J.B. Orlin, "Network Flows: Theory,
    Algorithms, and Applications", Prentice-Hall, 1993.

    .. [2] R. Tarjan, "Data Structures and Network Algorithms", Society for
    Industrial and Applied Mathematics, 1983.
    """
    inf = float('inf')
    graph_dict = graph.graph_dict

    if s not in graph_dict:
        raise ValueError("dijkstra: {} is not a vertex of graph.".format(s))

    if t not in graph_dict:
        raise ValueError("dijkstra: {} is not a vertex of graph.".format(t))

    # Search forward from `s`.
    forward_search = dijkstra_iterator(graph, s)
    forward_estimate = {}
    forward_path = {}

    # Search backward from`t`.
    backward_search = dijkstra_iterator(graph, t)
    backward_estimate = {}
    backward_path = {}

    # We will alternate taking one step in each direction. If our paths
    # cross, then we terminate. Otherwise, there is no path between the
    # vertices.
    directions = (
        (forward_estimate, backward_estimate, forward_path, forward_search),
        (backward_estimate, forward_estimate, backward_path, backward_search),
    )

    try:
        for estimate, other, path, search in itertools.cycle(directions):
            v, d, p = next(search)
            estimate[v] = d
            path[v] = p

            if v in other:
                break

    except StopIteration:
        return inf

    best_len = inf
    best_path = []

    for u in forward_estimate:
        for v in graph_dict[u]:
            if v not in backward_estimate:
                continue

            d = forward_estimate[u] + graph_dict[u][v] + backward_estimate[v]
            if d < best_len:
                best_len = d
                best_path = forward_path[u] + backward_path[v][1:][::-1] + [t]

    return best_len, best_path
