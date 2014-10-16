from heapq import heappush, heappop
from rosemary.data_structures.heaps import PairingHeap
from collections import deque, defaultdict
from itertools import cycle, count


################################################################################
# Algorithms for the single-source shortest path problem.
################################################################################


def dijkstra(graph, s):
    """
    Returns a shortest path tree of graph.

    Given a weighted graph with nonnegative edge weights and vertex s, this
    method returns a shortest path tree of graph rooted at s.

    Input:
        * graph: Graph
            Weighted graph with nonnegative edge weights.

        * s: vertex of graph
            Root vertex of shortest path tree.

    Output:
        * (distance, previous): tuple
            * distance: dict
                Dict containing the length of the shortest path from s to each
                vertex of graph.

            * previous: dict
                Predecessor dict of the shortest path tree. Use this to trace
                back the path from each vertex to the root.

    Examples:
        >>> G = Graph()
        >>> G.add_edges([('a', 'b', 2), ('a', 'c', 1), ('b', 'c', 1),
                         ('b', 'd', 2), ('b', 'e', 3), ('c', 'e', 4),
                         ('d', 'e', 2)])
        >>> disjkstra(G, 'a')
        ({'a': 0, 'b': 2, 'c': 1, 'd': 4, 'e': 5},
         {'a': None, 'b': 'a', 'c': 'a', 'd': 'b', 'e': 'c'})

    Details:
        This method uses Dijkstra's algorithm to compute the shortest-path
        spanning tree of the graph. Using the Python heapq module, we are able
        to obtain a O(|E|*log(|V|)) runtime.

        Dijkstra's algorithm is very similar to the minimum spanning tree
        algorithm of Prim. The major difference between the two is in the
        objective function: in the shortest path problem minimize path length
        from the root, while in the minimum spanning tree problem we minimize
        total sum of edge lengths.
    """
    graph_dict = graph.graph_dict
    estimate = {u: float('inf') for u in graph_dict}
    previous = {u: None for u in graph_dict}

    to_visit = [(0, s)]
    estimate[s] = 0
    visited = set()
    add_to_visited = visited.add

    while to_visit:
        (w, u) = heappop(to_visit)

        if u in visited:
            continue
        add_to_visited(u)

        for v in graph_dict[u]:
            new_estimate = w + graph_dict[u][v]
            if new_estimate < estimate[v]:
                estimate[v] = new_estimate
                previous[v] = u
                # Since Python heaps don't support the decrease_key operation,
                # we insert the new vertex with the new priority into the heap.
                # This doesn't change the asymptotic runtime, but could lead to
                # a larger heap.
                heappush(to_visit, (new_estimate, v))

    return estimate, previous


def dijkstra2(graph, s):
    """
    Returns a shortest path tree of graph.

    Given a weighted graph with nonnegative edge weights and vertex s, this
    method returns a shortest path tree of graph rooted at s.

    Input:
        * graph: Graph
            Weighted graph with nonnegative edge weights.

        * s: vertex of graph
            Root vertex of shortest path tree.

    Output:
        * (distance, previous): tuple
            * distance: dict
                Dict containing the length of the shortest path from s to each
                vertex of graph.

            * previous: dict
                Predecessor dict of the shortest path tree. Use this to trace
                back the path from each vertex to the root.

    Examples:
        >>> G = Graph()
        >>> G.add_edges([('a', 'b', 2), ('a', 'c', 1), ('b', 'c', 1),
                         ('b', 'd', 2), ('b', 'e', 3), ('c', 'e', 4),
                         ('d', 'e', 2)])
        >>> disjkstra(G, 'a')
        ({'a': 0, 'b': 2, 'c': 1, 'd': 4, 'e': 5},
         {'a': None, 'b': 'a', 'c': 'a', 'd': 'b', 'e': 'c'})

    Details:
        This method uses Dijkstra's algorithm to compute the shortest-path
        spanning tree of the graph. This method differs in that it avoids the
        Python heapq module in favor of a PairingHeap data structure.
    """
    graph_dict = graph.graph_dict
    inf = float('inf')

    estimate = {v: inf for v in graph_dict}
    estimate[s] = 0

    previous = {s: None}
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
            if d < estimate[v]:
                estimate[v] = d
                previous[v] = u
                if v in lookup:
                    decrease_key(v, d)
                else:
                    insert(d, v)

    return estimate, previous


def bellman_ford(graph, s):
    """
    Returns a shortest path tree of graph.

    Given a weighted graph (negative weight edges allowed) and vertex s, this
    method returns a shortest path tree of graph rooted at s.

    Input:
        * graph: Graph
            Weighted graph. Negative edges are allowed.

        * s: vertex of graph
            Root vertex of shortest path tree.

    Output:
        * (distance, previous): tuple
            * distance: dict
                Dict containing the length of the shortest path from s to each
                vertex of graph.

            * previous: dict
                Predecessor dict of the shortest path tree. Use this to trace
                back the path from each vertex to the root.

    Examples:
        >>> G = Graph()
        >>> G.add_edges([('a', 'b', 2), ('a', 'c', 1), ('b', 'c', 1),
                         ('b', 'd', 2), ('b', 'e', 3), ('c', 'e', 4),
                         ('d', 'e', 2)])
        >>> disjkstra(G, 'a')
        ({'a': 0, 'b': 2, 'c': 1, 'd': 4, 'e': 5},
         {'a': None, 'b': 'a', 'c': 'a', 'd': 'b', 'e': 'c'})

    Details:
        This method uses the Bellman-Ford algorithm to solve the general
        single-source shortest path problem in time O(|V|*|E|). As opposed to
        Dijkstra's algorithm which uses a shortest-first approach to scanning
        the vertices, this method uses a breadth-first approach.

        While negative edge weights are allowed, the notion of a shortest path
        only makes sense if there are no negative cycles in the graph. By
        counting the number of times that each vertex is removed from the queue,
        we are able to detect the presence of negative cycles, raising a
        ValueError in this case.

        Our implementation follows that given in section IV.7.3 of the book
        "Data Structures and Algorithms 2 - Graph Algorithms and
        NP-Completeness" by Mehlhorn. See also "Data Structures and Network
        Algorithms" by Tarjan for more details.
    """
    inf = float('inf')
    graph_dict = graph.graph_dict
    n = len(graph_dict)

    cost = {v: inf for v in graph_dict}
    cost[s] = 0

    previous = {v: None for v in graph_dict}
    count = {v: 0 for v in graph_dict}

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
            raise ValueError('bellman_ford: negative cycle detected.')

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


def dijkstra_iterator(graph, s):
    graph_dict = graph.graph_dict

    to_visit = [(0, s, [0])]
    visited = set()
    add_to_visited = visited.add
    inf = float('inf')

    estimate = {u: inf for u in graph_dict}
    estimate[s] = 0

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


def dijkstra_bidirectional(graph, s, t):
    inf = float('inf')

    forward_search = dijkstra_iterator(graph, s)
    forward_estimate = {}
    forward_path = {}

    backward_search = dijkstra_iterator(graph, t)
    backward_estimate = {}
    backward_path = {}

    directions = (
        (forward_estimate, backward_estimate, forward_path, forward_search),
        (backward_estimate, forward_estimate, backward_path, backward_search),
    )

    try:
        for estimate, other, path, search in cycle(directions):
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
        for v in graph[u]:
            if v not in backward_estimate:
                continue

            d = forward_estimate[u] + graph[u][v] + backward_estimate[v]
            if d < best_len:
                best_len = d
                best_path = forward_path[u] + backward_path[v][1:][::-1] + [t]

    return best_len, best_path


def dijkstra_buckets(graph, s):
    buckets = defaultdict(list)
    buckets[0] = [s]

    graph_dict = graph.graph_dict
    inf = float('inf')

    estimate = {v: inf for v in graph_dict}
    estimate[s] = 0

    previous = {s: None}
    min_weight = 0

    visited = set()
    add_to_visited = visited.add

    while buckets:
        for w in count(min_weight):
            if w in buckets:
                min_weight = w
                break

        vertices = buckets[min_weight]

        for u in vertices:
            if u in visited:
                continue
            add_to_visited(u)

            for v in graph_dict[u]:
                d = min_weight + graph_dict[u][v]
                if d < estimate[v]:
                    estimate[v] = d
                    previous[v] = u
                    buckets[d] += [v]

        del buckets[min_weight]

    return estimate, previous
