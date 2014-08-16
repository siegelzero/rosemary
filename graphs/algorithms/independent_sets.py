import random

################################################################################
# Exact algorithms for maximum indepdendent sets
################################################################################


def branch_and_bound(graph):
    """
    Finds a maximum indepdendent set of graph.

    Input:
        * graph: Graph

    Ouput:
        * (size, mis): tuple
            * size: int
                Size of the maximum indepdendent set.
            * mis: list
                Vertices of a maximum indepdendent set.

    Details:
        The algorithm used is a straighforward backtracking branch and bound
        method.
    """
    # Sort the vertices in order of decreasing degree for pruning.
    vertices = sorted(graph.vertices(), key=graph.degree, reverse=True)
    num_vertices = len(vertices)

    neighbors = {}
    for u in graph:
        neighbors[u] = graph.neighbors(u)

    best = [0, []]
    used = []
    push = used.append
    pop = used.pop

    def backtrack(allowed, size, idx):
        if size > best[0]:
            best[0] = size
            best[1] = list(used)

        last_i = idx
        discard = allowed.discard
        while allowed and size + len(allowed) > best[0]:
            for i in xrange(last_i, num_vertices):
                if vertices[i] in allowed:
                    u = vertices[i]
                    last_i = i + 1
                    break

            discard(u)
            push(u)
            backtrack(allowed - neighbors[u],
                      size + 1, i)
            pop()

    backtrack(set(vertices), 0, 0)

    return tuple(best)


def greedy(graph):
    """
    Returns an independent set of graph.

    Input:
        * graph: Graph

    Output:
        * indepdendent_set: list
            A list of the vertices of the independent set.

    Details:
        The algorithm used is the minimum-degree greedy algorithm. This
        algorithm incrementally constructs an independent set by selecting some
        vertex of minimum degree, removing it and its neighbors from the graph,
        and iterating on the remaining graph until empty. This simple algorithm
        achieves a performance ratio of (d + 2) / 3 for approximating maximum
        independent sets in graphs with degree bounded by d.

        See the paper "Greed is Good: Approximating Independent Sets in Sparse
        and Bounded-Degree Graphs" by Halldorsson and Radhakrishnan for details.
    """
    new_graph = graph.copy()
    new_vertex_set = new_graph.vertex_set()
    indepdendent_set = []

    while new_vertex_set:
        # Find the minimum degree in the graph.
        min_degree = len(new_vertex_set)
        degree = {}
        for u in new_vertex_set:
            deg = new_graph.degree(u)
            degree[u] = deg

            if deg < min_degree:
                min_degree = deg

        # Choose a random vertex from the graph with minimum degree, and insert
        # it into our independent set.
        v = random.choice([u for u in new_vertex_set if degree[u] == min_degree])
        indepdendent_set.append(v)

        # Remove this vertex and its neighbors
        neighbors = new_graph.neighbors(v)
        new_graph.delete_vertex(v)
        new_graph.delete_vertices(neighbors)

        new_vertex_set = new_graph.vertex_set()

    indepdendent_set.sort()

    return indepdendent_set


################################################################################
# Algorithms for maximal indepdendent sets
################################################################################


def maximal_independent_sets(graph):
    vertices = graph.vertices()
    num_vertices = len(vertices)

    neighbors = {}
    for u in graph:
        neighbors[u] = graph.neighbors(u)
        neighbors[u].add(u)

    def backtrack(used, allowed, idx):
        if not allowed:
            yield used

        for i in xrange(idx, num_vertices):
            u = vertices[i]

            if u not in allowed:
                continue

            for sol in backtrack(used + [u], allowed - neighbors[u], i + 1):
                yield sol

    for sol in backtrack([], set(vertices), 0):
        yield sol


def culberson(graph):
    def cutoff(n):
        if n >= 500:
            return 8
        elif n >= 200:
            return 6
        elif n >= 10:
            return 3
        else:
            return 10

    degree = {}
    for u in graph:
        degree[u] = graph.degree(u)

    def indset(graph):
        vertices = graph.vertices(order='degree')

        if not vertices:
            return (0, 0, [])

        max_size = 0
        max_deg = 0
        max_set = []

        vertices[0], vertices[-1] = vertices[-1], vertices[0]
        max_iterations = cutoff(len(vertices))

        for (i, v) in enumerate(vertices):
            if i >= max_iterations:
                break

            deg = degree[v]
            subgraph = graph.remove_neighborhood(v)
            (new_size, new_deg, new_set) = indset(subgraph)
            new_deg += deg

            if new_deg > max_deg:
                max_size = new_size + 1
                max_deg = new_deg
                max_set = new_set + [v]

        return (max_size, max_deg, max_set)

    (size, deg, mis) = indset(graph)
    return (size, mis)
