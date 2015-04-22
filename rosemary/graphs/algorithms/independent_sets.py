import random
import rosemary.graphs.algorithms.cliques


################################################################################
# Exact algorithms for maximum independent sets
################################################################################


def branch_and_bound(graph):
    """
    Finds a maximum independent set of graph.

    Input:
        * graph: Graph

    Output:
        * (size, mis): tuple
            * size: int
                Size of the maximum independent set.
            * mis: list
                Vertices of a maximum independent set.

    Details:
        The algorithm used is a straightforward backtracking branch and bound
        method. We keep the vertices in nonincreasing degree order, so that the
        maximum independent set returned is actually one of maximal total
        degree.
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
            backtrack(allowed - neighbors[u], size + 1, i)
            pop()

    backtrack(set(vertices), 0, 0)

    return tuple(best)


def greedy(graph):
    """
    Returns an independent set of graph.

    Input:
        * graph: Graph

    Output:
        * independent_set: list
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
    independent_set = []

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
        independent_set.append(v)

        # Remove this vertex and its neighbors
        neighbors = new_graph.neighbors(v)
        new_graph.delete_vertex(v)
        new_graph.delete_vertices(neighbors)

        new_vertex_set = new_graph.vertex_set()

    independent_set.sort()

    return len(independent_set), independent_set


def culberson(graph):
    """
    Returns a large independent set of graph.

    Input:
        * graph: Graph

    Output:
        * (size, mis): tuple
            * size: int
                Size of the independent set.
            * mis: list
                Vertices of the independent set.

    Details:
        The algorithms used here is a truncated backtracking algorithm, enhanced
        by a few heuristics. See the paper "Iterated Greedy Graph Coloring and the Difficulty
        Landscape" by Culberson for details.
    """
    def cutoff(n):
        """
        Determines the number of vertices to try at each level of recursion,
        based on the size of the graph.
        """
        if n >= 500:
            return 8
        elif n >= 200:
            return 6
        elif n >= 10:
            return 3
        else:
            return 10

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
            subgraph_vertices = graph.remove_closed_neighborhood(v)
            (new_size, new_deg, new_set) = indset(subgraph_vertices)
            new_deg += deg

            if new_deg > max_deg:
                max_size = new_size + 1
                max_deg = new_deg
                max_set = new_set + [v]

        return (max_size, max_deg, max_set)

    degree = {}
    for u in graph:
        degree[u] = graph.degree(u)

    (size, deg, mis) = indset(graph)
    return (size, mis)


def maximum_independent_set(graph):
    """
    Finds a maximum independent set of graph.

    In a simple, undirected graph G, an independent set is a set of pairwise
    nonadjacent vertices. An independent set is a maximum independent set if
    there are no larger independent sets of G.

    Input:
        * graph: Graph

    Output:
        * (size, mis): tuple
            * size: int
                Size of the maximum independent set.
            * mis: list
                Vertices of a maximum independent set.

    Details:
        To compute the maximum independent set, we find the maximum clique of
        the complement graph. See the function maximum_clique for algorithmic
        details.
    """
    complement = graph.complement()
    mis = rosemary.graphs.algorithms.cliques.maximum_clique(complement)
    return mis


################################################################################
# Algorithms for maximal independent sets
################################################################################


def maximal_independent_sets(graph):
    """
    Returns an iterator over all maximal independent sets of graph.

    In a simple, undirected graph G, an independent set is a set of pairwise
    nonadjacent vertices. An independent set is a maximal independent set if it
    is not contained in any larger independent sets of G.

    Input:
        * graph: Graph

    Output:
        * maximal_independent_sets: iterator
            Each independent set is given by a list of vertices of graph.

    Details:
        To compute the maximal independent sets, we find the maximal cliques of
        the complement graph. See the function maximal_cliques for algorithmic
        details.
    """
    complement = graph.complement()
    return rosemary.graphs.algorithms.cliques.maximal_cliques(complement)
