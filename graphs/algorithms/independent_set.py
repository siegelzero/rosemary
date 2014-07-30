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

    def backtrack(used, allowed, size, idx):
        if size > best[0]:
            best[0] = size
            best[1] = used

        last_i = idx
        while allowed:
            if size + len(allowed) <= best[0]:
                return

            for i in xrange(last_i, num_vertices):
                if vertices[i] in allowed:
                    u = vertices[i]
                    last_i = i + 1
                    break

            allowed.discard(u)
            backtrack(used + [u],
                      allowed - neighbors[u],
                      size + 1, i)

    backtrack([], set(vertices), 0, 0)

    return tuple(best)


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
