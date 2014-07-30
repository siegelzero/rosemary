import rosemary.graphs.algorithms.coloring


def maximal_cliques(graph):
    vertices = graph.vertex_set()

    neighbors = {}
    for v in vertices:
        neighbors[v] = set(graph[v].keys())

    greater_vertices = {}
    for v in vertices:
        greater_vertices[v] = {u for u in vertices if u > v}

    def backtrack(partial, last, choices, N, size):
        if size > 0:
            last_neighbors = neighbors[last]
            greater_than_last = greater_vertices[last]
            N = N & last_neighbors
            choices = (choices & last_neighbors) & greater_than_last

        if not N:
            yield partial

        for e in choices:
            for clique in backtrack(partial + [e], e, choices, N, size + 1):
                yield clique

    for clique in backtrack([], 0, vertices, vertices, 0):
        yield clique


def pardalos(graph):
    """
    Finds a maximum clique of graph.

    Input:
        * graph: Graph

    Ouput:
        * (size, clique): tuple
            * size: int
                Size of the maximum clique of graph.
            * clique: list
                Vertices of a maximum clique of graph.

    Details:
        This method is based on the algorithm by Carraghan and Pardalos from the
        paper "An Exact Algorithm for the Maximum Clique Problem". We follow the
        outline in the paper "A Fast Algorithm for the Maximum Clique Problem"
        by Ostergard.
    """
    vertices = graph.vertices(order='induced')
    num_vertices = len(vertices)
    max_clique = [0, []]

    neighbors = {}
    for u in graph:
        neighbors[u] = graph.neighbors(u)

    def backtrack(candidates, used, size):
        if not candidates:
            if size > max_clique[0]:
                max_clique[0] = size
                max_clique[1] = used
            return

        last_i = 0
        while candidates:
            if size + len(candidates) <= max_clique[0]:
                return

            for i in xrange(last_i, num_vertices):
                if vertices[i] in candidates:
                    u = vertices[i]
                    last_i = i + 1
                    break

            candidates.discard(u)
            backtrack(candidates & neighbors[u], used + [u], size + 1)

    backtrack(set(vertices), [], 0)
    return max_clique


def ostergard(graph):
    """
    Finds a maximum clique of graph.

    Input:
        * graph: Graph

    Ouput:
        * (size, clique): tuple
            * size: int
                Size of the maximum clique of graph.
            * clique: list
                Vertices of a maximum clique of graph.

    Details:
        This function follows the algorithm outlined in the paper "A Fast
        Algorithm for the Maximum Clique Problem" by Ostgergard.
    """
    neighbors = {}
    for v in graph:
        neighbors[v] = graph.neighbors(v)

    # We order the vertices by color class. In each class, the vertices are
    # ordered by their degree in the graph.
    vertices = []
    (num_colors, coloring) = rosemary.graphs.algorithms.coloring.greedy_sequential(graph)
    for color_class in coloring:
        vertices.extend(sorted(color_class, key=graph.degree))

    num_vertices = len(vertices)
    max_clique = [0, []]

    def backtrack(candidates, used, size):
        if not candidates:
            if size > max_clique[0]:
                max_clique[0] = size
                max_clique[1] = used
                found[0] = True
            return

        last_i = 0
        while candidates:
            if size + len(candidates) <= max_clique[0]:
                return

            # Find the candidate vertex v_i of least index in our ordering.
            for i in xrange(last_i, num_vertices):
                if vertices[i] in candidates:
                    u = vertices[i]
                    last_i = i + 1
                    break

            if size + largest[i] <= max_clique[0]:
                return

            candidates.discard(u)
            backtrack(candidates & neighbors[u], used + [u], size + 1)

            if found[0]:
                return

    # Each S[i] is the set of vertices {v_i, v_{i + 1}, ..., v_{n - 1}}.
    # largest[i] is the size of the maximum clique in S[i].
    S = [0]*num_vertices
    largest = [0]*num_vertices

    for i in xrange(num_vertices):
        S[i] = set(vertices[i:])

    for i in xrange(num_vertices - 1, -1, -1):
        found = [False]
        vi = vertices[i]
        backtrack(S[i] & neighbors[vi], [vi], 1)
        largest[i] = max_clique[0]

    return max_clique


def kreher(graph):
    vertices = graph.vertex_set()

    neighbors = {}
    for u in vertices:
        neighbors[u] = graph.neighbors(u)

    greater_vertices = {}
    for u in vertices:
        greater_vertices[u] = {v for v in vertices if v > u}

    max_clique = [0, []]

    def backtrack(partial, last, choices, size):
        if size > max_clique[0]:
            max_clique[0] = size
            max_clique[1] = partial

        if size == 0:
            choices = vertices
        else:
            choices = (choices & greater_vertices[last]) & neighbors[last]

        M = size + len(choices)

        for x in choices:
            if M <= max_clique[0]:
                break
            backtrack(partial + [x], x, choices, size + 1)

    backtrack([], 0, vertices, 0)

    return max_clique
