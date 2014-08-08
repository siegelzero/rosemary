import rosemary.graphs.algorithms.coloring


################################################################################
# Algorithms for enumerating all maximal cliques.
################################################################################


def bron_kerbosch(graph):
    """
    Returns an iterator over all maximal cliques of graph.

    Input:
        * graph: Graph

    Output:
        * cliques: iterator

    Details:
        The algorithm used in this method is the Bron-Kerbosch algorithm with
        pivoting. See Algorithm 2 in the paper "Enumerating All Connected
        Maximal Common Subgraphs in Two Graphs" by I. Koch for reference. We
        unroll the recursion to avoid recursion depth limitations.
    """
    vertices = graph.vertex_set()

    neighbors = {}
    for u in vertices:
        neighbors[u] = graph.neighbors(u)

    stack = [([], vertices, set())]
    pop = stack.pop
    push = stack.append

    while stack:
        (partial, allowed, seen) = pop()

        if not allowed and not seen:
            yield partial
        elif allowed:
            # Find a vertex u adjacent to the most allowed vertices.
            max_deg = -1
            for v in allowed:
                deg = len(neighbors[v] & allowed)
                if deg > max_deg:
                    max_deg = deg
                    u = v

            for v in allowed - neighbors[u]:
                allowed.remove(v)
                push((partial + [v],
                      allowed & neighbors[v],
                      seen & neighbors[v]))
                seen.add(v)


def bron_kerbosch_no_pivots(graph):
    """
    Returns an iterator over all maximal cliques of graph.

    Input:
        * graph: Graph

    Output:
        * cliques: iterator

    Details:
        The algorithm used in this method is the standard Bron-Kerbosch
        algorithm with no pivoting. See Algorithm 1 in the paper "Enumerating
        All Connected Maximal Common Subgraphs in Two Graphs" by I. Koch for
        reference. This algorithm is included mainly for reference, as the
        algorithm with pivoting is typically superior.
    """
    vertices = graph.vertex_set()

    neighbors = {}
    for u in vertices:
        neighbors[u] = graph.neighbors(u)

    def backtrack(used, allowed, forbidden):
        if not allowed and not forbidden:
            yield used
        else:
            while allowed:
                u = allowed.pop()
                new_allowed = allowed & neighbors[u]
                new_forbidden = forbidden & neighbors[u]
                for clique in backtrack(used + [u], new_allowed, new_forbidden):
                    yield clique
                forbidden.add(u)

    return backtrack([], vertices, set())


def bron_kerbosch_no_pivots_binary(graph):
    """
    Returns an iterator over all maximal cliques of graph.

    Input:
        * graph: Graph

    Output:
        * cliques: iterator

    Details:
        The algorithm used in this method is the standard Bron-Kerbosch
        algorithm with no pivoting. See Algorithm 1 in the paper "Enumerating
        All Connected Maximal Common Subgraphs in Two Graphs" by I. Koch for
        reference. This algorithm is included mainly for reference, as the
        algorithm with pivoting is typically superior.
    """
    vertex_list = graph.vertices()
    vertices = 2**(len(vertex_list)) - 1

    i2v = {}
    v2i = {}

    idx = 1
    for u in vertex_list:
        v2i[u] = idx
        i2v[idx] = u
        idx *= 2

    neighbors = {}
    for u in vertex_list:
        ui = v2i[u]
        neighbors[ui] = 0
        for v in graph.neighbors(u):
            neighbors[ui] |= v2i[v]

    def backtrack(used, allowed, forbidden):
        if not allowed and not forbidden:
            yield used
        else:
            while allowed:
                u = allowed & (-allowed)
                allowed &= (allowed - 1)
                new_allowed = allowed & neighbors[u]
                new_forbidden = forbidden & neighbors[u]

                for clique in backtrack(used | u, new_allowed, new_forbidden):
                    yield clique

                forbidden |= u

    return backtrack(0, vertices, 0)


def bron_kerbosch_binary(graph):
    """
    Returns an iterator over all maximal cliques of graph.

    Input:
        * graph: Graph

    Output:
        * cliques: iterator

    Details:
        The algorithm used in this method is the standard Bron-Kerbosch
        algorithm with no pivoting. See Algorithm 1 in the paper "Enumerating
        All Connected Maximal Common Subgraphs in Two Graphs" by I. Koch for
        reference. This algorithm is included mainly for reference, as the
        algorithm with pivoting is typically superior.
    """
    vertex_list = graph.vertices()
    vertices = 2**(len(vertex_list)) - 1

    i2v = {}
    v2i = {}

    idx = 1
    for u in vertex_list:
        v2i[u] = idx
        i2v[idx] = u
        idx *= 2

    neighbors = {}
    for u in vertex_list:
        ui = v2i[u]
        neighbors[ui] = 0
        for v in graph.neighbors(u):
            neighbors[ui] |= v2i[v]

    stack = [(0, vertices, 0)]
    pop = stack.pop
    push = stack.append

    while stack:
        (partial, allowed, seen) = pop()

        if not allowed and not seen:
            yield partial
        elif allowed:
            # Find a vertex u adjacent to the most allowed vertices.
            max_deg = -1
            t_allowed = allowed
            while t_allowed:
                v = t_allowed & (-t_allowed)
                k = allowed & neighbors[v]

                deg = 0
                while k:
                    k &= k - 1
                    deg += 1

                if deg > max_deg:
                    max_deg = deg
                    u = v

                t_allowed &= t_allowed - 1

            diff = (allowed ^ neighbors[u]) & allowed
            while diff:
                v = diff & (-diff)
                nbrs = neighbors[v]
                allowed &= (~v)
                push((partial | v, allowed & nbrs, seen & nbrs))
                seen |= v
                diff &= diff - 1


maximal_cliques = bron_kerbosch


################################################################################
# Algorithms for enumerating all maximal cliques.
################################################################################


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
        by Ostergard. Experimental evidence suggests that this algorithm is
        superior to the Ostargard algorithm for graphs with high edge density.
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
        while candidates and size + len(candidates) > max_clique[0]:
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
        Algorithm for the Maximum Clique Problem" by Ostergard. Experimental
        evidence suggests that this algorithm is superior to the Pardalos
        algorithm for graphs with low edge density.
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
        while candidates and size + len(candidates) > max_clique[0]:
            # Find the candidate vertex v_i of least index in our ordering.
            for i in xrange(last_i, num_vertices):
                if vertices[i] in candidates:
                    u = vertices[i]
                    last_i = i + 1
                    break

            if size + largest[i] <= max_clique[0]:
                break

            candidates.discard(u)
            backtrack(candidates & neighbors[u], used + [u], size + 1)

            if found[0]:
                break

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


maximum_clique = ostergard


################################################################################
# Algorithms for maximum weight cliques
################################################################################

def max_weight(graph):
    neighbors = {}
    for v in graph:
        neighbors[v] = graph.neighbors(v)

    num_vertices = graph.num_vertices()
    vertices = sorted(graph.vertex_set(), key=lambda u: sum(neighbors[u]))
    max_clique = [0, []]

    def backtrack(candidates, used, weight):
        if not candidates:
            if weight > max_clique[0]:
                max_clique[0] = weight
                max_clique[1] = used
            return

        last_i = 0
        while candidates and weight + sum(candidates) > max_clique[0]:
            # Find the candidate vertex v_i of least index in our ordering.
            for i in xrange(last_i, num_vertices):
                if vertices[i] in candidates:
                    u = vertices[i]
                    last_i = i + 1
                    break

            if weight + largest[i] <= max_clique[0]:
                break

            candidates.discard(u)
            backtrack(candidates & neighbors[u], used + [u], weight + u)

    # Each S[i] is the set of vertices {v_i, v_{i + 1}, ..., v_{n - 1}}.
    # largest[i] is the size of the maximum clique in S[i].
    S = [0]*num_vertices
    largest = [0]*num_vertices

    for i in xrange(num_vertices):
        S[i] = set(vertices[i:])

    for i in xrange(num_vertices - 1, -1, -1):
        vi = vertices[i]
        backtrack(S[i] & neighbors[vi], [vi], vi)
        largest[i] = max_clique[0]

    return max_clique
