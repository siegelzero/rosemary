from collections import defaultdict

import rosemary.graphs.algorithms.coloring


################################################################################
# Algorithms for generating all maximal cliques.
################################################################################


def bron_kerbosch(graph):
    """
    Returns an iterator over all maximal cliques of graph.

    Input:
        * graph: Graph

    Output:
        * cliques: iterator
            Each clique is given by a list of vertices of graph.

    Examples:
        >>> G = rosemary.graphs.graphs.random_graph(10, 0.5)
        >>> list(bron_kerbosch(G))
        [[0, 8, 9], [0, 8, 3], [0, 8, 4], [0, 9, 5], [1, 9, 5], [1, 3], [1, 6],
         [1, 7], [2, 5], [2, 6], [2, 7], [4, 8, 6], [7, 8]]
        >>> G = rosemary.graphs.graphs.coprime_pairs_graph(10)
        >>> list(bron_kerbosch(G))
        [[1, 2, 9, 5, 7], [1, 2, 3, 5, 7], [1, 3, 8, 5, 7], [1, 3, 10, 7],
         [1, 3, 4, 5, 7], [1, 4, 9, 5, 7], [1, 5, 8, 9, 7],
         [1, 5, 6, 7], [1, 7, 9, 10]]

    Details:
        This function uses the standard Bron-Kerbosch algorithm, with no
        pivoting. See the paper "Enumerating All Connected Maximal Common
        Subgraphs in Two Graphs" by Koch for details.
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


def bron_kerbosch_binary(graph):
    """
    Returns an iterator over all maximal cliques of graph.

    Input:
        * graph: Graph

    Output:
        * cliques: iterator
            Each clique is given by a list of vertices of graph.

    Examples:
        >>> G = random_graph(10, 0.5)
        >>> list(bron_kerbosch_binary(G))
        [[0, 4], [0, 6], [0, 7], [1, 2, 3], [1, 2, 4],
        [1, 3, 6], [2, 4, 8], [3, 5], [4, 5, 8, 9], [6, 9]]
        >>> G = coprime_pairs_graph(10)
        >>> list(bron_kerbosch_binary(G))
        [[1, 2, 3, 5, 7], [1, 2, 5, 7, 9], [1, 3, 4, 5, 7],
        [1, 3, 5, 7, 8], [1, 3, 7, 10], [1, 4, 5, 7, 9],
        [1, 5, 6, 7], [1, 5, 7, 8, 9], [1, 7, 9, 10]]

    Details:
        This function uses the standard Bron-Kerbosch algorithm for enumerating
        all maximal cliques of the graph. This version uses bit operations
        instead of set operations, so it is quite a bit faster than the other BK
        implementation offered, yet still considerably slower than the algorithm
        by Tomita, et al. This version is included mainly for reference.

        One feature of this algorithm is that the cliques are generated in
        lexicographic order.
    """
    vertex_list = graph.vertices()
    vertices = 2**(len(vertex_list)) - 1

    # Dicts for translating between bit form and vertex form.
    int_to_vertex_map = {}
    vertex_to_int_map = {}

    i = 1
    for u in vertex_list:
        vertex_to_int_map[u] = i
        int_to_vertex_map[i] = u
        i <<= 1

    neighbors = {}
    for u in vertex_list:
        i = vertex_to_int_map[u]
        neighbors[i] = 0
        for v in graph.neighbors(u):
            neighbors[i] |= vertex_to_int_map[v]

    def backtrack(used, allowed, forbidden):
        if not allowed and not forbidden:
            clique = []
            while used:
                # u is the least set bit
                u = used & (-used)
                clique.append(int_to_vertex_map[u])
                # clear the least set bit
                used &= used - 1
            yield clique
        else:
            while allowed:
                # least set bit
                u = allowed & (-allowed)
                # clear the least set bit
                allowed &= (allowed - 1)
                new_allowed = allowed & neighbors[u]
                new_forbidden = forbidden & neighbors[u]

                for clique in backtrack(used | u, new_allowed, new_forbidden):
                    yield clique

                forbidden |= u

    return backtrack(0, vertices, 0)


def tomita(graph):
    """
    Returns an iterator over all maximal cliques of graph.

    Input:
        * graph: Graph

    Output:
        * cliques: iterator
            Each clique is given by a list of vertices of graph.

    Examples:
        >>> G = random_graph(10, 0.5)
        >>> list(tomita(G))
        [[7, 0], [6, 9], [6, 0], [4, 8, 2], [4, 8, 9, 5],
        [4, 1, 2], [4, 0], [3, 5], [3, 1, 6], [3, 1, 2]]
        >>> G = coprime_pairs_graph(10)
        [[1, 7, 5, 6], [1, 7, 5, 3, 4], [1, 7, 5, 3, 2],
        [1, 7, 5, 3, 8], [1, 7, 5, 9, 4], [1, 7, 5, 9, 2],
        [1, 7, 5, 9, 8], [1, 7, 10, 3], [1, 7, 10, 9]]

    Details:
        This function uses a variation of the Bron-Kerbosch algorithm, as
        outlined in the paper "The Worst-Case Time Complexity for Generating All
        Maximal Cliques and Computational Experiments" by Tomita, Tanaka, and
        Takahashi. The algorithm has worst-case time complexity O(3^(n/3)) for
        an n-vertex graph, which is optimal as a function of n, since there
        exist n-vertex graphs with as many maximal cliques.
    """
    vertices = graph.vertex_set()

    neighbors = {}
    for u in vertices:
        neighbors[u] = graph.neighbors(u)

    stack = [(set(vertices), set(vertices), [])]
    push = stack.append
    pop = stack.pop

    while stack:
        (subgraph, candidates, used) = pop()

        if not subgraph:
            yield used
        else:
            # Find a vertex u adjacent to the most candidate vertices.
            max_deg = -1
            for v in subgraph:
                deg = len(neighbors[v] & candidates)
                if deg > max_deg:
                    max_deg = deg
                    u = v

            for q in candidates - neighbors[u]:
                nbrs = neighbors[q]
                push((subgraph & nbrs, candidates & nbrs, used + [q]))
                candidates.discard(q)


def maximal_cliques(graph, algorithm='tomita'):
    """
    Returns an iterator over all maximal cliques of graph.

    In a simple, undirected graph G, a clique is a complete subgraph of G; i.e.
    a subgraph in which any two vertices are adjacent. A clique is maximal if it
    is not contained in any larger cliques.

    Input:
        * graph: Graph

        * algorithm: string (default='tomita')
            The algorithm to use for finding the maximal cliques of graph.
            Currently, the supported algorithms are 'bron_kerbosch' for the
            Bron-Kerbosch algorithm, 'bron_kerbosch_binary' for the
            Bron-Kerbosch algorithm implemented using bitsets, and 'tomita' for
            the algorithm due to Tomita, et al.

    Output:
        * cliques: iterator
            Each clique is given by a list of vertices of graph.

    Examples:
        >>> G = random_graph(10, 0.5)
        >>> maximal_cliques(G)
        [[7, 0], [6, 9], [6, 0], [4, 8, 2], [4, 8, 9, 5],
        [4, 1, 2], [4, 0], [3, 5], [3, 1, 6], [3, 1, 2]]
        >>> G = coprime_pairs_graph(10)
        [[1, 7, 5, 6], [1, 7, 5, 3, 4], [1, 7, 5, 3, 2],
        [1, 7, 5, 3, 8], [1, 7, 5, 9, 4], [1, 7, 5, 9, 2],
        [1, 7, 5, 9, 8], [1, 7, 10, 3], [1, 7, 10, 9]]
    """
    if algorithm == 'tomita':
        return tomita(graph)
    if algorithm == 'bron_kerbosch':
        return bron_kerbosch(graph)
    elif algorithm == 'bron_kerbosch_binary':
        return bron_kerbosch_binary(graph)


################################################################################
# Algorithms for finding maximum cliques.
################################################################################


def pardalos(graph):
    """
    Returns a maximum clique of graph.

    Input:
        * graph: Graph

    Ouput:
        * (size, clique): tuple
            * size: int
                Size of the maximum clique.

            * clique: list
                Vertices of a maximum clique of graph.

    Examples:
        >>> G = coprime_pairs_graph(10)
        >>> pardalos(G)
        (5, [1, 2, 3, 5, 7])
        >>> G = random_graph(30, 0.5)
        (6, [0, 11, 20, 23, 26, 28])

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

    used = []
    push = used.append
    pop = used.pop

    def backtrack(candidates, size):
        if not candidates:
            if size > max_clique[0]:
                max_clique[0] = size
                max_clique[1] = list(used)
            return

        last_i = 0
        num_candidates = len(candidates)
        while candidates and size + num_candidates > max_clique[0]:
            for i in xrange(last_i, num_vertices):
                if vertices[i] in candidates:
                    u = vertices[i]
                    last_i = i + 1
                    break

            candidates.discard(u)
            num_candidates -= 1
            push(u)
            backtrack(candidates & neighbors[u], size + 1)
            pop()

    backtrack(set(vertices), 0)
    return tuple(max_clique)


def pardalos_binary(graph):
    """
    Returns a maximum clique of graph.

    Input:
        * graph: Graph

    Ouput:
        * (size, clique): tuple
            * size: int
                Size of the maximum clique.

            * clique: list
                Vertices of a maximum clique of graph.

    Examples:
        >>> G = coprime_pairs_graph(10)
        >>> pardalos(G)
        (5, [1, 2, 3, 5, 7])
        >>> G = random_graph(30, 0.5)
        (6, [0, 11, 20, 23, 26, 28])

    Details:
        This method is based on the algorithm by Carraghan and Pardalos from the
        paper "An Exact Algorithm for the Maximum Clique Problem". We follow the
        outline in the paper "A Fast Algorithm for the Maximum Clique Problem"
        by Ostergard.
    """
    vertices = graph.vertices(order='induced')
    num_vertices = len(vertices)
    all_vertices = 2**num_vertices - 1

    # Dicts for translating between bit form and vertex form.
    int_to_vertex_map = {}
    vertex_to_int_map = {}

    i = 1
    for u in vertices:
        vertex_to_int_map[u] = i
        int_to_vertex_map[i] = u
        i <<= 1

    neighbors = {}
    for u in vertices:
        i = vertex_to_int_map[u]
        neighbors[i] = 0
        for v in graph.neighbors(u):
            neighbors[i] |= vertex_to_int_map[v]

    max_clique = [0, 0]

    def backtrack(used, candidates, size):
        if not candidates:
            if size > max_clique[0]:
                max_clique[0] = size
                max_clique[1] = used
            return

        num_candidates = bin(candidates).count('1')
        while candidates and size + num_candidates > max_clique[0]:
            u = candidates & (-candidates)
            candidates &= (candidates - 1)
            num_candidates -= 1
            backtrack(used | u, candidates & neighbors[u], size + 1)

    backtrack(0, all_vertices, 0)

    bits = max_clique[1]
    clique = []
    for i in xrange(num_vertices):
        if bits & (1 << i):
            clique.append(vertices[i])

    return max_clique[0], clique


def ostergard(graph):
    """
    Returns a maximum clique of graph.

    Input:
        * graph: Graph

    Ouput:
        * (size, clique): tuple
            * size: int
                Size of the maximum clique.

            * clique: list
                Vertices of a maximum clique of graph.
    Examples:
        >>> G = coprime_pairs_graph(10)
        >>> ostergard(G)
        (5, [1, 5, 7, 9, 8])
        >>> G = random_graph(30, 0.5)
        >>> ostergard(G)
        (6, [0, 11, 20, 23, 26, 28])

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
        vertices.extend(sorted(color_class, key=graph.degree, reverse=True))

    vertices = vertices[::-1]
    num_vertices = len(vertices)
    max_clique = [0, []]

    used = []
    push = used.append
    pop = used.pop

    def backtrack(candidates, size):
        if not candidates:
            if size > max_clique[0]:
                max_clique[0] = size
                max_clique[1] = list(used)
                found[0] = True
            return

        last_i = 0
        num_candidates = len(candidates)
        while candidates and size + num_candidates > max_clique[0]:
            # Find the candidate vertex v_i of least index in our ordering.
            for i in xrange(last_i, num_vertices):
                if vertices[i] in candidates:
                    u = vertices[i]
                    last_i = i + 1
                    break

            if size + largest[i] <= max_clique[0]:
                return

            candidates.discard(u)
            num_candidates -= 1
            push(u)
            backtrack(candidates & neighbors[u], size + 1)
            pop()

            if found[0]:
                return

    # largest[i] is the size of the maximum clique in S[i].
    largest = [0]*num_vertices

    # S[i] is the set of vertices {v_i, v_{i + 1}, ..., v_{n - 1}}.
    S = [0]*num_vertices
    for i in xrange(num_vertices):
        S[i] = set(vertices[i:])

    for i in xrange(num_vertices - 1, -1, -1):
        found = [False]
        u = vertices[i]
        push(u)
        backtrack(S[i] & neighbors[u], 1)
        pop()
        largest[i] = max_clique[0]

    return tuple(max_clique)


def ostergard_binary(graph):
    """
    Returns a maximum clique of graph.

    Input:
        * graph: Graph

    Ouput:
        * (size, clique): tuple
            * size: int
                Size of the maximum clique.

            * clique: list
                Vertices of a maximum clique of graph.
    Examples:
        >>> G = coprime_pairs_graph(10)
        >>> ostergard(G)
        (5, [1, 5, 7, 9, 8])
        >>> G = random_graph(30, 0.5)
        >>> ostergard(G)
        (6, [0, 11, 20, 23, 26, 28])

    Details:
        This function follows the algorithm outlined in the paper "A Fast
        Algorithm for the Maximum Clique Problem" by Ostergard. Experimental
        evidence suggests that this algorithm is superior to the Pardalos
        algorithm for graphs with low edge density.
    """
    # We order the vertices by color class. In each class, the vertices are
    # ordered by their degree in the graph.
    vertices = []
    (num_colors, coloring) = rosemary.graphs.algorithms.coloring.greedy_sequential(graph)
    for color_class in coloring:
        vertices.extend(sorted(color_class, key=graph.degree, reverse=True))

    vertices = vertices[::-1]
    num_vertices = len(vertices)
    all_vertices = 2**num_vertices - 1

    # Dicts for translating between bit form and vertex form.
    int_to_vertex_map = {}
    vertex_to_int_map = {}

    i = 1
    for u in vertices:
        vertex_to_int_map[u] = i
        int_to_vertex_map[i] = u
        i <<= 1

    neighbors = {}
    for u in vertices:
        i = vertex_to_int_map[u]
        neighbors[i] = 0
        for v in graph.neighbors(u):
            neighbors[i] |= vertex_to_int_map[v]

    max_clique = [0, 0]

    def backtrack(used, candidates, size):
        if not candidates:
            if size > max_clique[0]:
                max_clique[0] = size
                max_clique[1] = used
                found[0] = True
            return

        num_candidates = bin(candidates).count('1')
        while candidates and size + num_candidates > max_clique[0]:
            u = candidates & (-candidates)

            if size + largest[u] <= max_clique[0]:
                return

            candidates &= (candidates - 1)
            num_candidates -= 1
            backtrack(used | u, candidates & neighbors[u], size + 1)

            if found[0]:
                return

    # largest[i] is the size of the maximum clique in S[i].
    largest = defaultdict(int)

    # S[i] is the set of vertices {v_i, v_{i + 1}, ..., v_{n - 1}}.
    S = [0]*num_vertices
    for i in xrange(num_vertices):
        S[i] = (all_vertices >> i) << i

    for i in xrange(num_vertices - 1, -1, -1):
        found = [False]
        u = 1 << i
        backtrack(u, S[i] & neighbors[u], 1)
        largest[u] = max_clique[0]

    bits = max_clique[1]
    clique = []
    for i in xrange(num_vertices):
        if bits & (1 << i):
            clique.append(vertices[i])

    return max_clique[0], clique


def maximum_clique(graph, algorithm='ostergard'):
    """
    Returns a maximum clique of graph.

    In a simple, undirected graph G, a clique is a complete subgraph of G; i.e.
    a subgraph in which any two vertices are adjacent. A clique of G is a
    maximum clique if there are no larger cliques in G.

    Input:
        * graph: Graph

        * algorithm: string (default='ostergard')
            The algorithm to use for finding the maximum clique. The supported
            values are 'ostergard' for the Ostergard algorithm, and 'pardalos'
            for the Pardalos algorithm.

    Ouput:
        * (size, clique): tuple
            * size: int
                Size of the maximum clique.

            * clique: list
                Vertices of a maximum clique of graph.

    Examples:
        >>> G = coprime_pairs_graph(10, algorithm='ostergard')
        >>> maximum_clique(G)
        (5, [1, 5, 7, 9, 8])
        >>> G = random_graph(30, 0.5)
        >>> maximum_clique(G, algorithm='pardalos')
        (6, [0, 11, 20, 23, 26, 28])
    """
    if algorithm == 'ostergard':
        return ostergard(graph)
    elif algorithm == 'pardalos':
        return pardalos(graph)


################################################################################
# Algorithms for maximum weight cliques
################################################################################


def maximum_weight_clique(graph, weight_map=None):
    """
    Returns a maximum weight clique of graph.

    Input:
        * graph: Graph

        * weight_map: dict (default=None)
            A map giving the weight for each vertex of graph. If weight_map is
            None, and graph has positive integer keys for vertices, then each
            vertex's key is used as its weight.

    Ouput:
        * (weight, clique): tuple
            * weight: int
                Weight of the maximum clique.

            * clique: list
                Vertices of a maximum weight clique of graph.

    Examples:
        >>> G = coprime_pairs_graph(30)
        >>> maximum_weight_clique(G)
        (193, [28, 27, 25, 13, 11, 29, 23, 19, 17, 1])
        >>> G = coprime_pairs_graph(500)
        >>> maximum_weight_clique(G)[0]
        5105

    Details:
        The algorithm used in this function is from the paper "A New Algorithm
        for the Maximum-Weight Clique Problem" by Ostergard.
    """
    # Use the key for each vertex as its weight if weight_map is None.
    if weight_map is None:
        weight_map = {}
        for u in graph:
            weight_map[u] = u

    neighbors = {}
    for v in graph:
        neighbors[v] = graph.neighbors(v)

    num_vertices = graph.num_vertices()
    vertices = sorted(graph.vertex_set(),
                      key=lambda u: sum(weight_map[v] for v in neighbors[u]))

    max_clique = [0, []]

    def backtrack(candidates, used, weight):
        if not candidates:
            if weight > max_clique[0]:
                max_clique[0] = weight
                max_clique[1] = used
            return

        last_i = 0
        candidates_weight = sum(weight_map[u] for u in candidates)

        while candidates:
            if weight + candidates_weight <= max_clique[0]:
                return

            # Find the candidate vertex v_i of least index in our ordering.
            for i in xrange(last_i, num_vertices):
                if vertices[i] in candidates:
                    u = vertices[i]
                    last_i = i + 1
                    break

            if weight + largest[i] <= max_clique[0]:
                return

            candidates.discard(u)
            candidates_weight -= weight_map[u]

            backtrack(candidates & neighbors[u],
                      used + [u],
                      weight + weight_map[u])

    # Each S[i] is the set of vertices {v_i, v_{i + 1}, ..., v_{n - 1}}.
    S = [0]*num_vertices
    for i in xrange(num_vertices):
        S[i] = set(vertices[i:])

    # largest[i] is the weight of the maximum weight clique in S[i].
    largest = [0]*num_vertices

    for i in xrange(num_vertices - 1, -1, -1):
        u = vertices[i]
        backtrack(S[i] & neighbors[u], [u], weight_map[u])
        largest[i] = max_clique[0]

    return tuple(max_clique)
