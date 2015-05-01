import rosemary.graphs.algorithms.independent_sets
import rosemary.graphs.algorithms.cliques


################################################################################
# Algorithms for approximate colorings
################################################################################


def dsatur(graph, classes=True):
    r"""Returns an approximation (upper bound) to the chromatic number of
    the graph, along with an approximate coloring.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph

    classes : bool, optional (default=True)
        If True, the algorithm returns a partition of the vertices into
        color classes. Otherwise, returns a dict `color_map` giving the
        color of each vertex.

    Returns
    -------
    (k, color_classes) : tuple
        The integer `k` is an upper bound to the chromatic number of the
        graph. The list `color_classes` is a partition of the vertices into
        color classes.

    Notes
    -----
    This is the classical DSATUR algorithm of Brelaz [1]. The algorithm
    colors a vertex, then selects the next vertex to be colored from the
    set which has the largest number of different colors assigned to its
    neighbors.  The chosen vertex is colored with the least possible color.
    This process is repeated until all vertices are colored.

    References
    ----------
    .. [1] D. Brelaz, "New Methods to Color the Vertices of a Graph",
    Communications of the ACM, Volume 22, Issue 4, April 1979, 251-256.

    Examples
    --------
    >>> graph = petersen_graph()
    >>> dsatur(graph)
    (3, [[0, 7, 8], [1, 3, 5, 9], [2, 4, 6]])
    >>> dsatur(graph, classes=False)
    (3, {0: 0, 1: 1, 2: 2, 3: 1, 4: 2, 5: 1, 6: 2, 7: 0, 8: 0, 9: 1})
    """
    vertices = graph.vertices()
    uncolored = set(vertices)
    color_map = {}

    # Pick a vertex of maximal degree and color it with color 0.
    u = graph.maximum_degree_vertex()
    color_map[u] = 0
    num_colors = 1
    uncolored.remove(u)

    while uncolored:
        # Choose a vertex of maximal saturation degree, breaking any ties
        # by choosing the vertex with the maximal number of uncolored
        # neighbors.
        saturation_degrees = []
        adjacent_colors = {}
        for u in uncolored:
            adjacent_colors[u] = set()
            num_uncolored_neighbors = 0
            for v in graph[u]:
                if v in color_map:
                    adjacent_colors[u].add(color_map[v])
                else:
                    num_uncolored_neighbors += 1
            num_adjacent_colors = len(adjacent_colors[u])
            saturation_degrees.append((num_adjacent_colors, num_uncolored_neighbors, u))

        (_, _, u) = max(saturation_degrees)
        uncolored.remove(u)
        neighbor_colors = adjacent_colors[u]

        # Color this vertex with the least possible color. If a new color
        # is used, increase the color count.
        for color in xrange(num_colors + 1):
            if color not in neighbor_colors:
                color_map[u] = color
                if color == num_colors:
                    num_colors += 1
                break

    if not classes:
        return (num_colors, color_map)
    else:
        color_classes = [[] for i in xrange(num_colors)]
        for (u, color) in color_map.iteritems():
            color_classes[color].append(u)

        return (num_colors, color_classes)


def greedy_sequential(graph, vertices=None, passes=1, classes=True):
    r"""Returns an approximation (upper bound) to the chromatic number of
    the graph, along with an approximate coloring.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph

    passes : int, optional (default=1)
        Optional number of passes to make. A standard is to make two
        passes; once forwards and once backwards.

    classes : bool, optional (default=True)
        If True, returns a partition of the vertices of graph into color
        classes. Otherwise, returns a color map of the vertices.

    vertices : list, optional (default=None)
        An optional list of the vertices in the order to be colored.

    Returns
    -------
    (k, color_classes) : tuple
        The integer `k` is an upper bound to the chromatic number of the
        graph. The list `color_classes` is a partition of the vertices into
        color classes.

    Notes
    -----
    This is the classical greedy sequential algorithm. The first vertex is
    put into the first color class. Thereafter, each vertex is colored with
    the minimum color such that no conflicts are created. See [2] for more
    details.

    A useful fact is that if we take any permutation of the vertices in
    which the vertices of each color class are adjacent in the permutation,
    then applying the greedy sequential algorithm will produce a coloring
    at least as good as the original. Because of this, it is often useful
    to iteratively apply the greedy sequential algorithm to refine a given
    coloring. See [1] for more details about this.

    References
    ----------
    .. [1] J.C. Culberson, "Iterated Greedy Graph Coloring and the
    Difficulty Landscape", Technical Report, University of Alberta, 1992.

    .. [2] D.S. Johnson, C.R. Aragon, L.A. McGeoch, C. Schevon,
    "Optimization by Simulated Annealing: An Experimental Evaluation; Part
    II, Graph Coloring and Number Partitioning", Operations Research,
    Volume 39, Issue 3, May-June 1991, 378-406.

    Examples
    --------
    >>> graph = petersen_graph()
    >>> greedy_sequential(graph)
    (3, [[0, 2, 8, 9], [4, 6, 7], [1, 3, 5]])
    >>> greedy_sequential(graph, classes=False)
    (3, {0: 0, 1: 2, 2: 0, 3: 2, 4: 1, 5: 2, 6: 1, 7: 1, 8: 0, 9: 0})
    >>> greedy_sequential(graph, passes=2)
    (3, [[1, 3, 5, 9], [4, 6, 7], [0, 2, 8]])
    """
    if vertices is None:
        vertices = graph.vertices(order='degree')[::-1]

    color_map = {}
    v = vertices[0]
    color_map[v] = 0
    num_colors = 1

    for i in xrange(1, len(vertices)):
        v = vertices[i]
        adjacent_colors = set([])
        for u in graph[v]:
            if u in color_map:
                adjacent_colors.add(color_map[u])

        for color in xrange(num_colors + 1):
            if color not in adjacent_colors:
                color_map[v] = color
                if color == num_colors:
                    num_colors += 1
                break

    if classes or passes > 1:
        color_classes = [[] for i in xrange(num_colors)]
        for (u, color) in color_map.iteritems():
            color_classes[color].append(u)

        if classes and passes == 1:
            return (num_colors, color_classes)

        elif passes > 1:
            reversed_vertices = []
            while color_classes:
                reversed_vertices.extend(color_classes.pop())

            return greedy_sequential(graph, vertices=reversed_vertices,
                                     passes=1, classes=classes)

    return (num_colors, color_map)


def maxis(graph, classes=True, color_limit=58, mis_limit=100, use_greedy=False, verbose=False):
    r"""Returns an approximation (upper bound) to the chromatic number of
    graph, along with an approximate coloring.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph

    classes : bool, optional (default=True)
        If True, returns a partition of the vertices of graph into color
        classes. Otherwise, returns a color map of the vertices.

    color_limit : int, optional (default=58)
        Limit on the size of the reduced graph for when a final exact
        coloring is performed.

    mis_limit : int, optional (default=100)
        Limit on the size of the reduced graph for when exact techniques
        are used for extracting maximum independent sets.

    use_greedy : bool, optional (default=False)
        If True, use a greedy algorithm for finding large independent sets.
        Otherwise, use a truncated backtracking approach.

    verbose : bool, optional (default=False)
        If True, extra information is printed to the terminal.

    Returns
    -------
    (k, color_classes) : tuple
        The integer `k` is an upper bound to the chromatic number of the
        graph. The list `color_classes` is a partition of the vertices into
        color classes.

    Notes
    -----
    The algorithm proceeds by repeatedly removing large independent sets
    from the graph, assigning the vertices in each set to a color class. By
    default, we use a truncated backtracking approach to find large
    independent sets, enhanced by some heuristics. Once a sufficient number
    of vertices have been removed, we use an exact coloring method on the
    remaining graph.

    This algorithm roughly follows the MAXIS method outlined in [3]. See
    also [1] and [2] for similar algorithms.

    References
    ----------
    .. [1] B. Bollobas, A. Thomason, "Random Graphs of Small Order", Annals
    of Discrete Mathematics 28, 1985.

    .. [2] J.C. Culberson, "Iterated Greedy Graph Coloring and the
    Difficulty Landscape", Technical Report, University of Alberta, 1992.

    .. [3] D.S. Johnson, C.R. Aragon, L.A. McGeoch, C. Schevon,
    "Optimization by Simulated Annealing: An Experimental Evaluation; Part
    II, Graph Coloring and Number Partitioning", Operations Research,
    Volume 39, Issue 3, May-June 1991, 378-406.

    Examples
    --------
    >>> graph = random_graph(30, 0.5)
    >>> maxis(graph)
    (7, [[7, 12, 18, 25], [5, 16, 17], [2, 22, 26, 27, 28],
    [4, 6, 9, 14, 20], [10, 11, 19, 29],
    [0, 3, 15, 21, 23, 24], [1, 8, 13]])
    >>> maxis(graph, classes=False)
    (7, {0: 5, 1: 6, 2: 2, 3: 5, 4: 3, 5: 1, 6: 3, 7: 0,
    8: 6, 9: 3, 10: 4, 11: 4, 12: 0, 13: 6, 14: 3, 15: 5,
    16: 1, 17: 1, 18: 0, 19: 4, 20: 3, 21: 5, 22: 2, 23: 5,
    24: 5, 25: 0, 26: 2, 27: 2, 28: 2, 29: 4})
    >>> maxis(graph, greedy=True)[0]
    7
    >>> maxis(graph, mis_limit=100)[0]
    7
    >>> maxis(graph, color_limit=20)[0]
    8
    """
    approximate = rosemary.graphs.algorithms.independent_sets.culberson
    exact = rosemary.graphs.algorithms.independent_sets.branch_and_bound
    greedy = rosemary.graphs.algorithms.independent_sets.greedy

    reduced_graph = graph.copy()
    vertices = graph.vertex_set()
    num_vertices = len(vertices)

    if verbose:
        print "Using color_limit={}, mis_limit={}".format(color_limit, mis_limit)
        print "Graph density: {}".format(graph.density())

    color_classes = []
    num_colors = 0

    while num_vertices > color_limit:
        if num_vertices > mis_limit:
            if use_greedy:
                if verbose:
                    print "Finding approximate (greedy) MIS..."
                (size, independent_set) = greedy(reduced_graph)

            else:
                if verbose:
                    print "Finding approximate MIS..."
                (size, independent_set) = approximate(reduced_graph)
        else:
            if verbose:
                print "Finding exact MIS..."
            (size, independent_set) = exact(reduced_graph)

        color_classes.append(independent_set)
        num_colors += 1

        reduced_graph.delete_vertices(independent_set)
        vertices = reduced_graph.vertex_set()
        num_vertices -= size

        if verbose:
            print "Color class of size {} found.".format(size),
            print "Remaining vertices: {}".format(num_vertices)

    if num_vertices:
        if verbose:
            print "Using exact coloring techniques for last {} vertices".format(num_vertices)
        num, colors = exact_coloring(reduced_graph, classes=True)
        num_colors += num
        color_classes += colors

    if classes:
        return num_colors, color_classes
    else:
        color_map = {}
        for (i, color_class) in enumerate(color_classes):
            for u in color_class:
                color_map[u] = i
        return num_colors, color_map


################################################################################
# Algorithms for exact colorings
################################################################################


def branch_and_bound(graph, classes=True):
    r"""Given a graph, this returns the chromatic number of the graph
    together with a vertex coloring.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph

    classes : bool, optional (default=True)
        If True, the algorithm returns a partition of the vertices into
        color classes. Otherwise, returns a dict `color_map` giving the
        color of each vertex.

    Returns
    -------
    (k, color_classes) : tuple
        The integer `k` is the chromatic number of the graph. The list
        `color_classes` is a partition of the vertices into color classes.

    Notes
    -----
    This follows the exhaustive search algorithm CHROM_NUM as outlined in
    [2]. This branch-and-bound algorithm implements most of the standard
    shortcuts in [1], such as vertex reordering and some simple heuristics.
    Due to the exponential nature of the problem, it is really only
    feasible to color graphs with <= 70 vertices with this method.

    References
    ----------
    .. [1] D. Brelaz, "New Methods to Color the Vertices of a Graph",
    Communications of the ACM, Volume 22, Issue 4, April 1979, 251-256.

    .. [2] D.S. Johnson, C.R. Aragon, L.A. McGeoch, C. Schevon,
    "Optimization by Simulated Annealing: An Experimental Evaluation; Part
    II, Graph Coloring and Number Partitioning", Operations Research,
    Volume 39, Issue 3, May-June 1991, 378-406.

    Examples
    --------
    >>> graph = petersen_graph()
    >>> branch_and_bound(graph)
    (3, [[0, 2, 8, 9], [4, 6, 7], [1, 3, 5]])
    >>> branch_and_bound(graph, classes=False)
    (3, {0: 1, 1: 3, 2: 1, 3: 3, 4: 2, 5: 3, 6: 2, 7: 2, 8: 1, 9: 1})
    """
    def color(best_number, num_colors):
        if len(uncolored) == 1:
            u = list(uncolored)[0]
            adjacent_colors = set(color_map[v] for v in neighbors[u] if v in color_map)

            # If there is any color j, 1 <= j <= num_colors, such that no
            # vertex adjacent to u has color j, then return num_colors.
            for j in xrange(1, num_colors + 1):
                if j not in adjacent_colors:
                    color_map[u] = j
                    if num_colors < best[0]:
                        best[0] = num_colors
                        best[1] = color_map.copy()
                    del color_map[u]
                    return num_colors

            if num_colors + 1 < best_number:
                color_map[u] = num_colors + 1
                if num_colors + 1 < best[0]:
                    best[0] = num_colors + 1
                    best[1] = color_map.copy()
                del color_map[u]
                return num_colors + 1

            return best_number
        else:
            # Choose the vertex u adjacent to the maximum number of
            # different colors, breaking ties in favor of vertices adjacent
            # to the maximum number of uncolored vertices.
            triples = []
            for u in uncolored:
                adjacent_colors = set()
                adjacent_colors_add = adjacent_colors.add
                num_uncolored_neighbors = 0
                for v in neighbors[u]:
                    if v in color_map:
                        adjacent_colors_add(color_map[v])
                    else:
                        num_uncolored_neighbors += 1
                num_adjacent_colors = len(adjacent_colors)
                triples.append((num_adjacent_colors, num_uncolored_neighbors, u))

            u = max(triples)[2]
            adjacent_colors = set(color_map[v] for v in neighbors[u] if v in color_map)
            num_adjacent_colors = len(adjacent_colors)

            if num_adjacent_colors == best_number - 1:
                return best_number

            for j in xrange(1, num_colors + 1):
                if j not in adjacent_colors:
                    color_map[u] = j
                    pop(u)
                    best_number = color(best_number, num_colors)
                    push(u)
                    del color_map[u]

            if num_colors < best_number - 1:
                color_map[u] = num_colors + 1
                pop(u)
                best_number = color(best_number, num_colors + 1)
                push(u)
                del color_map[u]

            return best_number

    vertices = graph.vertices()
    best = [len(vertices), []]

    # color_map and uncolored are global data structure used by the
    # recursive method.
    color_map = {}
    uncolored = set(vertices)

    # For speed.
    pop = uncolored.remove
    push = uncolored.add
    neighbors = {}
    for u in graph:
        neighbors[u] = graph.neighbors(u)

    # Main entry to the recursive function.
    color(len(vertices), 0)
    num_colors, color_map = best

    if not classes:
        return (num_colors, color_map)
    else:
        color_classes = [[] for _ in xrange(num_colors)]

        for v in color_map:
            color = color_map[v] - 1
            color_classes[color].append(v)

        return (num_colors, color_classes)


def korman(graph, classes=True):
    r"""Given a graph, this returns the chromatic number of the graph
    together with a vertex coloring.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph

    classes : bool, optional (default=True)
        If True, the algorithm returns a partition of the vertices into
        color classes. Otherwise, returns a dict `color_map` giving the
        color of each vertex.

    Returns
    -------
    (k, color_classes) : tuple
        The integer `k` is the chromatic number of the graph. The list
        `color_classes` is a partition of the vertices into color classes.

    Notes
    -----
    This method follows the algorithm in [3], as outlined in [2]. This
    algorithm is a generalization of the method in [1], with several
    optimizations applied. These optimizations include dynamic reordering
    of the vertices, using an approximate coloring to get an upper bound fo
    the chromatic number, and computing an initial clique to get a lower
    bound. These methods are all outlined in [2], as well as in [4]. We
    don't use any local search methods as in [2].

    For most cases, this algorithm is more efficient than the
    `branch_and_bound` exact algorithm, but we're still limited to graphs
    of about 70 vertices or less.

    References
    ----------
    .. [1] D. Brelaz, "New Methods to Color the Vertices of a Graph",
    Communications of the ACM, Volume 22, Issue 4, April 1979, 251-256.

    .. [2] N. Dubois, D. de Werra, "EPCOT: An Efficient Procedure for
    Coloring Optimally with Tabu Search", Computers and Mathematics with
    Applications, Vol 25, 1993.

    .. [3] S.M. Korman, "The Graph-Colouring Problem" in "Combinatorial
    Optimization" by N. Christofides, A. Mingozzi, P. Toth, and C. Sandi,
    Eds., Wiley, 1979.

    .. [4] M. Kubale, B. Jackowski, "A Generalized Implicit Enumeration
    Algorithm for Graph Coloring", Communications of the ACM, Volume 28,
    Number 4, April 1985.

    Examples
    --------
    >>> graph = petersen_graph()
    >>> korman(graph)
    (3, [[0, 7, 8], [1, 3, 5, 9], [2, 4, 6]])
    >>> korman(graph, classes=False)
    (3, {0: 1, 1: 3, 2: 1, 3: 3, 4: 2, 5: 3, 6: 2, 7: 2, 8: 1, 9: 1})
    """
    vertices = graph.vertices()
    num_vertices = len(vertices)
    graph_dict = graph.graph_dict

    degree = {u: len(graph_dict[u]) for u in vertices}
    color = {u: 0 for u in vertices}
    uncolored = set(vertices)

    best_coloring = None
    best, initial_coloring = dsatur(graph, classes=False)
    clique_size, clique = rosemary.graphs.algorithms.cliques.maximum_clique(graph)

    order = {}
    for (i, u) in enumerate(clique):
        color[u] = i + 1
        order[i] = u
        uncolored.discard(u)

    if clique_size != num_vertices:
        idx = clique_size
        find_next_s = True

        while idx >= clique_size:
            if find_next_s:
                max_adj = 0
                max_deg = 0
                for u in uncolored:
                    adjacent_colors = {color[v] for v in graph_dict[u]}
                    num_adjacent = len(adjacent_colors)

                    if num_adjacent == max_adj:
                        if degree[u] > max_deg:
                            max_deg = degree[u]
                            s = u
                    elif num_adjacent > max_adj:
                        max_adj = num_adjacent
                        max_deg = degree[u]
                        s = u

                order[idx] = s
                find_next_s = False

            adjacent_colors = {color[v] for v in graph_dict[s]}
            for i in xrange(color[s] + 1, num_vertices + 1):
                if i not in adjacent_colors:
                    s_color = i
                    break

            if s_color < best and s_color <= max([color[order[i]] + 1 for i in xrange(idx)]):
                color[s] = s_color
                uncolored.discard(s)
                if idx == num_vertices - 1:
                    best = max([color[u] for u in vertices])
                    best_coloring = color.copy()

                    idx = 0
                    while color[order[idx]] < best:
                        idx += 1

                    for i in xrange(idx, num_vertices):
                        color[order[i]] = 0
                        uncolored.add(order[i])

                    idx -= 1
                    s = order[idx]
                else:
                    idx += 1
                    find_next_s = True
            else:
                color[s] = 0
                uncolored.add(s)
                idx -= 1
                s = order[idx]

    if best_coloring is None:
        # color map from DSATUR has colors indexed starting at 0.
        # We need to fix this to make it consistent with the color map
        # found in the korman algorithm.
        best_coloring = {v: c + 1 for (v, c) in initial_coloring.iteritems()}

    if not classes:
        return (best, best_coloring)
    else:
        color_classes = [[] for _ in xrange(best)]

        for v in best_coloring:
            color_classes[best_coloring[v] - 1].append(v)

        return (best, color_classes)


def exact_coloring(graph, classes=True):
    r"""Given a graph, this returns the chromatic number of the graph
    together with a vertex coloring.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph

    classes : bool, optional (default=True)
        If True, the algorithm returns a partition of the vertices into
        color classes. Otherwise, returns a dict `color_map` giving the
        color of each vertex.

    Returns
    -------
    (k, color_classes) : tuple
        The integer `k` is the chromatic number of the graph. The list
        `color_classes` is a partition of the vertices into color classes.

    Notes
    -----
    Given a graph G with vertex set V, an exact coloring of G is a
    partition of V into a minimum number of of color classes where no two
    vertices u and v can be in the same color class if there is an edge
    between them. The minimum number of color classes for G is called the
    chromatic number of G.

    Determining whether the nodes of a graph can be colored with a given
    number k of colors is known to be NP-complete, so it is unlikely that
    efficient exact algorithms exist for graph coloring.

    References
    ----------
    .. [1] D. Brelaz, "New Methods to Color the Vertices of a Graph",
    Communications of the ACM, Volume 22, Issue 4, April 1979, 251-256.

    .. [2] J.C. Culberson, "Iterated Greedy Graph Coloring and the
    Difficulty Landscape", Technical Report, University of Alberta, 1992.

    .. [3] N. Dubois, D. de Werra, "EPCOT: An Efficient Procedure for
    Coloring Optimally with Tabu Search", Computers and Mathematics with
    Applications, Vol 25, 1993.

    .. [4] D.S. Johnson, C.R. Aragon, L.A. McGeoch, C. Schevon,
    "Optimization by Simulated Annealing: An Experimental Evaluation; Part
    II, Graph Coloring and Number Partitioning", Operations Research,
    Volume 39, Issue 3, May-June 1991, 378-406.

    .. [5] S.M. Korman, "The Graph-Colouring Problem" in "Combinatorial
    Optimization" by N. Christofides, A. Mingozzi, P. Toth, and C. Sandi,
    Eds., Wiley, 1979.

    .. [6] M. Kubale, B. Jackowski, "A Generalized Implicit Enumeration
    Algorithm for Graph Coloring", Communications of the ACM, Volume 28,
    Number 4, April 1985.

    Examples
    --------
    >>> graph = petersen_graph()
    >>> exact_coloring(graph)
    (3, [[0, 7, 8], [1, 3, 5, 9], [2, 4, 6]])
    >>> exact_coloring(graph, classes=False)
    (3, {0: 1, 1: 3, 2: 1, 3: 3, 4: 2, 5: 3, 6: 2, 7: 2, 8: 1, 9: 1})
    """
    return korman(graph, classes)
