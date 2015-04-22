import rosemary.graphs.algorithms.independent_sets

################################################################################
# Algorithms for approximate colorings
################################################################################


def dsatur(graph, classes=True):
    """
    Returns an approximation (upper bound) to the chromatic number of the graph,
    along with an approximate coloring.

    Input:
        * graph: Graph
            Graph to color.

        * classes: bool (default=True)
            If True, the algorithm returns a partition of the vertices into
            color classes. Otherwise, returns a color map.

    Output:
        * (k, color_classes): tuple
            * k: int
                An upper bound for the chromatic number of the graph.

            * color_classes: list
                A partition of the vertices into color classes.

    Examples:
        >>> graph = petersen_graph()
        >>> dsatur(graph)
        (3, [[0, 7, 8], [1, 3, 5, 9], [2, 4, 6]])
        >>> dsatur(graph, classes=False)
        (3, {0: 0, 1: 1, 2: 2, 3: 1, 4: 2, 5: 1, 6: 2, 7: 0, 8: 0, 9: 1})

    Details:
        This is the classical DSATUR algorithm of Brelaz. The algorithm colors a
        vertex, then selects the next vertex to be colored from the set which
        has the largest number of different colors assigned to its neighbors.
        The chosen vertex is colored with the least possible color.

        See the paper "New Methods to Color the Vertices of a Graph" by Brelaz
        for the original reference.
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
        # Choose a vertex of maximal saturation degree, breaking any ties by
        # choosing the vertex with the maximal number of uncolored neighbors.
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

        # Choose a vertex of maximal saturation degree.
        (_, _, u) = max(saturation_degrees)
        uncolored.remove(u)
        neighbor_colors = adjacent_colors[u]

        # Color this vertex with the least possible color. If a new color is
        # used, increase the color count.
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


def greedy_sequential(graph, **kwargs):
    """
    Returns an approximation (upper bound) to the chromatic number of the graph,
    along with an approximate coloring.

    Input:
        * graph: Graph
            Graph to color.

        * kwargs:
            Possible keywords:
            * vertices: list (default=None)
                An optional ordering of the vertices.

            * passes: int (default=1)
                Optional number of passes to make. A standard is to make two passes;
                once forwards and once backwards.

            * classes: bool (default=True)
                If True, returns a partition of the vertices of graph into color
                classes. Otherwise, returns a color map of the vertices.

    Output:
        * (k, color_classes): tuple
            * k: int
                An upper bound for the chromatic number of the graph.

            * color_classes: list
                A partition of the vertices into color classes.

    Examples:
        >>> graph = petersen_graph()
        >>> greedy_sequential(graph)
        (3, [[0, 2, 8, 9], [4, 6, 7], [1, 3, 5]])
        >>> greedy_sequential(graph, classes=False)
        (3, {0: 0, 1: 2, 2: 0, 3: 2, 4: 1, 5: 2, 6: 1, 7: 1, 8: 0, 9: 0})
        >>> greedy_sequential(graph, passes=2)
        (3, [[1, 3, 5, 9], [4, 6, 7], [0, 2, 8]])

    Details:
        This is the classical greedy sequential algorithm. The first vertex is
        put into the first color class. Thereafter, each vertex is colored with
        the minimum color such that no conflicts are created.

        A useful fact is that if we take any permutation of the vertices in
        which the vertices of each color class are adjacent in the permutation,
        then applying the greedy sequential algorithm will produce a coloring at
        least as good as the original. Because of this, it is often useful to
        iteratively apply the greedy sequential algorithm to refine a given
        coloring. See "Iterated Greedy Graph Coloring and the Difficulty
        Landscape" by Culberson for details about this.
    """
    vertices = kwargs.get('vertices', None)
    classes = kwargs.get('classes', True)
    passes = kwargs.get('passes', 1)

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

    if classes or passes == 2:
        color_classes = [[] for i in xrange(num_colors)]
        for (u, color) in color_map.iteritems():
            color_classes[color].append(u)

        if classes and passes == 1:
            return (num_colors, color_classes)

        elif passes == 2:
            reversed_vertices = []
            while color_classes:
                reversed_vertices.extend(color_classes.pop())

            return greedy_sequential(graph, vertices=reversed_vertices,
                                     passes=1, classes=classes)

    return (num_colors, color_map)


def maxis(graph, **kwargs):
    """
    Returns an approximation (upper bound) to the chromatic number of graph,
    along with an approximate coloring.

    Input:
        * graph: Graph
            Graph to color.

        * kwargs
            Possible keywords:
            * verbose: bool (default=False)
                If True, extra information is printed to the terminal.

            * classes: bool (default=True)
                If True, returns a partition of the vertices of graph into color
                classes. Otherwise, returns a color map of the vertices.

            * color_limit: int (default=58)
                Upper limit on the size of the reduced graph for when a final
                exact coloring is performed.

            * mis_limit: int (default=0)
                Upper limit on the size of the reduced graph for when exact
                techniques are used for extracting maximum independent sets.

            * greedy: bool (default=False)
                If True, use a greedy algorithm for finding large independent
                sets. Otherwise, use a truncated backtracking approach.

    Output:
        * (k, color_classes): tuple
            * k: int
                Approximation to the chromatic number of graph.

            * color_classes: list
                A partition of the vertices of graph into color classes.

    Examples:
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

    Details:
        The algorithm proceeds by repeatedly removing large independent sets
        from the graph, assigning the vertices in each set to a color class. By
        default, we use a truncated backtracking approach to find large
        independent sets, enhanced by some heuristics. Once a sufficient number
        of vertices have been removed, we use an exact coloring method on the
        remaining graph.

        This algorithm roughly follows the MAXIS method outlined in the paper
        "Iterated Greedy Graph Coloring and the Difficulty Landscape" by
        Culberson.
    """
    approximate = rosemary.graphs.algorithms.independent_sets.culberson
    exact = rosemary.graphs.algorithms.independent_sets.branch_and_bound
    greedy = rosemary.graphs.algorithms.independent_sets.greedy

    reduced_graph = graph.copy()
    vertices = graph.vertex_set()
    num_vertices = len(vertices)

    verbose = kwargs.get('verbose', False)
    classes = kwargs.get('classes', True)
    color_limit = kwargs.get('color_limit', 58)
    mis_limit = kwargs.get('mis_limit', 0)
    use_greedy = kwargs.get('greedy', False)

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
        num, colors = branch_and_bound(reduced_graph, classes=True)
    else:
        num, colors = 0, 0

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
    """
    Given a graph, this returns the chromatic number of the graph together with
    a vertex coloring.

    Input:
        * graph: Graph
            Graph to color.

        * classes: bool (default=True)
            If True, the algorithm returns a partition of the vertices into
            color classes. Otherwise, returns a color map.

    Output:
        * (k, color_classes) (tuple)
            * k: int
                The chromatic number of the graph

            * color_classes: list
                A partition of the vertices of graph into color classes.

    Examples:
        >>> graph = petersen_graph()
        >>> branch_and_bound(graph)
        (3, [[0, 2, 8, 9], [4, 6, 7], [1, 3, 5]])
        >>> branch_and_bound(graph, classes=False)
        (3, {0: 1, 1: 3, 2: 1, 3: 3, 4: 2, 5: 3, 6: 2, 7: 2, 8: 1, 9: 1})

    Details:
        This follows the algorithm CHROM_NUM as outlined in the paper
        "Optimization by Simulated Annealing: An Experimental Evaluation; Part
        II, Graph Coloring and Number Partitioning" By Johnson, Aragon, McGeoch,
        and Schevon.
    """
    def color(best_number, num_colors):
        if len(uncolored) == 1:
            u = list(uncolored)[0]
            adjacent_colors = set(color_map[v] for v in neighbors[u] if v in color_map)

            # If there is any color j, 1 <= j <= num_colors, such that no vertex
            # adjacent to u has color j, then return num_colors.
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

            color_map[u] = best_number
            if best_number < best[0]:
                best[0] = best_number
                best[1] = color_map.copy()
            del color_map[u]
            return best_number

        else:
            # Choose the vertex u adjacent to the maximum number of different
            # colors, breaking ties in favor of vertices adjacent to the maximum
            # number of uncolored vertices.
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

    # color_map and uncolored are global data structure used by the recursive
    # method.
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
