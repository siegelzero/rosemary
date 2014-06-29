import random

################################################################################
# Algorithms for approximate colorings
################################################################################

def dsatur(graph):
    """
    Returns an upper bound to the chromatic number of the graph, along with an
    approximate coloring.

    Input:
        * graph: (Graph)

    Output:
        * (num_colors, color_map)
            * num_colors: (int)
                This is an upper bound for the chromatic number of the graph.
            * color_map: (dict)
                This dict gives color info for the vertices of the graph. The
                keys are vertices, and the value corresponding to each key is
                the color in the approximate coloring.

    Details:
        This is the DSATUR algorithm outlined in "New Methods to Color the
        Vertices of a Graph" by Brelaz. This algorithm is exact for bipartite
        graphs.
    """
    vertices = graph.vertices()
    uncolored = set(vertices)
    color_map = {}

    # Arrange the vertices by decreasing order of degrees.
    degrees = []
    for v in vertices:
        v_degree = len(graph[v])
        degrees.append((v_degree, v))

    # Pick a vertex of maximal degree and color it with color 0.
    (_, u) = max(degrees)
    color_map[u] = 0
    num_colors = 1
    uncolored.remove(u)

    while uncolored:
        # Compute the saturation degree of each uncolored vertex.
        # The saturation degree of a vertex is the number of different colors to
        # which it is adjacent. We will choose a vertex of maximal saturation
        # degree, breaking any ties by choosing the vertex with the maximal
        # number of uncolored neighbors.
        saturation_degrees = []
        adjacent_colors = {}
        for u in uncolored:
            adjacent_colors[u] = set([])
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

    return (num_colors, color_map)


def greedy_sequential(graph):
    vertices = graph.vertices()
    colors = {}

    v = vertices[0]
    colors[v] = 0
    num_colors = 1

    for i in xrange(1, len(vertices)):
        v = vertices[i]
        adjacent_colors = set([])
        for u in graph[v]:
            if u in colors:
                adjacent_colors.add(colors[u])

        for color in xrange(num_colors + 1):
            if color not in adjacent_colors:
                colors[v] = color
                if color == num_colors:
                    num_colors += 1
                break

    return (num_colors, colors)

def xrlf(graph):
    vertices = graph.vertices()
    R = set(vertices)
    k = 0
    exact_limit = 60
    trial_num = [40]
    cand_num = 50 
    set_lim = [90]

    def find_subset(subgraph, vertices, W, C, F):
        """
        Returns a subset W' <= W that maximizes the size of the set
        {(u, v) in H : u in W' and v in U - (C u W')}
        """
        best = [0, set()]
        list_W = list(W)
        len_W = len(list_W)
        #print len_W

        def backtrack(used, neighbors, diff, last_index):
            ext = sum(1 for (u, v, _) in F if (u in used and v in diff) or (v in used and u in diff))

            if ext > best[0]:
                best[0] = ext
                best[1] = used

            for i in xrange(last_index + 1, len_W):
                v = list_W[i]
                if v in neighbors:
                    continue

                t_used = used.union([v])
                t_neighbors = neighbors.union([u for u in subgraph[v]])
                t_diff = diff.difference([v])

                backtrack(t_used, t_neighbors, t_diff, i)

        diff = vertices.difference(C)
        backtrack(set(), set(), diff, 0)
        return set(best[1])


    def ind_set(H):
        U_list = H.vertices()
        U = set(U_list)
        F = H.edges()

        best = -1
        CC = set()
        C0 = set()
        Dmin = min(len(H[v]) for v in U)
        Dmax = max(len(H[v]) for v in U)

        if trial_num == 1 and len(U) > set_lim[0]:
            vmax = [v for v in U if len(H[v]) == Dmax][0]
            C0.add(vmax)

        if min(trial_num[0], set_lim[0] + Dmin) >= len(U):
            set_lim[0] = len(U)
            trial_num[0] = 1

        for k in xrange(trial_num[0]):
            if C0:
                C = set(C0)
                X = set(u for u in U if u in H[vmax])
            elif len(U) > set_lim[0]:
                v = random.choice(U_list)
                C = set([v])
                X = set(u for u in U if v in H[u])
            else:
                C = set()
                X = set()

            W = U.difference(C.union(X))
            kill = False

            while len(W) > 0:
                if kill is True:
                    break

                if len(W) <= set_lim[0]:
                    WW = find_subset(H, U, W, C, F)
                    C.update(WW)

                    diff = U.difference(C)
                    new = sum(1 for (u, v, _) in F if (u in WW and v in diff) or (v in WW and u in diff))
                    if new > best:
                        CC = C
                        best = new

                    kill = True
                    continue

                bestdegree = -1
                candidates = random.sample(W, min(len(W), cand_num))
                for u in candidates:
                    num = sum(1 for (w, v, _) in F if w == u and v in X)
                    if num > bestdegree:
                        bestdegree = num
                        cand = u

                C.add(cand)
                X.update(v for v in W if cand in H[v])
                W = W.difference(X)
                W.discard(cand)

        return CC

    color_map = {}
    while len(R) > exact_limit:
        induced = graph.induced_subgraph(R)
        i_set = ind_set(induced)

        for v in i_set:
            color_map[v] = k

        R = R.difference(i_set)
        k += 1

    left = 0
    if R:
        induced = graph.induced_subgraph(R)
        left, exact = branch_and_bound(induced)

        for v in exact:
            exact[v] += k

        color_map.update(exact)

    return k + left, color_map


################################################################################
# Algorithms for exact colorings
################################################################################

def branch_and_bound(graph):
    """
    Given a graph, this returns the chromatic number of the graph together with
    a vertex coloring.

    Input:
        * graph (Graph)
    
    Output:
        * (k, color_map) (tuple)
            * k: int
                The chromatic number of the graph
            * color_map (dict)
                A dict with vertices as keys and associated colors as
                corresponding values.

    Details:
        This follows the algorithm CHROM_NUM as outlined in the paper
        "Optimization by Simulated Annealing: An Experimental Evaluation; Part
        II, Graph Coloring and Number Partitioning" By Johnson, Aragon, McGeoch,
        and Schevon.
    """
    def color(uncolored, color_map, best_number, num_colors):
        """
        Returns the least number of colors required to color the uncolored
        vertices of the parent graph.

        Input:
            * uncolored: (set)
                A set of the uncolored vertices of the graph.

            * color_map: (dict)
                A dict with the colored vertices of the graph as keys, and the
                corresponding colors as values.

            * best_number: (int)
                The number of colors in the best legal coloring seen so far.

            * num_colors: (int)
                The number of colors used in the coloring in color_map

        Output:
            * best_number:
                The least number of colors required to color the uncolored
                vertices of the graph.
        """
        new_color_map = color_map.copy()
        if len(uncolored) == 1:
            u = uncolored.pop()
            adjacent_colors = set(color_map[v] for v in graph[u] if v in color_map)

            # If there is any color j, 1 <= j <= num_colors, such that no vertex
            # adjacent to u has color j, then return num_colors.
            for j in xrange(1, num_colors + 1):
                if j not in adjacent_colors:
                    new_color_map[u] = j
                    if num_colors < best[0]:
                        best[0] = num_colors
                        best[1] = new_color_map
                    return num_colors

            if num_colors + 1 < best_number:
                new_color_map[u] = num_colors + 1
                if num_colors + 1 < best[0]:
                    best[0] = num_colors + 1
                    best[1] = new_color_map
                return num_colors + 1

            new_color_map[u] = best_number
            if best_number < best[0]:
                best[0] = best_number
                best[1] = new_color_map
            return best_number

        else:
            # Compute the saturation degree of each vertex.
            # We want the vertex u adjacent to colored vertices with the maximum
            # number of different colors. We break ties in favor of vertices
            # that are adjacent to the most as yet uncolored vertices.
            triples = []
            for u in uncolored:
                adjacent_colors = set()
                adjacent_colors_add = adjacent_colors.add
                num_uncolored_neighbors = 0
                for v in graph[u]:
                    if v in color_map:
                        adjacent_colors_add(color_map[v])
                    else:
                        num_uncolored_neighbors += 1
                num_adjacent_colors = len(adjacent_colors)
                triples.append((num_adjacent_colors, num_uncolored_neighbors, u))

            u = max(triples)[2]
            neighbor_colors = set(color_map[v] for v in graph[u] if v in color_map)
            num_adjacent_colors = len(neighbor_colors)

            # if u is adjacent to best_number - 1 colors, return best_number
            if num_adjacent_colors == best_number - 1:
                return best_number

            uncolored_difference = uncolored.difference
            for j in xrange(1, num_colors + 1):
                if j not in neighbor_colors:
                    new_color_map[u] = j
                    best_number = color(uncolored_difference([u]),
                            new_color_map, best_number, num_colors)

            if num_colors < best_number - 1:
                new_color_map[u] = num_colors + 1
                best_number = color(uncolored_difference([u]), new_color_map,
                        best_number, num_colors + 1)

            return best_number

    vertices = graph.vertices()
    best = [len(vertices), []]
    color(set(vertices), {}, len(vertices), 0)

    # The vertices are assigned colors >= 1, so we normalize to colors >= 0.
    color_map = best[1]
    for v in color_map:
        color_map[v] -= 1

    return best
