#!/usr/bin/python2

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


def hcd(graph):
    def pull_colors():
        # get the vertex with highest priority
        (_, v) = max((priority[u], u) for u in vvertices)

        adjacent_colors = set()
        for u in graph[v]:
            adjacent_colors.add(color[u])

        # assign to v the lowest admissible color
        old_color = color[v]
        for c in xrange(1, ub[0] + 1):
            if c not in adjacent_colors:
                color[v] = c
                priority[v] = c
                break

        if old_color == c:
            vvertices.remove(v)

            if len(vvertices) == 0:
                ub[0] = max(color[u] for u in vertices)
                push_colors()

    def push_colors():
        num_pushes[0] += 1
        modified = False
        for v in vertices:
            old_color = color[v]
            adjacent_colors = set()
            for u in graph[v]:
                adjacent_colors.add(color[u])

            for c in xrange(ub[0], old_color - 1, -1):
                if c not in adjacent_colors:
                    color[v] = c
                    priority[v] = 1.0 / c
                    if c != old_color:
                        modified = True
                    break

        if modified is False:
            ub[0] = max(color[u] for u in vertices)
            pop_colors()
        else:
            vvertices.update(vertices)

    def pop_colors():
        for v in vertices:
            priority[v] = color[v]

        for v in vertices:
            if color[v] == 1:
                color[v] = ub[0] + 1
        
        for v in vertices:
            #if color[v] == ub[0] + 1:
            #    priority[v] = 1.0 / color[v]
            priority[v] = 1.0 / color[v]
        
        vvertices.update(vertices)


    # Initialize
    priority = {}
    color = {}
    vertices = graph.vertices()
    vvertices = set(vertices)
    ub = [len(vertices)]
    num_pushes = [0]

    for i in xrange(len(vertices)):
        v = vertices[i]
        color[v] = i + 1
        priority[v] = i + 1


    while True:
        pull_colors()
        if num_pushes[0] > 1000:
            return len(set(color[v] for v in color)), color


def xrlf(graph):
    vertices = graph.vertices()
    R = set(vertices)
    k = 0
    exact_limit = 50
    trial_num = [30]
    cand_num = 50 
    set_lim = [20]

    def find_subset(H, U, W, C, F):
        """
        Returns a subset W' <= W that maximizes the size of the set
        {(u, v) in H : u in W' and v in U - (C u W')}
        """
        best = [0, []]
        set_U = set(U)
        list_W = list(W)

        def backtrack(used, last_index):
            diff = set_U.difference(C.union(used))
            ext = set((u, v) for (u, v, _) in F if (u in used and v in diff) or (v in used and u in diff))
            if len(ext) > best[0]:
                best[0] = len(ext)
                best[1] = used

            for i in xrange(last_index + 1, len(W)):
                v = list_W[i]
                skip = False
                for u in used:
                    if v in H[u]:
                        skip = True
                        break
                
                if skip is not True:
                    backtrack(used + [v], i)

        backtrack([], 0)
        return set(best[1])


    def ind_set(H):
        U = H.vertices()
        F = H.edges()

        if not F:
            return U

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
                v = random.choice(U)
                C = set([v])
                X = set(u for u in U if v in H[u])
            else:
                C = set()
                X = set()

            W = set(U).difference(C.union(X))
            kill = False

            while len(W) > 0:
                if kill is True:
                    break

                if len(W) <= set_lim[0]:
                    WW = find_subset(H, U, W, C, F)

                    C.update(WW)

                    for u in C:
                        for v in C:
                            if u in H[v]:
                                print u, v

                    diff = set(U).difference(C)
                    new = set((u, v) for (u, v, _) in F if (u in WW and v in diff) or (v in WW and u in diff))

                    if len(new) > best:
                        CC = set(C)
                        best = len(new)

                    kill = True
                    continue

                bestdegree = -1
                for i in xrange(cand_num):
                    u = random.choice(list(W))
                    s = set((w, v) for (w, v, _) in F if w == u and v in X)
                    if len(s) > bestdegree:
                        bestdegree = len(s)
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

        for u in i_set:
            for v in i_set:
                if u in graph[v]:
                    print "Error: ({}, {}) adjacent".format(u, v)

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
            L = []
            for u in uncolored:
                u_color = set([])
                u_uncolored = 0
                for v in graph[u]:
                    if v in color_map:
                        u_color.add(color_map[v])
                    else:
                        u_uncolored += 1
                L.append((len(u_color), u_uncolored, u))

            #L.sort(reverse=True)
            #u = L[0][2]
            u = max(L)[2]

            adjacent_colors = set(color_map[v] for v in graph[u] if v in color_map)
            missing_color = set(range(1, best_number + 1)) - adjacent_colors
            if len(adjacent_colors) == best_number - 1:
                new_color_map[u] = list(missing_color)[0]
                return best_number

            for j in xrange(1, num_colors + 1):
                if j not in adjacent_colors:
                    new_color_map[u] = j
                    best_number = color(uncolored - set([u]), new_color_map,
                            best_number, num_colors)

            if num_colors < best_number - 1:
                new_color_map[u] = num_colors + 1
                best_number = color(uncolored - set([u]), new_color_map,
                        best_number, num_colors + 1)

            return best_number

    vertices = graph.vertices()
    best = [len(vertices), []]
    color(set(vertices), {}, len(vertices), 0)

    color_map = best[1]
    for v in color_map:
        color_map[v] -= 1

    return best

