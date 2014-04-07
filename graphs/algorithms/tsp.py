import rosemary.graphs.graphs
import copy
import itertools
import math
import random

def euclidean_tsp_graph(points):
    graph = rosemary.graphs.graphs.Graph()
    num_points = len(points)

    for i in xrange(num_points):
        (x1, y1) = points[i]
        for j in xrange(i):
            (x2, y2) = points[j]
            dist = ((x2 - x1)**2 + (y2 - y1)**2)**(0.5)
            graph.add_edge(i, j, dist)

    return graph

def cost(path, graph, num_vertices):
    total = 0
    u = path[0]
    for i in xrange(1, num_vertices):
        v = path[i]
        total += graph[u][v]
        u = v
    return total + graph[u][0]

def depth_first(points):
    graph = euclidean_tsp_graph(points)
    vertices = graph.vertices()
    num_vertices = len(vertices)
    best = [float('inf'), []]

    def backtrack(path, l, last_choices):
        if l == num_vertices:
            total = cost(path, graph, num_vertices)
            if total < best[0]:
                best[0] = total
                best[1] = list(path)
        else:
            if l == 0:
                choices = [0]
            elif l == 1:
                choices = range(1, num_vertices)
            else:
                choices = [e for e in last_choices if e != path[-1]]

            for c in choices:
                backtrack(path + [c], l + 1, choices)

    backtrack([], 0, [])
    return best


def heldman_karp(points):
    graph = euclidean_tsp_graph(points)
    num_points = len(points)
    vertices = graph.vertices()
    others = vertices[1:]
    inf = float('inf')

    C = {}
    C[(0,), 0] = (0, (0,))

    for s in xrange(1, num_points):
        for S in itertools.combinations(others, s):
            subset = (0,) + S
            C[subset, 0] = (inf, [])

            for j in S:
                without = tuple(e for e in subset if e != j)
                mm = inf
                for i in subset:
                    if i == j:
                        continue
                    c, p = C[without, i]
                    c += graph[i][j]
                    if c < mm:
                        mm = c
                        path = p + (j,)
                C[subset, j] = (c, path)

    mm = inf
    for j in others:
        (c, p) = C[tuple(vertices), j]
        c += graph[0][j]
        if c < mm:
            mm = c
            path = p

    return (mm, path)


def branch_and_bound(points):
    graph = euclidean_tsp_graph(points)
    num_vertices = len(points)

    vertices = graph.vertices()
    inf = float('inf')
    best = [inf, []]
    cost_matrix = [[0]*num_vertices for _ in xrange(num_vertices)]

    for u in vertices:
        for v in vertices:
            if u == v:
                cost_matrix[u][v] = inf
            else:
                cost_matrix[u][v] = graph[u][v]

    def reduce(matrix):
        val = 0
        M = copy.deepcopy(matrix)
        num_vertices = len(matrix)

        for i in xrange(num_vertices):
            m = M[i][0]
            for j in xrange(1, num_vertices):
                if M[i][j] < m:
                    m = M[i][j]
            for j in xrange(num_vertices):
                M[i][j] -= m
            val += m

        for j in xrange(num_vertices):
            m = M[0][j]
            for i in xrange(1, num_vertices):
                if M[i][j] < m:
                    m = M[i][j]
            for i in xrange(num_vertices):
                M[i][j] -= m
            val += m

        return val

    def reduce_bound(X):
        m = len(X)
        if m == num_vertices:
            return cost(X, graph, num_vertices)

        other_vertices = [e for e in vertices if e not in X]
        d = len(other_vertices) + 1
        MM = [[0]*d for _ in xrange(d)]

        MM[0][0] = inf
        j = 1


        for y in other_vertices:
            MM[0][j] = cost_matrix[X[-1]][y]
            j += 1

        i = 1
        for x in other_vertices:
            MM[i][0] = cost_matrix[x][X[0]]
            i += 1

        i = 1
        for x in other_vertices:
            j = 1
            for y in other_vertices:
                MM[i][j] = cost_matrix[x][y]
                j += 1
            i += 1
        
        ans = reduce(MM)
        for i in xrange(1, m):
            ans += cost_matrix[X[i - 1]][X[i]]

        return ans

    def backtrack(path, l, last_choices):
        if l == num_vertices:
            total = cost(path, graph, num_vertices)
            if total < best[0]:
                best[0] = total
                best[1] = list(path)
        else:
            if l == 0:
                choices = [0]
            elif l == 1:
                choices = range(1, num_vertices)
            else:
                choices = [e for e in last_choices if e != path[-1]]

            L = []
            for x in choices:
                next_choice = x
                next_bound = reduce_bound(path + [x])
                L.append((next_bound, next_choice))

            L.sort()
            for (bound, choice) in L:
                if bound >= best[0]:
                    return
                backtrack(path + [choice], l + 1, choices)


    backtrack([], 0, [])
    return best

def simulated_annealing(points):
    graph = euclidean_tsp_graph(points)
    num_vertices = len(points)

    inf = float('inf')
    best_value = inf
    old_value = inf
    best_value, best_sol = farthest_insertion(points)
    count = 0


    while True:
        if best_value == old_value:
            count += 1
            if count == 10000:
                break
        else:
            print best_value
            old_value = best_value
            count = 0

        sol = best_sol
        cmax = 1000
        temp = 1000
        alpha = 0.99999

        for i in xrange(cmax):
            j = random.randint(2, num_vertices - 1)
            k = random.randint(1, j)

            new_sol = list(sol)
            new_sol[j], new_sol[k] = new_sol[k], new_sol[j]

            cost_new = cost(new_sol, graph, num_vertices)
            cost_old = cost(sol, graph, num_vertices)

            if cost_new < cost_old:
                sol = list(new_sol)
                if cost_new < best_value:
                    best_value = cost_new
                    best_sol = list(sol)
            else:
                r = random.random()
                diff = cost_new - cost_old
                if r < math.exp(diff / temp):
                    sol = list(new_sol)

            temp *= alpha

    return best_value, best_sol


def farthest_insertion(points, graph=None):
    if graph is None:
        graph = euclidean_tsp_graph(points)

    used = set([0])
    left = set(graph.vertices()).difference([0])

    edges = set([(0, 0)])
    tweight = 0
    
    dist = {}
    for v in left:
        dist[v] = graph[0][v]

    while left:
        max_dist, f = max((dist[v], v) for v in left)
        c = {}
        for (u, v) in edges:
            c[u, v] = graph[u][f] + graph[f][v]
            if (u, v) != (0, 0):
                c[u, v] -= graph[u][v]

        min_val, (t, h) = min((c[u, v], (u, v)) for (u, v) in edges)

        edges.update([(t, f), (f, h)])
        edges.discard((t, h))
        used.add(f)
        left.discard(f)
        tweight += c[t, h]

        for x in left:
            dist[x] = min(dist[x], graph[f][x])

    next_node = dict(edges)
    path = [0]
    next_hop = next_node[0]
    while next_hop != 0:
        path.append(next_hop)
        next_hop = next_node[next_hop]

    return tweight, path

