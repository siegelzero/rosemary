import math
import random
#import numpy

def euclidean_distance_matrix(points):
    num_points = len(points)
    distance_matrix = [[0]*num_points for _ in xrange(num_points)]
    #distance_matrix = numpy.zeros((num_points, num_points))
    inf = float('inf')

    for i in xrange(num_points):
        distance_matrix[i][i] = inf
        (x1, y1) = points[i]
        for j in xrange(i):
            (x2, y2) = points[j]
            dist = ((x2 - x1)**2 + (y2 - y1)**2)**(0.5)
            distance_matrix[i][j] = dist
            distance_matrix[j][i] = dist

    return distance_matrix

def cost(path, distance_matrix, num_vertices):
    total = distance_matrix[path[-1]][path[0]]
    for i in xrange(num_vertices - 1):
        total += distance_matrix[path[i]][path[i + 1]]
    return total

def depth_first(points, distance_matrix=None):
    if distance_matrix is None:
        distance_matrix = euclidean_distance_matrix(points)

    num_vertices = len(points)
    best = [float('inf'), []]

    def backtrack(path, l, last_choices):
        if l == num_vertices:
            total = cost(path, distance_matrix, num_vertices)
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


def branch_and_bound(points, distance_matrix=None):
    if distance_matrix is None:
        distance_matrix = euclidean_distance_matrix(points)

    num_vertices = len(points)
    vertices = range(num_vertices)
    inf = float('inf')
    best = [inf, []]

    def reduce_matrix(matrix):
        val = 0
        num_vertices = len(matrix)

        for i in xrange(num_vertices):
            m = min(matrix[i])
            for j in xrange(num_vertices):
                matrix[i][j] -= m
            val += m

        for j in xrange(num_vertices):
            m = matrix[0][j]
            for i in xrange(1, num_vertices):
                if matrix[i][j] < m:
                    m = matrix[i][j]
            for i in xrange(num_vertices):
                matrix[i][j] -= m
            val += m

        return val

    def reduce_bound(X):
        m = len(X)
        if m == num_vertices:
            return cost(X, distance_matrix, num_vertices)

        other_vertices = [e for e in vertices if e not in X]
        d = len(other_vertices) + 1
        MM = [[0]*d for _ in xrange(d)]

        MM[0][0] = inf
        j = 1


        for y in other_vertices:
            MM[0][j] = distance_matrix[X[-1]][y]
            j += 1

        i = 1
        for x in other_vertices:
            MM[i][0] = distance_matrix[x][X[0]]
            i += 1

        i = 1
        for x in other_vertices:
            j = 1
            for y in other_vertices:
                MM[i][j] = distance_matrix[x][y]
                j += 1
            i += 1

        ans = reduce_matrix(MM)
        for i in xrange(1, m):
            ans += distance_matrix[X[i - 1]][X[i]]

        return ans

    def backtrack(path, l, last_choices):
        if l == num_vertices:
            total = cost(path, distance_matrix, num_vertices)
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


def farthest_insertion(points, distance_matrix=None, **kwargs):
    if distance_matrix is None:
        distance_matrix = euclidean_distance_matrix(points)

    if 'edges' in kwargs:
        edges_only = True
    else:
        edges_only = False

    used = set([0])
    left = set(range(1, len(points)))

    edges = set([(0, 0)])
    tweight = 0

    dist = {}
    for v in left:
        dist[v] = distance_matrix[0][v]

    while left:
        max_dist, f = max((dist[v], v) for v in left)
        c = {}
        for (u, v) in edges:
            c[u, v] = distance_matrix[u][f] + distance_matrix[f][v]
            if (u, v) != (0, 0):
                c[u, v] -= distance_matrix[u][v]

        min_val, (t, h) = min((c[u, v], (u, v)) for (u, v) in edges)

        edges.update([(t, f), (f, h)])
        edges.discard((t, h))
        used.add(f)
        left.discard(f)
        tweight += c[t, h]

        for x in left:
            dist[x] = min(dist[x], distance_matrix[f][x])

    if edges_only:
        return edges

    next_node = dict(edges)
    path = [0]
    next_hop = next_node[0]
    while next_hop != 0:
        path.append(next_hop)
        next_hop = next_node[next_hop]

    return tweight, path


def random_two_opt(path):
    num_points = len(path)

    a = random.randint(0, num_points - 1)

    while True:
        b = random.randint(0, num_points - 1)
        if a - b in (-1, 0, 1):
            continue
        elif (a, b) in ((0, num_points - 1), (num_points - 1, 0)):
            continue
        else:
            break

    if a > b:
        a, b = b, a

    perm = list(path)
    perm[a:b + 1] = perm[a:b + 1][::-1]
    return perm


def steepest_ascent_two_opt(points, distance_matrix=None):
    if distance_matrix is None:
        distance_matrix = euclidean_distance_matrix(points)

    num_points = len(points)
    best_cost, best_sol = farthest_insertion(points)

    def G(X, i, j):
        xi = X[i]
        xii = X[i + 1] if i + 1 < num_points else X[0]
        xj = X[j]
        xjj = X[j + 1] if j + 1 < num_points else X[0]

        total = distance_matrix[xi][xii]
        total += distance_matrix[xj][xjj]
        total -= distance_matrix[xii][xjj]
        total -= distance_matrix[xi][xj]
        return total

    done = False
    while not done:
        done = True
        g0 = 0
        for i in xrange(num_points):
            for j in xrange(i + 2, num_points):
                g = G(best_sol, i, j)
                if g > g0:
                    g0 = g
                    i0 = i
                    j0 = j
        if g0 > 0.0000001:
            new_best = []
            for i in xrange(i0 + 1):
                new_best.append(best_sol[i])
            for i in xrange(j0, i0, -1):
                new_best.append(best_sol[i])
            for i in xrange(j0 + 1, num_points):
                new_best.append(best_sol[i])

            best_sol = new_best
            best_cost -= g0
            done = False

    return best_cost, best_sol


def two_opt(points, distance_matrix=None):
    if distance_matrix is None:
        distance_matrix = euclidean_distance_matrix(points)

    num_points = len(points)

    #best_cost, best_sol = farthest_insertion(points)

    best_sol = range(num_points)
    #random.shuffle(best_sol)
    best_cost = cost(best_sol, distance_matrix, num_points)

    done = False
    while not done:
        done = True
        g0 = 0
        for i in xrange(num_points):
            xi = best_sol[i]
            xii = best_sol[(i + 1) % num_points]

            for j in xrange(i + 2, num_points):
                xj = best_sol[j]
                xjj = best_sol[(j + 1) % num_points]

                gain = distance_matrix[xi][xii]
                gain += distance_matrix[xj][xjj]
                gain -= distance_matrix[xii][xjj]
                gain -= distance_matrix[xi][xj]

                if gain > g0:
                    g0 = gain
                    i0 = i
                    j0 = j

        if g0 > 0.0000001:
            best_sol[i0 + 1:j0 + 1] = best_sol[i0 + 1:j0 + 1][::-1]
            best_cost -= g0
            done = False

    return best_cost, best_sol

def simulated_annealing_two_opt(points, distance_matrix=None):
    if distance_matrix is None:
        distance_matrix = euclidean_distance_matrix(points)

    num_vertices = len(points)

    inf = float('inf')
    best_value = inf
    old_value = inf
    #best_value, best_sol = farthest_insertion(points)
    best_sol = range(num_vertices)
    random.shuffle(best_sol)
    best_value = cost(best_sol, distance_matrix, num_vertices)
    count = 0


    while True:
        if abs(best_value - old_value) < 0.000001:
            count += 1
            if count == 1000:
                break
        else:
            print best_value
            old_value = best_value
            count = 0

        sol = best_sol
        cmax = 10000
        temp = 1000
        alpha = 0.99

        for i in xrange(cmax):
            if temp < 0.0001:
                break
            k = random.randint(2, num_vertices - 1)
            j = random.randint(0, k - 2)

            new_sol = list(sol)
            new_sol[j:k + 1] = new_sol[j:k + 1][::-1]

            cost_new = cost(new_sol, distance_matrix, num_vertices)
            cost_old = cost(sol, distance_matrix, num_vertices)

            if cost_new < cost_old:
                sol = new_sol
                if cost_new < best_value:
                    best_value = cost_new
                    best_sol = sol
            else:
                r = random.random()
                diff = cost_new - cost_old
                if math.log(r) < diff/temp:
                    sol = new_sol

            temp *= alpha

    return best_value, best_sol
