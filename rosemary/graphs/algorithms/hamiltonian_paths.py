import itertools
import random


###########################################################################
# Algorithms for Hamiltonian Paths
###########################################################################


def hamiltonian_paths_con(graph):
    """
    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Weighted graph with nonnegative edge weights.
    """
    used = []
    unused = graph.vertex_set()

    def extend(last):
        if not unused:
            yield list(used)
        else:
            if is_connected2(graph, unused):
                for v in graph.neighbors(last) & unused:
                    used.append(v)
                    unused.discard(v)

                    for foo in extend(v):
                        yield foo

                    used.pop()
                    unused.add(v)

    vertices = graph.vertex_set()

    for v in vertices:
        used.append(v)
        unused.discard(v)

        for foo in extend(v):
            yield foo

        used.pop()
        unused.add(v)


def hamiltonian_paths(graph):
    """
    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Weighted graph with nonnegative edge weights.
    """
    used = []
    unused = graph.vertex_set()

    def extend(last):
        if not unused:
            yield list(used)
        else:
            for v in graph.neighbors(last) & unused:
                used.append(v)
                unused.discard(v)

                for foo in extend(v):
                    yield foo

                used.pop()
                unused.add(v)

    vertices = graph.vertex_set()

    for v in vertices:
        used.append(v)
        unused.discard(v)

        for foo in extend(v):
            yield foo

        used.pop()
        unused.add(v)


def number_of_hamiltonian_paths(graph):
    """
    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Weighted graph with nonnegative edge weights.
    """
    def extend(last, left):
        if not left:
            return 1
        else:
            total = 0
            for u in graph.neighbors(last) & left:
                total += extend(u, left - {u})
            return total

    vertices = graph.vertex_set()
    ss = 0
    for v in vertices:
        ss += extend(v, vertices - {v})
    return ss


def number_of_hamiltonian_paths_bits(graph):
    """
    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Weighted graph with nonnegative edge weights.
    """
    vertices = graph.vertices()
    random.shuffle(vertices)
    vertex_to_bit = {}
    for (i, v) in enumerate(vertices):
        vertex_to_bit[v] = 2**i

    neighbors = {}
    for (i, u) in enumerate(vertices):
        neighbors[2**i] = 0
        for v in graph.neighbors(u):
            neighbors[2**i] |= vertex_to_bit[v]

    n = graph.num_vertices()
    all_bits = 2**n - 1

    def extend(bits, cache={}):
        if bits not in cache:
            remaining = bits % (2**n)

            if remaining == 0:
                return 1
            else:
                last = bits >> n
                total = 0
                candidates = remaining & neighbors[last]

                while candidates:
                    u = candidates & (-candidates)
                    total += extend((u << n) + (remaining - u))
                    candidates -= u

                cache[(last << n) + remaining] = total

        return cache[bits]

    count = 0
    for i in xrange(n):
        count += extend(2**(n + i) + (all_bits - 2**i))
    return count


def number_of_hamiltonian_paths_bits_con(graph):
    """
    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Weighted graph with nonnegative edge weights.
    """
    vertices = graph.vertices()
    random.shuffle(vertices)

    vertex_to_bit = {}
    for (i, v) in enumerate(vertices):
        vertex_to_bit[v] = 2**i

    neighbors = {}
    for (i, u) in enumerate(vertices):
        neighbors[2**i] = 0
        for v in graph.neighbors(u):
            neighbors[2**i] |= vertex_to_bit[v]

    n = graph.num_vertices()
    all_bits = 2**n - 1

    def extend(bits, cache={}):
        if bits not in cache:
            remaining = bits % (2**n)

            if remaining == 0:
                return 1
            else:
                # Check connectivity:
                # Traverse from the least remaining bit until all vertices
                # connected to it have been visited.
                visited = 0
                start = remaining & -remaining
                stack = start

                while stack:
                    # Pop the least bit off the stack, and mark as visited.
                    u = stack & -stack
                    stack -= u
                    visited |= u

                    # This is the set difference
                    # (neighbors[u] \cap remaining) \ visited
                    unvisited_neighbors = (neighbors[u] & remaining) & ~visited
                    while unvisited_neighbors:
                        v = unvisited_neighbors & -unvisited_neighbors
                        stack |= v
                        unvisited_neighbors -= v

                # If we haven't visited all remaining vertices, then the
                # graph is not connected.
                if visited != remaining:
                    return 0

                last = bits >> n
                total = 0
                candidates = remaining & neighbors[last]

                while candidates:
                    u = candidates & -candidates
                    total += extend((u << n) + (remaining - u))
                    candidates -= u

                cache[(last << n) + remaining] = total

        return cache[bits]

    count = 0
    for i in xrange(n):
        count += extend(2**(n + i) + (all_bits - 2**i))
    return count


def walks(graph, a, n, allowed):
    """
    number of walks in graph of length n, starting at vertex a, only
    traversing the allowed vertices
    """
    neighbors = {u: graph.neighbors(u) & allowed for u in allowed}

    num_walks = {u: 0 for u in allowed}
    for u in neighbors[a]:
        num_walks[u] += 1

    for i in xrange(n - 1):
        num_walks2 = {u: 0 for u in allowed}
        for u in allowed:
            for v in neighbors[u]:
                num_walks2[u] += num_walks[v]
        num_walks = num_walks2.copy()

    return num_walks


def number_of_hamiltonian_paths5(graph):
    """
    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Weighted graph with nonnegative edge weights.
    """
    total = 0
    vertices = graph.vertex_set()
    n = len(vertices)

    for s in xrange(n + 1):
        for allowed in itertools.combinations(vertices, s):
            for a in allowed:
                num_walks = walks(graph, a, n - 1, set(allowed))
                for u in allowed:
                    total += (-1)**(n - s)*num_walks[u]
    return total


def number_of_hamiltonian_paths5_con(graph):
    """
    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Weighted graph with nonnegative edge weights.
    """
    vertices = graph.vertices()
    n = len(vertices)

    def extend(partial, idx):
        if is_connected2(graph, graph.vertex_set() - set(partial)):
            yield partial

            for i in xrange(idx, n):
                for foo in extend(partial + [vertices[i]], i + 1):
                    yield foo

    total = 0
    for allowed in extend([], 0):
        for a in allowed:
            num_walks = walks(graph, a, n - 1, set(allowed))
            for u in allowed:
                total += (-1)**(n - len(allowed))*num_walks[u]

    return total


def is_connected2(graph, allowed):
    if not allowed:
        return True
    u = list(allowed)[0]
    visited = set()
    stack = [u]
    nbrs = {u: graph.neighbors(u) & allowed for u in allowed}

    while stack:
        u = stack.pop()
        if u in visited:
            continue

        for v in nbrs[u]:
            if v in visited:
                continue
            stack.append(v)

        visited.add(u)

    return visited == allowed
