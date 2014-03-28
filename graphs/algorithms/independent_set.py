def maximum_independent_set(graph):
    vertices = graph.vertices()
    num_vertices = len(vertices)
    best = [0, []]

    def backtrack(used, neighbors, size, idx):
        if size > best[0]:
            best[0] = size
            best[1] = list(used)
            print best

        for i in xrange(idx + 1, num_vertices):
            u = vertices[i]

            if size + (num_vertices - i) <= best[0]:
                continue

            if u not in neighbors:
                new_neighbors = neighbors.union(graph[u].keys())
                backtrack(used + [u], new_neighbors, size + 1, i)

    backtrack([], set(), 0, -1)

    return best[1]

