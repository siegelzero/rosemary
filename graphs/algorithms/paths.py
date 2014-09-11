from heapq import heappush, heappop


def dijkstra(graph, s, t=None):
    graph_dict = graph.graph_dict
    estimate = {u: float('inf') for u in graph_dict}
    previous = {u: None for u in graph_dict}

    to_visit = [(0, s)]
    estimate[s] = 0
    visited = set()

    while to_visit:
        (w, u) = heappop(to_visit)

        if u in visited:
            continue
        visited.add(u)

        if u == t:
            break

        for v in graph_dict[u]:
            d = w + graph_dict[u][v]
            if d < estimate[v]:
                estimate[v] = d
                previous[v] = u
            heappush(to_visit, (estimate[v], v))

    return estimate, previous


def boroujerdi(graph, s, bad=None):
    if bad is None:
        bad = []

    estimate = {u: float('inf') for u in graph}
    previous = {u: None for u in graph}

    to_visit = [(0, s)]
    visited = set()
    estimate[s] = 0

    while to_visit:
        (_, vj) = heappop(to_visit)

        if vj in visited:
            continue
        visited.add(vj)

        vi = previous[vj]

        for vk in graph[vj]:
            if (vi, vj, vk) in bad:
                continue

            d = estimate[vj] + graph[vj][vk]
            if d < estimate[vk]:
                estimate[vk] = d
                previous[vk] = vj
            heappush(to_visit, (estimate[vk], vk))

    return estimate, previous
