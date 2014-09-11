from heapq import heappush, heappop
from rosemary.data_structures.heaps import PairingHeap


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
                heappush(to_visit, (d, v))

    path = []
    while t is not None:
        path.append(t)
        t = previous[t]

    return path[::-1]


def dijkstra2(graph, s, t=None):
    graph_dict = graph.graph_dict
    inf = float('inf')

    estimate = {s: 0}
    previous = {s: None}
    heap = PairingHeap()

    for v in graph_dict:
        if v != s:
            estimate[v] = inf
            previous[v] = None
        heap.insert(estimate[v], v)

    delete_min = heap.delete_min
    decrease_key = heap.decrease_key

    while True:
        node = delete_min()
        if node is None:
            break

        u = node.value
        w = node.key

        if u == t:
            break

        for v in graph_dict[u]:
            d = w + graph_dict[u][v]
            if d < estimate[v]:
                estimate[v] = d
                previous[v] = u
                decrease_key(v, d)

    path = []
    while t is not None:
        path.append(t)
        t = previous[t]

    return path[::-1]


def dijkstra3(graph, s, t=None):
    graph_dict = graph.graph_dict
    inf = float('inf')

    estimate = {v: inf for v in graph_dict}
    estimate[s] = 0

    previous = {s: None}
    heap = PairingHeap()
    heap.insert(0, s)

    delete_min = heap.delete_min
    decrease_key = heap.decrease_key
    insert = heap.insert
    lookup = heap.lookup

    while True:
        node = delete_min()
        if node is None:
            break

        u = node.value
        w = node.key

        if u == t:
            break

        for v in graph_dict[u]:
            d = w + graph_dict[u][v]
            if d < estimate[v]:
                estimate[v] = d
                previous[v] = u
                if v in lookup:
                    decrease_key(v, d)
                else:
                    insert(d, v)

    path = []
    while t is not None:
        path.append(t)
        t = previous[t]

    return path[::-1]
