from collections import deque


def breadth_first_search(graph, start, max_depth=None):
    """
    Returns an iterator over the vertices of graph in a breadth-first ordering.

    Input:
        * graph: Graph

        * start: vertex
            Vertex from which to start the traversal.

        * max_depth: int (default=None)
            Maximum depth to traverse from the start vertex.

    Ouput:
        * vertices: generator
            Iterator over the vertices of graph.
    """
    # Starting vertex has depth 0.
    stack = deque([(start, 0)])
    visited = set()

    while stack:
        (u, depth) = stack.popleft()
        if max_depth is None or depth < max_depth:
            for v in graph.graph_dict[u]:
                if v in visited:
                    continue
                stack.append((v, depth + 1))
                visited.add(v)
                yield v


def depth_first_search(graph, start, max_depth=None):
    """
    Returns an iterator over the vertices of graph in a depth-first ordering.

    Input:
        * graph: Graph

        * start: vertex
            Vertex from which to start the traversal.

        * max_depth: int (default=None)
            Maximum depth to traverse from the start vertex.

    Ouput:
        * vertices: generator
            Iterator over the vertices of graph.
    """
    # Starting vertex has depth 0.
    stack = [(start, 0)]
    visited = set()

    while stack:
        (u, depth) = stack.pop()
        if u in visited:
            continue
        yield u
        visited.add(u)
        if max_depth is None or depth < max_depth:
            stack.extend((v, depth + 1) for v in graph.graph_dict[u] if v not in visited)
