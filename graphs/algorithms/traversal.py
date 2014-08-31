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
    pop = stack.popleft
    append = stack.append

    visited = set()
    add = visited.add

    while stack:
        (u, depth) = pop()
        if max_depth is None or depth < max_depth:
            for v in graph.graph_dict[u]:
                if v in visited:
                    continue
                append((v, depth + 1))
                add(v)
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
    pop = stack.pop
    extend = stack.extend

    visited = set()
    add = visited.add

    while stack:
        (u, depth) = pop()
        if u in visited:
            continue
        yield u

        add(u)
        if max_depth is None or depth < max_depth:
            extend((v, depth + 1) for v in graph.graph_dict[u] if v not in visited)
