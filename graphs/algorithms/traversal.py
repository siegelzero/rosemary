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


def breadth_first_search_tree(graph, start):
    """
    Returns the breadth-first search tree of the graph.

    Given a graph and a starting vertex, this method returns the breadth-first
    search tree of graph, rooted at the vertex start.

    Input:
        * graph: Graph

        * start: vertex of graph

    Output:
        * previous: dict
            Predecessor dict of the search tree.

    Details:
        >>> graph = Graph()
        >>> graph.add_edges([('a', 'b'), ('a', 's'), ('b', 'c'), ('c', 's'),
                             ('d', 'e'), ('d', 's'), ('e', 'd'), ('e', 's')])
        >>> breadth_first_search_tree(graph)
        {'a': 's', 'b': 'a', 'c': 's', 'd': 's', 'e': 's', 's': None}
    """
    graph_dict = graph.graph_dict
    stack = deque([start])
    pop = stack.popleft
    append = stack.append
    previous = {start: None}

    while stack:
        u = pop()
        for v in graph_dict[u]:
            if v in previous:
                continue
            previous[v] = u
            append(v)

    return previous


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
