from rosemary.graphs.algorithms.traversal import breadth_first_search


def connected_components(graph):
    """
    Returns a list of the vertices of each connected component of graph.

    Given a graph, this method returns a list containing lists of the vertices
    of each connected component of the graph.

    Input:
        * graph: Graph

    Output:
        * components: list

    Examples:
        >>> G = Graph()
        >>> G.add_vertices(['a', 'b', 'c', 'd', 'e', 'f',
                            'g', 'h', 'i', 'j', 'k', 'l'])
        >>> G.add_edges([('a', 'b'), ('a', 'e'), ('e', 'i'), ('e', 'j'),
                         ('i', 'j'), ('c', 'd'), ('c', 'g'), ('c', 'h'),
                         ('d', 'h'), ('g', 'h'), ('g', 'k'), ('h', 'k'),
                         ('h', 'l')])
        >>> connected_components(G)
        [['a', 'b', 'e', 'i', 'j'], ['c', 'h', 'd', 'g', 'k', 'l'], ['f']]

    Details:
        This method uses a breadth-first-search approach to finding the
        connected components of the given graph. We begin at a start vertex, and
        use visit all vertices connected to the start vertex, storing each
        vertex along the way. Once the search is complete, we have a connected
        component of the graph. We then proceed to find an unvisited vertex, and
        start the search again from this new start, repeating this until all
        vertices are visited. Since this algorithm is basically a
        breadth-first-search of the entire graph, we see that the runtime of
        this method is O(|V| + |E|).
    """
    visited = set()
    components = []

    for u in graph.graph_dict:
        if u in visited:
            continue
        neighbors = list(breadth_first_search(graph, u))
        components.append(neighbors)
        visited.update(neighbors)

    return components


def connected_component(graph, u):
    """
    Returns a list of the vertices of graph in the same connected component as
    the vertex u

    Input:
        * graph: Graph

        * u: vertex of graph

    Output:
        * component: list
            A list of vertices of graph.

    Examples:
        >>> G = Graph()
        >>> G.add_vertices(['a', 'b', 'c', 'd', 'e', 'f',
                            'g', 'h', 'i', 'j', 'k', 'l'])
        >>> G.add_edges([('a', 'b'), ('a', 'e'), ('e', 'i'), ('e', 'j'),
                         ('i', 'j'), ('c', 'd'), ('c', 'g'), ('c', 'h'),
                         ('d', 'h'), ('g', 'h'), ('g', 'k'), ('h', 'k'),
                         ('h', 'l')])
        >>> connected_component(G, 'a')
        ['a', 'b', 'e', 'i', 'j']

    Details:
        This method uses performs a breadth-first-search beginning at vertex u
        to visit all vertices connected to u. We return a list of these
        vertices.  Since we are performing a BFS on the graph, the complexity of
        this algorithm is O(|V| + |E|).
    """
    connected_vertices = list(breadth_first_search(graph, u))
    return connected_vertices


def is_connected(graph):
    r"""Returns True if the graph is connected and False otherwise.
    """
    vertex_set = graph.vertex_set()

    if len(vertex_set) == 0:
        return True

    u = list(vertex_set)[0]
    connected_vertices = set(breadth_first_search(graph, u))
    return connected_vertices == vertex_set
