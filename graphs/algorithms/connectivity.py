from rosemary.graphs.algorithms.traversal import depth_first_search


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
        >>> G = random_graph(10, 0.4)
        >>> connected_components(G)
        [[0, 1, 2, 9, 6], [3, 8, 7], [4], [5]]

    Details:
        This method uses a standard depth-first-search approach to finding the
        connected components of the given graph. We begin at a start vertex, and
        use DFS to visit all vertices connected to the start vertex, storing
        each vertex along the way. Once the search is complete, we have a
        connected component of the graph. We then proceed to find an unvisited
        vertex, and start the DFS again from this new start, repeating this
        until all vertices are visited. Since this algorithm is basically a
        depth-first-search of the entire graph, we see that the runtime of this
        method is O(|V| + |E|).
    """
    visited = set()
    components = []

    for u in graph.graph_dict:
        if u in visited:
            continue
        neighbors = list(depth_first_search(graph, u))
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
        >>> G = random_graph(10, 0.4)
        >>> connected_components(G)
        [[0, 1, 2, 9, 6], [3, 8, 7], [4], [5]]
        >>> connected_component(G, 0)
        [0, 1, 2, 9, 6]

    Details:
        This method uses performs a depth-first-search beginning at vertex u to
        visit all vertices connected to u. We return a list of these vertices.
        Since we are performing a DFS on the graph, the complexity of this
        algorithm is O(|V| + |E|).
    """
    connected_vertices = list(depth_first_search(graph, u))
    return connected_vertices
