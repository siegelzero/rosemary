from rosemary.graphs.algorithms.traversal import depth_first_search


def connected_components(graph):
    """
    Returns a list of the vertices of each connected component of graph.

    Input:
        * graph: Graph

    Output:
        * components: list
            A list of lists, each sublist corresponding to the vertices of a
            connected component of the graph.
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
    """
    connected_vertices = list(depth_first_search(graph, u))
    return connected_vertices
