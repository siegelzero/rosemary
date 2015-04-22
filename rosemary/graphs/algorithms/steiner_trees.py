import itertools

from rosemary.graphs.graphs import Graph
from rosemary.graphs.algorithms.spanning_trees import prim
from rosemary.graphs.algorithms.paths import dijkstra


inf = float('inf')


###########################################################################
# Approximation Algorithms
###########################################################################


def distance_network_heuristic(graph, terminals, improve=False):
    r"""Returns an approximate minimal Steiner tree connecting `terminals`
    in `graph`.

    Given a connected, undirected graph `graph` with positive edge weights
    and a subset of the vertices `terminals`, this method finds a subgraph
    with near-minimal total edge cost among all connected subgraphs that
    contain these vertices.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph

    terminals : list or set or tuple
        A collection of vertices of `graph`.

    Returns
    -------
    (weight, edges) : tuple
        The number `weight` is the weight of the Steiner tree. The list
        `edges` is a list of the edges in the Steiner tree.

    Notes
    -----
    The Steiner problem in graphs asks to find the minimal tree connecting
    the terminal vertices. This problem is NP-complete and has proven to be
    extremely difficult even for small problem instances. Given this, it is
    of practical importance to have heuristics for finding low-cost Steiner
    trees.

    This method uses the heuristic algorithm given in [1], which produces a
    tree spanning the terminals with total cost <= 2*(1 - 1/t)*MST, where t
    is the number of terminal vertices and MST is the weight of a minimal
    steiner tree. The implementation given runs in time O(t*|E|*log(|V|)),
    where |E| is the number of edges of `graph` and |V| is the number of
    vertices.

    References
    ----------
    .. [1] L. Kou, G. Markowsky, L. Berman, "A Fast Algorithm for Steiner
    Trees", Acta Informatica, Volume 15, Issue 2, June 1981, 141-145.

    .. [2] P. Winter, "Steiner Problem in Networks: A Survey", Networks,
    Volume 17, Issue 2, Summer 1987, 129-167.

    .. [3] H.J. Proemel, A. Steger, "The Steiner Tree Problem - A Tour
    through Graphs, Algorithms, and Complexity", Vieweg, 2002.
    """
    # Create the distance network induced by the terminal set.
    terminal_set = set(terminals)
    distance_network = Graph()

    # shortest_prev[u] holds the predecessor dict for the shortest path
    # tree rooted at u.
    shortest_prev = {}
    shortest_dist = {}

    for u in terminals:
        u_dist, u_prev = dijkstra(graph, u)
        shortest_dist[u] = u_dist
        shortest_prev[u] = u_prev

    # For each pair (u, v) of terminal vertices, add an edge with weight
    # equal to the length of the shortest u, v path.
    distance_network.add_edges([(u, v, shortest_dist[u][v]) for (u, v) in itertools.combinations(terminals, 2)])

    # Determine the minimum spanning tree of the distance network.
    _, mst_edges = prim(distance_network, edge_list=True)
    subnetwork = Graph()

    # Construct a subnetwork of the graph by replacing each edge in the
    # minimum spanning tree by the corresponding minimum cost path.
    for (u, v) in mst_edges:
        a, b = shortest_prev[u][v], v
        while a is not None:
            subnetwork.add_edge(a, b, graph[a][b])
            a, b = shortest_prev[u][a], a

    # Determine the minimum spanning tree of the subnetwork.
    __, subnetwork_mst_edges = prim(subnetwork, edge_list=True)

    if not improve:
        tree_edges = _delete_nonterminal_leaves(subnetwork_mst_edges, terminal_set)
        tree_weight = sum(graph[u][v] for (u, v) in tree_edges)
    else:
        tree_weight, tree_edges = _improve(graph, subnetwork_mst_edges, terminals)

    return (tree_weight, tree_edges)


def minimum_cost_paths_heuristic(graph, terminals):
    r"""Returns an approximate minimal Steiner tree connecting `terminals`
    in `graph`.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph

    terminals : list or set or tuple

    Returns
    -------
    (weight, edges) : tuple
        The number `weight` is the weight of the Steiner tree. The list
        `edges` is a list of the edges in the Steiner tree.

    Notes
    -----
    This method uses the heuristic algorithm given in [1].
    """
    # shortest_prev[u] holds the predecessor dict for the shortest path
    # tree rooted at u.
    shortest_prev = {}
    shortest_dist = {}

    for u in terminals:
        u_dist, u_prev = dijkstra(graph, u)
        shortest_dist[u] = u_dist
        shortest_prev[u] = u_prev

    best_weight = inf

    for root in terminals:
        tree_edges = []
        tree_vertices = set()
        terminal_set = set(terminals)
        terminal_set.discard(root)
        tree_vertices.add(root)

        while terminal_set:
            min_dist = inf
            for u in terminal_set - tree_vertices:
                for v in tree_vertices:
                    if shortest_dist[u][v] < min_dist:
                        min_dist = shortest_dist[u][v]
                        min_vertices = (u, v)

            u, v = min_vertices
            a, b = shortest_prev[u][v], v
            tree_vertices.add(v)
            terminal_set.discard(u)

            while a is not None:
                if a < b:
                    tree_edges.append((a, b))
                else:
                    tree_edges.append((b, a))

                tree_vertices.add(a)
                a, b = shortest_prev[u][a], a

        tree_weight = sum(graph[u][v] for (u, v) in tree_edges)
        if tree_weight < best_weight:
            best_weight = tree_weight
            best_edges = tree_edges

    return best_weight, best_edges


def _delete_nonterminal_leaves(edges, terminals):
    tree = Graph()

    for (u, v) in edges:
        tree.add_edge(u, v)

    # Now delete all leaves that are not terminals
    while True:
        leaves = tree.leaves()
        changed = False

        for leaf in leaves:
            if leaf not in terminals:
                tree.delete_vertex(leaf)
                changed = True

        if not changed:
            break

    return tree.edges()


def _improve(graph, edges, terminals):
    r"""
    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
    """
    tree_vertices = set()
    for (u, v) in edges:
        tree_vertices.add(u)
        tree_vertices.add(v)

    induced = graph.induced_subgraph(tree_vertices)
    __, mst_edges = prim(induced)

    improved = _delete_nonterminal_leaves(mst_edges, terminals)
    weight = sum(graph[u][v] for (u, v) in improved)

    return weight, improved
