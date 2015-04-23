import itertools

from rosemary.graphs.graphs import Graph
from rosemary.graphs.algorithms.spanning_trees import prim
from rosemary.graphs.algorithms.paths import dijkstra


inf = float('inf')


###########################################################################
# Approximation Algorithms
###########################################################################


def distance_network_heuristic(graph, terminals):
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
    Steiner tree. We also apply the improvement procedure given in [3].

    The implementation given runs in time O(t*|E|*log(|V|)), where |E| is
    the number of edges of `graph` and |V| is the number of vertices.

    References
    ----------
    .. [1] L. Kou, G. Markowsky, L. Berman, "A Fast Algorithm for Steiner
    Trees", Acta Informatica, Volume 15, Issue 2, June 1981, 141-145.

    .. [2] P. Winter, "Steiner Problem in Networks: A Survey", Networks,
    Volume 17, Issue 2, Summer 1987, 129-167.

    .. [3] S. Voss, "Steiner's Problem in Graphs: Heuristic Methods",
    Discrete Applied Mathematics 40, 1992, 45-72.

    .. [4] H.J. Proemel, A. Steger, "The Steiner Tree Problem - A Tour
    through Graphs, Algorithms, and Complexity", Vieweg, 2002.

    Examples
    --------
    >>> G = Graph()
    >>> G.add_edges([('u1', 'u2', 1), ('u1', 'v1', 2), ('u1', 'v2', 2),
                     ('u2', 'v3', 4), ('u2', 'v4', 3), ('u2', 'v5', 5),
                     ('u3', 'v1', 2), ('u3', 'v3', 8), ('u4', 'v2', 8),
                     ('u4', 'v5', 8), ('v1', 'v2', 8), ('v1', 'v3', 9),
                     ('v2', 'v5', 5), ('v3', 'v4', 8)])
    >>> distance_network_heuristic(G, ['v1', 'v2', 'v3', 'v4', 'v5'])
    (17, [('u1', 'u2'), ('u1', 'v1'), ('u1', 'v2'), ('u2', 'v3'),
          ('u2', 'v4'), ('v2', 'v5')])
    """
    # Create the distance network induced by the terminal set.
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
    _, subnetwork_mst_edges = prim(subnetwork, edge_list=True)
    tree_weight, tree_edges = _improve(graph, subnetwork_mst_edges, terminals)

    return (tree_weight, tree_edges)


def minimum_cost_paths_heuristic(graph, terminals):
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
    Steiner tree. We use the modification given in [2], using each terminal
    as a root, and keeping the best tree. We also apply the improvement
    procedure given in [2].

    The implementation given runs in time O(t**2*|E|*log(|V|)), where |E|
    is the number of edges of `graph` and |V| is the number of vertices.

    References
    ----------
    .. [1] H. Takahashi, A. Matsuyama, "An Approximate Solution for the
    Steiner Problem in Graphs", Math. Japonica 24, No. 6, 1980, 573-577.

    .. [2] S. Voss, "Steiner's Problem in Graphs: Heuristic Methods",
    Discrete Applied Mathematics 40, 1992, 45-72.

    .. [3] P. Winter, "Steiner Problem in Networks: A Survey", Networks,
    Volume 17, Issue 2, Summer 1987, 129-167.

    Examples
    --------
    >>> G = Graph()
    >>> G.add_edges([('u1', 'u2', 1), ('u1', 'v1', 2), ('u1', 'v2', 2),
                     ('u2', 'v3', 4), ('u2', 'v4', 3), ('u2', 'v5', 5),
                     ('u3', 'v1', 2), ('u3', 'v3', 8), ('u4', 'v2', 8),
                     ('u4', 'v5', 8), ('v1', 'v2', 8), ('v1', 'v3', 9),
                     ('v2', 'v5', 5), ('v3', 'v4', 8)])
    >>> distance_network_heuristic(G, ['v1', 'v2', 'v3', 'v4', 'v5'])
    (17, [('u1', 'u2'), ('u1', 'v1'), ('u1', 'v2'), ('u2', 'v3'),
          ('u2', 'v4'), ('v2', 'v5')])
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

        # Repeat until the tree contains all terminals.
        while terminal_set:
            min_dist = inf
            # Find vertices u and v with minimal distance with v in the
            # tree and u a terminal not in the tree.
            for u in terminal_set - tree_vertices:
                for v in tree_vertices:
                    if shortest_dist[u][v] < min_dist:
                        min_dist = shortest_dist[u][v]
                        min_vertices = (u, v)

            u, v = min_vertices
            a, b = shortest_prev[u][v], v

            tree_vertices.add(u)
            terminal_set.discard(u)

            # Add the edges (a, b) in the shortest u, v path to the tree.
            while a is not None:
                if a < b:
                    tree_edges.append((a, b))
                else:
                    tree_edges.append((b, a))

                tree_vertices.add(a)
                a, b = shortest_prev[u][a], a

        # Apply the improvement procedure to the tree.
        tree_weight, tree_edges = _improve(graph, tree_edges, terminals)

        if tree_weight < best_weight:
            best_weight = tree_weight
            best_edges = tree_edges

    return best_weight, best_edges


def _delete_nonterminal_leaves(edges, terminals):
    r"""Prunes all non-terminal leaves from the tree.

    Given the edges of a tree and a set of terminal vertices, this method
    removes all non-terminal leaves and returns the new tree edges.

    Parameters
    ----------
    edges : list
        A list of the edges of the tree.

    terminals : set
        A set of the terminal vertices.

    Returns
    -------
    tree_edges : list
        A list of the edges of the pruned tree.

    Notes
    -----
    This method deletes all non-terminal leaves from the tree. This process
    is repeated until all leaves are terminals. There is room to improve
    this algorithm.
    """
    # Form a graph from the edge list.
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
    r"""Attempts to improve the quality of a tree.

    Given a graph, the edges of a subtree, and a list of terminal vertices,
    this method attempts to improve the quality of the tree.

    Parameters
    ----------
    graph : rosemary.graphs.graphs.Graph
        Parent graph.

    edges : list
        A list of the edges of a tree in `graph`.

    terminals : list
        A list of the terminal vertices.

    Returns
    -------
    (weight, edges) : tuple
        The number `weight` is the weight of the new tree. The list `edges`
        contains the edges of the new tree.

    Notes
    -----
    This method implements the improvement procedure outlined in [1]. Given
    a candidate tree with edges in `edges`, we construct the subgraph of
    `graph` induced by the vertices of the tree. Next, we construct a
    minimum spanning tree of this induced subgraph. Finally, we delete all
    non-terminal leaves from this new tree.

    References
    ----------
    .. [1] S. Voss, "Steiner's Problem in Graphs: Heuristic Methods",
    Discrete Applied Mathematics 40, 1992, 45-72.
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
