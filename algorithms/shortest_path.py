from rosemary.utilities import cached_function

"""
Store a weighted graph as a dict G, where 'for v in G' runs through the
vertices of G, and G[v] is a dict W, where the keys of W are the neighbors of
v, and W[v] is the weight of the edge between them.
"""

################################################################################

def bellman_ford(G, s):
    """
    bellman_ford(G, s):
    Given a weighted graph G and a node s, this returns D, P, where D is the
    length of the shortest path from s to each other node in G, and P is a dict
    containing the nodes in the shortest path

    Examples:
    >>> M = [[131, 673],
    >>>      [201, 96]]
    >>> W = { (0,0): {(0,1):673, (1,0):201}, (0,1): {(1,1):96}, (1,0): {(1,1):96} }
    >>> D, P = bellman_ford(W, (0, 0))
    >>> D
    {(0, 0): 0, (0, 1): 673, (1, 0): 201, (1, 1): 297}
    >>> P
    {(0, 1): (0, 0), (1, 0): (0, 0), (1, 1): (1, 0)}
    >>> D[(1, 1)] # this is the shortest path length from (0,0) to (1,1)
    297
    >>> 201 + 96
    297
    """
    inf = float('inf')

    def relax(W, u, v, D, P):
        d = D.get(u, inf) + W[u][v]
        if d < D.get(v, inf):
            D[v], P[v] = d, u
            return True

    D, P = {s:0}, {}
    for rnd in G:
        changed = False
        for u in G:
            for v in G[u]:
                if relax(G, u, v, D, P):
                    changed = True
        if not changed:
            break
    else:
        raise ValueError('negative cycle')

    return D, P

################################################################################

def shortest_path_dag(G, s, t):
    """
    shortest_path_dag(G, s, t):
    This finds the length of the shortest path from node s to node t in the
    directed acyclic graph G.
    """
    @cached_function
    def d(u):
        if u == t:
            return 0
        return min(G[u][v] + d(v) for v in G[u])
    return d(s)

################################################################################
