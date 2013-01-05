import heapq

"""
Store a weighted graph as a dict G, where 'for v in G' runs through the
vertices of G, and G[v] is a dict W, where the keys of W are the neighbors of
v, and W[v] is the weight of the edge between them.
"""

################################################################################

def kruskal(G):
    """
    kruskal(G):
    Given a weighted graph G, this algorithm returns a minimal spanning tree
    for G using Kruskal's algorithm. The algorithm builds a spanning tree by
    successively inserting edges from G in order of increasing cost. We
    insert each edge e as long as it does not create a cycle when added to
    the edges we've already inserted. If, on the other hand, inserting e
    would result in a cycle, then we simply discard e and continue.
    """
    def find(parent, u):
        """
        This subroutine returns a representative element from C that contains
        u. Using this, we can determine whether two vertices u and v belong to
        the same tree by testing whether find(C, u) = find(C, v).
        """
        if C[u] != u:
            C[u] = find(C, C[u])
        return C[u]

    def union(C, R, u, v):
        """
        This subroutine
        """
        u = find(C, u)
        v = find(C, v)
        if R[u] > R[v]:
            C[v] = u
        else:
            C[u] = v
        if R[u] == R[v]:
            R[v] += 1

    E = [(G[u][v], u, v) for u in G for v in G[u]]
    T = set()
    C = {u:u for u in G}
    R = {u:0 for u in G}
    
    for (_, u, v) in sorted(E):
        if find(C, u) != find(C, v):
            T.add((u, v))
            union(C, R, u, v)

    return T

################################################################################

def prim(G, s):
    P = {}
    Q = [(0, None, s)]

    while Q:
        (_, p, u) = heapq.heappop(Q)
        if u in P:
            continue
        P[u] = p
        for (v, w) in G[u].items():
            heapq.heappush(Q, (w, u, v))
    return P

################################################################################

def prob107():
    f = open('/home/brownkenny/network.txt')
    L = [ e.replace('-', '0') for e in f ]
    M = [ [ int(e) for e in l.split(',') ] for l in L ]

    G = {}
    for i in xrange(len(M)):
        W = {}
        for j in xrange(len(M)):
            if M[i][j] > 0:
                W[j] = M[i][j]
        G[i] = W.copy()

    T = kruskal(G)
    min_wt = sum(G[u][v] for (u, v) in T)
    orig_wt = sum([ sum(e) for e in M ]) // 2

    print "Original Weight:", orig_wt
    print "Min Weight:", min_wt
    print "Savings:", orig_wt - min_wt
