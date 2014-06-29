def maximal_cliques(G):
    vertices = G.vertex_set()

    A = {}
    for v in vertices:
        A[v] = {u for u in vertices if u in G[v]}

    B = {}
    for v in vertices:
        B[v] = {u for u in vertices if u > v}

    def backtrack(X, C, N, l):
        x = X[l - 1]

        if l == 0:
            print []
        else:
            print X[:l]

        if l == 0:
            Nl = vertices
        else:
            Nl = N.intersection(A[x])

        if not Nl:
            print X[:l]

        if l == 0:
            Cl = vertices
        else:
            Cl = (A[x].intersection(B[x])).intersection(C)

        for e in C:
            Xl = X[:]
            Xl[l] = e
            backtrack(Xl, Cl, Nl, l + 1)

    backtrack([0]*len(vertices), {}, {}, 0)

