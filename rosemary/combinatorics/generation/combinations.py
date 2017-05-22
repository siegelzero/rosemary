def combinations(n, t):
    """
    generates all combinations of 0, 1, ..., n - 1 of size t.
    """
    c = [0] + [j - 1 for j in range(1, t + 1)] + [n, 0]

    while True:
        yield c[t:0:-1]
        j = 1
        while c[j] + 1 == c[j + 1]:
            c[j] = j - 1
            j += 1
        if j > t:
            return
        c[j] += 1
