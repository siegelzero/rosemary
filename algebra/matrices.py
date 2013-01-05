################################################################################
# Functions associated to matrices
################################################################################

def matrix_add(M, N):
    nrows = len(M)
    ncols = len(M[0])

    if len(N) != nrows or len(N[0]) != ncols:
        raise ValueError("Invalid Dimensions")

    S = M[:]
    for i in xrange(nrows):
        for j in xrange(ncols):
            S[i][j] += N[i][j]
    return S

def matrix_mul(M, N):
    nrows = len(M)
    ncols = len(M[0])

    if len(N) != nrows or len(N[0]) != ncols:
        raise ValueError("Invalid Dimensions")

 
