import collections

################################################################################

def three_term_progressions(A):
    """
    three_term_progression(A):
    Given an array of integers A, this determines if A has a three term
    arithmetic progression. If True, this returns their indices in A.
    Otherwise, returns False.
    """
    n = len(A)
    L = []
    for j in range(1, n):
        # we start at A[j] and try to find an element larger than A[j] and
        # and element smaller than A[j] whose average is A[j]
        i = j - 1
        k = j + 1
        while i >= 0 and k < n:
            if A[i] + A[k] < 2 * A[j]:
                k += 1
            elif A[i] + A[k] > 2 * A[j]:
                i -= 1
            else:
                L += [(i, j, k)]
                break
                #return (i, j, k)

    #return False
    return L

################################################################################

def longest_ap(A, length = None):
    # L[i, j] will store the maximum length of an arithmetic progression
    # whose first two terms are A[i] and A[j]
    L = collections.defaultdict(int)
    n = len(A)
    longest = 2
    for j in range(n - 2, 0, -1):
        i = j - 1
        k = j + 1
        while i >= 0 and k <= n - 1:
            if A[i] + A[k] < 2 * A[j]:
                k += 1
            elif A[i] + A[k] > 2 * A[j]:
                L[i, j] = 2
                i -= 1
            else:
                L[i, j] = L[j, k] + 1
                longest = max(longest, L[i, j])
                i -= 1
                k += 1
        """
        this next block of code isn't necessary to find the longest AP, but it's
        needed if we wish to compute the lengths of all APs.
        """
        while i >= 0:
            L[i, j] = 2
            i -= 1

    if length is None:
        length = longest

    progs = []

    for (a, b) in L:
        if L[(a, b)] == length:
            diff = A[b] - A[a]
            progs.append([ A[a] + diff * i for i in range(length) ])

    return progs

################################################################################
