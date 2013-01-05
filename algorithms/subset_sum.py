# Subset sum algorithms

################################################################################

def subset_sum_recursive(L, target):
    """
    subset_sum_recursive(L, s):
    This finds the largest sum of elements of L that is <= s.
    """
    if sum(L) < target:
        return sum(L)

    max_arr = [0, []]
    A = sorted(L)
    ll = len(A)

    def backtrack(used, partial_sum, last):
        for part in A:
            if part >= last:
                break
            t_sum = partial_sum + part
            if t_sum <= target:
                t_used = used + [ part ]
                if t_sum >= max_arr[0]:
                    max_arr[0] = t_sum
                    max_arr[1] = t_used
                if t_sum == target:
                    return 0
                else:
                    backtrack(t_used, t_sum, part)

    backtrack([], 0, A[-1] + 1)
    return max_arr

################################################################################

def subset_sum_middle(L, target):
    """
    subset_sum_middle(L, s):
    This finds the largest sum of elements of L that is <= s.
    The algorithm used is a meet-in-the-middle greedy algorithm.
    """
    S = sorted(L)
    ll = len(S)
    P1 = S[:ll // 2]
    P2 = S[ll // 2:]

    A = set([0])
    for a in P1:
        A |= set(a + e for e in A if a + e <= target)

    B = set([0])
    for b in P2:
        B |= set(b + e for e in B if b + e <= target)

    sums1 = sorted(A)
    sums2 = sorted(B, reverse=True)

    l1, l2 = len(sums1), len(sums2)
    i = j = best = 0

    while i < l1 and j < l2:
        current = sums1[i] + sums2[j]
        if current > target:
            j += 1
        elif current < target:
            i += 1
            if current > best:
                best = current
        else:
            return current

    return best

################################################################################

def subset_sum_dynamic(L, target):
    """
    subset_sum_dynamic(L, target):
    This finds the largest sum of elements of L that is <= s.
    The algorithm is a standard dynamic programming algorithm.
    """
    M = {}
    n = len(L)
    for w in xrange(target + 1):
        M[0, w] = 0

    for i in xrange(1, n):
        for w in xrange(0, target + 1):
            if w < L[i]:
                M[i, w] = M[i - 1, w]
            else:
                M[i, w] = max(M[i - 1, w], L[i] + M[i - 1, w - L[i]])

    return M[n - 1, target]

