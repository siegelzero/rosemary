# Algorithms for solving the knapsack problem.

def exhaustive(items, capacity):
    """
    Returns optimal solution to the knapsack problem.

    Input:
        * items: (list)
            Each item in the list is a pair (v, w) of value and weight.

        * capacity: (int)
            Maximum capacity of the knapsack.
    """
    num_items = len(items)
    best = [0, []]

    def backtrack(idx, current_value, current_weight, used):
        if current_weight <= capacity:
            if current_value > best[0]:
                best[0] = current_value
                best[1] = used

            for i in xrange(idx + 1, num_items):
                v, w = items[i]
                new_value = current_value + v
                new_weight = current_weight + w
                backtrack(i, new_value, new_weight, used + [(v, w)])

    backtrack(-1, 0, 0, [])

    return best


def dynamic_programming(items, capacity):
    n = len(items)
    M = [[0]*(n + 1) for _ in xrange(capacity + 1)]
    keep = set()

    for i in xrange(n):
        v, w = items[i]
        for k in xrange(capacity + 1):
            if k < w:
                M[k][i + 1] = M[k][i]
            else:
                if M[k][i] > v + M[k - w][i]:
                    M[k][i + 1] = M[k][i]
                else:
                    M[k][i + 1] = v + M[k - w][i]
                    keep.add((k, i + 1))

    used = []
    w = capacity
    for i in xrange(n, 0, -1):
        if (w, i) in keep:
            used.append(items[i - 1])
            w -= items[i - 1][1]

    best = [M[capacity][n], used]
    return best


def branch_and_bound_dfs(items, capacity):
    unit_cost = [(float(v)/float(w), v, w, i) for (i, (v, w)) in enumerate(items)]
    unit_cost.sort(reverse=True)
    num_items = len(items)
    best = [0, []]

    def optimistic_bound(idx, target):
        if idx == num_items:
            return 0
        value = 0
        weight = 0
        triples = (unit_cost[i] for i in xrange(idx, num_items))
        for (u, v, w, _) in triples:
            if weight + w > target:
                break
            weight += w
            value += v
        return value + u*(target - weight)

    def backtrack(idx, current_value, current_weight, used):
        if current_weight <= capacity:
            if current_value > best[0]:
                best[0] = current_value
                best[1] = used

            if idx == num_items:
                choices = []
            else:
                (_, v, w, i) = unit_cost[idx]
                if current_weight + w <= capacity:
                    choices = [1, 0]
                else:
                    choices = [0]

            bound = current_value + optimistic_bound(idx, capacity - current_weight)

            for choice in choices:
                if bound <= best[0]:
                    return
                new_value = current_value + choice*v
                new_weight = current_weight + choice*w
                new_used = used + [(v, w, i)] if choice else used
                backtrack(idx + 1, new_value, new_weight, new_used)

    backtrack(0, 0, 0, [])

    return best

