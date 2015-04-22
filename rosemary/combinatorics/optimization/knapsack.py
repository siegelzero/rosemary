"""
This module contains algorithms to find exact and approximate solutions to the
optimization knapsack problem:

Given values v_0, v_1, ..., v_{n - 1}, weights w_0, w_1, ..., w_{n - 1}, and a
maximum capacity M, find a 0,1 n-tuple [x_0, x_1, ..., x_{n - 1}] such that
P = \sum_{i = 0}^{n - 1} v_i x_i is maximized, subject to the capacity
constraint \sum_{i = 0}^{n - 1} w_i x_i <= M.
"""

import math
import random


def depth_first(items, capacity):
    """
    Returns the optimal solution to the given knapsack problem.

    This function performs a depth-first search on the solution space to find
    the optimum solution to the given knapsack problem. Only elementary pruning
    is done to the solution tree, so this algorithm is only effective for small
    problem instances.

    Input:
        * items: list
            Each item in the list is a pair (v, w) of value and weight.

        * capacity: int
            Maximum capacity of the knapsack.

    Output:
        * (best_value, best_items): tuple
            * best_value: int
                The best value found.

            * best_items: list
                The items giving this value.

    Examples:
        >>> items = [(135, 70), (139, 73), (149, 77), (150, 80), (156, 82),\
                (163, 87), (173, 90), (184, 94), (192, 98), (201, 106),\
                (210, 110), (214, 113), (221, 115), (229, 118), (240, 120)]
        >>> capacity = 750
        >>> depth_first(items, capacity)
        (1458, [1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1])
    """
    num_items = len(items)
    best = [0, []]

    def backtrack(idx, current_value, current_weight, used):
        if current_weight <= capacity:
            if idx == num_items:
                if current_value > best[0]:
                    best[0] = current_value
                    best[1] = used
            else:
                (v, w) = items[idx]
                backtrack(idx + 1, current_value + v, current_weight + w, used + [(v, w)])
                backtrack(idx + 1, current_value, current_weight, used)

    backtrack(0, 0, 0, [])

    return tuple(best)


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

    # We sort so that the output agrees with the other functions
    best = (M[capacity][n], sorted(used))
    return best


def branch_and_bound_depth_first(items, capacity):
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


def simulated_annealing(items, capacity):
    """
    Returns an approximate solution to the knapsack problem.

    This function uses a simulated annealing strategy to approximate the optimal
    value to the given knapsack problem.

    Input:
        * items: list
            Each item in the list is a pair (v, w) of value and weight.

        * capacity: int
            Maximum capacity of the knapsack.

    Output:
        * (best_value, best_items): tuple
            * best_value: int
                The best value found.

            * best_items: list
                The items giving this value.

    Examples:
        >>> items = [(135, 70), (139, 73), (149, 77), (150, 80), (156, 82),\
                (163, 87), (173, 90), (184, 94), (192, 98), (201, 106),\
                (210, 110), (214, 113), (221, 115), (229, 118), (240, 120)]
        >>> capacity = 750
        >>> simulated_annealing(items, capacity)
        (1458, [1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1])

    Details:
        This algorithm is based on Algorithm 5.20 from "Combinatorial
        Algorithms: Generation, Enumeration, and Search" by Kreher and Stinson.
    """
    cmax = 100000
    temp = 1000
    alpha = 0.99999
    n = len(items)

    current_items = [0]*n
    current_weight = 0
    best_items = list(current_items)
    best_value = 0

    def P(X):
        total = 0
        for i in xrange(n):
            if X[i]:
                total += items[i][0]
        return total

    for i in xrange(cmax):
        # Let Y be a random neighbor of the current configuration.
        j = random.randint(0, n - 1)
        neighbor = list(current_items)
        neighbor[j] = 1 - neighbor[j]

        # If this new neighbor is not feasible, continue
        if neighbor[j] and current_weight + items[j][1] > capacity:
            continue

        if neighbor[j]:
            # Here, we know that P(neighbor) > P(current_items)
            current_items = neighbor
            current_weight += items[j][1]
            current_value = P(current_items)

            if current_value > best_value:
                best_items = list(current_items)
                best_value = current_value
        else:
            r = random.random()
            diff = P(neighbor) - P(current_items)
            if r < math.exp(diff / temp):
                current_items = neighbor
                current_weight -= items[j][1]

        temp *= alpha

    return best_value, best_items
