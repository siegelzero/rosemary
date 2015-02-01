from rosemary.number_theory.continued_fractions import QuadraticIrrational
from rosemary.number_theory.core import integer_sqrt
from rosemary.number_theory.zn_arithmetic import sqrts_mod_n

import sys
import itertools


def pell_fundamental_solution(D, n=1):
    if D <= 0 or n not in (1, -1):
        raise ValueError("pell_fundamental_solution: Must have D > 0 and n in (-1, 1).")

    for (x, y) in QuadraticIrrational(D).convergents(D):
        if x*x - D*y*y == n:
            yield (x, y)


def _pell_general(D, N, one_solution=False):
    """
    test with D = 1121311213, N = 11889485036288588
    """
    if D < 0 or N < 0:
        raise ValueError("_pell_general: Must have D > 0 and N > 0")

    # Find fundamental solution to Pell
    for (t, u) in QuadraticIrrational(D).convergents():
        if t*t - D*u*u == 1:
            break

    # Now find fundamental solution to general Pell
    fundamental_solutions = []
    residues = sqrts_mod_n(N, D)
    L = N*(t - 1)

    for r in residues:
        for x in itertools.count(r, D):
            w = (x*x - N)//D

            if w < 0:
                continue
            elif 2*D*w > L:
                break

            y = integer_sqrt(w)

            if y*y == w:
                if one_solution:
                    yield (x, y)
                    return

                fundamental_solutions.append((x, y))
                (r, s) = (-x, y)

                if (x*r - D*y*s) % N != 0 or (x*s - y*r) % N != 0:
                    fundamental_solutions.append((r, s))

    minimal_positive_solutions = []

    for (x, y) in fundamental_solutions:
        if x < 0:
            (r, s) = (-x, -y)
            p = r*t + s*u*D
            q = r*u + s*t
            minimal_positive_solutions.append((p, q))
        else:
            minimal_positive_solutions.append((x, y))

    minimal_positive_solutions.sort()

    if not minimal_positive_solutions:
        return

    while True:
        for i in xrange(len(minimal_positive_solutions)):
            x, y = minimal_positive_solutions[i]
            yield (x, y)
            minimal_positive_solutions[i] = (x*t + y*u*D, x*u + y*t)


def pivotal(m):
    for (x, y) in _pell_general(4*m*(m + 1), 4*m**2*(m + 1)**3):
        if y % (m + 1) == 0 and x % (2*m*(m + 1)) == 0:
            k2 = y//(m + 1) + m
            n2 = x//(2*m*(m + 1)) - m - 1

            if k2 % 2 == 0 and n2 % 2 == 0 and k2 >= 2 and n2 >= k2:
                yield (k2//2, n2//2)


def solve(N):
    sols = set()
    old = 0

    for m in itertools.count(1):
        not_changed2 = 0
        for (k, n) in pivotal(m):
            if k <= N:
                if k not in sols:
                    sols.add(k)
            else:
                not_changed2 += 1

            if not_changed2 >= 10:
                break

        new = sum(sols)
        if new != old:
            print new
            old = new


if __name__ == "__main__":
    N = int(sys.argv[1])
    print solve(N)
