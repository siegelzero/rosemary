from rosemary.number_theory.classification import is_square
from rosemary.number_theory.continued_fractions import QuadraticIrrational
from rosemary.number_theory.core import integer_sqrt
from rosemary.number_theory.zn_arithmetic import sqrts_mod_n

import sys
import itertools


def pell_fundamental_solution(D, n=1):
    """Returns the fundamental solution to the Pell equation x**2 - D*y**2
    == n, where n in (-1, 1).

    Given D > 0 not a square, and n in (-1, 1), this method returns the
    fundamental solution to the Pell equation described above. The
    fundamental solution (x, y) is the one with least positive value of x,
    and correspondingly the least positive value of y.

    Input:
        * D: int (D > 0)

        * n: int (n in (-1, 1)).

    Returns:
        * (x, y): tuple

    Raises:
        * ValueError: If D <= 0, D is perfect square, n not in (-1, 1), or if no
        solutions exist.

    Examples:
        >>> pell_fundamental_solution(61)
        (1766319049, 226153980)
        >>> 1766319049**2 - 61*226153980**2
        1
        >>> pell_fundamental_solution(17, -1)
        (4, 1)
        >>> 4**2 - 17*1**2
        -1
        >>> pell_fundamental_solution(15, -1)
        Traceback (most recent call last):
        ...
        ValueError: pell_fundamental_solution: Solution nonexistent.
        >>> pell_fundamental_solution(15, -2)
        Traceback (most recent call last):
        ...
        ValueError: pell_fundamental_solution: Must have D > 0 not a perfect square and n in (-1, 1).

    Details:
        For D > 0 not a perfect square, the equation x**2 - D*y**2 == 1
        always has solutions, while the equation x**2 - D*y**2 == -1 only
        has solutions when the continued fraction expansion of sqrt(D) has
        odd period length.

        See Corollary 5.7 of "Fundamental Number Theory with Applications"
        by Mollin for details. See also the article "Simple Continued
        Fraction Solutions for Diophantine Equations" by Mollin.
    """
    if D <= 0 or is_square(D) or n not in (1, -1):
        raise ValueError("pell_fundamental_solution: Must have D > 0 not a perfect square and n in (-1, 1).")

    contfrac = QuadraticIrrational(D)

    # No solution for n == -1 if the period length is even.
    if n == -1 and contfrac.period_length % 2 == 0:
        raise ValueError("pell_fundamental_solution: Solution nonexistent.")

    # Otherwise, solutions always exist for D not a perfect square.
    for (x, y) in contfrac.convergents():
        if x*x - D*y*y == n:
            return (x, y)


def _mollin(D, N):
    positive_roots = sqrts_mod_n(D, abs(N))

    roots = []
    print positive_roots
    for r in positive_roots:
        for e in [-r, r]:
            if -abs(N) < 2*e <= abs(N):
                roots.append(e)

    return roots


def _pell_general(D, N, one_solution=False):
    """
    test with D = 1121311213, N = 11889485036288588
    """
    # Find fundamental solution to Pell
    for (t, u) in QuadraticIrrational(D).convergents():
        if t*t - D*u*u == 1:
            break

    # Now find fundamental solution to general Pell
    fundamental_solutions = []
    residues = sqrts_mod_n(N, D)
    L = N*(t - 1)
    print residues

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

    print minimal_positive_solutions

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
