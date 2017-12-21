from rosemary.number_theory.classification import is_square
from rosemary.number_theory.continued_fractions import QuadraticIrrational
from rosemary.number_theory.zn_arithmetic import sqrts_mod_n


class PellEquation(object):
    def __init__(self, D, N):
        self.D = D
        self.N = N

        if D <= 0 or is_square(D):
            raise ValueError("PellEquation: Must have D > 0 not a perfect square.")

    def __repr__(self):
        return "Pell Equation x^2 - {}*y^2 = {}".format(self.D, self.N)

    def resolvent_solution(self):
        r"""Returns the fundamental solution to the Pell equation
        x**2 - D*y**2 == n, where n in (-1, 1).

        Given D > 0 not a square, and n in (-1, 1), this method returns the
        fundamental solution to the Pell equation described above. The
        fundamental solution (x, y) is the one with least positive value of x,
        and correspondingly the least positive value of y.

        Parameters
        ----------
        D : int (D > 0)

        n : int (n == 1 or n == -1)

        Returns
        -------
        (x, y): tuple

        Raises
        ------
        ValueError : If D <= 0, D is perfect square, n not in (-1, 1), or if no
        solution exists.

        Notes
        -----
        If D is a positive integer that is not a perfect square, then the
        equation x**2 - D*y**2 == 1 has a solution in positive integers.

        On the other hand, of the continued fraction expansion of sqrt(D) has
        even period length, then there are no solutions to the equation
        x**2 - D*y**2 == -1. Otherwise, the equation has a solution.

        See Theorem 3.2.1 in [1] for the positive case and Theorem 3.6.1 in [1]
        for the negative case.

        References
        ----------
        .. [1] T. Andreescu, D. Andrica, "Quadratic Diophantine Equations",
        Springer-Verlag, New York, 2015.

        Examples
        --------
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


    def pell_fundamental_solution(D, n=1):
        r"""Returns the fundamental solution to the Pell equation
        x**2 - D*y**2 == n, where n in (-1, 1).

        Given D > 0 not a square, and n in (-1, 1), this method returns the
        fundamental solution to the Pell equation described above. The
        fundamental solution (x, y) is the one with least positive value of x,
        and correspondingly the least positive value of y.

        Parameters
        ----------
        D : int (D > 0)

        n : int (n == 1 or n == -1)

        Returns
        -------
        (x, y): tuple

        Raises
        ------
        ValueError : If D <= 0, D is perfect square, n not in (-1, 1), or if no
        solution exists.

        Notes
        -----
        If D is a positive integer that is not a perfect square, then the
        equation x**2 - D*y**2 == 1 has a solution in positive integers.

        On the other hand, of the continued fraction expansion of sqrt(D) has
        even period length, then there are no solutions to the equation
        x**2 - D*y**2 == -1. Otherwise, the equation has a solution.

        See Theorem 3.2.1 in [1] for the positive case and Theorem 3.6.1 in [1]
        for the negative case.

        References
        ----------
        .. [1] T. Andreescu, D. Andrica, "Quadratic Diophantine Equations",
        Springer-Verlag, New York, 2015.

        Examples
        --------
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


def pell_small(D, N):
    if D <= 0 or is_square(D) or D*D <= N:
        raise ValueError("pell_small: Must have D > 0 not a perfect square and N < sqrt(D)")

    contfrac = QuadraticIrrational(D)

    # Otherwise, solutions always exist for D not a perfect square.
    for (x, y) in contfrac.convergents():
        print (x, y)
        if x*x - D*y*y == N:
            yield (x, y)


def _pell_small_minimal_positive_solutions(D, N):
    if D <= 0 or is_square(D) or D*D <= N:
        raise ValueError("pell_small: Must have D > 0 not a perfect square and N < sqrt(D)")

    contfrac = QuadraticIrrational(D)
    terms = []

    # Otherwise, solutions always exist for D not a perfect square.
    for (x, y) in contfrac.convergents():
        if x*x - D*y*y == 1:
            break
        terms.append((x,  y))

    sols = []
    for (x, y) in terms:
        k = x*x - D*y*y

        if N % k != 0:
            continue

        f = is_square(N // k)
        if f is not False:
            sols.append((f*x, f*y))

    sols.sort()
    return sols


def pell_minimal_positive_solutions(D, N):
    r"""Returns the minimal positive solutions to the Pell equation
    x**2 - D*y**2 == N.

    Given D > 0 not a square, and N != 0, this method returns a list of the
    minimal positive solutions in each class for the Pell equation
    described above.

    Parameters
    ----------
    D : int (D > 0)

    N : int (N != 0)

    Returns
    -------
    L : list
        List of the minimal positive solutions of the equation. The list is
        empty if there are no solutions.

    Notes
    -----

    References
    ----------
    .. [1] T. Andreescu, D. Andrica, "Quadratic Diophantine Equations",
    Springer-Verlag, New York, 2015.

    .. [2] R.A. Mollin, "Fundamental Number Theory with Applications",
    Chapman & Hall/CRC, 2008.

    Examples
    --------
    >>> pell_minimal_positive_solutions(13, 27)
    [(12, 3), (40, 11), (220, 61), (768, 213)]
    >>> pell_minimal_positive_solutions(
    """
    (u, v) = pell_fundamental_solution(D, 1)

    # If N == 1, there is only one class of solutions.
    if N == 1:
        return [(u, v)]

    if N > 0:
        # See Theorem 7.1 in [2]
        B = int(((u + 1)*N / 2.0)**(0.5)) + 1
    else:
        # See Theorem 7.2 in [2]
        B = int(((u - 1)*(-N) / 2.0)**(0.5)) + 1

    # We can proceed in two different ways at this point:
    # 1) We can check all values of y in the range described in Theorems
    # 7.1 and 7.2, finding the corresponding valid values of x for each.
    #
    # 2) By solving a quadratic congruence, we get congruences that x must
    # satisfy. By iterating through these progressions, we are able to
    # check the entire range faster than in 1.
    #
    # We use the second method here.

    residues = sqrts_mod_n(N, D)
    fundamental_solutions = []

    # Find the fundamental solutions for each class.
    for r in residues:
        for x in xrange(r, B, D):
            y = is_square((x*x - N) // D)
            if y > 0 and y is not False:
                if (-x*x - D*y*y) % N == 0 and 2*x*y % N == 0:
                    fundamental_solutions.append((x, y))
                else:
                    fundamental_solutions.append((x, y))
                    fundamental_solutions.append((-x, y))

    # Find the minimal positive solutions for each class.
    minimal_positive_solutions = []
    for (x, y) in fundamental_solutions:
        if x < 0:
            if N > 0:
                (p, q) = (-x, -y)
            else:
                (p, q) = (x, y)
            p, q = p*u + q*v*D, p*v + q*u
        else:
            (p, q) = (x, y)
        minimal_positive_solutions.append((p, q))

    minimal_positive_solutions.sort()

    return minimal_positive_solutions


def pell_solution_generator(D, N, B):
    (u, v) = pell_fundamental_solution(D, 1)
    minimal = pell_minimal_positive_solutions(D, N)

    for (x, y) in minimal:
        (p, q) = (x, y)
        while p <= B and q <= B:
            yield (p, q)
            p, q = p*u + q*v*D, p*v + q*u


def farey(D, N):
    n1, d1 = 0, 1
    n2, d2 = 1, 0

    while True:
        a = n1 + n2
        b = d1 + d2
        print (a, b)

        t = a*a - D*b*b

        if N == t:
            return (a, b)

        if t > 0:
            n2, d2 = a, b
        else:
            n1, d1 = a, b
