################################################################################
# Algorithms related to Continued Fractions.
################################################################################

from rosemary.number_theory.core import integer_sqrt


def quadratic_expansion(d):
    """Returns the continued fraction expansion of sqrt(d).

    Input:
        * d: int (d > 0)

    Returns:
        * terms: list
            A list of the terms in the continued fraction expansion.

    Examples:
        >>> quadratic_expansion(5)
        [2, 4]
        >>> quadratic_expansion(31)
        [5, 1, 1, 3, 5, 3, 1, 1, 10]

    Details:
        If d > 0 is not a perfect square, then the continued fraction expansion
        of sqrt(d) is periodic of the form [a0; a1, a2, ..., a2, a1, 2*a0],
        where a0 = floor(sqrt(d)), the periodic part beginning after the first
        term. See Theorem 5.10 of "Fundamental Number Theory with Applications" by
        Mollin for details.
    """
    r = integer_sqrt(d)

    if r*r == d:
        return [r]

    Q = 1
    P = 0
    terms = []

    while True:
        q = (P + r)//Q

        terms.append(q)
        # The only terms equal to 2*floor(sqrt(d)) are the last terms in a
        # period, so we can terminate once this occurs.
        if q == 2*r:
            return terms

        P = q*Q - P
        Q = (d - P*P)//Q


def quadratic_convergents(d, number=0):
    """Returns the continued fraction convergents to sqrt(d).

    Input:
        * d: int (d > 0)

    Returns::q


    Details:
        See Theorem 1.12 of "Fundamental Number Theory with Applications" by
        Mollin for details.
    """
    expansion = quadratic_expansion(d)

    a0, a1 = 0, 1
    b0, b1 = 1, 0

    for q in expansion:
        a2 = q*a1 + a0
        b2 = q*b1 + b0

        yield (a2, b2)

        a0, a1 = a1, a2
        b0, b1 = b1, b2
