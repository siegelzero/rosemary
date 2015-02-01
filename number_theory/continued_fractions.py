import itertools

from rosemary.number_theory.core import integer_sqrt


class QuadraticIrrational(object):
    """A real number a is a quadratic irrational if it can be written in the
    form (p + sqrt(d))/q, for some integers p, q != 0, and d > 0, where d is not
    a perfect square.
    """
    def __init__(self, d, p=0, q=1):
        """Initializes a new quadratic irrational (p + sqrt(d))/q.
        """
        self.p = p
        self.q = q
        self.d = d

        if p == 0 and q == 1:
            self._initialize_sqrt(d)
        else:
            self._initialize_quadratic(d, p, q)

    def __repr__(self):
        """String representation of self.
        """
        if self.p == 0:
            if self.q == 1:
                s = "sqrt({})".format(self.d)
            else:
                s = "sqrt({})/{}".format(self.d, self.q)
        else:
            if self.q == 1:
                s = "{} + sqrt({})".format(self.q, self.d)
            else:
                s = "({} + sqrt({}))/{}".format(self.p, self.d, self.q)

        return "Continued Fraction Expansion of {}".format(s)

    def _initialize_quadratic(self, d, p, q):
        r = integer_sqrt(d)

        if r*r == d:
            raise ValueError("d cannot be a perfect square")

        Q = q
        P = p
        pq_pairs = {(p, q): 0}
        partial_quotients = []

        for i in itertools.count(1):
            q = (P + r)//Q
            P = q*Q - P
            Q = (d - P*P)//Q
            partial_quotients.append(q)

            if (P, Q) in pq_pairs:
                j = pq_pairs[(P, Q)]
                k = i
                break
            else:
                pq_pairs[(P, Q)] = i

        self.period_length = k - j
        self.pre_period = partial_quotients[0:j]
        self.fundamental_period = partial_quotients[j:]

    def _initialize_sqrt(self, d):
        """If d > 0 is not a perfect square, then the continued fraction
        expansion of sqrt(d) is periodic of the form [a0; a1, a2, ..., a2, a1,
        2*a0], where a0 = floor(sqrt(d)), the periodic part beginning after the
        first term. See Theorem 5.10 of "Fundamental Number Theory with
        Applications" by Mollin for details.
        """
        r = integer_sqrt(d)
        if r*r == d:
            raise ValueError("d cannot be a perfect square")

        Q = 1
        P = 0
        terms = []

        while True:
            q = (P + r)//Q
            terms.append(q)
            if q == 2*r:
                break
            P = q*Q - P
            Q = (d - P*P)//Q

        self.period_length = len(terms) - 1
        self.pre_period = [r]
        self.fundamental_period = terms[1:]

    def convergents(self):
        """Returns the continued fraction convergents to self.

        Returns:
            * X: generator
                The values output from this generator are pairs (a, b) representing
                the convergents a/b to sqrt(d).

        Examples:
            >>> quadratic_convergents(2, 6)
            [(1, 1), (3, 2), (7, 5), (17, 12), (41, 29), (99, 70), (239, 169)]

        Details:
            See Theorem 1.12 of "Fundamental Number Theory with Applications" by
            Mollin for details.
        """
        a0, a1 = 0, 1
        b0, b1 = 1, 0

        terms = itertools.chain(self.pre_period,
                                itertools.cycle(self.fundamental_period))

        for q in terms:
            a2 = q*a1 + a0
            b2 = q*b1 + b0
            a0, a1 = a1, a2
            b0, b1 = b1, b2
            yield (a2, b2)
