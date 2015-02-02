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
        """Initializes the continued fraction expansion of (p + sqrt(d))/q, for
        (p, q) != (0, 1).
        """
        r = integer_sqrt(d)

        if r*r == d:
            raise ValueError("QuadraticIrrational: d cannot be a perfect square.")

        # A quadratic irrational can be written in the form (P + sqrt(D))/Q,
        # for some integers P, Q != 0, and D not a perfect square, with
        # (D - P*P) = 0 (mod D). If we don't have this last divisibility
        # condition, then we modify P, Q, and D as appropriate. See Theorem 9.23
        # of "The Theory of Numbers - A Text and Source Book of Problems" by
        # Adler and Coury for details.
        if (d - p*p) % q == 0:
            Q = q
            P = p
            D = d
        else:
            abs_q = abs(q)
            Q = abs_q*q
            P = abs_q*p
            D = q*q*d
            r = integer_sqrt(D)

        pq_pairs = {(p, q): 0}
        partial_quotients = []

        for i in itertools.count(1):
            a = (P + r)//Q
            P = a*Q - P
            Q = (D - P*P)//Q
            partial_quotients.append(a)

            # Encountering a (P, Q) pair which we have already seen means that
            # we have started a new period, and so we can break out of the loop.
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
        """Initializes the continued fraction expansion of sqrt(d).
        """
        r = integer_sqrt(d)
        if r*r == d:
            raise ValueError("QuadraticIrrational: d cannot be a perfect square.")

        Q = 1
        P = 0
        terms = []

        while True:
            q = (P + r)//Q
            terms.append(q)
            # The continued fraction expansion of sqrt(d) is periodic of the
            # form [a0; a1, a2, ..., a2, a1, 2*a0], so we can terminate when we
            # encounter a partial quotient equal to 2*a0. See Theorem 5.10 of
            # "Fundamental Number with Applications" by Mollin for details.
            if q == 2*r:
                break
            P = q*Q - P
            Q = (d - P*P)//Q

        self.period_length = len(terms) - 1
        self.pre_period = [r]
        self.fundamental_period = terms[1:]

    def partial_quotients(self):
        """Returns the partial quotients of self.

        Returns:
            * X: generator
                The values output from this generator are the partial quotients
                in the continued fraction expansion of self.

        Examples:
            >>> X = QuadraticIrrational(13).partial_quotients()
            >>> [X.next() for _ in xrange(10)]
            [3, 1, 1, 1, 1, 6, 1, 1, 1, 1]
        """
        terms = itertools.chain(self.pre_period,
                                itertools.cycle(self.fundamental_period))
        for t in terms:
            yield t

    def convergents(self):
        """Returns the continued fraction convergents of self.

        Returns:
            * X: generator
                The values output from this generator are pairs (a, b) representing
                the convergents a/b to sqrt(d).

        Examples:
            >>> X = QuadraticIrrational(2).convergents()
            >>> [X.next() for _ in xrange(7)]
            [(1, 1), (3, 2), (7, 5), (17, 12), (41, 29), (99, 70), (239, 169)]

        Details:
            See Theorem 1.12 of "Fundamental Number Theory with Applications" by
            Mollin for details.
        """
        a0, a1 = 0, 1
        b0, b1 = 1, 0

        for q in self.partial_quotients():
            a2 = q*a1 + a0
            b2 = q*b1 + b0
            a0, a1 = a1, a2
            b0, b1 = b1, b2
            yield (a2, b2)
