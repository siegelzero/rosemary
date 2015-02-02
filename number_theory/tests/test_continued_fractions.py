import unittest
import itertools

from rosemary.number_theory.continued_fractions import (
    QuadraticIrrational,
)


class TestCore(unittest.TestCase):
    def test_pure_sqrt1(self):
        # Example 5.10 from "Fundamental Number Theory with Applications" by Mollin.
        cf = QuadraticIrrational(425)
        self.assertEqual(str(cf), 'Continued Fraction Expansion of sqrt(425)')
        self.assertEqual((cf.d, cf.p, cf.q), (425, 0, 1))
        self.assertEqual(cf.pre_period, [20])
        self.assertEqual(cf.fundamental_period, [1, 1, 1, 1, 1, 1, 40])
        self.assertEqual(list(itertools.islice(cf.convergents(), 6)),
                         [(20, 1), (21, 1), (41, 2), (62, 3), (103, 5), (165, 8)])

    def test_quadratic_irrational1(self):
        # Example 5.8 from "Fundamental Number Theory with Applications" by Mollin.
        cf = QuadraticIrrational(22, 4, 6)
        self.assertEqual(str(cf), 'Continued Fraction Expansion of (4 + sqrt(22))/6')
        self.assertEqual((cf.d, cf.p, cf.q), (22, 4, 6))
        self.assertEqual(cf.fundamental_period, [1, 2, 4, 2, 1, 8])
        self.assertEqual(cf.pre_period, [])
        self.assertEqual(cf.period_length, 6)
        self.assertEqual(list(itertools.islice(cf.convergents(), 6)),
                         [(1, 1), (3, 2), (13, 9), (29, 20), (42, 29), (365, 252)])

    def test_pure_sqrt2(self):
        self.assertRaisesRegexp(ValueError,
                                "QuadraticIrrational: d cannot be a perfect square.",
                                QuadraticIrrational, 100)

    def test_quadratic_irrational2(self):
        self.assertRaisesRegexp(ValueError,
                                "QuadraticIrrational: d cannot be a perfect square.",
                                QuadraticIrrational, 100, 2, 3)

    def test_quadratic_irrational3(self):
        """Test string representation."""
        cf = QuadraticIrrational(10, 1)
        self.assertEqual(str(cf), 'Continued Fraction Expansion of 1 + sqrt(10)')
        pq = cf.partial_quotients()
        self.assertEqual([pq.next() for _ in xrange(10)], [4, 6, 6, 6, 6, 6, 6, 6, 6, 6])

    def test_quadratic_irrational4(self):
        """Test string representation."""
        cf = QuadraticIrrational(7, 0, 2)
        self.assertEqual(str(cf), 'Continued Fraction Expansion of sqrt(7)/2')

    def test_quadratic_irrational5(self):
        # This example is one where (D - P*P) is not divisible by Q, and so they
        # must be modified.
        cf = QuadraticIrrational(10, 1, 2)
        self.assertEqual(str(cf), 'Continued Fraction Expansion of (1 + sqrt(10))/2')
        self.assertEqual((cf.d, cf.p, cf.q), (10, 1, 2))
        self.assertEqual(cf.pre_period, [2])
        self.assertEqual(cf.fundamental_period, [12, 3])
        self.assertEqual(cf.period_length, 2)
        self.assertEqual(list(itertools.islice(cf.convergents(), 6)),
                         [(2, 1), (25, 12), (77, 37), (949, 456), (2924, 1405), (36037, 17316)])
