import random
import unittest

from rosemary.number_theory.classification import (
    is_power,
    is_square,
    is_squarefree,
)


class TestClassification(unittest.TestCase):
    def test_is_squarefree(self):
        return

    def test_is_power(self):
        self.assertRaisesRegexp(ValueError, "is_power: Must have k >= 1.", is_power, 10, 0)
        self.assertRaisesRegexp(ValueError, "is_power: Must have n >= 1.", is_power, 0, 3)
        self.assertRaisesRegexp(ValueError, "is_power: Must have n >= 1.", is_power, [(-1, 1), (2, 2)])

        for k in [2, 3, 5, 7, 11]:
            self.assertEqual(is_power(k), False)

        for k in xrange(10):
            a = random.randint(10, 100)
            b = random.randint(10, 100)
            n = a**b
            x, y = is_power(n)

            self.assertTrue(x**y == n)

        for k in xrange(10):
            a = random.randint(10, 100)
            b = random.randint(10, 100)
            n = a**b
            x, y = is_power(n, b)

            self.assertTrue(x**y == n)

        values = [
            4, 8, 9, 16, 25, 27, 32, 36, 49, 64, 81, 100, 121, 125, 128, 144, 169, 196, 216, 225, 243, 256, 289, 324,
            343, 361, 400, 441, 484, 512, 529, 576, 625, 676, 729, 784, 841, 900, 961, 1000
        ]

        powers = [e for e in xrange(1, 1001) if is_power(e)]

        self.assertEqual(values, powers)

        self.assertEqual(is_power(10, 1), (10, 1))
        self.assertEqual(is_power(125, 3), (5, 3))
        self.assertEqual(is_power(125, 4), False)

        self.assertEqual(is_power([(2, 3), (5, 2)], 3), False)
        self.assertEqual(is_power([(2, 4), (3, 2)]), (12, 2))
        self.assertEqual(is_power([(2, 4), (3, 3)]), False)
        self.assertEqual(is_power([(3, 4), (5, 4)]), (15, 4))
        self.assertEqual(is_power([(3, 4), (5, 4)], 2), (225, 2))

    def test_is_square(self):
        for k in [-2, -1, 2, 3, 5, 7, 11, 2545]:
            self.assertEqual(is_square(k), False)

        for (k, r) in [(9, 3), (16, 4), (289, 17), (130321, 361)]:
            self.assertEqual(is_square(k), r)

        self.assertEqual(is_square([(2, 4), (5, 2)]), 20)
        self.assertEqual(is_square([(2, 4), (5, 3)]), False)

        values = [
            0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225, 256, 289, 324, 361, 400, 441, 484, 529,
            576, 625, 676, 729, 784, 841, 900, 961
        ]

        squares = [a for a in xrange(1000) if is_square(a)]

        self.assertEqual(squares, values)


