import unittest
from collections import defaultdict

from rosemary.number_theory.arithmetic_functions.functions import euler_phi

from rosemary.number_theory.arithmetic_functions.inverse_functions import (
    euler_phi_inverse,
)


class TestCore(unittest.TestCase):
    def test_euler_phi_inverse(self):
        self.assertRaisesRegexp(ValueError, "euler_phi_inverse: Must have n > 0.", euler_phi_inverse, -1)
        self.assertEqual(euler_phi_inverse([(2, 2), (5, 2)]), [101, 125, 202, 250])
        values = defaultdict(list)

        for k in xrange(1, 1001):
            phi = euler_phi(k)
            values[phi].append(k)

        for v in xrange(1, 201):
            inv = euler_phi_inverse(v)
            self.assertEqual(inv, values[v])


if __name__ == "__main__":
    unittest.main()
