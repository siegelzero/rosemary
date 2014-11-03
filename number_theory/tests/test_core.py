import unittest
import random

from rosemary.number_theory.core import (
    bit_count,
    chinese,
    chinese_preconditioned,
    crt_preconditioning_data,
    ext_gcd,
    gcd,
    gcd_list,
    integer_log,
    integer_nth_root,
    integer_sqrt,
    inverse_mod,
    is_power,
    is_square,
    jacobi_symbol,
    lcm,
    lcm_list,
    power_mod,
    valuation,
)


class TestCore(unittest.TestCase):
    def test_bit_count(self):
        for k in xrange(10):
            self.assertTrue(bit_count(-k) == bit_count(k) == bin(k).count('1'))

    def test_chinese(self):
        solution = chinese([(1, 5), (2, 7)])
        self.assertTrue(solution % 5 == 1)
        self.assertTrue(solution % 7 == 2)

        solution = chinese([(0, 4), (1, 25)])
        self.assertEqual(solution, 76)

        solution = chinese([(2, 3), (3, 5), (2, 7)])
        self.assertEqual(solution, 23)

        self.assertRaisesRegexp(ValueError, 'chinese: Moduli must be coprime.', chinese, [(3, 8), (5, 10)])

    def test_crt_preconditioning_data(self):
        moduli = [3, 5, 7]
        preconditioning_data = crt_preconditioning_data(moduli)
        self.assertEqual(preconditioning_data, (3, [3, 5, 7], [1, 3, 15], [1, 2, 1], 105))

    def test_chinese_preconditioned(self):
        moduli = [3, 5, 7]
        preconditioning_data = crt_preconditioning_data(moduli)

        system = [(2, 3), (3, 5), (2, 7)]
        solution = chinese_preconditioned(system, preconditioning_data)
        self.assertEqual(solution, 23)

        system = [(1, 3), (2, 5), (4, 7)]
        solution = chinese_preconditioned(system, preconditioning_data)
        self.assertEqual(solution, 67)

    def test_ext_gcd(self):
        self.assertEqual(ext_gcd(5, 7), (3, -2, 1))
        self.assertEqual(ext_gcd(12, 10), (1, -1, 2))
        self.assertEqual(ext_gcd(12, -18), (-1, -1, 6))
        self.assertEqual(ext_gcd(-3, -1), (0, -1, 1))

        for i in xrange(10):
            x = random.randint(10, 100)
            y = random.randint(10, 100)
            (a, b, d) = ext_gcd(x, y)
            self.assertEqual(a*x + b*y, d)

    def test_gcd(self):
        self.assertEqual(gcd(20, 10), 10)
        self.assertEqual(gcd(10, 20), 10)

        self.assertEqual(gcd(-10, 20), 10)
        self.assertEqual(gcd(-20, 10), 10)

        self.assertEqual(gcd(-10, -10), 10)

        self.assertEqual(gcd(1, 11213), 1)
        self.assertEqual(gcd(11213, 1), 1)

        self.assertEqual(gcd(0, 10), 10)
        self.assertEqual(gcd(10, 0), 10)
        self.assertEqual(gcd(0, -10), 10)

    def test_gcd_list(self):
        self.assertEqual(gcd_list([2, 4, 6, 8]), 2)
        self.assertEqual(gcd_list([1, 4, 6, 8]), 1)
        self.assertEqual(gcd_list([-2, 4, 6, 8]), 2)
        self.assertEqual(gcd_list([0, 4, 6, 8]), 2)

    def test_integer_log(self):
        self.assertEqual(integer_log(2, 100), 6)
        self.assertEqual(integer_log(5, 30), 2)
        self.assertRaisesRegexp(ValueError, "integer_log: Must have b >= 2.", integer_log, 1, 40)
        self.assertRaisesRegexp(ValueError, "integer_log: Must have n >= 1.", integer_log, 3, 0)

    def test_integer_nth_root(self):
        self.assertRaisesRegexp(ValueError, "integer_nth_root: Must have n >= 1", integer_nth_root, 0, 10)
        self.assertRaisesRegexp(ValueError, "integer_nth_root: Must have m >= 0", integer_nth_root, 3, -1)

        for i in xrange(10):
            a = random.randint(100, 1000)
            b = random.randint(100, 1000)
            m = a**b

            n = random.randint(2, 200)

            r = integer_nth_root(n, m)
            self.assertTrue(r**n <= m < (r + 1)**n)


if __name__ == "__main__":
    unittest.main()
