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
        for k in xrange(20):
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

        for i in xrange(20):
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

        values = [
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
            1, 2, 1, 2, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 1, 2, 1, 4, 1, 2, 1, 4, 1, 2, 1, 4,
            1, 2, 1, 4, 1, 2, 1, 4, 1, 1, 1, 1, 5, 1, 1, 1, 1, 5, 1, 1, 1, 1, 5, 1, 1, 1, 1, 5, 1, 2, 3, 2, 1, 6, 1, 2,
            3, 2, 1, 6, 1, 2, 3, 2, 1, 6, 1, 2, 1, 1, 1, 1, 1, 1, 7, 1, 1, 1, 1, 1, 1, 7, 1, 1, 1, 1, 1, 1, 1, 2, 1, 4,
            1, 2, 1, 8, 1, 2, 1, 4, 1, 2, 1, 8, 1, 2, 1, 4, 1, 1, 3, 1, 1, 3, 1, 1, 9, 1, 1, 3, 1, 1, 3, 1, 1, 9, 1, 1,
            1, 2, 1, 2, 5, 2, 1, 2, 1, 10, 1, 2, 1, 2, 5, 2, 1, 2, 1, 10, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 11, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 2, 3, 4, 1, 6, 1, 4, 3, 2, 1, 12, 1, 2, 3, 4, 1, 6, 1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 13, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 7, 2, 1, 2, 1, 2, 1, 14, 1, 2, 1, 2, 1, 2, 1, 1, 3, 1, 5, 3,
            1, 1, 3, 5, 1, 3, 1, 1, 15, 1, 1, 3, 1, 5, 1, 2, 1, 4, 1, 2, 1, 8, 1, 2, 1, 4, 1, 2, 1, 16, 1, 2, 1, 4, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 17, 1, 1, 1, 1, 2, 3, 2, 1, 6, 1, 2, 9, 2, 1, 6, 1, 2, 3, 2, 1,
            18, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 19, 1, 1, 2, 1, 4, 5, 2, 1, 4, 1, 10, 1, 4,
            1, 2, 5, 4, 1, 2, 1, 20
        ]

        i = 0
        for a in xrange(1, 21):
            for b in xrange(1, 21):
                self.assertEqual(gcd(a, b), values[i])
                i += 1

    def test_gcd_list(self):
        self.assertEqual(gcd_list([2, 4, 6, 8]), 2)
        self.assertEqual(gcd_list([1, 4, 6, 8]), 1)
        self.assertEqual(gcd_list([-2, 4, 6, 8]), 2)
        self.assertEqual(gcd_list([0, 4, 6, 8]), 2)
        self.assertEqual(gcd_list([2]), 2)

    def test_integer_log(self):
        self.assertRaisesRegexp(ValueError, "integer_log: Must have b >= 2.", integer_log, 1, 40)
        self.assertRaisesRegexp(ValueError, "integer_log: Must have n >= 1.", integer_log, 3, 0)

        self.assertEqual(integer_log(2, 100), 6)
        self.assertEqual(integer_log(5, 30), 2)
        self.assertEqual(integer_log(5, 2), 0)
        self.assertEqual(integer_log(5, 125), 3)

        for i in xrange(20):
            b = random.randint(2, 20)
            n = random.randint(10**6, 10**7)
            k = integer_log(b, n)
            self.assertTrue(b**k <= n < b**(k + 1))

    def test_integer_nth_root(self):
        self.assertRaisesRegexp(ValueError, "integer_nth_root: Must have n >= 1", integer_nth_root, 0, 10)
        self.assertRaisesRegexp(ValueError, "integer_nth_root: Must have m >= 0", integer_nth_root, 3, -1)

        self.assertEqual(integer_nth_root(1, 10), 10)
        self.assertEqual(integer_nth_root(3, 1), 1)

        for i in xrange(100):
            a = random.randint(100, 1000)
            b = random.randint(100, 1000)
            m = a**b

            n = random.randint(2, 200)

            r = integer_nth_root(n, m)
            self.assertTrue(r**n <= m < (r + 1)**n)

    def test_integer_sqrt(self):
        self.assertRaisesRegexp(ValueError, "integer_sqrt: Must have n >= 0.", integer_sqrt, -1)

        for i in xrange(100):
            a = random.randint(100, 1000)
            b = random.randint(100, 1000)
            n = a**b
            r = integer_sqrt(n)

            self.assertTrue(r**2 <= n < (r + 1)**2)

    def test_inverse_mod(self):
        self.assertRaisesRegexp(ValueError, "inverse_mod: Integers must be relatively prime.", inverse_mod, 2, 4)
        self.assertRaisesRegexp(ValueError, "inverse_mod: Must have m >= 2.", inverse_mod, 10, 1)

        self.assertEqual(inverse_mod(5, 17), 7)
        self.assertEqual(inverse_mod(2, 9), 5)

        for k in xrange(1, 101):
            a = inverse_mod(k, 101)
            self.assertEqual(a*k % 101, 1)

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
        self.assertRaisesRegexp(ValueError, "is_square: Must have n >= 0.", is_square, -2)
        self.assertRaisesRegexp(ValueError, "is_square: Must have n >= 0.", is_square, [(-1, 1), (2, 2)])

        for k in [2, 3, 5, 7, 11]:
            self.assertEqual(is_square(k), False)

        for (k, r) in [(9, 3), (16, 4), (289, 17), (130321, 361)]:
            self.assertEqual(is_square(k), r)

        self.assertEqual(is_square([(2, 4), (5, 2)]), 20)
        self.assertEqual(is_square([(2, 4), (5, 3)]), False)

        values = [
            1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225, 256, 289, 324, 361, 400, 441, 484, 529, 576,
            625, 676, 729, 784, 841, 900, 961
        ]

        squares = [a for a in xrange(1, 1000) if is_square(a)]

        self.assertEqual(squares, values)

    def test_jacobi_symbol(self):
        self.assertRaisesRegexp(ValueError, "jacobi_symbol: Must have m odd.", jacobi_symbol, 19, 4)

        values = [
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, -1, 0, -1, 1, 1, -1, 0, -1, 1, 1, 0, 1, -1, 0, -1, -1, 0, -1, 1, 1, 1, -1, 0, 1, -1,
            -1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, -1, 1, -1, 1, 0,
            -1, -1, 1, -1, 1, 1, 1, 0, -1, -1, -1, 1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 1, 1, -1, -1, 1, -1, 0, -1, 1,
            -1, 1, -1, 1, 0, 1, 1, 1, -1, -1, -1, 1, 0, 0, 1, 0, 1, -1, 0, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1,
            -1, -1, 1, -1, 1, 1, 0, 1, 1, 0, -1, 1, 0, -1, -1, 0, 1, -1, 1, 1, 1, -1, 1, -1, -1, 1, 1, 0, 1, -1, 0, -1,
            1, 1, -1, 0, -1, 1
        ]

        i = 0
        for a in xrange(1, 21):
            for b in xrange(1, 21, 2):
                self.assertEqual(jacobi_symbol(a, b), values[i])
                i += 1

    def test_lcm(self):
        self.assertTrue(lcm(10, 0) == lcm(0, 10) == 0)
        self.assertEqual(lcm(0, 0), 0)

        for a in xrange(10):
            for b in xrange(10):
                self.assertTrue(lcm(a, b) == lcm(b, a) == lcm(-a, b) == lcm(a, -b))

        values = [
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 2, 2, 6, 4, 10, 6, 14, 8, 18, 10, 22,
            12, 26, 14, 30, 16, 34, 18, 38, 20, 3, 6, 3, 12, 15, 6, 21, 24, 9, 30, 33, 12, 39, 42, 15, 48, 51, 18, 57,
            60, 4, 4, 12, 4, 20, 12, 28, 8, 36, 20, 44, 12, 52, 28, 60, 16, 68, 36, 76, 20, 5, 10, 15, 20, 5, 30, 35,
            40, 45, 10, 55, 60, 65, 70, 15, 80, 85, 90, 95, 20, 6, 6, 6, 12, 30, 6, 42, 24, 18, 30, 66, 12, 78, 42, 30,
            48, 102, 18, 114, 60, 7, 14, 21, 28, 35, 42, 7, 56, 63, 70, 77, 84, 91, 14, 105, 112, 119, 126, 133, 140, 8,
            8, 24, 8, 40, 24, 56, 8, 72, 40, 88, 24, 104, 56, 120, 16, 136, 72, 152, 40, 9, 18, 9, 36, 45, 18, 63, 72,
            9, 90, 99, 36, 117, 126, 45, 144, 153, 18, 171, 180, 10, 10, 30, 20, 10, 30, 70, 40, 90, 10, 110, 60, 130,
            70, 30, 80, 170, 90, 190, 20, 11, 22, 33, 44, 55, 66, 77, 88, 99, 110, 11, 132, 143, 154, 165, 176, 187,
            198, 209, 220, 12, 12, 12, 12, 60, 12, 84, 24, 36, 60, 132, 12, 156, 84, 60, 48, 204, 36, 228, 60, 13, 26,
            39, 52, 65, 78, 91, 104, 117, 130, 143, 156, 13, 182, 195, 208, 221, 234, 247, 260, 14, 14, 42, 28, 70, 42,
            14, 56, 126, 70, 154, 84, 182, 14, 210, 112, 238, 126, 266, 140, 15, 30, 15, 60, 15, 30, 105, 120, 45, 30,
            165, 60, 195, 210, 15, 240, 255, 90, 285, 60, 16, 16, 48, 16, 80, 48, 112, 16, 144, 80, 176, 48, 208, 112,
            240, 16, 272, 144, 304, 80, 17, 34, 51, 68, 85, 102, 119, 136, 153, 170, 187, 204, 221, 238, 255, 272, 17,
            306, 323, 340, 18, 18, 18, 36, 90, 18, 126, 72, 18, 90, 198, 36, 234, 126, 90, 144, 306, 18, 342, 180, 19,
            38, 57, 76, 95, 114, 133, 152, 171, 190, 209, 228, 247, 266, 285, 304, 323, 342, 19, 380, 20, 20, 60, 20,
            20, 60, 140, 40, 180, 20, 220, 60, 260, 140, 60, 80, 340, 180, 380, 20
        ]

        i = 0
        for a in xrange(1, 21):
            for b in xrange(1, 21):
                self.assertTrue(lcm(a, b) == lcm(b, a) == values[i])
                i += 1

    def test_lcm_list(self):
        self.assertEqual(lcm_list(range(1, 31)), 2329089562800)

    def test_power_mod(self):
        self.assertRaisesRegexp(ValueError, "power_mod: Must have a >= 0.", power_mod, -2, 3, 10)
        self.assertRaisesRegexp(ValueError, "power_mod: Must have k >= 0.", power_mod, 2, -1, 10)
        self.assertRaisesRegexp(ValueError, "power_mod: Must have m >= 1.", power_mod, 2, 3, 0)

        self.assertEqual(power_mod(2, 45, 17), 15)
        self.assertEqual(power_mod(3, 100, 101), 1)
        self.assertEqual(power_mod(3, 100, 3), 0)
        self.assertEqual(power_mod(3, 0, 19), 1)

        for i in xrange(20):
            a = random.randint(2, 100)
            b = random.randint(2, 100)
            m = random.randint(2, 100)

            self.assertEqual(power_mod(a, b, m), pow(a, b, m))

    def test_valuation(self):
        values = [
            1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 0,
            2, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 4, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 3, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0,
            0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0,
            2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 1, 0, 1, 3, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 1, 0, 6, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0,
            1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 1, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 1, 0, 0, 0,
            4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 2, 0, 0, 1, 2, 0, 2, 0, 0
        ]

        i = 0
        for a in xrange(2, 101):
            for p in [2, 3, 5, 7, 11]:
                self.assertEqual(valuation(p, a), values[i])
                i += 1

if __name__ == "__main__":
    unittest.main()
