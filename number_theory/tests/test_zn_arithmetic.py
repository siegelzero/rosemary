import unittest

from collections import defaultdict

from rosemary.number_theory.sieves import primes

from rosemary.number_theory.zn_arithmetic import (
    discrete_log,
    fibonacci_primitive_roots,
    is_primitive_root,
    linear_congruence,
    nth_roots_of_minus1_mod_p,
    nth_roots_of_unity_mod_p,
    primitive_root,
    primitive_roots,
    quadratic_congruence,
    sqrts_mod_n,
    sqrts_mod_p,
)


class TestCore(unittest.TestCase):
    def setUp(self):
        roots = {
            2: [1],
            3: [2],
            4: [3],
            5: [2, 3],
            6: [5],
            7: [3, 5],
            9: [2, 5],
            10: [3, 7],
            11: [2, 6, 7, 8],
            13: [2, 6, 7, 11],
            14: [3, 5],
            17: [3, 5, 6, 7, 10, 11, 12, 14],
            18: [5, 11],
            19: [2, 3, 10, 13, 14, 15],
            22: [7, 13, 17, 19],
            23: [5, 7, 10, 11, 14, 15, 17, 19, 20, 21],
            25: [2, 3, 8, 12, 13, 17, 22, 23],
            26: [7, 11, 15, 19],
            27: [2, 5, 11, 14, 20, 23],
            29: [2, 3, 8, 10, 11, 14, 15, 18, 19, 21, 26, 27],
            31: [3, 11, 12, 13, 17, 21, 22, 24],
            34: [3, 5, 7, 11, 23, 27, 29, 31],
            37: [2, 5, 13, 15, 17, 18, 19, 20, 22, 24, 32, 35],
            38: [3, 13, 15, 21, 29, 33],
            41: [6, 7, 11, 12, 13, 15, 17, 19, 22, 24, 26, 28, 29, 30, 34, 35],
            43: [3, 5, 12, 18, 19, 20, 26, 28, 29, 30, 33, 34],
            46: [5, 7, 11, 15, 17, 19, 21, 33, 37, 43],
            47: [5, 10, 11, 13, 15, 19, 20, 22, 23, 26, 29, 30, 31, 33, 35, 38, 39, 40, 41, 43, 44, 45],
            49: [3, 5, 10, 12, 17, 24, 26, 33, 38, 40, 45, 47],
            50: [3, 13, 17, 23, 27, 33, 37, 47],
            53: [2, 3, 5, 8, 12, 14, 18, 19, 20, 21, 22, 26, 27, 31, 32, 33, 34, 35, 39, 41, 45, 48, 50, 51],
            54: [5, 11, 23, 29, 41, 47],
            58: [3, 11, 15, 19, 21, 27, 31, 37, 39, 43, 47, 55],
            59: [2, 6, 8, 10, 11, 13, 14, 18, 23, 24, 30, 31, 32, 33, 34, 37, 38, 39, 40, 42, 43, 44, 47, 50, 52, 54,
                 55, 56],
            61: [2, 6, 7, 10, 17, 18, 26, 30, 31, 35, 43, 44, 51, 54, 55, 59],
            62: [3, 11, 13, 17, 21, 43, 53, 55],
            67: [2, 7, 11, 12, 13, 18, 20, 28, 31, 32, 34, 41, 44, 46, 48, 50, 51, 57, 61, 63],
            71: [7, 11, 13, 21, 22, 28, 31, 33, 35, 42, 44, 47, 52, 53, 55, 56, 59, 61, 62, 63, 65, 67, 68, 69],
            73: [5, 11, 13, 14, 15, 20, 26, 28, 29, 31, 33, 34, 39, 40, 42, 44, 45, 47, 53, 58, 59, 60, 62, 68],
            74: [5, 13, 15, 17, 19, 35, 39, 55, 57, 59, 61, 69],
            79: [3, 6, 7, 28, 29, 30, 34, 35, 37, 39, 43, 47, 48, 53, 54, 59, 60, 63, 66, 68, 70, 74, 75, 77],
            81: [2, 5, 11, 14, 20, 23, 29, 32, 38, 41, 47, 50, 56, 59, 65, 68, 74, 77],
            82: [7, 11, 13, 15, 17, 19, 29, 35, 47, 53, 63, 65, 67, 69, 71, 75],
            83: [2, 5, 6, 8, 13, 14, 15, 18, 19, 20, 22, 24, 32, 34, 35, 39, 42, 43, 45, 46, 47, 50, 52, 53, 54, 55, 56,
                 57, 58, 60, 62, 66, 67, 71, 72, 73, 74, 76, 79, 80],
            86: [3, 5, 19, 29, 33, 55, 61, 63, 69, 71, 73, 77],
            89: [3, 6, 7, 13, 14, 15, 19, 23, 24, 26, 27, 28, 29, 30, 31, 33, 35, 38, 41, 43, 46, 48, 51, 54, 56, 58,
                 59, 60, 61, 62, 63, 65, 66, 70, 74, 75, 76, 82, 83, 86],
            94: [5, 11, 13, 15, 19, 23, 29, 31, 33, 35, 39, 41, 43, 45, 57, 67, 69, 73, 77, 85, 87, 91],
            97: [5, 7, 10, 13, 14, 15, 17, 21, 23, 26, 29, 37, 38, 39, 40, 41, 56, 57, 58, 59, 60, 68, 71, 74, 76, 80,
                 82, 83, 84, 87, 90, 92],
            98: [3, 5, 17, 33, 45, 47, 59, 61, 73, 75, 87, 89],
        }

        self.roots = roots

    def test_discrete_log(self):
        self.assertEqual(discrete_log(2, 6, 19), 14)
        self.assertEqual(discrete_log(59, 67, 113), 11)
        self.assertEqual(discrete_log(5, 3, 2017), 1030)

    def test_fibonacci_primitive_roots(self):
        values = [
            (5, [3]),
            (11, [8]),
            (19, [15]),
            (31, [13]),
            (41, [7, 35]),
            (59, [34]),
            (61, [18, 44]),
            (71, [63]),
            (79, [30]),
            (109, [11, 99]),
            (131, [120]),
            (149, [41, 109]),
            (179, [105]),
            (191, [89]),
        ]

        computed = []
        for p in primes(200):
            roots = fibonacci_primitive_roots(p)
            if roots:
                computed.append((p, roots))

        self.assertEqual(values, computed)

    def test_is_primitive_root(self):
        for n in xrange(2, 100):
            roots = [a for a in xrange(1, n) if is_primitive_root(a, n)]
            if n not in self.roots:
                self.assertEqual(roots, [])
            else:
                self.assertEqual(roots, self.roots[n])

        p = 61
        phi_divisors = [2, 3, 5]
        roots = [a for a in xrange(1, p) if is_primitive_root(a, p, phi_divisors=phi_divisors)]
        self.assertEqual(roots, self.roots[61])

        roots = [a for a in xrange(1, p) if is_primitive_root(a, p, phi_divisors=phi_divisors, phi=p - 1)]
        self.assertEqual(roots, self.roots[p])

        roots = [a for a in xrange(1, p) if is_primitive_root(a, p, phi=p - 1)]
        self.assertEqual(roots, self.roots[p])

        self.assertRaisesRegexp(ValueError, 'is_primitive_root: n must be >= 2.', is_primitive_root, 3, -1)

    def test_linear_congruence(self):
        self.assertEqual(linear_congruence(10, 6, 12), [3, 9])
        self.assertEqual(linear_congruence(12, 9, 15), [2, 7, 12])
        self.assertEqual(linear_congruence(10, 3, 12), [])
        self.assertRaisesRegexp(ValueError, 'linear_congruence: Must have n >= 2.', linear_congruence, 2, 3, 0)

    def test_nth_roots_of_minus1_mod_p(self):
        for (p, n) in [(101, 4), (17, 8), (19, 2)]:
            values = [a for a in xrange(p) if a**n % p == p - 1]
            computed = nth_roots_of_minus1_mod_p(n, p)
            self.assertEqual(values, computed)

    def test_nth_roots_of_unity_mod_p(self):
        for (p, n) in [(101, 4), (17, 8), (19, 2)]:
            values = [a for a in xrange(p) if a**n % p == 1]
            computed = nth_roots_of_unity_mod_p(n, p)
            self.assertEqual(values, computed)

    def test_primitive_root(self):
        for n in xrange(2, 100):
            try:
                g = primitive_root(n)
            except:
                continue
            self.assertEqual(g, self.roots[n][0])

    def test_quadratic_congruence(self):
        return

    def test_sqrts_mod_n(self):
        for n in xrange(2, 200):
            values = defaultdict(list)
            for a in xrange(n):
                values[a*a % n].append(a)

            for a in xrange(n):
                self.assertEqual(values[a], sqrts_mod_n(a, n))

    def test_sqrts_mod_p(self):
        self.assertRaisesRegexp(ValueError, "sqrts_mod_p: Must have p >= 2 be prime.", sqrts_mod_p, 10, 1)

        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
            values = defaultdict(list)
            for a in xrange(p):
                values[a*a % p].append(a)

            for a in xrange(p):
                self.assertEqual(values[a], sqrts_mod_p(a, p))

        # Exercise 2.32 of Crandall and Pomerance
        pairs = {
            (3615, 2**16 + 1): [367, 65170],
            (552512556430486016984082237, 2**89 - 1): [1000000000000000000L, 618970018642690137449562111L],
        }

        for (a, p) in pairs:
            self.assertEqual(sqrts_mod_p(a, p), pairs[a, p])
