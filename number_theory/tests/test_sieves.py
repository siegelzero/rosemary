import unittest

from rosemary.number_theory.sieves import (
    chartres,
    eratosthenes1,
    eratosthenes2,
    eratosthenes3,
    factored_xrange,
    luo,
    moebius_xrange,
    prime_xrange,
    primes,
    primes_first_n,
)


class TestCore(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.primes_to_1000 = [
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103,
            107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223,
            227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
            349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
            467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
            613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743,
            751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883,
            887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997
        ]

        self.sum_of_primes_to_200000 = 1709600813
        self.sum_of_primes_to_2000000 = 142913828922

        self.primes_first_200 = [
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103,
            107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223,
            227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
            349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
            467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
            613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743,
            751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883,
            887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031,
            1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153,
            1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223
        ]

    def test_chartres(self):
        values = chartres(1000)
        self.assertEqual(values, self.primes_to_1000)

        values = chartres(0)
        self.assertEqual(values, [])

        values = chartres(-1)
        self.assertEqual(values, [])

        value = sum(chartres(2*10**5))
        self.assertEqual(value, self.sum_of_primes_to_200000)

    def test_eratosthenes1(self):
        values = eratosthenes1(1000)
        self.assertEqual(values, self.primes_to_1000)

        values = eratosthenes1(0)
        self.assertEqual(values, [])

        values = eratosthenes1(-1)
        self.assertEqual(values, [])

        value = sum(eratosthenes1(2*10**5))
        self.assertEqual(value, self.sum_of_primes_to_200000)

    def test_eratosthenes2(self):
        values = eratosthenes2(1000)
        self.assertEqual(values, self.primes_to_1000)

        values = eratosthenes2(0)
        self.assertEqual(values, [])

        values = eratosthenes2(-1)
        self.assertEqual(values, [])

        value = sum(eratosthenes2(2*10**5))
        self.assertEqual(value, self.sum_of_primes_to_200000)

    def test_eratosthenes3(self):
        values = eratosthenes3(1000)
        self.assertEqual(values, self.primes_to_1000)

        values = eratosthenes3(0)
        self.assertEqual(values, [])

        values = eratosthenes3(-1)
        self.assertEqual(values, [])

        value = sum(eratosthenes3(2*10**5))
        self.assertEqual(value, self.sum_of_primes_to_200000)

    def test_factored_xrange(self):
        values = factored_xrange(10000)
        for (n, n_fac) in values:
            prod = 1
            for (p, e) in n_fac:
                prod *= p**e
            self.assertEqual(n, prod)

    def test_luo(self):
        values = luo(1000)
        self.assertEqual(values, self.primes_to_1000)

        values = luo(0)
        self.assertEqual(values, [])

        values = luo(-1)
        self.assertEqual(values, [])

        value = sum(luo(2*10**5))
        self.assertEqual(value, self.sum_of_primes_to_200000)

    def test_moebius_xrange(self):
        values = [
            (1, 1), (2, -1), (3, -1), (4, 0), (5, -1), (6, 1), (7, -1), (8, 0), (9, 0), (10, 1), (11, -1), (12, 0),
            (13, -1), (14, 1), (15, 1), (16, 0), (17, -1), (18, 0), (19, -1), (20, 0), (21, 1), (22, 1), (23, -1),
            (24, 0), (25, 0), (26, 1), (27, 0), (28, 0), (29, -1), (30, -1), (31, -1), (32, 0), (33, 1), (34, 1),
            (35, 1), (36, 0), (37, -1), (38, 1), (39, 1), (40, 0), (41, -1), (42, -1), (43, -1), (44, 0), (45, 0),
            (46, 1), (47, -1), (48, 0), (49, 0), (50, 0), (51, 1), (52, 0), (53, -1), (54, 0), (55, 1), (56, 0),
            (57, 1), (58, 1), (59, -1), (60, 0), (61, -1), (62, 1), (63, 0), (64, 0), (65, 1), (66, -1), (67, -1),
            (68, 0), (69, 1), (70, -1), (71, -1), (72, 0), (73, -1), (74, 1), (75, 0), (76, 0), (77, 1), (78, -1),
            (79, -1), (80, 0), (81, 0), (82, 1), (83, -1), (84, 0), (85, 1), (86, 1), (87, 1), (88, 0), (89, -1),
            (90, 0), (91, 1), (92, 0), (93, 1), (94, 1), (95, 1), (96, 0), (97, -1), (98, 0), (99, 0),
        ]

        self.assertEqual(values, list(moebius_xrange(100)))

    def test_prime_xrange(self):
        primes = [
            1000000000039, 1000000000061, 1000000000063, 1000000000091, 1000000000121, 1000000000163, 1000000000169,
            1000000000177, 1000000000189, 1000000000193, 1000000000211, 1000000000271, 1000000000303, 1000000000331,
            1000000000333, 1000000000339, 1000000000459, 1000000000471, 1000000000537, 1000000000543, 1000000000547,
            1000000000561, 1000000000609, 1000000000661, 1000000000669, 1000000000721, 1000000000751, 1000000000787,
            1000000000789, 1000000000799, 1000000000841, 1000000000903, 1000000000921, 1000000000931, 1000000000933,
            1000000000949, 1000000000997,
        ]

        values = prime_xrange(10**12, 10**12 + 10**3)
        self.assertEqual(list(values), primes)

        values = prime_xrange(10**3)
        self.assertEqual(list(values), self.primes_to_1000)

        self.assertEqual(sum(prime_xrange(2*10**6)), self.sum_of_primes_to_2000000)

    def test_primes(self):
        values = primes(1000)
        self.assertEqual(values, self.primes_to_1000)

        values = primes(0)
        self.assertEqual(values, [])

        values = primes(-1)
        self.assertEqual(values, [])

        value = sum(primes(2*10**6))
        self.assertEqual(value, self.sum_of_primes_to_2000000)

    def test_primes_first_n(self):
        values = primes_first_n(200)
        self.assertEqual(values, self.primes_first_200)

        values = primes_first_n(0)
        self.assertEqual(values, [])

        values = primes_first_n(-1)
        self.assertEqual(values, [])
