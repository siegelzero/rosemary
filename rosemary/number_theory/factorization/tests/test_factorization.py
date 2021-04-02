import unittest

from rosemary.number_theory.factorization.factorization import (
    divisors,
    factor,
    factor_back,
    factored_range,
    prime_divisors,
    xdivisors,
)


class TestCore(unittest.TestCase):
    def test_divisors(self):
        return

    def test_factor(self):
        factorizations = {
            5**77 - 1: [
                (2, 2),
                (19531, 1),
                (12207031, 1),
                (527093491, 1),
                (8090594434231, 1),
                (162715052426691233701, 1)
            ],
            5**93 + 1: [
                (2, 1),
                (3, 2),
                (7, 1),
                (1303, 1),
                (21207101, 1),
                (258065887, 1),
                (28086211607, 1),
                (75005167927, 1),
                (53345671490722200466369, 1)
            ],
            1303**6*11213**4: [
                (1303, 6),
                (11213, 4)
            ],
            2**91 - 1: [
                (127, 1),
                (911, 1),
                (8191, 1),
                (112901153, 1),
                (23140471537, 1)
            ],
            10**22 + 1: [
                (89, 1),
                (101, 1),
                (1052788969, 1),
                (1056689261, 1)
            ]
        }

        for (n, n_fac) in factorizations.items():
            self.assertEqual(n_fac, factor(n))

    def test_factored_range(self):
        values = factored_range(10000)
        for (n, n_fac) in values:
            prod = 1
            for (p, e) in n_fac:
                prod *= p**e
            self.assertEqual(n, prod)

        self.assertEqual(list(factored_range(0, 100)), list(factored_range(100)))
        self.assertEqual(list(factored_range(0)), [])

    def test_factor_back(self):
        return

    def test_prime_divisors(self):
        return

    def test_xdivisors(self):
        return
