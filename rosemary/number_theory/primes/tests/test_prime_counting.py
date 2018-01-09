import unittest

from rosemary.number_theory.primes.prime_counting import (
    legendre,
    lehmer,
    lmo,
    lmo_bit,
    meissel,
    pi_table,
    prime_count,
    prime_sum,
    prime_sum2,
)


class TestCore(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.prime_counts = {
            1000: 168,
            1000000: 78498,
            12345678: 809227,
            123456789: 7027260,
        }

        self.prime_sums = {
            200000: 1709600813,
            2000000: 142913828922,
        }

    def test_legendre(self):
        for (n, pi_n) in self.prime_counts.iteritems():
            value = legendre(n)
            self.assertEqual(value, pi_n)

    def test_lehmer(self):
        for (n, pi_n) in self.prime_counts.iteritems():
            value = lehmer(n)
            self.assertEqual(value, pi_n)

    def test_lmo(self):
        for (n, pi_n) in self.prime_counts.iteritems():
            value = lmo(n)
            self.assertEqual(value, pi_n)

    def test_lmo_bit(self):
        for (n, pi_n) in self.prime_counts.iteritems():
            value = lmo_bit(n)
            self.assertEqual(value, pi_n)

    def test_meissel(self):
        for (n, pi_n) in self.prime_counts.iteritems():
            value = meissel(n)
            self.assertEqual(value, pi_n)
