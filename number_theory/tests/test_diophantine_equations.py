import unittest

from rosemary.number_theory.diophantine_equations import (
    pell_fundamental_solution,
)


class TestCore(unittest.TestCase):
    def test_pell_fundamental_solution(self):
        self.assertEqual(pell_fundamental_solution(61), (1766319049, 226153980))
        self.assertEqual(pell_fundamental_solution(17, -1), (4, 1))

        self.assertRaisesRegexp(ValueError,
                                "pell_fundamental_solution: Solution nonexistent.",
                                pell_fundamental_solution, 15, -1)

        self.assertRaisesRegexp(ValueError,
                                "pell_fundamental_solution: Must have D > 0 not a perfect square and n*",
                                pell_fundamental_solution, 15, -2)
