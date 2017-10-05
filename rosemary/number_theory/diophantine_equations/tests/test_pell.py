import unittest

from rosemary.number_theory.diophantine_equations.pell import (
    pell_fundamental_solution,
    pell_minimal_positive_solutions,
)


class TestPell(unittest.TestCase):
    def test_pell_fundamental_solution(self):
        self.assertEqual(pell_fundamental_solution(61), (1766319049, 226153980))
        self.assertEqual(pell_fundamental_solution(17, -1), (4, 1))
        self.assertEqual(pell_fundamental_solution(13, 1), (649, 180))
        self.assertEqual(pell_fundamental_solution(13, -1), (18, 5))

        self.assertRaisesRegexp(ValueError,
                                "pell_fundamental_solution: Solution nonexistent.",
                                pell_fundamental_solution, 15, -1)

        self.assertRaisesRegexp(ValueError,
                                "pell_fundamental_solution: Must have D > 0 not a perfect square and n*",
                                pell_fundamental_solution, 15, -2)

        self.assertRaisesRegexp(ValueError,
                                "pell_fundamental_solution: Must have D > 0 not a perfect square and n*",
                                pell_fundamental_solution, 16, 1)

    def test_pell_minimal_positive_solutions(self):
        sols = [
            (11, 1), (15, 3), (24, 6), (41, 11), (80, 22), (141, 39), (249, 69),
            (440, 122), (869, 241), (1536, 426), (2715, 753), (4799, 1331)
        ]
        self.assertEqual(pell_minimal_positive_solutions(13, 108), sols)

        sols = [
            (13, 1), (10663, 851), (579160, 46222), (483790960, 38610722),
            (26277068347, 2097138361), (21950079635497, 1751807067011)
        ]
        self.assertEqual(pell_minimal_positive_solutions(157, 12), sols)

        sols = [(12, 3), (40, 11), (220, 61), (768, 213)]
        self.assertEqual(pell_minimal_positive_solutions(13, 27), sols)
