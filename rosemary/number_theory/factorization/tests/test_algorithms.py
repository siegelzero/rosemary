import unittest

from rosemary.number_theory.factorization.algorithms import (
    _cfrac_aq_pairs,
    _cfrac_multiplier,
    _z2_gaussian_elimination,
    cfrac,
    fermat,
    lehman,
    pollard_p_minus_1,
    pollard_rho,
    pollard_rho_brent,
    smooth_factor,
    trial_division,
)


class TestCore(unittest.TestCase):
    def test__cfrac_aq_pairs(self):
        n = 13290059
        # This table comes from Table 1 in the Morrison/Brillhart paper.
        table = {
            1: (3645, 4034),
            2: (3646, 3257),
            3: (7291, 1555),
            4: (32810, 1321),
            5: (171341, 2050),
            10: (6700527, 1333),
            22: (5235158, 4633),
            23: (1914221, 226),
            26: (11455708, 3286),
            31: (1895246, 5650),
            40: (3213960, 4558),
            52: (2467124, 25)
        }

        for (i, a, q) in _cfrac_aq_pairs(n):
            if i in table:
                self.assertEqual(table[i], (a, q))
            if i > 52:
                break

    def test__cfrac_multiplier(self):
        self.assertEquals(_cfrac_multiplier(5**77 - 1), 781)
        self.assertEquals(_cfrac_multiplier(2**128 + 1), 273)

    def test__z2_gaussian_elimination(self):
        vectors = [
            [1, 1, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 1, 0, 0],
            [0, 0, 1, 0, 0, 1, 0],
            [1, 1, 0, 0, 0, 1, 0],
            [0, 1, 0, 1, 0, 0, 1],
            [1, 1, 0, 0, 0, 1, 0],
            [0, 1, 0, 0, 1, 0, 1]
        ]
        results = _z2_gaussian_elimination(vectors)
        expected = [[1, 0, 1, 1, 0, 0, 0], [1, 0, 1, 0, 0, 1, 0], [0, 1, 0, 0, 1, 0, 1]]
        self.assertEqual(results, expected)

    def test_cfrac(self):
        factorizations = {
            12007001: 4001,
            1112470797641561909L: 1052788969L,
            2175282241519502424792841: 513741730823L,
        }
        for (n, p) in factorizations.items():
            self.assertEqual(cfrac(n), p)

    def test_fermat(self):
        factorizations = {
            100: 2,
            1112470797641561909L: 1052788969L,
        }
        for (n, p) in factorizations.items():
            self.assertEqual(fermat(n), p)

    def test_lehman(self):
        factorizations = {
            100: 2,
            1112470797641561909L: 1056689261L,
        }
        for (n, p) in factorizations.items():
            self.assertEqual(lehman(n), p)

    def test_pollard_p_minus_1(self):
        factorizations = {
            1112470797641561909L: 1056689261L,
            615028785027771754392800840284305067659311379685188908444467L: 75005167927L,
        }
        for (n, p) in factorizations.items():
            self.assertEqual(pollard_p_minus_1(n), p)

    def test_pollard_rho(self):
        factorizations = {
            414876509374946224887319096986555104270201L: set([711628063L]),
            615028785027771754392800840284305067659311379685188908444467L: set([21207101L, 258065887L])
        }
        for (n, p_set) in factorizations.items():
            self.assertIn(pollard_rho(n), p_set)

    def test_pollard_rho_brent(self):
        factorizations = {
            414876509374946224887319096986555104270201L: set([711628063L]),
            615028785027771754392800840284305067659311379685188908444467L: set([21207101L, 258065887L])
        }
        for (n, p_set) in factorizations.items():
            self.assertIn(pollard_rho_brent(n), p_set)

    def test_smooth_factor(self):
        self.assertEqual(smooth_factor(100, [2, 5]), [2, 2])
        self.assertEqual(smooth_factor(100, [2, 3, 5, 7]), [2, 0, 2, 0])
        self.assertRaisesRegexp(ValueError,
                                "smooth_factor: n does not factor over the given factor base",
                                smooth_factor, 100, [3, 5])

    def test_trial_division(self):
        self.assertEqual(trial_division(100), 2)
        self.assertEqual(trial_division(10000004400000259), 100000007)
        self.assertAlmostEqual(trial_division(10000004400000259, 100), 10000004400000259)
