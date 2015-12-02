import unittest

from rosemary.number_theory.factorization import (
    _cfrac_aq_pairs,
    _cfrac_multiplier,
    _z2_gaussian_elimination,
    cfrac,
    divisors,
    factor,
    factor_back,
    fermat,
    is_squarefree,
    lehman,
    pollard_p_minus_1,
    pollard_rho,
    pollard_rho_brent,
    prime_divisors,
    smooth_factor,
    trial_division,
    xdivisors,
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

    def test_divisors(self):
        return

    def test_factor(self):
        factorizations = {
            5**77 - 1: [
                (2, 2),
                (19531, 1),
                (12207031L, 1),
                (527093491L, 1),
                (8090594434231L, 1),
                (162715052426691233701L, 1)
            ],
            5**93 + 1: [
                (2, 1),
                (3, 2),
                (7, 1),
                (1303, 1),
                (21207101L, 1),
                (258065887L, 1),
                (28086211607L, 1),
                (75005167927L, 1),
                (53345671490722200466369L, 1)
            ],
            1303**6*11213**4: [
                (1303, 6),
                (11213, 4)
            ],
            2**91 - 1: [
                (127, 1),
                (911, 1),
                (8191, 1),
                (112901153L, 1),
                (23140471537L, 1)
            ],
            10**22 + 1: [
                (89, 1),
                (101, 1),
                (1052788969L, 1),
                (1056689261L, 1)
            ]
        }

        for (n, n_fac) in factorizations.items():
            self.assertEqual(n_fac, factor(n))

    def test_factor_back(self):
        return

    def test_fermat(self):
        factorizations = {
            100: 2,
            1112470797641561909L: 1052788969L,
        }
        for (n, p) in factorizations.items():
            self.assertEqual(fermat(n), p)

    def test_is_squarefree(self):
        return

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
            414876509374946224887319096986555104270201L: 711628063L,
            615028785027771754392800840284305067659311379685188908444467L: 21207101L,
        }
        for (n, p) in factorizations.items():
            self.assertEqual(pollard_rho(n), p)

    def test_pollard_rho_brent(self):
        factorizations = {
            414876509374946224887319096986555104270201L: 711628063L,
            615028785027771754392800840284305067659311379685188908444467L: 21207101L,
        }
        for (n, p) in factorizations.items():
            self.assertEqual(pollard_rho(n), p)

    def test_prime_divisors(self):
        return

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

    def test_xdivisors(self):
        return
