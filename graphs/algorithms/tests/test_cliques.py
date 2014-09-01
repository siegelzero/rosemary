import unittest

from rosemary.graphs.graphs import (
    Graph,
    petersen_graph,
    random_graph,
    coprime_pairs_graph,
)

from rosemary.graphs.algorithms.cliques import (
    bron_kerbosch,
    bron_kerbosch_binary,
    tomita,
    maximal_cliques,
    pardalos,
    ostergard,
    maximum_clique,
    maximum_weight_clique,
)


class TestCliques(unittest.TestCase):
    def setUp(self):
        self.petersen_graph = petersen_graph()
        self.random_graph10_5 = random_graph(10, 0.5)
        self.random_graph100_3 = random_graph(100, 0.3)

        self.petersen_cliques = [
            [0, 1], [0, 4], [0, 5], [1, 2], [1, 6],
            [2, 3], [2, 7], [3, 4], [3, 8], [4, 9],
            [5, 7], [5, 8], [6, 8], [6, 9], [7, 9],
        ]

        self.random_graph_edges = [
            (0, 2), (0, 3), (0, 5), (0, 7),
            (0, 9), (1, 2), (1, 3), (1, 5),
            (1, 6), (1, 7), (1, 9), (2, 3),
            (2, 5), (2, 6), (2, 7), (2, 8),
            (3, 5), (3, 7), (4, 6), (4, 7),
            (5, 6), (5, 7), (5, 9), (6, 7),
            (6, 8), (6, 9), (7, 9),
        ]

        self.random_graph = Graph()
        self.random_graph.add_edges(self.random_graph_edges)
        self.random_graph_clique = [0, 2, 3, 5, 7]

        self.coprime_graph10 = coprime_pairs_graph(10)
        self.coprime_graph10_weight = 30
        self.coprime_graph10_clique = [1, 5, 7, 8, 9]

        self.coprime_graph100 = coprime_pairs_graph(100)
        self.coprime_graph100_weight = 1356
        self.coprime_graph100_clique = [
            1, 17, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
            67, 71, 73, 79, 81, 83, 88, 89, 91, 95, 97,
        ]

    def test_maximal_cliques(self):
        print "Testing Bron-Kerbosh algorithm"
        cliques = [sorted(e) for e in bron_kerbosch(self.petersen_graph)]
        self.assertEqual(sorted(cliques), self.petersen_cliques)

        print "Testing binary Bron-Kerbosh algorithm"
        cliques = [sorted(e) for e in bron_kerbosch_binary(self.petersen_graph)]
        self.assertEqual(sorted(cliques), self.petersen_cliques)

        print "Testing Tomita algorithm"
        cliques = [sorted(e) for e in tomita(self.petersen_graph)]
        self.assertEqual(sorted(cliques), self.petersen_cliques)

        print "Testing generic maximal cliques algorithm"
        for algorithm in ['tomita', 'bron_kerbosch', 'bron_kerbosch_binary']:
            cliques = [sorted(e) for e in maximal_cliques(self.petersen_graph, algorithm=algorithm)]
            self.assertEqual(sorted(cliques), self.petersen_cliques)

        print "Testing with random graphs"
        cliques1 = list(bron_kerbosch(self.random_graph10_5))
        cliques2 = list(bron_kerbosch_binary(self.random_graph10_5))
        cliques3 = list(tomita(self.random_graph10_5))

        self.assertEqual(len(cliques1), len(cliques2))
        self.assertEqual(len(cliques2), len(cliques3))

        cliques1 = list(bron_kerbosch(self.random_graph100_3))
        cliques2 = list(bron_kerbosch_binary(self.random_graph100_3))
        cliques3 = list(tomita(self.random_graph100_3))

        self.assertEqual(len(cliques1), len(cliques2))
        self.assertEqual(len(cliques2), len(cliques3))

    def test_maximum_clique(self):
        print "\nTesting Ostergard algorithm"
        size, clique = ostergard(self.random_graph)
        self.assertEqual(size, 5)
        self.assertEqual(sorted(clique), self.random_graph_clique)

        print "Testing Pardalos algorithm"
        size, clique = pardalos(self.random_graph)
        self.assertEqual(size, 5)
        self.assertEqual(sorted(clique), self.random_graph_clique)

        print "Testing generic maximum-clique algorithm"
        for algorithm in ['pardalos', 'ostergard']:
            size, clique = maximum_clique(self.random_graph, algorithm=algorithm)
            self.assertEqual(size, 5)
            self.assertEqual(sorted(clique), self.random_graph_clique)

    def test_maximum_weight_clique(self):
        print "\nTesting Maximum-weight clique algorithm"
        size, clique = maximum_weight_clique(self.coprime_graph10)
        self.assertEqual(size, self.coprime_graph10_weight)
        self.assertEqual(sorted(clique), self.coprime_graph10_clique)

        size, clique = maximum_weight_clique(self.coprime_graph100)
        self.assertEqual(size, self.coprime_graph100_weight)
        self.assertEqual(sorted(clique), self.coprime_graph100_clique)


if __name__ == "__main__":
    unittest.main()
