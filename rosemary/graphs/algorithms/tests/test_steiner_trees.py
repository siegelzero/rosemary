import unittest

from rosemary.graphs.graphs import Graph
from rosemary.graphs.algorithms.steiner_trees import (
    distance_network_heuristic,
    minimum_cost_paths_heuristic,
)


class TestSteinerTrees(unittest.TestCase):
    def setUp(self):
        # Example in Section 5.1.5 in Dasgupta
        graph = Graph()
        graph.add_edges([('u1', 'u2', 1), ('u1', 'v1', 2), ('u1', 'v2', 2),
                         ('u2', 'v3', 4), ('u2', 'v4', 3), ('u2', 'v5', 5),
                         ('u3', 'v1', 2), ('u3', 'v3', 8), ('u4', 'v2', 8),
                         ('u4', 'v5', 8), ('v1', 'v2', 8), ('v1', 'v3', 9),
                         ('v2', 'v5', 5), ('v3', 'v4', 8)])

        self.graph1 = graph
        self.terminals1 = ['v1', 'v2', 'v3', 'v4', 'v5']
        self.weight1 = 17
        self.edges1 = [('u1', 'u2'), ('u1', 'v1'), ('u1', 'v2'), ('u2', 'v3'), ('u2', 'v4'), ('u2', 'v5')]

        graph = Graph()
        graph.add_edges([(0, 1, 1), (0, 4, 2), (1, 2, 3), (1, 4, 1), (1, 5, 2),
                         (2, 3, 2), (2, 6, 6), (3, 6, 3), (3, 7, 1), (4, 5, 3),
                         (4, 8, 3), (4, 9, 4), (5, 6, 1), (5, 9, 4), (5, 10, 3),
                         (6, 7, 4), (6, 10, 2), (7, 10, 1), (7, 11, 2), (8, 9, 2),
                         (8, 12, 1), (9, 10, 5), (9, 13, 6), (10, 11, 3), (10, 13, 1),
                         (10, 14, 1), (10, 15, 4), (11, 15, 1), (12, 13, 3),
                         (13, 14, 4), (14, 15, 6)])

        self.graph2 = graph
        self.terminals2 = [0, 3, 6, 9, 12, 15]
        self.weight2 = 18
        self.edges2 = [(0, 1), (1, 5), (3, 6), (3, 7), (5, 6), (5, 9), (7, 11), (8, 9), (8, 12), (11, 15)]

        graph = Graph()
        edges = [
            (0, 10, 0.18711581985575576), (1, 9, 0.8355759080627565), (1, 16, 0.7333895522013495),
            (2, 8, 0.2179642234669944), (2, 13, 0.8912036215549415), (2, 19, 0.6431373258117827),
            (2, 24, 0.04038712999015037), (3, 6, 0.12324233297257259), (4, 6, 0.8716054380512974),
            (4, 14, 0.8394000635595052), (5, 13, 0.5808432084428543), (5, 15, 0.3817889531355405),
            (6, 28, 0.7104386125476283), (7, 17, 0.8196386513589949), (7, 22, 0.5415599110939858),
            (7, 26, 0.48022044306517697), (10, 13, 0.2289107839123704), (10, 27, 0.6078612772657713),
            (10, 28, 0.18320583236837396), (11, 14, 0.21670116683074658), (11, 23, 0.814622874490316),
            (11, 24, 0.2741012085352842), (12, 14, 0.020838710732781984), (12, 26, 0.32305190820208796),
            (13, 16, 0.7798216089004076), (14, 24, 0.42823907465971456), (14, 29, 0.8005872492974788),
            (15, 24, 0.5714901649169496), (16, 23, 0.3254825323952686), (17, 19, 0.4371620863350263),
            (17, 20, 0.7305916551727157), (17, 26, 0.20057807143450757), (18, 25, 0.99192719132646),
            (19, 25, 0.4630494211413302), (20, 29, 0.02443067673210797), (21, 22, 0.19217126655249206),
            (21, 29, 0.5188924055605231), (22, 25, 0.7632549163886935), (24, 28, 0.8954891075360596)
        ]
        graph.add_edges(edges)
        self.graph3 = graph
        self.terminals3 = range(5)
        self.weight3 = 4.708933602364698
        self.edges3 = [(0, 10), (1, 16), (2, 13), (3, 6), (4, 6), (6, 28), (10, 13), (10, 28), (13, 16)]

    def test_distance_network_heuristic(self):
        w, edges = distance_network_heuristic(self.graph1, self.terminals1)
        self.assertEqual(w, self.weight1)
        self.assertEqual(edges, self.edges1)

        w, edges = distance_network_heuristic(self.graph2, self.terminals2)
        self.assertEqual(w, self.weight2)
        self.assertEqual(edges, self.edges2)

        w, edges = distance_network_heuristic(self.graph3, self.terminals3)
        self.assertEqual(w, self.weight3)
        self.assertEqual(edges, self.edges3)

    def test_minimum_cost_paths_heuristic(self):
        w, edges = minimum_cost_paths_heuristic(self.graph1, self.terminals1)
        self.assertEqual(w, self.weight1)
        self.assertEqual(edges, self.edges1)

        w, edges = minimum_cost_paths_heuristic(self.graph2, self.terminals2)
        self.assertEqual(w, self.weight2)
        self.assertEqual(edges, self.edges2)

        w, edges = minimum_cost_paths_heuristic(self.graph3, self.terminals3)
        self.assertEqual(w, self.weight3)
        self.assertEqual(edges, self.edges3)

if __name__ == "__main__":
    unittest.main()
