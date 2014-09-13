import unittest

from rosemary.graphs.graphs import Graph
from rosemary.graphs.algorithms.spanning_trees import (
    cheriton_tarjan,
    kruskal,
    prim,
)


class TestSpanningTrees(unittest.TestCase):
    def setUp(self):
        # Example in Section 5.1.5 in Dasgupta
        graph = Graph()
        graph.add_edges([('a', 'b', 5), ('a', 'c', 6), ('a', 'd', 4),
                         ('b', 'c', 1), ('b', 'd', 2), ('c', 'd', 2),
                         ('c', 'e', 5), ('c', 'f', 3), ('d', 'f', 4),
                         ('e', 'f', 4)])
        self.dasgupta_graph = graph
        self.dasgupta_weight = 14
        self.dasgupta_edges = [('a', 'd'), ('b', 'c'), ('b', 'd'), ('c', 'f'), ('e', 'f')]
        self.dasgupta_dict = {'a': None, 'b': 'd', 'c': 'b', 'd': 'a', 'e': 'f', 'f': 'c'}

    def test_prim(self):
        w, edges = prim(self.dasgupta_graph, 'a')
        self.assertEqual(w, self.dasgupta_weight)
        self.assertEqual(edges, self.dasgupta_edges)

        w, pred = prim(self.dasgupta_graph, 'a', edge_list=False)
        self.assertEqual(w, self.dasgupta_weight)
        self.assertEqual(pred, self.dasgupta_dict)

    def test_kruskal(self):
        w, edges = kruskal(self.dasgupta_graph)
        self.assertEqual(w, self.dasgupta_weight)
        self.assertEqual(edges, self.dasgupta_edges)

    def test_cheriton_tarjan(self):
        w, edges = cheriton_tarjan(self.dasgupta_graph)
        self.assertEqual(w, self.dasgupta_weight)
        self.assertEqual(edges, self.dasgupta_edges)

if __name__ == "__main__":
    unittest.main()
