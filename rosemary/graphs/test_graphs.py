import unittest

from rosemary.graphs.graphs import (
    Graph,
    petersen_graph,
    random_graph,
    coprime_pairs_graph,
)


class TestGraphs(unittest.TestCase):
    def test_vertices(self):
        graph = random_graph(10, 0.5)

        self.assertTrue(graph.vertices() == range(10))
        self.assertTrue(graph.vertex_set() == set(range(10)))

        self.assertTrue(graph.num_vertices() == 10)
        self.assertTrue(len(graph.graph_dict.keys()) == 10)

        new_graph = graph.remove_vertex(0)
        self.assertTrue(new_graph.num_vertices() == 9)
        self.assertTrue(len(new_graph.graph_dict.keys()) == 9)

        self.assertTrue(graph.num_vertices() == 10)
        self.assertTrue(len(graph.graph_dict.keys()) == 10)

        new_graph = graph.remove_vertices([0, 1, 2])
        self.assertTrue(new_graph.num_vertices() == 7)
        self.assertTrue(len(new_graph.graph_dict.keys()) == 7)

        graph.delete_vertex(0)
        self.assertTrue(graph.num_vertices() == 9)
        self.assertTrue(len(graph.graph_dict.keys()) == 9)


if __name__ == "__main__":
    unittest.main()
