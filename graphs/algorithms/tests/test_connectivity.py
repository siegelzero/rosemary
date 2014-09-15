import unittest

from rosemary.graphs.graphs import Graph
from rosemary.graphs.algorithms.connectivity import (
    connected_component,
    connected_components,
)


class TestSpanningTrees(unittest.TestCase):
    def setUp(self):
        # Graph in Section 3.2.2 of Dasgupta
        self.graph = Graph()
        self.graph.add_vertices(['a', 'b', 'c', 'd', 'e', 'f',
                                 'g', 'h', 'i', 'j', 'k', 'l'])

        self.graph.add_edges([('a', 'b'), ('a', 'e'), ('e', 'i'), ('e', 'j'),
                              ('i', 'j'), ('c', 'd'), ('c', 'g'), ('c', 'h'),
                              ('d', 'h'), ('g', 'h'), ('g', 'k'), ('h', 'k'),
                              ('h', 'l')])

        self.graph_components = [
            ['a', 'b', 'e', 'i', 'j'],
            ['c', 'h', 'd', 'g', 'k', 'l'],
            ['f']
        ]

        self.graph_component = ['a', 'b', 'e', 'i', 'j']

    def test_connected_component(self):
        component = connected_component(self.graph, 'a')
        self.assertEqual(component, self.graph_component)

    def test_connected_components(self):
        components = connected_components(self.graph)
        self.assertEqual(components, self.graph_components)


if __name__ == "__main__":
    unittest.main()
