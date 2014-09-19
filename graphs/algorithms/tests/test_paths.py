import unittest

from rosemary.graphs.graphs import Graph
from rosemary.graphs.algorithms.paths import (
    bellman_ford,
    dijkstra,
    dijkstra2,
)


class TestShortestPaths(unittest.TestCase):
    def setUp(self):
        graph = Graph()
        graph.add_edges([('a', 'b', 2), ('a', 'c', 1), ('b', 'c', 1),
                         ('b', 'd', 2), ('b', 'e', 3), ('c', 'e', 4),
                         ('d', 'e', 2)])
        self.graph = graph

        self.distance = {'a': 0, 'b': 2, 'c': 1, 'd': 4, 'e': 5}
        self.previous = {'a': None, 'b': 'a', 'c': 'a', 'd': 'b', 'e': 'c'}

    def test_bellman_ford(self):
        distance, previous = bellman_ford(self.graph, 'a')
        self.assertEqual(distance, self.distance)
        self.assertEqual(previous, self.previous)

    def test_dijkstra(self):
        distance, previous = dijkstra(self.graph, 'a')
        self.assertEqual(distance, self.distance)
        self.assertEqual(previous, self.previous)

    def test_dijkstra2(self):
        distance, previous = dijkstra2(self.graph, 'a')
        self.assertEqual(distance, self.distance)
        self.assertEqual(previous, self.previous)


if __name__ == "__main__":
    unittest.main()
