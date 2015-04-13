import unittest

from rosemary.graphs.graphs import Graph
from rosemary.graphs.algorithms.paths import (
    bellman_ford,
    dijkstra,
    dijkstra_bidirectional,
    dijkstra_buckets,
    dijkstra_iterator,
    dijkstra_pairing_heap,
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

        graph = Graph()
        edges = [('a', 'd', 4), ('a', 'e', 1), ('a', 'h', 10), ('b', 'c', 2), ('b', 'f', 1), ('c', 'f', 3),
                 ('d', 'h', 1), ('e', 'f', 3), ('e', 'h', 5), ('f', 'g', 7), ('f', 'i', 1), ('g', 'j', 1),
                 ('h', 'i', 9), ('i', 'j', 2)]

        graph.add_edges(edges)

        self.graph2 = graph
        self.distance2 = {'a': 4, 'c': 11, 'b': 9, 'e': 5, 'd': 0, 'g': 12, 'f': 8, 'i': 9, 'h': 1, 'j': 11}
        self.previous2 = {'a': 'd', 'c': 'f', 'b': 'f', 'e': 'a', 'd': None,
                          'g': 'j', 'f': 'e', 'i': 'f', 'h': 'd', 'j': 'i'}

    def test_bellman_ford(self):
        distance, previous = bellman_ford(self.graph, 'a')
        self.assertEqual(distance, self.distance)
        self.assertEqual(previous, self.previous)

        self.assertRaisesRegexp(ValueError,
                                "bellman_ford: 0 is not a vertex of graph.",
                                bellman_ford, self.graph, 0)

        self.graph['a']['b'] = -10
        self.assertRaisesRegexp(ValueError,
                                'bellman_ford: Negative cycle detected.',
                                bellman_ford, self.graph, 'a')

        distance, previous = bellman_ford(self.graph2, 'd')
        self.assertEqual(distance, self.distance2)
        self.assertEqual(previous, self.previous2)

    def test_dijkstra(self):
        distance, previous = dijkstra(self.graph, 'a')
        self.assertEqual(distance, self.distance)
        self.assertEqual(previous, self.previous)

        self.assertRaisesRegexp(ValueError,
                                "dijkstra: 0 is not a vertex of graph.",
                                dijkstra, self.graph, 0)

        distance, previous = dijkstra(self.graph2, 'd')
        self.assertEqual(distance, self.distance2)
        self.assertEqual(previous, self.previous2)

    def test_dijkstra_bidirectional(self):
        self.assertRaisesRegexp(ValueError,
                                "dijkstra_bidirectional: 0 is not a vertex of graph.",
                                dijkstra_bidirectional, self.graph, 0, 'a')

        self.assertRaisesRegexp(ValueError,
                                "dijkstra_bidirectional: 2 is not a vertex of graph.",
                                dijkstra_bidirectional, self.graph, 'a', 2)

        dist, path = dijkstra_bidirectional(self.graph, 'a', 'e')
        ddist, dprev = dijkstra(self.graph, 'a')
        self.assertEqual(dist, ddist['e'])

        self.graph.add_vertex('f')
        self.assertEqual(dijkstra_bidirectional(self.graph, 'a', 'f'), float('inf'))

    def test_dijkstra_buckets(self):
        distance, previous = dijkstra_buckets(self.graph, 'a')
        self.assertEqual(distance, self.distance)
        self.assertEqual(previous, self.previous)

        self.assertRaisesRegexp(ValueError,
                                "dijkstra_buckets: 0 is not a vertex of graph.",
                                dijkstra_buckets, self.graph, 0)

        distance, previous = dijkstra_buckets(self.graph2, 'd')
        self.assertEqual(distance, self.distance2)
        self.assertEqual(previous, self.previous2)

        self.graph['a']['b'] = 1.4
        self.assertRaisesRegexp(ValueError,
                                "dijkstra_buckets: Weights must be integral.",
                                dijkstra_buckets, self.graph, 'a')

    def test_dijkstra_iterator(self):
        for graph in (self.graph, self.graph2):
            paths = list(dijkstra_iterator(graph, 'a'))
            ddist, dprev = dijkstra(graph, 'a')
            for (node, dist, path) in paths:
                self.assertEqual(dist, ddist[node])

    def test_dijkstra_pairing_heap(self):
        distance, previous = dijkstra_pairing_heap(self.graph, 'a')
        self.assertEqual(distance, self.distance)
        self.assertEqual(previous, self.previous)

        distance, previous = dijkstra_pairing_heap(self.graph2, 'd')
        self.assertEqual(distance, self.distance2)
        self.assertEqual(previous, self.previous2)


if __name__ == "__main__":
    unittest.main()
