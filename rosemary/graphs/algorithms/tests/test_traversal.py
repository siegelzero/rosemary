import unittest

from rosemary.graphs.graphs import Graph
from rosemary.graphs.algorithms.traversal import (
    breadth_first_search,
    breadth_first_search_tree,
    depth_first_search,
    depth_first_search_tree
)


class TestSpanningTrees(unittest.TestCase):
    def setUp(self):
        self.graph = Graph()
        self.graph.add_edges([
            ('a', 'b'), ('a', 's'), ('b', 'c'), ('c', 's'),
            ('d', 'e'), ('d', 's'), ('e', 'd'), ('e', 's')
        ])

    def test_breadth_first_search(self):
        vertices = list(breadth_first_search(self.graph, 's'))
        self.assertEqual(vertices, ['s', 'a', 'c', 'e', 'd', 'b'])

        vertices = list(breadth_first_search(self.graph, 's', 1))
        self.assertEqual(vertices, ['s', 'a', 'c', 'e', 'd'])

    def test_breadth_first_search_tree(self):
        pred = breadth_first_search_tree(self.graph, 's')
        ans = {'a': 's', 'b': 'a', 'c': 's', 'd': 's', 'e': 's', 's': None}
        self.assertEqual(pred, ans)

    def test_depth_first_search(self):
        vertices = list(depth_first_search(self.graph, 's'))
        self.assertEqual(vertices, ['s', 'd', 'e', 'c', 'b', 'a'])

        vertices = list(depth_first_search(self.graph, 's', 1))
        self.assertEqual(vertices, ['s', 'd', 'e', 'c', 'a'])

    def test_depth_first_search_tree(self):
        pred = depth_first_search_tree(self.graph, 's')
        ans = {'a': 'b', 'b': 'c', 'c': 's', 'd': 's', 'e': 'd', 's': None}
        self.assertEqual(pred, ans)

if __name__ == "__main__":
    unittest.main()
