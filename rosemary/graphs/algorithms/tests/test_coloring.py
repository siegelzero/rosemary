import unittest

from rosemary.graphs.graphs import (
    coprime_pairs_graph,
    petersen_graph,
    random_graph
)

from rosemary.graphs.algorithms.coloring import (
    branch_and_bound,
    dsatur,
    greedy_sequential,
    maxis,
)


class TestColoring(unittest.TestCase):
    def setUp(self):
        self.petersen_graph = petersen_graph()
        self.petersen_colors = 3

        self.coprime_graph30 = coprime_pairs_graph(30)
        self.coprime_graph30_colors = 11

    def test_branch_and_bound(self):
        print "Testing branch and bound coloring algorithm"
        (num_colors, classes) = branch_and_bound(self.petersen_graph)
        self.assertEqual(num_colors, self.petersen_colors)
        self.assertTrue(is_valid_color_classes(self.petersen_graph, classes))

        (num_colors, classes) = branch_and_bound(self.coprime_graph30)
        self.assertEqual(num_colors, self.coprime_graph30_colors)
        self.assertTrue(is_valid_color_classes(self.coprime_graph30, classes))

        for i in xrange(10):
            graph = random_graph(10, 0.8)
            _, classes = branch_and_bound(graph)
            self.assertTrue(is_valid_color_classes(graph, classes))

    def test_greedy_sequential(self):
        print "Testing greedy sequential coloring algorithm"
        (num_colors, classes) = greedy_sequential(self.petersen_graph)
        self.assertTrue(num_colors <= self.petersen_colors)
        self.assertTrue(is_valid_color_classes(self.petersen_graph, classes))

        (num_colors, classes) = greedy_sequential(self.coprime_graph30, passes=2)
        self.assertTrue(num_colors <= self.coprime_graph30_colors)
        self.assertTrue(is_valid_color_classes(self.coprime_graph30, classes))

        for i in xrange(10):
            graph = random_graph(10, 0.8)
            _, classes = branch_and_bound(graph)
            self.assertTrue(is_valid_color_classes(graph, classes))

    def test_dsatur(self):
        print "Testing DSATUR coloring algorithm"
        (num_colors, classes) = dsatur(self.petersen_graph)
        self.assertTrue(num_colors <= self.petersen_colors)
        self.assertTrue(is_valid_color_classes(self.petersen_graph, classes))

        (num_colors, classes) = dsatur(self.coprime_graph30)
        self.assertTrue(num_colors <= self.coprime_graph30_colors)
        self.assertTrue(is_valid_color_classes(self.coprime_graph30, classes))

        for i in xrange(10):
            graph = random_graph(10, 0.8)
            _, classes = branch_and_bound(graph)
            self.assertTrue(is_valid_color_classes(graph, classes))

    def test_maxis(self):
        print "Testing MAXIS coloring algorithm"
        graph = random_graph(50, 0.5)
        (num_colors, classes) = maxis(graph, color_limit=30)
        self.assertTrue(is_valid_color_classes(graph, classes))

        graph = random_graph(50, 0.5)
        (num_colors, classes) = maxis(graph, mis_limit=50, color_limit=30)
        self.assertTrue(is_valid_color_classes(graph, classes))

        graph = random_graph(50, 0.5)
        (num_colors, classes) = maxis(graph, greedy=True, color_limit=30)
        self.assertTrue(is_valid_color_classes(graph, classes))

        for i in xrange(10):
            graph = random_graph(10, 0.8)
            _, classes = maxis(graph)
            self.assertTrue(is_valid_color_classes(graph, classes))


def is_valid_color_classes(graph, classes):
    for color_class in classes:
        for i in xrange(len(color_class)):
            u = color_class[i]
            for j in xrange(i):
                v = color_class[j]
                if u in graph[v]:
                    return False

    return True


if __name__ == "__main__":
    unittest.main()
