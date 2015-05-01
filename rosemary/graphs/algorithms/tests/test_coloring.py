import unittest

from rosemary.graphs.graphs import (
    coprime_pairs_graph,
    petersen_graph,
    random_graph,
    Graph,
)

from rosemary.graphs.algorithms.coloring import (
    branch_and_bound,
    dsatur,
    greedy_sequential,
    korman,
    maxis,
)


class TestColoring(unittest.TestCase):
    def setUp(self):
        self.petersen_graph = petersen_graph()
        self.petersen_colors = 3

        self.coprime_graph30 = coprime_pairs_graph(30)
        self.coprime_graph30_colors = 11

        graph = Graph()
        graph.add_edges([(0, 2), (0, 4), (0, 5), (0, 8), (1, 3), (1, 4), (2, 6), (2, 8), (3, 4), (3, 7), (4, 6), (4, 7),
                         (4, 9), (5, 8)])
        self.small_random = graph
        self.small_colors = 3
        self.small_classes = [[4, 8], [0, 3, 6, 9], [1, 2, 5, 7]]

        graph = Graph()
        graph.add_edges([(0, 4), (0, 6), (1, 2), (1, 4), (1, 8), (1, 9), (2, 3), (2, 4), (3, 5), (3, 7), (4, 6), (4, 8),
                         (4, 9), (5, 7), (5, 9), (6, 8), (7, 8)])
        self.small_random2 = graph
        self.small_colors2 = 4
        self.small_classes2 = [[4, 7], [0, 2, 8, 9], [1, 5, 6], [3]]

        graph = Graph()
        graph.add_edges([(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 7), (0, 8), (0, 9), (0, 10), (0, 11), (0, 12),
                         (0, 13), (0, 14), (0, 15), (0, 16), (0, 17), (0, 18), (0, 19), (1, 2), (1, 3), (1, 4), (1, 5),
                         (1, 6), (1, 7), (1, 8), (1, 9), (1, 10), (1, 12), (1, 14), (1, 15), (1, 17), (1, 18), (1, 19),
                         (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (2, 9), (2, 11), (2, 12), (2, 13), (2, 14), (2, 15),
                         (2, 17), (2, 18), (2, 19), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (3, 10), (3, 11),
                         (3, 12), (3, 14), (3, 15), (3, 16), (3, 17), (3, 18), (4, 5), (4, 6), (4, 7), (4, 9), (4, 10),
                         (4, 11), (4, 13), (4, 14), (4, 15), (4, 16), (4, 17), (4, 18), (4, 19), (5, 6), (5, 7), (5, 9),
                         (5, 10), (5, 13), (5, 14), (5, 15), (5, 17), (5, 19), (6, 8), (6, 9), (6, 10), (6, 11),
                         (6, 12), (6, 13), (6, 14), (6, 15), (6, 16), (6, 17), (6, 19), (7, 8), (7, 9), (7, 12),
                         (7, 15), (7, 16), (7, 17), (7, 19), (8, 9), (8, 10), (8, 11), (8, 12), (8, 13), (8, 14),
                         (8, 16), (8, 17), (8, 18), (8, 19), (9, 10), (9, 11), (9, 14), (9, 17), (9, 18), (9, 19),
                         (10, 11), (10, 12), (10, 13), (10, 15), (10, 16), (10, 17), (10, 18), (10, 19), (11, 12),
                         (11, 13), (11, 14), (11, 15), (11, 16), (11, 17), (11, 18), (11, 19), (12, 13), (12, 14),
                         (12, 15), (12, 16), (12, 17), (12, 18), (12, 19), (13, 14), (13, 15), (13, 16), (13, 17),
                         (13, 18), (13, 19), (14, 16), (14, 17), (14, 18), (15, 17), (15, 19), (16, 18), (16, 19),
                         (17, 18), (18, 19)])
        self.big_random = graph
        self.big_colors = 10
        self.big_classes = [[0, 6], [17, 19], [4, 8], [5, 11], [7, 10, 14], [2, 3], [15, 18], [9, 12], [1, 16], [13]]

    def test_korman(self):
        (num_colors, classes) = korman(self.petersen_graph)
        self.assertEqual(num_colors, self.petersen_colors)
        self.assertTrue(is_valid_color_classes(self.petersen_graph, classes))

        (num_colors, classes) = korman(self.small_random)
        self.assertEqual(num_colors, self.small_colors)
        self.assertEqual(classes, self.small_classes)

        (num_colors, classes) = korman(self.small_random2)
        self.assertEqual(num_colors, self.small_colors2)
        self.assertEqual(classes, self.small_classes2)

        (num_colors, classes) = korman(self.big_random)
        self.assertEqual(num_colors, self.big_colors)
        self.assertEqual(classes, self.big_classes)

        (num_colors, classes) = korman(self.coprime_graph30)
        self.assertEqual(num_colors, self.coprime_graph30_colors)
        self.assertTrue(is_valid_color_classes(self.coprime_graph30, classes))

        (num_colors, color_map) = korman(self.coprime_graph30, classes=False)
        self.assertEqual(color_map, {1: 7, 2: 1, 3: 2, 4: 1, 5: 3, 6: 1, 7: 4, 8: 1, 9: 2, 10: 1, 11: 5, 12: 1, 13: 6,
                                     14: 1, 15: 2, 16: 1, 17: 8, 18: 1, 19: 9, 20: 1, 21: 2, 22: 1, 23: 10, 24: 1, 25:
                                     3, 26: 1, 27: 2, 28: 1, 29: 11, 30: 1})

        for i in xrange(10):
            graph = random_graph(20, 0.5)
            k_num, classes = korman(graph)
            b_num, classes = branch_and_bound(graph)
            self.assertEqual(b_num, k_num)
            self.assertTrue(is_valid_color_classes(graph, classes))

    test_exact = test_korman

    def test_branch_and_bound(self):
        print "Testing branch and bound coloring algorithm"
        (num_colors, classes) = branch_and_bound(self.petersen_graph)
        self.assertEqual(num_colors, self.petersen_colors)
        self.assertTrue(is_valid_color_classes(self.petersen_graph, classes))

        (num_colors, classes) = branch_and_bound(self.coprime_graph30)
        self.assertEqual(num_colors, self.coprime_graph30_colors)
        self.assertTrue(is_valid_color_classes(self.coprime_graph30, classes))

        (num_colors, color_map) = branch_and_bound(self.coprime_graph30, classes=False)
        self.assertEqual(color_map, {1: 5, 2: 11, 3: 10, 4: 11, 5: 9, 6: 10, 7: 8, 8: 11, 9: 10, 10: 9, 11: 7, 12: 10,
                                     13: 6, 14: 8, 15: 9, 16: 11, 17: 4, 18: 10, 19: 3, 20: 9, 21: 8, 22: 7, 23: 2, 24:
                                     10, 25: 9, 26: 6, 27: 10,
                                     28: 8, 29: 1, 30: 9})

        (num_colors, classes) = branch_and_bound(self.small_random)
        self.assertEqual(num_colors, self.small_colors)
        self.assertEqual(classes, self.small_classes)

        (num_colors, classes) = branch_and_bound(self.small_random2)
        self.assertEqual(num_colors, self.small_colors2)
        self.assertEqual(classes, self.small_classes2)

        (num_colors, classes) = branch_and_bound(self.big_random)
        self.assertEqual(num_colors, self.big_colors)
        self.assertEqual(classes, self.big_classes)

        graph = random_graph(30, 0.5)
        (num_colors, classes) = branch_and_bound(graph)
        self.assertTrue(is_valid_color_classes(graph, classes))

        for i in xrange(10):
            graph = random_graph(10, 0.8)
            _, classes = branch_and_bound(graph)
            self.assertTrue(is_valid_color_classes(graph, classes))

    def test_greedy_sequential(self):
        print "Testing greedy sequential coloring algorithm"
        (num_colors, classes) = greedy_sequential(self.petersen_graph)
        self.assertTrue(num_colors == self.petersen_colors)
        self.assertTrue(is_valid_color_classes(self.petersen_graph, classes))

        (num_colors, classes) = greedy_sequential(self.petersen_graph, classes=False)
        self.assertTrue(num_colors == self.petersen_colors)
        self.assertTrue(classes, {0: 0, 1: 2, 2: 0, 3: 2, 4: 1, 5: 2, 6: 1, 7: 1, 8: 0, 9: 0})

        (num_colors, classes) = greedy_sequential(self.coprime_graph30, passes=2)
        self.assertTrue(num_colors <= self.coprime_graph30_colors)
        self.assertTrue(is_valid_color_classes(self.coprime_graph30, classes))

        for i in xrange(10):
            graph = random_graph(10, 0.8)
            _, classes = greedy_sequential(graph)
            self.assertTrue(is_valid_color_classes(graph, classes))

    def test_dsatur(self):
        print "Testing DSATUR coloring algorithm"
        (num_colors, classes) = dsatur(self.petersen_graph)
        self.assertTrue(num_colors <= self.petersen_colors)
        self.assertTrue(is_valid_color_classes(self.petersen_graph, classes))

        (num_colors, color_map) = dsatur(self.petersen_graph, classes=False)
        self.assertEqual(num_colors, self.petersen_colors)
        self.assertEqual(color_map, {0: 0, 1: 1, 2: 2, 3: 1, 4: 2, 5: 1, 6: 2, 7: 0, 8: 0, 9: 1})

        (num_colors, classes) = dsatur(self.coprime_graph30)
        self.assertTrue(num_colors <= self.coprime_graph30_colors)
        self.assertTrue(is_valid_color_classes(self.coprime_graph30, classes))

        for i in xrange(10):
            graph = random_graph(10, 0.8)
            _, classes = dsatur(graph)
            self.assertTrue(is_valid_color_classes(graph, classes))

    def test_maxis(self):
        print "Testing MAXIS coloring algorithm"

        (num_colors, color_map) = maxis(self.petersen_graph, classes=False)
        self.assertEqual(num_colors, self.petersen_colors)
        self.assertEqual(color_map, {0: 0, 1: 1, 2: 2, 3: 1, 4: 2, 5: 1, 6: 2, 7: 0, 8: 0, 9: 1})

        graph = random_graph(50, 0.5)
        (num_colors, classes) = maxis(graph, color_limit=30)
        self.assertTrue(is_valid_color_classes(graph, classes))

        graph = random_graph(60, 0.5)
        (num_colors, classes) = maxis(graph, use_greedy=False, color_limit=30, mis_limit=30, verbose=True)
        self.assertTrue(is_valid_color_classes(graph, classes))

        graph = random_graph(50, 0.5)
        (num_colors, classes) = maxis(graph, mis_limit=50, color_limit=30)
        self.assertTrue(is_valid_color_classes(graph, classes))

        graph = random_graph(50, 0.5)
        (num_colors, classes) = maxis(graph, use_greedy=True, color_limit=30, mis_limit=40, verbose=True)
        self.assertTrue(is_valid_color_classes(graph, classes))

        graph = random_graph(50, 0.5)
        (num_colors, classes) = maxis(graph, use_greedy=True, color_limit=0, mis_limit=40, verbose=True)
        self.assertTrue(is_valid_color_classes(graph, classes))

        graph = random_graph(50, 0.5)
        (num_colors, classes) = maxis(graph, use_greedy=True, color_limit=30)
        self.assertTrue(is_valid_color_classes(graph, classes))

        for i in xrange(10):
            graph = random_graph(10, 0.8)
            _, classes = maxis(graph, verbose=True)
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
