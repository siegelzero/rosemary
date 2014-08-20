import rosemary.algebra.matrices.matrices
import rosemary.combinatorics.enumeration
import rosemary.graphs.algorithms.traversal

from rosemary.number_theory.core import gcd

import copy
import itertools
import random


class Graph(object):
    """
    Class of simple graphs.
    """
    def __init__(self):
        """
        Initializes a new Graph object.
        """
        self.graph_dict = {}

    def __repr__(self):
        """
        Representation of Graph object.
        """
        num_vertices = len(self.graph_dict)
        return "Graph on {} vertices".format(num_vertices)

    def __getitem__(self, key):
        """
        Gets the item from the graph dict corresponding to key.
        """
        return self.graph_dict[key]

    def __iter__(self):
        """
        Iterator over the keys of the grap dict of self.
        """
        return self.graph_dict.iterkeys()

    def copy(self):
        """
        Returns a copy of self.
        """
        new_graph = Graph()
        new_graph.graph_dict = copy.deepcopy(self.graph_dict)
        return new_graph

    def num_vertices(self):
        """
        Returns the number of vertices of self.

        Input:
            * self: Graph

        Output:
            * num_vertices: int

        Examples:
            >>> G = Graph()
            >>> G.num_vertices()
            0
            >>> G.add_vertices([0, 1, 2])
            >>> G.num_vertices()
            3
        """
        return len(self.graph_dict.keys())

    def num_edges(self):
        """
        Returns the number of edges of self.
        """
        return len(self.edge_set())

    def vertices(self, **kwargs):
        """
        Returns a sorted list of the vertices of self.
        """
        vertices = self.graph_dict.keys()
        if 'order' in kwargs:
            order = kwargs['order']

            if order == 'degree':
                vertex_list = sorted(vertices, key=self.degree)
            elif order == 'induced':
                num_vertices = self.num_vertices()
                new_graph = self.copy()
                vertex_list = []

                while len(vertex_list) < num_vertices:
                    vertex = new_graph.min_degree_vertex()
                    vertex_list.append(vertex)
                    new_graph.delete_vertex(vertex)

            else:
                vertex_list = sorted(vertices)
        else:
            vertex_list = sorted(vertices)

        return vertex_list

    def vertex_set(self):
        """
        Returns a set containing the vertices of self.
        """
        vertex_set = set(self.graph_dict.keys())
        return vertex_set

    def edges(self):
        """
        Returns a sorted list of the edges of self.
        """
        edge_list = sorted(self.edge_set())
        return edge_list

    def edge_set(self):
        """
        Returns a set containing the vertices of self.
        """
        edges_seen = set()
        graph_dict = self.graph_dict

        for u in graph_dict:
            for v in graph_dict[u]:
                weight = graph_dict[u][v]
                triple = (min(u, v), max(u, v), weight)
                edges_seen.add(triple)

        return edges_seen

    def max_degree_vertex(self):
        """
        Returns a vertx of self of maximal degree.
        """
        vertices = self.vertices(order='degree')

        if len(vertices) == 0:
            raise ValueError('min_degree_vertex: Graph has no vertices.')

        return vertices[-1]

    def min_degree_vertex(self):
        """
        Returns a vertx of self of minimal degree.
        """
        vertices = self.vertices(order='degree')

        if len(vertices) == 0:
            raise ValueError('min_degree_vertex: Graph has no vertices.')

        return vertices[0]

    def degree(self, u):
        """
        Returns the degree of the vertex; i.e. the number of vertices adjacent
        to vertex.
        """
        degree = len(self.neighbors(u))
        return degree

    def total_degree(self, vertices):
        """
        Returns the sum of the degrees of all vertices in the list.
        """
        total = sum(self.degree(u) for u in vertices)
        return total

    def density(self):
        """
        Returns the edge density of self.
        """
        num_vertices = self.num_vertices()
        num_edges = self.num_edges()
        total_possible = rosemary.combinatorics.enumeration.binomial(num_vertices, 2)
        density = num_edges / (1.0*total_possible)
        return density

    def neighbors(self, u):
        """
        Returns a set of the neighbors of u; i.e. the vertices adjacent to u.
        """
        neighbors = set(self.graph_dict[u].keys())
        return neighbors

    def neighborhood(self, u):
        """
        Returns the neighbors of u along with u.
        """
        neighbors = self.neighbors(u)
        neighbors.add(u)
        return neighbors

    ############################################################
    # Methods for adding / deleting vertices
    ############################################################

    def add_vertex(self, u):
        """
        Adds the vertex u to self, modifying self. If u in already a vertex of
        self, nothing is done.

        Input:
            * self: Graph
            * u: vertex of self

        Examples:
            >>> G = Graph()
            >>> G.vertices()
            []
            >>> G.add_vertex(0)
            >>> G.vertices()
            [0]
        """
        if u not in self.graph_dict:
            self.graph_dict[u] = {}

    def add_vertices(self, vertex_list):
        """
        Adds the vertices in vertex_list to self, modifying self.. Any vertices
        which are already in self are skipped.

        Input:
            * self: Graph
            * vertex_list: iterable (list, set, tuple)
                List of vertices to add to self.

        Examples:
            >>> G = Graph()
            >>> G.add_vertices([0, 1, 2])
            >>> G.vertices()
            [0, 1, 2]
        """
        graph_dict = self.graph_dict

        for u in vertex_list:
            if u not in graph_dict:
                graph_dict[u] = {}

    def delete_vertex(self, u):
        """
        Deletes the vertex u along with any edges adjacent to u, modifying self.
        If u is not in self, then nothing is done.
        """
        graph_dict = self.graph_dict
        if u in graph_dict:
            del graph_dict[u]
            for v in graph_dict:
                if u in graph_dict[v]:
                    del graph_dict[v][u]

    def delete_vertices(self, vertices):
        """
        Deletes vertices from self, along with any adjacent edges. Any vertices
        which are not in self are skipped.
        """
        graph_dict = self.graph_dict

        for u in vertices:
            if u in graph_dict:
                del graph_dict[u]
                for v in graph_dict:
                    if u in graph_dict[v]:
                        del graph_dict[v][u]

    def remove_vertex(self, u):
        """
        Returns a copy of self with vertex u removed.
        """
        new_graph = self.copy()
        new_graph.delete_vertex(u)
        return new_graph

    def remove_vertices(self, vertices):
        """
        Returns a copy of self with vertex u removed.
        """
        new_graph = self.copy()
        new_graph.delete_vertices(vertices)
        return new_graph

    def add_edge(self, u, v=None, weight=1):
        """
        Adds the edge (u, v) to self.

        Given vertices u, v, and the optional weight, this adds the edge (u, v)
        to self, modifying self.

        Input:
            * self: Graph
            * u: vertex of self
            * v: vertex of self
            * weight (default=1)

        Examples:
            >>> G = Graph()
            >>> G
            Graph on 0 vertices
            >>> G.add_edge(0, 1)
            >>> G.add_edge(0, 2, 2)
            >>> G.add_edge((1, 3))
            >>> G.add_edge((1, 3, 10))
        """
        if v is None:
            if len(u) == 3:
                (u, v, weight) = u
            else:
                (u, v) = u

        pairs = [(u, v), (v, u)]

        for (a, b) in pairs:
            if a not in self.graph_dict:
                self.graph_dict[a] = {b: weight}
            else:
                self.graph_dict[a][b] = weight

    def add_edges(self, edge_list):
        """
        Adds the edges in edge_list to self, modifying self.

        Input:
            * self: Graph

            * edge_list: list
                The edges in edge_list can either be tuples (u, v) or (u, v,
                weight). If no weight is given, the default is weight=1.

        Examples:
            >>> G = Graph()
            >>> G.add_edges([(0, 1), (0, 2), (1, 2)])
            >>> G
            Graph on 3 vertices
        """
        for edge in edge_list:
            self.add_edge(*edge)

    def delete_edge(self, u, v=None):
        """
        Deletes edge from self. If edge is not in self, nothing is done.
        """
        if v is None:
            if len(u) == 2:
                (u, v) = u
            elif len(u) == 3:
                (u, v, _) = u
            else:
                raise ValueError("delete_edge: Invalid edge passed.")

        graph_dict = self.graph_dict

        if v in graph_dict:
            if u in graph_dict[v]:
                del graph_dict[v][u]

        if u in graph_dict:
            if v in graph_dict[u]:
                del graph_dict[u][v]

    def remove_edge(self, u, v=None):
        """
        Deletes edge from self. If edge is not in self, nothing is done.
        """
        if v is None:
            if len(u) == 2:
                (u, v) = u
            elif len(u) == 3:
                (u, v, _) = u
            else:
                raise ValueError("delete_edge: Invalid edge passed.")

        new_graph_dict = copy.deepcopy(self.graph_dict)

        if v in new_graph_dict:
            if u in new_graph_dict[v]:
                del new_graph_dict[v][u]

        if u in new_graph_dict:
            if v in new_graph_dict[u]:
                del new_graph_dict[u][v]

        new_graph = Graph()
        new_graph.graph_dict = new_graph_dict
        return new_graph

    def remove_neighborhood(self, u):
        """
        Returns a copy of self with the neighborhood of u removed; i.e. u along
        with all vertices adjacent to u.
        """
        neighborhood = self.neighborhood(u)
        new_graph = self.remove_vertices(neighborhood)
        return new_graph

    def adjacency_matrix(self):
        """
        Returns the adjacency matrix of self.

        The adjacency matrix of a graph G with vertices v1, v2, ..., vn is the
        n x n matrix A in which entry A[i][j] is the number of edges in G with
        endpoints (v_i, v_j).

        Input:
            * self: Graph

        Output:
            * mat: MatrixZZ

        Details:
            The adjacency matrix is determined by vertex ordering. Here, we use
            lex order of the vertices. For undirected graphs, the adjacency
            matrix is symmetric, and for simple graphs the matrix has entries 0
            or 1, with 0s on the diagonal.

        Examples:
            >>> G = Graph()
            >>> G.add_edges([(1, 2), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5), (4, 6)])
            >>> G.adjacency_matrix()
            [0 1 0 0 1 0]
            [1 0 1 0 1 0]
            [0 1 0 1 0 0]
            [0 0 1 0 1 1]
            [1 1 0 1 0 0]
            [0 0 0 1 0 0]
        """
        vertices = self.vertices()
        num_vertices = len(vertices)
        mat = rosemary.algebra.matrices.matrices.MatrixZZ(num_vertices)

        for (i, u) in enumerate(vertices):
            for (j, v) in enumerate(vertices):
                if v in self.graph_dict[u]:
                    mat[i][j] = 1
        return mat

    def degree_matrix(self):
        """
        Returns the degree matrix of self.

        The degree matrix of a graph G is the diagonal matrix of degrees of the
        vertices of G.

        Input:
            * self: Graph

        Output:
            * mat: MatrixZZ

        Examples:
            >>> G = Graph()
            >>> G.add_edges([(1, 2), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5), (4, 6)])
            >>> G.degree_matrix()
            [2 0 0 0 0 0]
            [0 3 0 0 0 0]
            [0 0 2 0 0 0]
            [0 0 0 3 0 0]
            [0 0 0 0 3 0]
            [0 0 0 0 0 1]
        """
        vertices = self.vertices()
        num_vertices = len(vertices)
        mat = rosemary.algebra.matrices.matrices.MatrixZZ(num_vertices)

        for (i, v) in enumerate(vertices):
            mat[i][i] = len(self[v])
        return mat

    def laplacian_matrix(self):
        """
        Returns the Laplacian matrix of self.

        The Laplacian matrix of a graph G is D - A, where D is the degree matrix
        of G, and A is the adjacency matrix.

        Input:
            * self: Graph

        Output:
            * mat: MatrixZZ

        Details:
            The Laplacian matrix can be used via the Matrix-Tree Theorem to
            compute the number of spanning trees of a matrix.

        Examples:
            >>> G = Graph()
            >>> G.add_edges([(1, 2), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5), (4, 6)])
            >>> G.laplcian_matrix()
            [ 2 -1  0  0 -1  0]
            [-1  3 -1  0 -1  0]
            [ 0 -1  2 -1  0  0]
            [ 0  0 -1  3 -1 -1]
            [-1 -1  0 -1  3  0]
            [ 0  0  0 -1  0  1]
        """
        degree_mat = self.degree_matrix()
        adjacency_mat = self.adjacency_matrix()
        return degree_mat - adjacency_mat

    def induced_subgraph(self, vertices):
        """
        Returns the subgraph of self induced by vertices.
        """
        induced = Graph()
        vertex_set = set(vertices)

        for v in vertices:
            induced.graph_dict[v] = {}

        for u in vertices:
            for v in (vertex_set & self.neighbors(u)):
                induced.add_edge(u, v)

        return induced

    def complement(self):
        """
        Returns the complement graph of self.

        The complement G' of a simple graph G is the simple graph with the same
        vertex set as G, where an edge uv appears in G' if and only if edge uv
        does not appear in G. Note that this concept isn't well defined for
        weighted graphs. In the case where self is weighted, we give all
        complementary edges weight 1.

        Input:
            * self: Graph

        Output:
            * complement: Graph
        """
        complement_graph = Graph()
        vertices = self.vertices()

        for u in vertices:
            complement_graph.add_vertex(u)

        for (i, u) in enumerate(vertices):
            for j in xrange(i):
                v = vertices[j]
                if u not in self[v]:
                    complement_graph.add_edge(u, v)

        return complement_graph


def complete_graph(n):
    """
    Returns the complete graph on vertices 0, 1, ..., n - 1.
    """
    Kn = Graph()

    for u in xrange(n):
        for v in xrange(u):
            Kn.add_edge(u, v)

    return Kn


def random_graph(n, density):
    vertices = range(n)
    graph = Graph()
    for edge in itertools.combinations(vertices, 2):
        r = random.uniform(0, 1)
        if r <= density:
            graph.add_edge(edge)
    return graph


def load_dimacs_graph(path):
    graph = Graph()
    f = open(path, 'r')

    for line in f:
        if line.startswith('e'):
            (u, v) = [int(e) for e in line.split(' ')[1:]]
            graph.add_edge(u, v)

    return graph


def load_graph(path):
    graph = Graph()
    f = open(path, 'r')

    for line in f:
        (u, v) = [int(e) for e in line.split(' ')]
        graph.add_edge(u, v)

    return graph


def coprime_pairs_graph(n):
    """
    Graph with vertices 1, 2, ..., n, with edges between coprime pairs.
    """
    graph = Graph()
    for i in xrange(1, n + 1):
        for j in xrange(1, i):
            if gcd(i, j) == 1:
                graph.add_edge(i, j)

    return graph
