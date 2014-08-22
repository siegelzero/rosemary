import rosemary.algebra.matrices.matrices

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

        Input:
            * self: Graph

        Output:
            * num_edges: int

        Examples:
            >>> G = Graph()
            >>> G.add_edge(0, 1)
            >>> G.num_edges()
            1
        """
        return len(self.edge_set())

    def vertices(self, **kwargs):
        """
        Returns a sorted list of the vertices of self.

        Input:
            * self: Graph

            * kwargs:
                Possible keywords:
                * order: string (default=None)
                    An optional ordering of the vertices. This parameter can
                    take the values 'degree', 'induced', or None.

                    If order is 'degree', the vertices are returned in
                    nondecreasing order of degree.

                    If order is 'induced', the vertices are returned in the
                    following order: The first vertex v0 has least degree in
                    self. Next, v1 has least degree in the graph we obtain from
                    deleting v0. We repeat in this fasion until all vertices are
                    listed.

                    If order is None, then the vertices are returned in
                    lexicographic order.

                * reverse: bool (default=False)
                    If True, then the vertices are given in reverse order.

        Output:
            * vertices: list
                An ordered list of the vertices of self.

        Examples:
            >>> G = random_graph(5, 0.5)
            >>> G.vertices()
            [0, 1, 2, 3, 4]
            >>> G.vertices(reverse=True)
            [4, 3, 2, 1, 0]
            >>> G.vertices(order='degree')
            [1, 2, 4, 3, 0]
            >>> G.vertices(order='induced')
            [1, 2, 3, 0, 4]
        """
        vertices = self.graph_dict.keys()
        order = kwargs.get('order', None)
        reverse = kwargs.get('reverse', False)

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

        if reverse:
            vertex_list = vertex_list[::-1]

        return vertex_list

    def vertex_set(self):
        """
        Returns a set containing the vertices of self.

        Input:
            * self: Graph

        Output:
            * vertex_set: set

        Examples:
            >>> G = random_graph(10, 0.5)
            >>> G.vertex_set()
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
        """
        vertex_set = set(self.graph_dict.keys())
        return vertex_set

    def edges(self, weights=False):
        """
        Returns a sorted list of the edges of self.

        Input:
            * self: Graph

            * weights: bool (default=False)
                If True, each edge is in the format (u, v, w), where u and v are
                the endpoints of the edge, and w is the weight (default=1).
                Otherwise, each edges is in the format (u, v)

        Output:
            * edge_list: list

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.edges()
            [(0, 1), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4)]
            >>> G.edges(weights=True)
            [(0, 1, 1), (0, 3, 1), (0, 4, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1)]
        """
        edge_list = sorted(self.edge_set(weights=weights))
        return edge_list

    def edge_set(self, weights=False):
        """
        Returns a set containing the edges of self.

        Input:
            * self: Graph

            * weights: bool (default=False)
                If True, each edge is in the format (u, v, w), where u and v are
                the endpoints of the edge, and w is the weight (default=1).
                Otherwise, each edges is in the format (u, v)

        Output:
            * edge_set: set

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.edge_set()
            {(0, 1), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4)}
            >>> G.edge_set(weights=True)
            {(0, 1, 1), (0, 3, 1), (0, 4, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1)}
        """
        edge_set = set()
        graph_dict = self.graph_dict

        for u in graph_dict:
            for v in graph_dict[u]:
                if weights:
                    weight = graph_dict[u][v]
                    triple = (min(u, v), max(u, v), weight)
                    edge_set.add(triple)
                else:
                    edge_set.add((min(u, v), max(u, v)))

        return edge_set

    def maximum_degree(self):
        """
        Returns the largest degree over all vertices of self.

        Input:
            * self: Graph

        Output:
            * deg: int

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.maximum_degree()
            4
            >>> G.maximum_degree_vertex()
            0
            >>> G.degree(0)
            4
        """
        vertices = self.vertex_set()

        if len(vertices) == 0:
            raise ValueError('maximum_degree: Graph has no vertices.')

        deg = max(self.degree(u) for u in vertices)
        return deg

    def maximum_degree_vertex(self):
        """
        Returns a vertex of self of maximum degree.

        Input:
            * self: Graph

        Output:
            * vertex: vertex of self

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.maximum_degree()
            4
            >>> G.maximum_degree_vertex()
            0
            >>> G.degree(0)
            4
        """
        vertices = self.vertex_set()

        if len(vertices) == 0:
            raise ValueError('maximum_degree_vertex: Graph has no vertices.')

        deg = max(vertices, key=self.degree)
        return deg

    def minimum_degree(self):
        """
        Returns the smallest degree over all vertices of self.

        Input:
            * self: Graph

        Output:
            * deg: int

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.minimum_degree()
            2
            >>> G.minimum_degree_vertex()
            1
            >>> G.degree(1)
            2
        """
        vertices = self.vertex_set()

        if len(vertices) == 0:
            raise ValueError('minimum_degree: Graph has no vertices.')

        deg = min(self.degree(u) for u in vertices)
        return deg

    def minimum_degree_vertex(self):
        """
        Returns a vertex of self of minimum degree.

        Input:
            * self: Graph

        Output:
            * vertex: vertex of self

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.minimum_degree()
            2
            >>> G.minimum_degree_vertex()
            1
            >>> G.degree(1)
            2
        """
        vertices = self.vertex_set()

        if len(vertices) == 0:
            raise ValueError('minimum_degree_vertex: Graph has no vertices.')

        return min(vertices, key=self.degree)

    def degree(self, u):
        """
        Returns the degree of u.

        The degree of a vertex u is the number of edges incident to u. In simple
        graphs, this is equivalent to the number of vertices adjacent to u.

        Input:
            * self: Graph

            * u: vertex of self

        Output:
            * degree: int

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.degree(3)
            2
        """
        degree = len(self.graph_dict.keys())
        return degree

    def total_degree(self, vertex_list=None):
        """
        Returns the total degree of self.

        The total degree of a graph is the sum of the degrees of all its
        vertices.

        Input:
            * self: Graph

            * vertex_list: iterable (list, set, tuple) (default=None)
                Optional list of vertices of self. If given, the sum of the
                degrees of all vertices in this list is returned.

        Output:
            * total: int

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.total_degree([0, 1, 2])
            7
        """
        if vertex_list is None:
            vertex_list = self.vertex_set()

        total = sum(self.degree(u) for u in vertex_list)
        return total

    def density(self):
        """
        Returns the density of self.

        For a simple graph with n vertices, the density is defined as the ratio
        of the number of edges in the graph over the maximum number of edges of
        an n-vertex graph, namely n*(n - 1) / 2.

        Input:
            * self: Graph

        Output:
            * density: float

        Examples:
            >>> G = random_graph(100, 0.3)
            >>> G.density()
            0.3058585858585859
            >>> H = G.complement()
            >>> H.density()
            0.6941414141414142
        """
        num_vertices = self.num_vertices()
        num_edges = self.num_edges()
        total_possible = num_vertices*(num_vertices - 1) / 2.0
        density = num_edges / total_possible
        return density

    def neighbors(self, u):
        """
        Returns the set of neighbors of u.

        The neighbors of u are the vertices adjacent to u, not including u
        itself.

        Input:
            * self: Graph

            * u: vertex of self

        Output:
            * neighbors: set

        Examples:
            >>> G = random_graph(20, 0.4)
            >>> G.neighbors(0)
            {3, 4, 6, 9, 14, 15, 16, 19
        """
        neighbors = set(self.graph_dict[u].keys())
        return neighbors

    neighborhood = neighbors

    def closed_neighborhood(self, u):
        """
        Returns the closed neighborhood of u.

        The closed neighborhood of u is defined to be the set of neighbors of u
        along with u itself.

        Input:
            * self: Graph

            * u: vertex of self

        Output:
            * closed_neighborhood: set

        Examples:
            >>> G = random_graph(20, 0.4)
            >>> G.neighborhood(0)
            {0, 3, 4, 6, 9, 14, 15, 16, 19}
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
        Adds the vertices in vertex_list to self, modifying self. Any vertices
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

        Input:
            * self: Graph

            * u: vertex of self

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.vertices()
            [0, 1, 2, 3, 4]
            >>> G.edges()
            [(0, 1), (0, 3), (0, 4), (1, 3), (1, 4), (2, 4), (3, 4)]
            >>> G.delete_vertex(3)
            >>> G.vertices()
            [0, 1, 2, 4]
            >>> G.edges()
            [(0, 1), (0, 4), (1, 4), (2, 4)]
        """
        graph_dict = self.graph_dict
        if u in graph_dict:
            del graph_dict[u]
            for v in graph_dict:
                if u in graph_dict[v]:
                    del graph_dict[v][u]

    def delete_vertices(self, vertex_list):
        """
        Deletes the vertices in vertex_list from self along with any adjacent
        edges, modifying self. Any vertices not present in self are skipped.

        Input:
            * self: Graph

            * vertex_list: iterable (list, set, tuple)
                List of vertices to delete from self.

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.vertices()
            [0, 1, 2, 3, 4]
            >>> G.delete_vertices([0, 1])
            >>> G.vertices()
            [2, 3, 4]
        """
        graph_dict = self.graph_dict

        for u in vertex_list:
            if u in graph_dict:
                del graph_dict[u]
                for v in graph_dict:
                    if u in graph_dict[v]:
                        del graph_dict[v][u]

    def remove_vertex(self, u):
        """
        Returns a copy of self with vertex u removed.

        Input:
            * self: Graph

            * u: vertex of self

        Output:
            * graph: Graph

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.vertices()
            [0, 1, 2, 3, 4]
            >>> G.edges()
            [(0, 1), (0, 2), (0, 4), (1, 2), (2, 3)]
            >>> H = G.remove_vertex(0)
            >>> H.vertices()
            [1, 2, 3, 4]
            >>> H.edges()
            [(1, 2), (2, 3)]
        """
        new_graph = self.copy()
        new_graph.delete_vertex(u)
        return new_graph

    def remove_vertices(self, vertex_list):
        """
        Returns a copy of self with vertex u removed.

        Input:
            * self: Graph

            * vertex_list: iterable (list, tuple, set)
                List of vertices to remove from self.

        Output:
            * new_graph: Graph

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.vertices()
            [0, 1, 2, 3, 4]
            >>> H = G.remove_vertices([0, 1])
            >>> H.vertices()
            [2, 3, 4]
            >>> H.edges()
            [(2, 3)]
        """
        new_graph = self.copy()
        new_graph.delete_vertices(vertex_list)
        return new_graph

    def add_edge(self, u, v=None, weight=1):
        """
        Adds the edge (u, v) to self.

        Given vertices u, v, and the optional weight, this adds the edge (u, v)
        to self, modifying self.

        Input:
            * self: Graph

            * u: vertex of self

            * v: vertex of self (default=None)

            * weight: number (default=1)

        Examples:
            >>> G = Graph()
            >>> G.add_edge(0, 1)
            >>> G.add_edge(0, 2, 2)
            >>> G.add_edge((1, 3))
            >>> G.add_edge((1, 3, 10))
            >>> G.vertices()
            [0, 1, 2, 3]
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
            >>> G.vertices()
            [0, 1, 2]
        """
        for edge in edge_list:
            self.add_edge(*edge)

    def delete_edge(self, u, v=None):
        """
        Deletes edge from self, modifyign self. If edge is not in self, nothing
        is done.

        Input:
            * self: Graph

            * u: vertex of self

            * v: vertex of self (default=None)

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.edges()
            [(0, 1), (0, 2), (0, 4), (1, 2), (2, 3)]
            >>> G.delete_edge((0, 1))
            >>> G.edges()
            [(0, 2), (0, 4), (1, 2), (2, 3)]
            >>> G.delete_edge(1, 2)
            >>> G.edges()
            [(0, 2), (0, 4),(2, 3)]
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
        Returns a copy of self with the edge (u, v) deleted. If (u, v) is not an
        edge of self, then nothing is deleted.

        Input:
            * self: Graph

            * u: vertex of self

            * v: vertex of self (default=None)

        Output:
            * new_graph: Graph

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.edges()
            [(0, 1), (0, 2), (0, 4), (1, 2), (2, 3)]
            >>> H = G.remove_edge(0, 1)
            >>> H.edges()
            [(0, 2), (0, 4), (1, 2), (2, 3)]
        """
        if v is None:
            if len(u) == 2:
                (u, v) = u
            elif len(u) == 3:
                (u, v, _) = u
            else:
                raise ValueError("remove_edge: Invalid edge passed.")

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
        Returns a copy of self with the neighborhood of u removed.

        Input:
            * self: Graph

            * u: vertex of self

        Output:
            * new_graph: Graph

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.vertices()
            [0, 1, 2, 3, 4]
            >>> G.edges()
            [(0, 1), (0, 2), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
            >>> H = G.remove_neighborhood(0)
            >>> H.vertices()
            [0, 3]
            >>> H.edges()
            []
        """
        neighborhood = self.neighborhood(u)
        new_graph = self.remove_vertices(neighborhood)
        return new_graph

    def remove_closed_neighborhood(self, u):
        """
        Returns a copy of self with the closed neighborhood of u removed.

        Input:
            * self: Graph

            * u: vertex of self

        Output:
            * new_graph: Graph

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.vertices()
            [0, 1, 2, 3, 4]
            >>> G.edges()
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4)]
            >>> H = G.remove_closed_neighborhood(0)
            >>> H.vertices()
            [4]
            >>> H.edges()
            []
        """
        closed_neighborhood = self.closed_neighborhood(u)
        new_graph = self.remove_vertices(closed_neighborhood)
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
            >>> G.laplacian_matrix()
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

    def induced_subgraph(self, vertex_list):
        """
        Returns the subgraph of self induced by the vertices in vertex_list.

        Given a subset S of the vertices of a graph G, the subgraph induced by S
        is the graph with vertex set S and edge set consisting precisely of the
        edges of G whose endpoints are in S.

        Input:
            * self: Graph

            * vertex_list: iterable (list, set, tuple)
                List of vertices of self.

        Output:
            * subgraph: Graph
                Subgraph of self induced by the given vertices.

        Examples:
            >>> G = petersen_graph()
            >>> H = G.induced_subgraph([0, 1, 2])
            >>> H.vertices()
            [0, 1, 2]
            >>> H.edges()
            [(0, 1), (1, 2)]
        """
        induced = Graph()
        vertex_set = set(vertex_list)

        for v in vertex_set:
            induced.graph_dict[v] = {}

        for u in vertex_set:
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

        Examples:
            >>> G = random_graph(5, 0.8)
            >>> G.edges()
            [(0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4)]
            >>> H = G.complement()
            >>> H.edges()
            [(0, 1), (0, 2), (3, 4)]
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

    The complete graph on n vertices is the graph containing an edge between all
    pairs of vertices.

    Input:
        * self: Graph

        * n: int
            Number of vertices

    Output:
        * graph: Graph
            Complete graph on n vertices.

    Examples:
        >>> G = complete_graph(4)
        >>> G.edges()
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        >>> G.num_edges()
        6
    """
    Kn = Graph()

    for u in xrange(n):
        for v in xrange(u):
            Kn.add_edge(u, v)

    return Kn


def petersen_graph():
    """
    Returns the Petersen graph.

    The Petersen graph has a vertex for each 2-element subset of a 5-element
    set, and two vertices are adjacent if and only if the corresponding
    2-element subsets are disjoint from each other.

    Examples:
        >>> G = petersen_graph()
        >>> G.num_edges()
        15
        >>> G.num_vertices()
        10
    """
    graph = Graph()
    edges = [(0, 1), (0, 4), (0, 5), (1, 2), (1, 6), (2, 3), (2, 7), (3, 4),
             (3, 8), (4, 9), (5, 7), (5, 8), (6, 8), (6, 9), (7, 9)]
    graph.add_edges(edges)
    return graph


def random_graph(n, d):
    """
    Returns a random graph on vertices 0, 1, ..., n - 1 with approximate
    density d.

    We construct a graph on n vertices, where the probability that two vertices
    will be adjacent is precisely d. Note that the exact density will is
    unlikely to be exactly d, but will be close (especially for large graphs).

    Input:
        * n: int
            Number of vertices.

        * d: float (0 <= d <= 1)
            Edge probability (approximate density)

    Output:
        * graph: Graph
            Graph with approximate density d.

    Examples:
        >>> G = random_graph(40, 0.3)
        >>> G.density()
        0.2987179487179487
    """
    vertices = range(n)
    graph = Graph()
    for edge in itertools.combinations(vertices, 2):
        r = random.uniform(0, 1)
        if r <= d:
            graph.add_edge(edge)
    return graph


def load_dimacs_graph(path):
    """
    Loads DIMACS graph from the given path.

    The DIMACS graphs are given as text files, with lines in the format "e u v"
    signifying an edge between vertices u and v.

    Input:
        * path: str
            Absolute path to the DIMACS graph file.

    Output:
        * graph: Graph

    Examples:
        >>> G = load_dimacs_graph('~/downloads/DIMACS/coloring/dsjc250.5.txt')
        >>> G.num_vertices()
        250
        >>> G.density()
        0.5033895582329317
    """
    graph = Graph()
    f = open(path, 'r')

    for line in f:
        if line.startswith('e'):
            (u, v) = [int(e) for e in line.split(' ')[1:]]
            graph.add_edge(u, v)

    return graph


def load_graph(path):
    """
    Loads graph from the given path.

    This method reads a file with lines of the form "u v" signifying an edge
    between vertices u and v.

    Input:
        * path: str
            Absolute path to the DIMACS graph file.

    Output:
        * graph: Graph

    Examples:
        >>> G = load_graph("r1000.5.txt")
        >>> G.num_vertices()
        1000
        >>> G.density()
        0.5008368368368369
    """
    graph = Graph()
    f = open(path, 'r')

    for line in f:
        (u, v) = [int(e) for e in line.split(' ')]
        graph.add_edge(u, v)

    return graph


def coprime_pairs_graph(n):
    """
    Graph with vertices 1, 2, ..., n, with edges between coprime pairs.

    Input:
        n: int

    Output:
        graph: Graph

    Examples:
        >>> G = coprime_pairs_graph(6)
        >>> G.num_vertices()
        6
        >>> G.edges()
        [(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 3), (2, 5), (3, 4), (3, 5), (4, 5), (5, 6)]
    """
    graph = Graph()
    for i in xrange(1, n + 1):
        for j in xrange(1, i):
            if gcd(i, j) == 1:
                graph.add_edge(i, j)

    return graph
