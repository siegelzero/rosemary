import numpy
import gmpy

def LLL(B):
    """
    Given a basis b_1, b_2, ..., b_n of a lattice (L, q), this algorithm
    transforms the vectors b_i so that when the algorithm terminates, the
    b_i form an LLL-reduced basis. In addition, the algorithm outputs a
    matrix H giving the coordinates of the LLL-reduced basis in terms of
    the initial basis.
    """
