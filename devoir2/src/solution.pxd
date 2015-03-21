"""declarations for solution to cvrp problem"""

cimport numpy as np
import numpy as np

cimport routes
from routes cimport Route


cdef class Solution:
    """a solution to the cvrp problem"""
    cdef public list routes
    cdef public double score

    cpdef Solution copy(Solution self)
    cpdef np.ndarray get_information(Solution self, np.ndarray distance_matrix, np.ndarray weights)
    cpdef double get_distance(Solution self, np.ndarray distance_matrix)

cdef inline bint richcmp_helper(int compare, int op):
    """Returns True/False for each compare operation given an op code.
    Compare should act similarly to Java's comparable interface"""
    if op == 2: # ==
        return compare == 0
    elif op == 3: # !=
        return compare != 0
    elif op == 0: # <
        return compare < 0
    elif op == 1: # <=
        return compare <= 0
    elif op == 4: # >
        return compare > 0
    elif op == 5: # >=
        return compare >= 0


cpdef list get_centroids(Solution sol, np.ndarray positions)

cpdef double get_angle(double x1, double y1, double x2, double y2)

cpdef list sort_routes_by_angle(Solution sol, np.ndarray positions)

