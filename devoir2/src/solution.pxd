"""declarations for solution to cvrp problem"""

cimport numpy as np

cdef class Solution:
    """a solution to the cvrp problem"""
    cdef public list routes
    cdef public double score

    cpdef Solution copy(Solution self)


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

cpdef list get_solution_information(Solution sol, np.ndarray distance_matrix, np.ndarray weights)
