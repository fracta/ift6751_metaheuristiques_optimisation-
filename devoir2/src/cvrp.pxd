
cimport numpy as np


cpdef list read_vrp(str file_problem_name)

cpdef np.ndarray get_distance_matrix(np.ndarray client_positions)

cdef class CVRPProblem:
    """data for the constrained vrp problem"""
    cdef readonly double vehicle_capacity
    cdef readonly np.ndarray positions
    cdef readonly np.ndarray weights
    cdef readonly str problem_name
    cdef readonly np.ndarray distance_matrix
    cdef readonly int num_clients
