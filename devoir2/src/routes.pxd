"""declarations for routes.pyx"""


cimport numpy as np

cdef class Route:
    cpdef public list nodes
    cpdef public weight
    cpdef copy(self)


cpdef tuple get_information(Route route, np.ndarray distance_matrix, np.ndarray weights)


cpdef two_opt(route, int ind1, int ind3)


cpdef steepest_improvement(route, np.ndarray distance_matrix)

