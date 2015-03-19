"""declarations for routes.pyx"""


cimport numpy as np

cdef class Route:
    cpdef public list nodes
    cpdef public weight
    cpdef copy(self)
    cpdef int remove_client_index(self, int index, np.ndarray weights)
    cpdef int remove_client(self, int client, np.ndarray weights)
    cpdef add_client(self, int index, int client, np.ndarray weights)
    cpdef bint is_equal(self, Route other)

cpdef tuple get_information(Route route, np.ndarray distance_matrix, np.ndarray weights)


cpdef two_opt(route, int ind1, int ind3)


cpdef steepest_improvement(route, np.ndarray distance_matrix)

