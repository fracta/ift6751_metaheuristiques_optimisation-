"""declarations for routes.pyx"""


cimport numpy as np

cdef class Route:
    cpdef public list nodes
    cpdef public weight
    cpdef Route copy(self)
    cpdef bint is_equal(self, Route other)
    cpdef int remove_client_index(self, int index, np.ndarray weights)
    cpdef int remove_client(self, int client, np.ndarray weights)
    cpdef add_client(self, int index, int client, np.ndarray weights)
    cpdef double get_distance(Route self, np.ndarray distance_matrix)
    cpdef double get_weight(Route self, np.ndarray weights)
