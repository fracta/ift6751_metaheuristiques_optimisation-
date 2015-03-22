"""declarations for the lambda interchange"""
cimport routes
from routes cimport Route

cimport solution
from solution cimport Solution

cimport numpy as np

cdef class Move:
    cdef public float value
    cdef public int client1
    cdef public int client2
    cdef public int r1_index
    cdef public int r2_index


cpdef double removal_cost(Route route, int index, np.ndarray distance_matrix)


cpdef double insertion_cost(Route route, int index, int client, np.ndarray distance_matrix)


cpdef tuple least_insertion_cost(int client, Route route, np.ndarray distance_matrix, np.ndarray weights, double vehicle_capacity)


cpdef Move transfer_to(Route route1, Route route2, np.ndarray distance_matrix, np.ndarray weights, double vehicle_capacity)


cpdef Move transfer_from(Route route1, Route route2, np.ndarray distance_matrix, np.ndarray weights, double vehicle_capacity)


cpdef Move best_client_interchange(Route route1, Route route2, np.ndarray distance_matrix, np.ndarray weights, double vehicle_capacity)


cpdef Move find_best_move(Route route1, Route route2, np.ndarray distance_matrix, np.ndarray weights, double vehicle_capacity)


cpdef apply_move(Route route1, Route route2, Move move, np.ndarray weights)


cdef class MovesMatrix:
    cdef np.ndarray matrix

    cpdef Move get(self, int i, int j)

    cpdef tuple min(self)

    cpdef update(self, Solution sol, int index1, int index2, np.ndarray distance_matrix, np.ndarray weights, double vehicle_capacity)


cpdef steepest_descent(Solution sol, np.ndarray distance_matrix, np.ndarray weights, double vehicle_capacity, int max_iteration)

