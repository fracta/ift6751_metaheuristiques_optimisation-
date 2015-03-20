"""route and other objects for CVRP optimization"""

import copy
import numpy as np

cdef class Route:
    """represents a route, sequence of integers"""
    def __init__(self, list nodes, weight):
        assert(nodes[0] == 0)
        assert(nodes[-1]== 0)
        assert(len(nodes) > 1), "depot to depot routes are allowed"
        for i in range(1, len(nodes)-1):
            assert(i != 0)
        self.nodes = nodes
        self.weight = weight
    def __str__(self):
        return str(self.nodes)
    def __repr__(self):
        return self.__str__()
    cpdef Route copy(self):
        return Route(copy.copy(self.nodes), self.weight)
    cpdef bint is_equal(self, Route other):
        return np.array_equal(self.nodes, other.nodes) and self.weight == other.weight

    cpdef int remove_client_index(self, int index, np.ndarray weights):
        """remove the client at specified index and update the weight"""
        cdef int client = self.nodes.pop(index)
        self.weight -= weights[client]
        return client

    cpdef int remove_client(self, int client, np.ndarray weights):
        """remove specified client and update the weight"""
        cdef int index=0
        for index, c in enumerate(self.nodes):
            if c == client:
                return self.remove_client_index(index, weights)
        raise LookupError("Specified client isn't in the route")

    cpdef add_client(self, int index, int client, np.ndarray weights):
        """add the client at index and update the weight"""
        self.weight += weights[client]
        self.nodes.insert(index, client)
        return

    cpdef double get_distance(Route self, np.ndarray distance_matrix):
        """get the distance of the route"""
        cdef double total_distance = 0.
        cdef int client1, client2, index
        for index in range(0, len(self.nodes)-1):
            client1 = self.nodes[index]
            client2 = self.nodes[index+1]
            total_distance += distance_matrix[client1, client2]
        return total_distance

    cpdef double get_weight(Route self, np.ndarray weights):
        """get the distance of the route"""
        cdef double total_weight = 0.
        cdef int client
        for client in self.nodes:
            total_weight += weights[client]
        return total_weight
