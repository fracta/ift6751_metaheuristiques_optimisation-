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
    cpdef copy(self):
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


cpdef tuple get_information(Route route,
                            np.ndarray distance_matrix,
                            np.ndarray weights):
    """calculate the distance and the capacity used by the route"""
    cdef double distance = 0.
    cdef double weight = 0
    for (index, node) in enumerate(route.nodes[:-1]):
        # calculate the distance from this node to the next
        distance += distance_matrix[node][route.nodes[index+1]]
        weight += weights[node]
    return (distance, weight)


cpdef two_opt(route, int ind1, int ind3):
    """2-opt procedure for local optimization"""
    assert(ind1 != ind3 and ind1 + 1 != ind3)
    assert(ind1 < ind3)
    rev = route.nodes[ind1+1:ind3+1]
    rev = rev[::-1]
    route.nodes[ind1+1:ind3+1] = rev
    return


cpdef steepest_improvement(route, np.ndarray distance_matrix):
    """route reorganization optimization, greedy local search
       as described in: Solving the Vehicle Routing Problem with Genetic Algorithms,
       Áslaug Sóley Bjarnadóttir"""
    if len(route.nodes) < 5:
        # 2 nodes routes are empty, 3 and 4 are automatically optimal
        return
    cdef int ind1, ind3, n1, n2, n3, n4
    cdef int best_ind1 = 0
    cdef int best_ind3 = 0
    cdef double savings = 0.
    cdef double proposed_savings = 0.
    while True:  # iterate until there isn't any better local choice (2-opt)
        savings = 0.
        for ind1 in range(0, len(route.nodes)-2):
            for ind3 in range(ind1+2, len(route.nodes)-1):
                n1 = route.nodes[ind1]
                n2 = route.nodes[ind1 + 1]
                n3 = route.nodes[ind3]
                n4 = route.nodes[ind3+1]
                actual = distance_matrix[n1][n2] + distance_matrix[n3][n4]
                proposed = distance_matrix[n1][n3] + distance_matrix[n2][n4]
                proposed_savings = actual - proposed
                if proposed_savings > savings:
                    best_ind1 = ind1
                    best_ind3 = ind3
                    savings = proposed_savings
        if savings > 0.:
            two_opt(route, best_ind1, best_ind3)
        else:
            return
    return
