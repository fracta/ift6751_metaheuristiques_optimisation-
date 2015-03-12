"""tabu solver for the CVRP, inspired by Taillard 1993"""

import numpy as np
cimport numpy as np

import copy

import cvrp


cdef class Route:
    """represents a route, sequence of integers.
       routes are separated by a chosen symbol
       in this case, the depot"""
    cpdef public list nodes

    def __init__(self, list nodes):
        assert(nodes[0] == 0)
        assert(nodes[-1]== 0)
        assert(len(nodes) > 1), "depot to depot routes are allowed"
        for i in range(1, len(nodes)-1):
            assert(i != 0)
        self.nodes = nodes

    def __getitem__(self, index):
        return self.nodes[index]

    def __len__(self):
        return len(self.nodes)
    def __str__(self):
        return str(self.nodes)
    def __repr__(self):
        return self.__str__()
    def pop(self, int i):
        self.nodes.pop(i)
        return
    def insert(self, int position, int client):
        self.nodes.insert(position, client)
        return
    def get_copy(self):
        return Route(self.nodes)


cpdef get_route_information(Route route,
                            np.ndarray distance_matrix,
                            np.ndarray weights):
    """calculate the distance and the capacity used by the route"""
    cdef double distance = 0.
    cdef double capacity_used = 0.
    for (index, node) in enumerate(route.nodes[:-1]):
        # calculate the distance from this node to the next
        distance += distance_matrix[node][route.nodes[index+1]]
        capacity_used += weights[node]
    return (distance, capacity_used)



###############################################################################
# LOCAL OPTIMIZATION METHOD FOR ROUTE

cpdef two_opt(route, int ind1, int ind3):
    """2-opt procedure for local optimization"""
    assert(ind1 != ind3 and ind1 + 1 != ind3)
    assert(ind1 < ind3)
    cdef np.ndarray rev = route.nodes[ind1+1:ind3+1]
    rev = rev[::-1]
    route.nodes[ind1+1:ind3+1] = rev
    return


cpdef steepest_improvement(route, np.ndarray distance_matrix):
    """route reorganization optimization, greedy local search
       as described in: Solving the Vehicle Routing Problem with Genetic Algorithms,
       Áslaug Sóley Bjarnadóttir"""
    if len(route) < 5:
        # 2 nodes routes are empty, 3 and 4 are automatically optimal
        return
    cdef int ind1, ind3, n1, n2, n3, n4
    cdef int best_ind1 = 0
    cdef int best_ind3 = 0
    cdef double savings = 0.
    cdef double proposed_savings = 0.
    while True:  # iterate until there isn't any better local choice (2-opt)
        savings = 0.
        for ind1 in range(0, len(route)-2):
            for ind3 in range(ind1+2, len(route)-1):
                n1 = route[ind1]
                n2 = route[ind1 + 1]
                n3 = route[ind3]
                n4 = route[ind3+1]
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



# neighborhood definition
# having chosen two tours A and B, we have the following
# (1, 1): client exchange
# (1, 0): client insertion, basically like CW savings


cpdef tuple (Route route1, Route route2):
    """go through every possible neighbors"""
    # find the 
    


cpdef tuple find_insert_cost_index(int client,
                                   Route route,
                                   np.ndarray distance_matrix):
    """figure out the best insertion cost and associated index"""
    cdef int best_index = 0
    cdef double best_cost = np.inf
    cdef double cost
    for i in range(len(route)-1):
        cost = (distance_matrix[route[i], client] +
                distance_matrix[client, route[i+1]])
        if cost < best_cost:
            best_index = i
            best_cost = cost
    return (best_cost, best_index)




