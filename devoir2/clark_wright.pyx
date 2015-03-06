

import numpy as np
cimport numpy as np

import cvrp


# generate high quality routes to feed in the genetic algorithm


###############################################################################
# ROUTE OBJECT

cdef class Route:
    """represents a route, sequence of integers.
       routes are separated by a chosen symbol
       in this case, the depot"""
    cpdef readonly list nodes
    cpdef readonly weight

    def __init__(self, list nodes, weight):
        assert(nodes[0] == 0)
        assert(nodes[-1]== 0)
        assert(len(nodes) > 1), "depot to depot routes are allowed"
        for i in range(1, len(nodes)-1):
            assert(i != 0)
        self.nodes = nodes
        self.weight = weight

    def __getitem__(self, index):
        return self.nodes[index]

    def __len__(self):
        return len(self.nodes)
    def __str__(self):
        return str(self.nodes)
    def __repr__(self):
        return self.__str__()


cpdef double distance_savings(Route route1, Route route2,
                              np.ndarray distance_matrix,
                              double route_capacity):
    """compute the Clark & Wright distance savings"""
    # s_i,j = c_i0 + c0j - \cij
    if route1 == None or route2 == None:
        return -np.inf
    elif route1.nodes == route2.nodes:  # any route with itself is not feasable
        return -np.inf
    elif route1.weight + route2.weight > route_capacity:  # any route that exceeds cap is not feasable
        return -np.inf
    else:
        return distance_matrix[route1[-2], 0] + distance_matrix[0, route2[-2]] - distance_matrix[route1[-2], route2[-2]]


cpdef Route merge_routes(Route route1, Route route2):
    """merge two roads, taking the first one and concatenating the second"""
    cdef list concatenated = route1[0:-1]
    concatenated.extend(route2[1:])
    return Route(concatenated, route1.weight + route2.weight)


cpdef calculate_savings(list routes,
                        np.ndarray distance_matrix,
                        np.ndarray savings_matrix,
                        double vehicle_capacity):
    """ """
    assert(savings_matrix.shape[0] == len(routes) and savings_matrix.shape[1] == len(routes)), "wrong dimensionality for distance matrix"
    for index1, route1 in enumerate(routes):
        for index2, route2 in enumerate(routes):
            savings_matrix[index1, index2] = distance_savings(route1, route2, distance_matrix, vehicle_capacity)
    return





cpdef list cw_parallel(routes, np.ndarray distance_matrix, np.ndarray weights, double vehicle_capacity):
    """solve the cvrp problem using the original clark & wright parallel heuristic"""

    # calculate all the savings!
    savings = np.zeros((len(routes), len(routes)), dtype=float)
    calculate_savings(routes, distance_matrix, savings, vehicle_capacity)

    # loop until no good savings left (max (savings) <= 0)
    valid_routes = [i for i in range(len(routes)) if routes[i] != None]
    while True:
        index1, index2 = np.unravel_index(savings.argmax(), savings.shape)
        # stop if there is still a valid saving possible
        if savings[index1, index2] <= 0:
            break
        # merge the routes
        routes[index1] = merge_routes(routes[index1], routes[index2])
        routes[index2] = None
        # set the savings implicating index2 to -infinity
        for i in range(savings.shape[0]):
            savings[i][index2] = -np.inf
            savings[index2][i] = -np.inf
        valid_routes.remove(index2)
        # recalculate all the savings implying the index1
        for i in valid_routes:
            savings[i][index1] = distance_savings(routes[i], routes[index1], distance_matrix, vehicle_capacity)
            savings[index1][i] = distance_savings(routes[index1], routes[i], distance_matrix, vehicle_capacity)

   # remove the non-routes (None) from the route list
    result = []
    for route_index in valid_routes:
        result.append(routes[route_index])
    return result


###############################################################################
# LOCAL OPTIMIZATION METHOD FOR ROUTE

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


###############################################################################
# INDIVIDUALS USED IN THE GENETIC ALGORITHM

cdef class Solution:
    """individuals upon which the evolution acts"""

    cdef readonly list routes
    cdef public double score

    def __init__(self, list routes, double score=-1):
        self.routes = routes
        self.score = score

    def __str__(self):
        return str(self.genes)
    def __repr__(self):
        return self.__str__()


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


cpdef get_score(list routes, distance_matrix, weights):
    res = []
    for route in routes:
        res.append(get_route_information(route, distance_matrix, weights))
    return res
    