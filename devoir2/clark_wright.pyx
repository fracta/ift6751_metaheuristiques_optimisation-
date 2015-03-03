

import numpy as np
cimport numpy as np

import cvrp


###############################################################################
# ROUTE OBJECT

cdef class Route:
    """represents a route, sequence of integers.
       routes are separated by a chosen symbol
       in this case, the depot"""
    cpdef readonly list nodes

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


cpdef double distance_savings(Route route1, Route route2, np.ndarray distance_matrix):
    """compute the Clark & Wright distance savings"""
    # s_i,j = c_i0 + c0j - \cij
    return distance_matrix[route1[-2], 0] + distance_matrix[0, route2[-2]] - distance_matrix[route1[-2], route2[-2]]


cpdef solve_cw_parallel(cvrp_problem):
    """solve the cvrp problem using the original clark & wright parallel heuristic"""
    
    routes = [Route([0, i, 0]) for i in range(1, cvrp_problem.num_clients)]

    # keep a list of savings
    memoized_savings = dict()
    for route in range(cvrp.num_clients):
        memoized_savings = dict()
    
    
    
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



###############################################################################
# ROUTES -> GENES, GENES -> ROUTES CONVERTERS

cpdef np.ndarray genes_to_routes(np.ndarray genes):
    """GENES -> ROUTES
       0 is used as separator between routes"""
    assert(genes[0] == 0)
    assert(genes[-1] == 0)
    cdef current_route = [0]
    cdef all_routes = []
    for client in genes[1:]:
        if client == 0:
            # end of the route
            current_route.append(0)
            all_routes.append(Route(np.array(current_route)))
            current_route = [0]
        else:
            current_route.append(client)
    return np.array(all_routes)


cpdef np.ndarray routes_to_genes(routes):
    """ROUTES -> GENES"""
    concatenated = np.array([0])
    for route in routes:
        concatenated = np.concatenate((concatenated, route[1:]))
    return concatenated



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



cpdef np.ndarray get_individual_information(Solution individual,
                                            np.ndarray distance_matrix,
                                            np.ndarray weights):
    """get both the capacity and the distance used by the route"""
    cdef np.ndarray information = np.zeros(len(individual.routes),
         dtype= ([("distance", np.float), ("capacity", np.float)]))

    for (index, route) in enumerate(individual.routes):
        information[index] = get_route_information(route, distance_matrix, weights)
    return information


cpdef optimize_routes(Solution individual,
                      np.ndarray distance_matrix):
    """optimize the routes using steepest improvement"""
    for route in individual.routes:
        steepest_improvement(route, distance_matrix)

    individual.genes = routes_to_genes(individual.routes)
    return
