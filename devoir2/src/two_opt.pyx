"""2-opt local search for TSP"""

cimport routes
from routes cimport Route

cimport numpy as np


cpdef two_opt(Route route, int ind1, int ind3):
    """2-opt procedure for local optimization"""
    assert(ind1 != ind3 and ind1 + 1 != ind3)
    assert(ind1 < ind3)
    rev = route.nodes[ind1+1:ind3+1]
    rev = rev[::-1]
    route.nodes[ind1+1:ind3+1] = rev
    return


cpdef steepest_improvement(Route route, np.ndarray distance_matrix):
    """route reorganization optimization, greedy local search
       from "Solving the Vehicle Routing Problem with Genetic Algorithms"
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
