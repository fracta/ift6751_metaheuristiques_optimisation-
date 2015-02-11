"""implementations from Solving the Vehicle Routing Problem with
Genetic Algorithms"""

import numpy as np
cimport numpy as np


# REPRESENTATION OF A SINGLE ROUTE (solution is multiple routes)
cdef class Route:
    """individuals upon which the evolution acts"""
    cdef np.ndarray nodes
    
    def __init__(self, np.ndarray nodes):
        assert(nodes[0] == 0)
        assert(nodes[-1]== 0)
        for i in range(1, len(nodes)-1):
            assert(i != 0)
        assert(len(nodes) > 2)
        self.nodes = nodes
    
    cpdef two_opt(self, int ind1, int ind3):
        """2-opt procedure for vertice exchange"""
        assert(ind1 != ind3 and ind1 + 1 != ind3)
        assert(ind1 < ind3)
        cdef np.ndarray rev = self.nodes[ind1+1:ind3+1]
        rev = rev[::-1]
        self.nodes[ind1+1:ind3+1] = rev
        return
    
    cpdef respects_capacity(self, cvrp_problem):
        """used to assert that the route created respects capacity limit"""
        cdef double total = 0.
        for client in self.nodes:
            total+=cvrp_problem.get_weights()[client]
        if total <= cvrp_problem.get_vehicule_capacity():
            return True
        else:
            return False
    
    cpdef get_distance(self, cvrp_problem):
        """calculate the distance used in the route"""
        cdef double distance = 0
        cdef np.ndarray distance_matrix = cvrp_problem.get_distance_matrix()
        cdef int i
        for i in range(len(self.nodes)-1):
            vertex = distance_matrix[self.nodes[i]][self.nodes[i+1]]
            distance += vertex
        return distance
    
    def get_nodes(self):
        return self.nodes
    def __getitem__(self, index):
        return self.nodes[index]
    def __len__(self):
        return len(self.nodes)
    def __str__(self):
        return str(self.nodes)
    def __repr__(self):
        return self.__str__()


# LOCAL SEARCH METHOD
cpdef steepest_improvement(Route route, np.ndarray distance_matrix):
    """route reorganization optimization, greedy local search"""
    if len(route) < 5:
        # 2 nodes are impossible, 3 and 4 are automatically optimal
        return
    cdef int ind1, ind3, n1, n2, n3, n4
    cdef double savings, proposed_savings
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
            route.two_opt(best_ind1, best_ind3)
        else:
            return
    return


def convex_hull(cities, positions):
    """find the convex hull of a TSP, basically a route"""
    


# CROSSOVER OPERATOR

