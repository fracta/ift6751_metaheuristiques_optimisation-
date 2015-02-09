"""implementations from Solving the Vehicle Routing Problem with 
Genetic Algorithms"""

import numpy as np
cimport numpy as np

cdef class Route:
    """individuals upon which the evolution acts"""
    cdef np.ndarray nodes
    
    def __init__(self, np.ndarray nodes):
        assert(nodes[0] == 0)
        assert(nodes[-1]== 0)
        for i in range(1, len(nodes)-1):
            assert(i != 0)
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





def steepest_improvement_debug(route, prob):
    """route reorganization optimization, greedy local search"""
    if len(route) < 4:
        return
    
    dist = prob.get_distance_matrix()
    distance = route.get_distance(prob)
    print("initial distance " + str(distance))
    while True:
        savings = 0.
        for ind1 in range(0, len(route)-2):
            for ind3 in range(ind1+2, len(route)-1):
                if (ind3 != ind1 + 1):
                    #print("{0} and {1}".format(ind1,ind3))
                    t1 = route[ind1]
                    t2 = route[ind1 + 1]
                    t3 = route[ind3]
                    t4 = route[ind3+1]
                    actual = dist[t1][t2] + dist[t3][t4]
                    proposed = dist[t1][t3] + dist[t2][t4]
                    if proposed < actual:
                        savings = actual - proposed
                        best_ind1 = ind1
                        best_ind3 = ind3

        if savings > 0.:
            print("best savings = {0} for {1} {2}".format(savings, best_ind1, best_ind3))
            print("{0}->{1}; {2}->{3}".format(best_ind1, best_ind1+1, best_ind3, best_ind3+1))
            route.two_opt(best_ind1, best_ind3)
            new_dist = route.get_distance(prob)
            assert( new_dist < distance)
            print("distance = " + str(new_dist))
            distance = new_dist
        else:
            break
    return


#def test_random_routes(nroutes, prob):
    #for i in range(nroutes):
        #print("iter " + str(i))
        #cur = np.zeros(7, dtype='i')
        #permut = np.random.permutation(np.arange(1, 6))
        #cur[1:6] = permut
        #route = simple_ga.Route(cur)
        #print(route)
        #steepest_improvement_debug(route, prob)
        #print("\n\n\n")

#def lambda_opt(graph, distance_matrix):
    #"""local search method for the VRP"""
    
    




#def steepest_improvement():
    #""" """
    