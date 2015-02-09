"""implementations from Solving the Vehicle Routing Problem with 
Genetic Algorithms"""

import numpy as np
cimport numpy as np

cdef class Route:
    """individuals upon which the evolution acts"""
    cdef np.ndarray nodes
    
    def __init__(self, np.ndarray nodes):
        self.nodes = nodes
    
    def two_opt(self, ind1, ind3):
        assert(ind1 != ind3 and ind1 + 1 != ind3 and ind1 != ind3 + 1)
        assert(ind1 < ind3)
        ind2 = ind1 + 1
        ind4 = ind3 + 1
        rev = self.nodes[ind1:ind4]
        rev = rev[::-1]
        self.nodes[ind1:ind4] = rev
        return
    
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


#def steepest_improvement(tour, dist):
    #"""tour reorganization optimization, goes through every possible 2-opt combinations"""
    #if len(tour) <= 2:
        #return
    
    #savings = 1.
    
    #while savings > 0:
        #for t1 in range(len(route)):
            #for t4 in range(t1, len(route)):
                #if (t1 != t4) and (t4 != t1 + 1) and (t4+1 != t1):
                    #t1 = tour[t1]
                    #t2 = tour[t1 + 1]
                    #t3 = tour[t4 + 1]
                    #t4 = tour[t4]
                    #diff = dist[t1][t2] + distance[t4][t3] - dist[t2][t3] - dist[t1][t4]
                    #if diff > savings:
                        #savings = diff
                        #t1_best = t1
                        #t4_best = t4
        #if savings > 0:
            #tour.rearrange(t1_best, t1_best + 1, t4_best+1, t4_best)
    #return




#def lambda_opt(graph, distance_matrix):
    #"""local search method for the VRP"""
    
    




#def steepest_improvement():
    #""" """
    