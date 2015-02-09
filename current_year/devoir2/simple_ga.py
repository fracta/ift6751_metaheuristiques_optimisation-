"""implementations from
Solving the Vehicle Routing Problem with Genetic Algorithms"""

import numpy as np

@cython.cclass
@cython.locals(nodes=np.ndarray)
class Route(object):
    """individuals upon which the evolution acts"""
    def __init__(self, nodes):
        self.nodes = nodes
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


@cython.ccall
@cython.locals(tour=Route, dist=np.ndarray, savings=cython.double)
def steepest_improvement(tour, dist):
    """tour reorganization optimization, goes through every possible 2-opt combinations"""
    savings = 1.
    
    while savings > 0:
        for t1 in range(len(route)):
            for t4 in range(len(route)):
                if (t1 != t4) and (t4 != t1 + 1) and (t4+1 != t1):
                    t1 = tour[t1]
                    t2 = tour[t1 + 1]
                    t3 = tour[t4 + 1]
                    t4 = tour[t4]
                    diff = dist[t1][t2] + distance[t4][t3] - dist[t2][t3] - dist[t1][t4]
                    if diff > savings:
                        savings = diff
                        t1_best = t1
                        t4_best = t4
        if savings > 0:
            rearrange(tour, t1_best, t1_best + 1, t4_best)



