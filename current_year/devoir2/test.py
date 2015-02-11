import pyximport
import numpy as np
pyximport.install(reload_support=True)

import simple_ga, cvrp

data=cvrp.read_vrp("lit_instances/vrp_50.txt")
prob = cvrp.CVRPProblem(data)

route1 = simple_ga.Route (np.array([0,1,2,3,4,5,0]))
route2 = simple_ga.Route(np.array([0,1,3,2,5,4,0]))
route3 = simple_ga.Route(np.array([0,1,2,3,0]))





def test_random_routes(nroutes, prob):
    for i in range(nroutes):
        print("iter " + str(i))
        cur = np.zeros(7, dtype='i')
        permut = np.random.permutation(np.arange(1, 6))
        cur[1:6] = permut
        route = simple_ga.Route(cur)
        print(route)
        steepest_improvement_debug(route, prob)
        print("\n\n\n")
        
        
def steepest_improvement_debug(route, prob):
    """route reorganization optimization, greedy local search"""
    if len(route) < 5:
        # 2 nodes are impossible, 3 and 4 are automatically optimal
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
            return
    return
