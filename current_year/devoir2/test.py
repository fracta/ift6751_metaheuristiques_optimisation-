
#------------------------------------------------------------------------------
# setup
%pylab
import pyximport
import numpy as np
pyximport.install(reload_support=True)

import ga_solver
import cvrp

data=cvrp.read_vrp("lit_instances/vrp_50.txt")
prob = cvrp.CVRPProblem(data)
#------------------------------------------------------------------------------




# graphics stuff
import matplotlib.pyplot as plt
import matplotlib

def show_routes(coordinates, routes):
    # plot the clients
    plt.scatter(coordinates[1:]['x'], coordinates[1:]['y']) # 
    plt.scatter(coordinates[0]['x'], coordinates[0]['y'], marker='s', color='r') # depot
    cmap = matplotlib.cm.Set1
    
    # plot the routes
    num_routes = len(routes)
    for (route_index, route) in enumerate(routes):
        path = np.zeros(len(route.get_nodes()-1), dtype=[("x", float), ("y", float)])
        for (i, e) in enumerate(route.get_nodes()):
            path[i] = coordinates[e]
        plt.plot(path['x'], path['y'], color=cmap(route_index / float(num_routes)))
    
    # show the plot
    plt.show()
    return



#------------------------------------------------------------------------------

num_vehicles = ga_solver.approx_num_vehicles(prob.get_weights(), prob.get_vehicle_capacity())
num_clients = prob.get_num_clients()

pop = ga_solver.initialize_population(1, num_clients, num_vehicles)
ind = pop[0]


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
