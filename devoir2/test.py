
#------------------------------------------------------------------------------
# general imports

import pyximport
pyximport.install(reload_support=True)

import numpy as np
import ga_solver
import cvrp


#------------------------------------------------------------------------------
# graphics setup

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
        path = np.zeros(len(route.nodes-1), dtype=[("x", float), ("y", float)])
        for (i, e) in enumerate(route.nodes):
            path[i] = coordinates[e]
        plt.plot(path['x'], path['y'], color=cmap(route_index / float(num_routes)))

    # show the plot
    plt.show()
    return


#------------------------------------------------------------------------------
# start the process

data = cvrp.read_vrp("lit_instances/vrp_75.txt")
prob = cvrp.CVRPProblem(data)



