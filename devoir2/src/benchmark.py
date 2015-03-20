"""graphic benchmark for the homework (comparison to best know solutions"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib


# GRAPHICS

def display(problem, solution, cmap=matplotlib.cm.Set1):
    """print a cvrp solution with route annotation"""

    # make local refs
    positions = problem.positions
    distance_matrix = problem.distance_matrix
    weights = problem.weights
    vehicle_capacity = problem.vehicle_capacity

    # get the distance and weight information
    routes = solution.routes
    D = [route.get_distance(distance_matrix) for route in routes]
    W = [route.get_weight(weights) for route in routes]

    total_distance = sum(D)
    maximum_distance = max(D)

    # start the layout of the figure
    plt.figure(figsize=(15, 10))
    plt.title("Solved with {} routes, total distance of {:.2f}".format(len(routes), total_distance))
    plt.scatter(positions[1:]['x'], positions[1:]['y']) # clients
    plt.scatter(positions[0]['x'],  positions[0]['y'], marker='s', color='r') # depot

    for index, tup in enumerate(positions[1:]):
        plt.annotate(tup["id"], (tup["x"], tup["y"]))

    # plot the routes
    for (route_index, route) in enumerate(routes):
        path = np.zeros(len(route.nodes), dtype=[("x", float), ("y", float)])
        for (i, e) in enumerate(route.nodes):
            path[i] = positions[e]
        d = D[route_index]
        w = W[route_index]

        plt.plot(path['x'], path['y'], label="{:10.2f}; {}".format(d, w))

    plt.legend()
    plt.show()
    return


# BENCHMARKING

best_known = {
"vrpnc1" : 524.61,
"vrpnc2" : 835.26,
"vrpnc3" : 826.14,
"vrpnc4" : 1028.42,
"vrpnc5" : 1291.29,
"vrpnc11" : 1042.11,
"vrpnc12" : 819.56,
}


def percentage_of_score(score, best_known):
    """returns the percentage of best known score"""
    return (score - best_known) / best_known
