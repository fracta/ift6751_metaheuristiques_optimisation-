"""graphic benchmark for the homework (comparison to best know solutions"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib


# GRAPHICS

def get_route_information(route, distance_matrix, weights):
    distance = 0.
    capacity_used = 0.
    for (index, node) in enumerate(route.nodes[:-1]):
        # calculate the distance from this node to the next
        distance += distance_matrix[node][route.nodes[index+1]]
        capacity_used += weights[node]
    return (distance, capacity_used)


def print_solution(client_positions, distance_matrix, weights, routes):
    """print a cvrp solution with route annotation"""
    routes_zipped = zip(routes, [get_route_information(i, distance_matrix, weights) for i in routes])
    routes_zipped = sorted(routes_zipped, key=lambda x: x[1][0])
    routes, route_information = zip(*routes_zipped) 
    total_dist = sum(map(lambda x:x[0], route_information))
    max_dist = max(map(lambda x:x[0], route_information))

    plt.figure(figsize=(15, 10))
    plt.title("Solved with {} routes, total distance of {:.2f}".format(len(routes), total_dist))
    plt.scatter(client_positions[1:]['x'], client_positions[1:]['y']) # clients
    plt.scatter(client_positions[0]['x'], client_positions[0]['y'], marker='s', color='r') # depot

    for index, tup in enumerate(client_positions[1:]):
        plt.annotate(tup["id"], (tup["x"], tup["y"]))

    # plot the routes
    cmap = matplotlib.cm.Set1
    for (route_index, route) in enumerate(routes):
        if route_information[route_index][0] == 0:
            continue
        path = np.zeros(len(route.nodes), dtype=[("x", float), ("y", float)])
        for (i, e) in enumerate(route.nodes):
            path[i] = client_positions[e]
        plt.plot(path['x'], path['y'], label="{:10.2f}; {}".format(route_information[route_index][0], int(route_information[route_index][1])))
    plt.legend()
    plt.show()
    return

def get_distance(solution, distance_matrix, weights):
    total_distance = 0.
    for route in solution:
        dist, cap = get_route_information(route, distance_matrix, weights)
        total_distance += dist
    return total_distance

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


