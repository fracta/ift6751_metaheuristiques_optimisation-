"""graphic benchmark for the homework (comparison to best know solutions"""


# GRAPHIC

import matplotlib.pyplot as plt
import matplotlib

def get_route_information(Route route, distance_matrix, weights):
    distance = 0.
    capacity_used = 0.
    for (index, node) in enumerate(route.nodes[:-1]):
        # calculate the distance from this node to the next
        distance += distance_matrix[node][route.nodes[index+1]]
        capacity_used += weights[node]
    return (distance, capacity_used)


def print_solution(client_positions, distance_matrix, weights, routes):
    """print a cvrp solution with route annotation"""
    route_information = [get_route_information(i, distance_matrix, weights) for i in routes]
    
    plt.figure(figsize=(15, 10))
    plt.title("total score = {0} with {1} routes".format(score, len(routes)))
    plt.scatter(client_positions[1:]['x'], client_positions[1:]['y']) # clients
    plt.scatter(client_positions[0]['x'], client_positions[0]['y'], marker='s', color='r') # depot
    
    # plot the routes
    cmap = matplotlib.cm.Blues
    num_routes = len(routes)
    for (route_index, route) in enumerate(routes):
        path = np.zeros(len(route.nodes-1), dtype=[("x", float), ("y", float)])
        for (i, e) in enumerate(route.nodes):
            path[i] = client_positions[e]
        plt.plot(path['x'], path['y'], color=cmap(route_index / float(num_routes)), label="{dist = {0}, cap = {1}".format(route_information[i][0], route_information[i][1]))

    plt.show()
    return

christofides = ["vrpnc"+str(i) for i [1, 2, 3, 4, 5, 11, 12]]

