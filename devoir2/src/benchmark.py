"""graphic benchmark for the homework (comparison to best know solutions"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import get_cmap


# GRAPHICS

def display(problem, solution, cmap=matplotlib.cm.Set1, figsize=(15,10), show_legend=True):
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
    fig = plt.figure(figsize=figsize)
    plt.title("{} routes, total distance {:.2f}".format(len(routes), total_distance))
    plt.scatter(positions[1:]['x'], positions[1:]['y']) # clients
    plt.scatter(positions[0]['x'],  positions[0]['y'], marker='s', color='r') # depot

    # plot the routes
    num_routes = float(len(routes))
    for (route_index, route) in enumerate(routes):
        path = np.zeros(len(route.nodes), dtype=[("x", float), ("y", float)])
        for (i, e) in enumerate(route.nodes):
            path[i] = positions[e]
        d = D[route_index]
        w = W[route_index]
        fract = route_index / num_routes
        plt.plot(path['x'], path['y'], label="{:10.2f}; {}".format(d, w), color=cmap(fract))

    if show_legend==True:
        plt.legend()
    for index, tup in enumerate(positions[1:]):
        plt.annotate(tup["id"], (tup["x"], tup["y"]))

    return fig



def display_cm_ax(ax, problem, solution, cmap=matplotlib.cm.Set1, figsize=(15,10), show_legend=True):
    """print a cvrp solution with route annotation"""

    # make local refs
    routes = [route for route in solution.routes if route.nodes != [0, 0]]
    positions = problem.positions
    distance_matrix = problem.distance_matrix
    weights = problem.weights
    vehicle_capacity = problem.vehicle_capacity

    # get the distance and weight information
    D = [route.get_distance(distance_matrix) for route in routes]
    W = [route.get_weight(weights) for route in routes]

    total_distance = sum(D)
    maximum_distance = max(D)

    # start the layout of the figure

    ax.set_title("{} routes, total distance {:.2f}".format(len(routes), total_distance))
    ax.scatter(positions[1:]['x'], positions[1:]['y']) # clients
    ax.scatter(positions[0]['x'],  positions[0]['y'], marker='s', color='r') # depot

    # plot the routes
    num_routes = float(len(routes))
    for (route_index, route) in enumerate(routes):
        path = np.zeros(len(route.nodes), dtype=[("x", float), ("y", float)])
        for (i, e) in enumerate(route.nodes):
            path[i] = positions[e]
        d = D[route_index]
        w = W[route_index]
        fract = route_index / num_routes
        ax.plot(path['x'], path['y'], label="{:10.2f}; {}".format(d, w), color=cmap(fract))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    return



def display_routes(ax, problem, routes, figsize=(15,10), show_legend=False):
    """print a cvrp solution with route annotation"""

    # make local refs
    positions = problem.positions
    distance_matrix = problem.distance_matrix
    weights = problem.weights
    vehicle_capacity = problem.vehicle_capacity

    # get the distance and weight information
    D = [route.get_distance(distance_matrix) for route in routes]
    W = [route.get_weight(weights) for route in routes]

    total_distance = sum(D)
    maximum_distance = max(D)

    # start the layout of the figure
    num_valid_routes = len([r for r in routes if r.nodes != [0, 0]])
    ax.set_title("{} routes, total distance {:.2f}".format(num_valid_routes, total_distance))
    ax.scatter(positions[1:]['x'], positions[1:]['y']) # clients
    ax.scatter(positions[0]['x'],  positions[0]['y'], marker='s', color='black') # depot

    # plot the routes
    num_routes = float(len(routes))
    for (route_index, route) in enumerate(routes):
        path = np.zeros(len(route.nodes), dtype=[("x", float), ("y", float)])
        for (i, e) in enumerate(route.nodes):
            path[i] = positions[e]
        d = D[route_index]
        w = W[route_index]
        fract = route_index / num_routes
        ax.plot(path['x'], path['y'], label="{:10.2f}; {}".format(d, w), color="black")
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    return


def display_routes_cm(problem, routes, cmap, figsize=(15,10), show_legend=True):
    """print a cvrp solution with route annotation"""

    # make local refs
    positions = problem.positions
    distance_matrix = problem.distance_matrix
    weights = problem.weights
    vehicle_capacity = problem.vehicle_capacity

    # get the distance and weight information
    D = [route.get_distance(distance_matrix) for route in routes]
    W = [route.get_weight(weights) for route in routes]

    total_distance = sum(D)
    maximum_distance = max(D)

    # start the layout of the figure
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    num_valid_routes = len([r for r in routes if r.nodes != [0, 0]])
    plt.title("{} routes, total distance {:.2f}".format(num_valid_routes, total_distance))
    ax.scatter(positions[1:]['x'], positions[1:]['y']) # clients
    ax.scatter(positions[0]['x'],  positions[0]['y'], marker='s', color='black') # depot

    # plot the routes
    num_routes = float(len(routes))
    for (route_index, route) in enumerate(routes):
        path = np.zeros(len(route.nodes), dtype=[("x", float), ("y", float)])
        for (i, e) in enumerate(route.nodes):
            path[i] = positions[e]
        d = D[route_index]
        w = W[route_index]
        fract = route_index / num_routes
        ax.plot(path['x'], path['y'], label="{:10.2f}; {}".format(d, w), color=cmap(fract), linewidth=4)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if show_legend==True:
        plt.legend()
    for index, tup in enumerate(positions[1:]):
        plt.annotate(tup["id"], (tup["x"], tup["y"]))
    return fig



def display_routes_col(problem, routes, col, figsize=(15,10), ax=None, show_legend=True, title=""):
    """print a cvrp solution with route annotation"""

    # make local refs
    positions = problem.positions
    distance_matrix = problem.distance_matrix
    weights = problem.weights
    vehicle_capacity = problem.vehicle_capacity

    # get the distance and weight information
    D = [route.get_distance(distance_matrix) for route in routes]
    W = [route.get_weight(weights) for route in routes]

    total_distance = sum(D)
    maximum_distance = max(D)

    # start the layout of the figure
    if ax == None:
        fig, ax = plt.subplots(1, figsize=figsize)

    num_valid_routes = len([r for r in routes if r.nodes != [0, 0]])
    ax.set_title(title)
    ax.scatter(positions[1:]['x'], positions[1:]['y']) # clients
    ax.scatter(positions[0]['x'],  positions[0]['y'], marker='s', color=col) # depot

    # plot the routes
    num_routes = float(len(routes))
    for (route_index, route) in enumerate(routes):
        path = np.zeros(len(route.nodes), dtype=[("x", float), ("y", float)])
        for (i, e) in enumerate(route.nodes):
            path[i] = positions[e]
        d = D[route_index]
        w = W[route_index]
        fract = route_index / num_routes
        ax.plot(path['x'], path['y'], label="{:10.2f}; {}".format(d, w), color=col, linewidth=4)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if show_legend==True:
        plt.legend()
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
