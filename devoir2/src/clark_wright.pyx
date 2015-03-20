# cython: profile=True
"""Clark & Wright savings for the CVRP, also used to feed into Tabu Search and Genetic Algorithm solvers"""

import numpy as np
cimport numpy as np

import copy
import progress_bar
import random

# get the route representation
cimport routes
from routes cimport Route

# two opt improvement
import two_opt
from two_opt import steepest_improvement

# get the solution representation
cimport solution
from solution import Solution


import timeit

cdef inline double distance_savings(Route route1, Route route2,
                                    np.ndarray distance_matrix,
                                    double route_capacity):
    """compute the Clark & Wright distance savings"""
    # s_ij = c_i0 + c_0j - c_ij
    if route1 == None or route2 == None:
        return -np.inf
    elif route1.weight + route2.weight > route_capacity:  # any route that exceeds cap is not feasable
        return -np.inf
    else:
        return distance_matrix[route1.nodes[-2], 0] + distance_matrix[0, route2.nodes[-2]] - distance_matrix[route1.nodes[-2], route2.nodes[-2]]


cpdef merge_routes(list routes,
                   list valid_routes,
                   int index1,
                   int index2,
                   np.ndarray savings,
                   np.ndarray distance_matrix,
                   double vehicle_capacity):
    """merge two roads, updating their value in the savings"""
    cdef list concatenated = routes[index1].nodes[0:-1]
    concatenated.extend(routes[index2].nodes[1:])
    # update the routes
    routes[index1] = Route(concatenated, routes[index1].weight + routes[index2].weight)
    routes[index2] = None
    cdef int i
    # set the savings for merged route to -infinity
    for i in range(savings.shape[0]):
        savings[i][index2] = -np.inf
        savings[index2][i] = -np.inf
    # remove the route merged from the group of valid routes
    valid_routes.remove(index2)
    # recalculate all the savings implying the index1
    for i in valid_routes:
        if i != index1:  # invalid merge
            savings[i][index1] = distance_savings(routes[i], routes[index1], distance_matrix, vehicle_capacity)
            savings[index1][i] = distance_savings(routes[index1], routes[i], distance_matrix, vehicle_capacity)
    return


cpdef np.ndarray calculate_savings(list routes,
                                   list valid_routes,
                                   np.ndarray distance_matrix,
                                   double vehicle_capacity):
    """calculate all the savings for merging"""
    cdef np.ndarray savings = np.empty((len(routes), len(routes)), dtype=float)
    savings.fill(-np.inf)
    cdef int index
    cdef int index2
    cdef Route route1
    cdef Route route2
    for index1 in valid_routes:
        route1 = routes[index1]
        for index2 in valid_routes:
            if index1 != index2:  # illegal merge
                route2 = routes[index2]
                savings[index1, index2] = distance_savings(route1, route2, distance_matrix, vehicle_capacity)
    return savings


cpdef list cw_parallel(list original_routes,
                       np.ndarray distance_matrix,
                       double vehicle_capacity):
    """solve the cvrp problem using the original clark & wright parallel heuristic"""
    # setup
    cdef list routes = copy.copy(original_routes)
    cdef list valid_routes = [i for i in range(len(routes)) if routes[i] != None]
    cdef np.ndarray savings = calculate_savings(routes, valid_routes, distance_matrix, vehicle_capacity)

    # loop until no good savings left (max (savings) <= 0)
    cdef int iter_index = 1
    while True:
        iter_index += 1
        index1, index2 = np.unravel_index(savings.argmax(), (savings.shape[0], savings.shape[1]))
        # stop if there is still a valid merge possible
        if savings[index1, index2] <= 0:
            break
        else:
            merge_routes(routes, valid_routes, index1, index2, savings, distance_matrix, vehicle_capacity)
    # remove the non-routes (None) from the route list
    result = []
    for route_index in valid_routes:
        result.append(routes[route_index])
    return result


cpdef tuple select_from_k(int k, np.ndarray savings):
    """select the index amongst the k best (or less if there isn't enough feasable)
       must make sure that the saving > 0"""
    # need to flatten the matrix and then apply the partition
    cdef np.ndarray partition = np.argpartition(np.matrix.flatten(savings), -k)
    cdef tuple indices = np.unravel_index(partition, (savings.shape[0], savings.shape[1]))
    # return a random choice of indices
    cdef int index, x_index, y_index
    cdef bint feasable = False
    cdef list adequate_indices = []
    # from the k largest savings, find the ones that are > 0
    cdef int l = len(indices[0])
    for index in range(l -k, l):
        x_index = indices[0][index]
        y_index = indices[1][index]
        if savings[x_index, y_index] > 0:
            feasable = True
            adequate_indices.append((x_index, y_index))
    # choose a random index from the feasable
    cdef tuple chosen_index
    if feasable == True:
        chosen_index = random.choice(adequate_indices)
    else:
        chosen_index = (indices[0][0], indices[1][0])
    return chosen_index


cpdef list cw_parallel_random(list original_routes,
                              np.ndarray distance_matrix,
                              double vehicle_capacity,
                              int k):
    """solve the cvrp problem using the original clark & wright parallel heuristic"""
    # calculate all the savings!
    cdef list routes = copy.copy(original_routes)
    cdef list valid_routes = [i for i in range(len(routes)) if routes[i] != None]
    cdef np.ndarray savings = calculate_savings(routes, valid_routes, distance_matrix, vehicle_capacity)

    # loop until no good savings left (max (savings) <= 0)
    iter_index = 1
    while True:
        indices = select_from_k(k, savings)
        index1, index2 = indices
        # stop if there is still a valid merge possible
        iter_index += 1
        if k == 1:
            assert(savings.max() == savings[index1, index2])
        if savings[index1, index2] <= 0:
            break
        else:
            merge_routes(routes, valid_routes, index1, index2, savings, distance_matrix, vehicle_capacity)
    result = []
    for route_index in valid_routes:
        result.append(routes[route_index])
    return result


cpdef list random_savings(prob,
                          double time_limit,  # in seconds
                          int k,
                          np.ndarray distance_matrix,
                          double vehicle_capacity):
    """run random choices of the k-best savings at each step, num times"""
    # as suggested in "A new enhancement of the Clark and Wright Savings.."
    assert(time_limit > 0)
    assert(k > 0)
    cdef list routes = [Route([0, i, 0], prob.weights[i]) for i in range(1, prob.num_clients+1)]
    cdef list solutions = []
    p = progress_bar.ProgressBar("Monte Carlo")
    start_time = timeit.default_timer()
    while True:
        current_time = timeit.default_timer()
        time_difference = current_time - start_time
        if time_limit < time_difference:
            break
        p.update(time_difference/time_limit)
        sol = cw_parallel_random(routes, distance_matrix, vehicle_capacity, k)
        for route in sol:
            steepest_improvement(route, distance_matrix)
        solutions.append(sol)
    p.clean()
    return solutions
