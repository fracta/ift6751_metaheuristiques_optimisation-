"""lambda-interchange with lambda=1, from Osman"""
cimport routes
from routes cimport Route

cimport solution
from solution cimport Solution

cimport numpy as np

import two_opt

import numpy as np

import sys


cdef class Move:
    def __init__(self, float value, int c1, int r2, int c2, int r1):
        self.value=value
        self.client1=c1
        self.client2=c2
        self.r1_index=r1
        self.r2_index=r2

    def __repr__(self):
        return str(self)
    def __str__(self):
        return "(c1: {1} -> r2: {2})\t\
                (c2: {3} -> r1: {4}) \t {0}".format(self.value,
                                                    self.client1,
                                                    self.r2_index,
                                                    self.client2,
                                                    self.r1_index)


# INSERTION / REMOVAL COST

cpdef double removal_cost(Route route,
                          int index,
                          np.ndarray distance_matrix):
    """compute the difference of distance by removal of a client"""
    assert(index != 0 and index != len(route.nodes))
    cdef double v1, v2, v3
    v1 = distance_matrix[route.nodes[index-1], route.nodes[index]]
    v2 = distance_matrix[route.nodes[index], route.nodes[index+1]]
    v3 = distance_matrix[route.nodes[index-1], route.nodes[index+1]]
    return v3 - (v1 + v2)


cpdef double insertion_cost(Route route,
                            int index, int client,
                            np.ndarray distance_matrix):
    """compute the difference of distance by insertion of a client"""
    cdef double v1, v2, v3,
    v1 = distance_matrix[route.nodes[index], client]
    v2 = distance_matrix[client, route.nodes[index+1]]
    v3 = distance_matrix[route.nodes[index], route.nodes[index+1]]
    return (v1 + v2) - v3


# BEST INSERTION

cpdef tuple least_insertion_cost(int client,
                                 Route route,
                                 np.ndarray distance_matrix,
                                 np.ndarray weights,
                                 double vehicle_capacity):
    cdef int insertion_point
    cdef int best_insertion_point = 0
    cdef double cost
    cdef double best_cost = np.inf

    # if not feasable (weight is too much), return infinite distance
    if (weights[client] + route.weight) > vehicle_capacity:
        return (best_cost, best_insertion_point)

    for insertion_point in range(len(route.nodes)-1):
        cost = insertion_cost(route, insertion_point, client, distance_matrix)

        if cost < best_cost:
            best_cost = cost
            best_insertion_point = insertion_point
    return (best_cost, best_insertion_point+1)


# BEST MOVES FROM, TO, INTER

cpdef Move transfer_to(Route route1, Route route2,
                       np.ndarray distance_matrix,
                       np.ndarray weights,
                       double vehicle_capacity):
    """c1 in route1 transfers to route2"""
    cdef int client
    cdef int insertion_point
    cdef double distance_difference

    cdef int best_client = 0
    cdef int best_insertion_point = 0
    cdef double best_distance_difference = np.inf

    for index in range(1, len(route1.nodes)-1):
        # assign the client
        client = route1.nodes[index]

        # the distance decrease when the current client is taken out of the route1
        removal_difference = removal_cost(route1, index, distance_matrix)

        # the distance increase when the current client is inserted in route2
        insertion_difference, insertion_point = least_insertion_cost(client,
                                                                     route2,
                                                                     distance_matrix,
                                                                     weights,
                                                                     vehicle_capacity)

        # assign
        distance_difference = removal_difference + insertion_difference

        # assign if new optimal found
        if distance_difference < best_distance_difference:
            best_client = client
            best_insertion_point = insertion_point
            best_distance_difference = distance_difference
    return Move(best_distance_difference, best_client, best_insertion_point, 0, 0)


cpdef Move transfer_from(Route route1, Route route2,
                                 np.ndarray distance_matrix,
                                 np.ndarray weights,
                                 double vehicle_capacity):
    """inverse of transfer to"""
    cdef Move move = transfer_to(route2, route1,
                                 distance_matrix,
                                 weights,
                                 vehicle_capacity)
    # invert the insertion point indices
    move.r1_index, move.r2_index = move.r2_index, move.r1_index
    # invert the client indices
    move.client1, move.client2 = move.client2, move.client1
    return move


cpdef Move best_client_interchange(Route route1, Route route2,
                                   np.ndarray distance_matrix,
                                   np.ndarray weights,
                                   double vehicle_capacity):
    """lambda interchange (1, 1), more costly than the simple shifts"""
    cdef int ind1, ind2, client1, client2
    cdef double removal_savings1, removal_savings2
    cdef double insertion_cost1, insertion_cost2

    # keep the best move seen yet
    cdef Move best_move = Move(np.inf, 0, 0, 0, 0)

    for ind1 in range(1, len(route1.nodes)-1):
        # remove the client in route 1
        removal_savings1 = removal_cost(route1, ind1, distance_matrix)
        client1 = route1.remove_client_index(ind1, weights)

        for ind2 in range(1, len(route2.nodes)-1):
            # remove the client in route 2
            removal_savings2 = removal_cost(route2, ind2, distance_matrix)
            client2 = route2.remove_client_index(ind2, weights)

            # calculate the savings now
            insertion_cost1, insertion_point2 = least_insertion_cost(client1, route2, distance_matrix, weights, vehicle_capacity)
            insertion_cost2, insertion_point1 = least_insertion_cost(client2, route1, distance_matrix, weights, vehicle_capacity)
            overall_dist_diff = insertion_cost1 + insertion_cost2 + removal_savings1 + removal_savings2

            # update best solution
            if overall_dist_diff < best_move.value:
                best_move = Move(overall_dist_diff, client1, insertion_point2, client2, insertion_point1)

            # add back the client in the second route at previous spot
            route2.add_client(ind2, client2, weights)

        # add back the client in the first route at previous spot
        route1.add_client(ind1, client1, weights)
    return best_move


cpdef Move find_best_move(Route route1, Route route2,
                          np.ndarray distance_matrix,
                          np.ndarray weights,
                          double vehicle_capacity):
    """return the best move between the two routes"""
    cdef Move removal, insertion, interchange, best
    # route 1 -> route 2
    removal = transfer_to(route1, route2,
                          distance_matrix,
                          weights,
                          vehicle_capacity)
    best = removal
    # route 2 -> route 1
    insertion  = transfer_from(route1, route2,
                               distance_matrix,
                               weights,
                               vehicle_capacity)
    if insertion.value < best.value:
        best = insertion
    # route 1 <-> route 2
    swap = best_client_interchange(route1, route2,
                                   distance_matrix,
                                   weights,
                                   vehicle_capacity)
    if swap.value < best.value:
        best = swap
    return best


cpdef apply_move(Route route1, Route route2,
                 Move move, np.ndarray weights):
    """apply the specified move to the solution (local exploration)"""
    # 4 cases: invalid, insertion, deletion, interchange
    cdef int client1, client2
    # invalid
    if (move.client1 == 0) and (move.client2 == 0):
        raise ValueError("Move is void")

    # interchange
    elif (move.client1 != 0) and (move.client2 != 0):
        client1 = route1.remove_client(move.client1, weights)
        client2 = route2.remove_client(move.client2, weights)
        route1.add_client(move.r1_index, client2, weights)
        route2.add_client(move.r2_index, client1, weights)

    # deletion
    elif (move.client1 != 0) and (move.client2 == 0):
        client1 = route1.remove_client(move.client1, weights)
        route2.add_client(move.r2_index, client1, weights)

    # insertion
    else:
        client2 = route2.remove_client(move.client2, weights)
        route1.add_client(move.r1_index, client2, weights)
    return


cdef class MovesMatrix:

    def __init__(self, Solution sol,
                 np.ndarray distance_matrix,
                 np.ndarray weights,
                 double vehicle_capacity):
        """fill the best moves upper triangular matrix"""
        cdef int num_routes = len(sol.routes)
        self.matrix = np.empty((num_routes, num_routes), dtype=Move)
        cdef Move default = Move(np.inf, 0, 0, 0, 0)
        self.matrix.fill(default)

        cdef Route route1, route2
        for index1 in range(0, num_routes-1):
            route1 = sol.routes[index1]
            for index2 in range((index1+1), num_routes):
                route2 = sol.routes[index2]
                self.matrix[index1, index2] = find_best_move(route1, route2,
                                                             distance_matrix,
                                                             weights,
                                                             vehicle_capacity)

    cpdef Move get(self, int i, int j):
        """get the move at index i, j"""
        assert(i < j), "illegal index"
        return self.matrix[i, j]

    cpdef tuple min(self):
        """search for the min value in upper triangular part of the matrix"""
        cdef int index1, index2
        cdef int best1 = 0
        cdef int best2 = 0
        cdef double best_score = np.inf
        for index1 in range(self.matrix.shape[0]-1):
            for index2 in range(index1+1, self.matrix.shape[1]):
                if self.matrix[index1, index2].value < best_score:
                    best_score = self.matrix[index1, index2].value
                    best1 = index1
                    best2 = index2
        return (best1, best2)

    cpdef update(self, Solution sol,
                 int index1, int index2,
                 np.ndarray distance_matrix,
                 np.ndarray weights,
                 double vehicle_capacity):
        """update all moves implying route 1 and route 2"""
        cdef Move move
        cdef Route route1, route2
        cdef int num_routes = len(sol.routes)
        cdef int i1, i2
        for i1 in range(0, num_routes-1):
            route1 = sol.routes[i1]
            for i2 in range((i1+1), num_routes):
                if i1==index1 or i1==index2 or i2==index1 or i2==index2:
                    route2 = sol.routes[i2]
                    move = find_best_move(route1, route2,
                                          distance_matrix,
                                          weights,
                                          vehicle_capacity)
                    self.matrix[i1, i2] = move
                    self.matrix[i2, i1] = move


cpdef steepest_descent(Solution sol, np.ndarray distance_matrix,
                     np.ndarray weights, double vehicle_capacity,
                     int max_iteration):
    """greedily make the move that improves most the value"""
    cdef MovesMatrix possible_moves = MovesMatrix(sol,
                                      distance_matrix,
                                      weights,
                                      vehicle_capacity)
    cdef int iteration = 0
    cdef int x,y
    cdef Move move

    while True:
        if iteration > max_iteration:
            break
        x, y = possible_moves.min()
        move = possible_moves.get(x, y)
        #print iteration
        #print move
        # check if stuck in a local minima
        if move.value >= 0:
            break

        # update the solution
        apply_move(sol.routes[x], sol.routes[y], move, weights)

        # improve the path of the concerned routes
        two_opt.steepest_improvement(sol.routes[x], distance_matrix)
        two_opt.steepest_improvement(sol.routes[y], distance_matrix)

        # update the moves matrix
        possible_moves.update(sol, x, y, distance_matrix, weights, vehicle_capacity)

        # update iteration count
        iteration += 1
    return sol
