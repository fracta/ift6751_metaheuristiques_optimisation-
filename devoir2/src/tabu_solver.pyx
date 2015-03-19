"""tabu solver for the CVRP, inspired by Taillard 1993"""


import numpy as np
cimport numpy as np

import clark_wright as cw

cimport cvrp
from cvrp cimport CVRPProblem

cimport routes
from routes cimport Route, steepest_improvement

cimport solution
from solution cimport Solution, get_solution_information

import progress_bar


cdef class Move:
    cdef public float value
    cdef public int client1
    cdef public int client2
    cdef public int r1_index
    cdef public int r2_index
    def __init__(self, float value, int c1, int r2, int c2, int r1):
        self.value=value
        self.client1=c1
        self.client2=c2
        self.r1_index=r1
        self.r2_index=r2

    def __repr__(self):
        return str(self)
    def __str__(self):
        return "(c1: {1} -> r2: {2})\t (c2: {3} -> r1: {4}) \t {0}".format(self.value,
                                                                           self.client1,
                                                                           self.r2_index,
                                                                           self.client2,
                                                                           self.r1_index)


cpdef Solution generate_initial_solution(CVRPProblem prob):
    """generate initial solution with Clark & Wright savings"""
    cdef list routes = [Route([0, i, 0], prob.weights[i])
                                         for i in range(1, prob.num_clients+1)]
    cdef Solution result = Solution(cw.cw_parallel(routes,
                                                   prob.distance_matrix,
                                                   prob.vehicle_capacity))
    # sort the routes by their angle to the depot
    solution.sort_routes_by_angle(result, prob.positions)
    return result


cpdef np.ndarray initial_calculate_savings(Solution initial_solution, CVRPProblem problem):
    """fill the savings matrix"""

    cdef int num_routes = len(initial_solution.routes)

    distance_matrix = problem.distance_matrix
    weights = problem.weights
    vehicle_capacity = problem.vehicle_capacity

    # best_moves matrix
    # stores the (value, client1 index, client2 index) of the best move at route best_move[i, j]
    # a client index value of 0 simbolizes no interchange (would move the depot...)
    cdef np.ndarray best_moves = np.empty((num_routes, num_routes), dtype=Move)
    cdef Move default = Move(np.inf, 0, 0, 0, 0)
    best_moves.fill(default)  # initialize

    # calculate triangular only (symmetric, never itself with itself)
    cdef Route route1, route2
    cdef Move insertion, deletion, swap, best

    for index1 in range(0, num_routes-1):
        route1 = initial_solution.routes[index1]
        for index2 in range((index1+1), num_routes):
            route2 = initial_solution.routes[index2]
            # figure out the best move
            best_moves[index1, index2] = best_move(route1, route2,
                                                   distance_matrix,
                                                   weights,
                                                   vehicle_capacity)
    return best_moves


cpdef Move best_move(Route route1, Route route2,
                     np.ndarray distance_matrix,
                     np.ndarray weights,
                     double vehicle_capacity):
    """return the best move between the two routes"""
    cdef Move removal, insertion, interchange, best
    # route 1 -> route 2
    removal = best_client_transfer(route1, route2, distance_matrix, weights, vehicle_capacity)
    best = removal
    # route 2 -> route 1
    insertion  = best_client_insertion(route1, route2, distance_matrix, weights, vehicle_capacity)
    if insertion.value < best.value:
        best = insertion
    # route 1 <-> route 2
    swap = best_client_interchange(route1, route2, distance_matrix, weights, vehicle_capacity)
    if swap.value < best.value:
        best = swap
    return best


cpdef tuple min_index_ut(np.ndarray matrix):
    """search for the minimal value on the upper triangular part of the matrix"""
    cdef int index1, index2
    cdef int best1 = 0
    cdef int best2 = 0
    cdef double best_score = np.inf
    for index1 in range(matrix.shape[0]-1):
        for index2 in range(index1+1, matrix.shape[1]):
            if matrix[index1, index2].value < best_score:
                best_score = matrix[index1, index2].value
                best1 = index1
                best2 = index2
    return (best1, best2)


cpdef tuple least_insertion_cost(int client,
                                 Route route,
                                 np.ndarray distance_matrix,
                                 np.ndarray weights,
                                 double vehicle_capacity):
    cdef double best_insertion_cost = np.inf
    cdef int best_insertion_point = 0

    cdef double v1, v2, v3, insertion_cost
    cdef int index

    # if not feasable (weight is too much), return infinite distance
    if (weights[client] + route.weight) > vehicle_capacity:
        return (best_insertion_cost, best_insertion_point)

    for index in range(len(route.nodes)-1):
        v1 = distance_matrix[route.nodes[index], client]
        v2 = distance_matrix[client, route.nodes[index+1]]
        v3 = distance_matrix[route.nodes[index], route.nodes[index+1]]
        insertion_cost = v1 + v2 - v3

        if insertion_cost < best_insertion_cost:
            best_insertion_cost = insertion_cost
            best_insertion_point = index
    return (best_insertion_cost, best_insertion_point+1)


cpdef Move best_client_transfer(Route route1, Route route2, distance_matrix,  weights, vehicle_capacity):
    """c1 in route1 transfers to route2. reverse the routes to get (1, 0)"""

    cdef double distance_difference
    cdef int current_client, index, insertion_point

    cdef double best_distance_difference = np.inf
    cdef int best_client = 0
    cdef int best_insertion_point = 0

    for index in range(1, len(route1.nodes)-1):
        current_client = route1.nodes[index]

        # the distance decrease when the current client is taken out of the route1
        # will be negative
        removal_difference = removal_cost(route1, index, distance_matrix)
        #assert(removal_difference <= 0)

        # the distance increase when the current client is inserted in route2
        # will be positive
        insertion_difference, insertion_point = least_insertion_cost(current_client,
                                                                     route2,
                                                                     distance_matrix,
                                                                     weights,
                                                                     vehicle_capacity)
        #assert(insertion_difference >= 0)

        # assign
        distance_difference = removal_difference + insertion_difference
        #print "({0} -> {1}) : {2} {3} {4}".format(index, insertion_point, distance_difference, removal_difference, insertion_difference)
        # assign if new optimal found
        if distance_difference < best_distance_difference:
            best_distance_difference = distance_difference
            best_client = current_client
            best_insertion_point = insertion_point
    return Move(best_distance_difference, best_client, best_insertion_point, 0, 0)


cpdef Move best_client_insertion(Route route1, Route route2, distance_matrix,  weights, vehicle_capacity):
    """reuse best_client_transfer by inverting the routes and inverting the indices at the end"""
    cdef Move move = best_client_transfer(route2, route1, distance_matrix, weights, vehicle_capacity)
    move.r1_index, move.r2_index = move.r2_index, move.r1_index  # invert the insertion point indices
    move.client1, move.client2 = move.client2, move.client1  # invert the client indices
    return move


cpdef double removal_cost(Route route,
                          int index,
                          np.ndarray distance_matrix):
    """compute the difference of distance by religating before and after"""
    assert(index != 0 and index != len(route.nodes))
    cdef double v1, v2, v3
    v1 = distance_matrix[route.nodes[index-1], route.nodes[index]]
    v2 = distance_matrix[route.nodes[index], route.nodes[index+1]]
    v3 = distance_matrix[route.nodes[index-1], route.nodes[index+1]]
    return v3 - (v1 + v2)


cpdef Move best_client_interchange(Route route1, Route route2, distance_matrix, weights, vehicle_capacity):
    """lambda interchange (1, 1), more costly than the simple shifts"""
    rt1 = route1.copy()
    rt2 = route2.copy()

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
    assert(route1.is_equal(rt1))
    assert(route2.is_equal(rt2))
    return best_move


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


cdef class TabuList:
    cdef np.ndarray matrix

    def __init__(self, int num_clients, int num_routes):
        """client 0 will be the null customer"""
        self.matrix = np.empty((num_clients + 1, num_routes), dtype=int)
        self.matrix.fill(np.iinfo(int).min)

    cpdef check_tabu(self, int i, int j, int current_iteration):
        """verify if the move is tabu"""
        if self.matrix[i, j] > current_iteration:
            return False
        else:
            return True

    cpdef set_tabu(self, int i, int j, int until_iteration):
        """set the move tabu until specified iteration"""
        self.matrix[i, j] = until_iteration


cpdef steepest_descent(CVRPProblem prob, Solution sol):
    """ """
    cdef Solution copied = sol.copy()
    cdef np.ndarray mat = initial_calculate_savings(copied, prob)
    cdef int x, y
    cdef double saving
    cdef Move mv
    cdef int iteration = 0
    while True:
        iteration += 1
        #print iteration
        x, y = min_index_ut(mat)
        mv = mat[x, y]
        #print "{0}, {1}: {2}".format(x, y, mv)
        if mv.value > 0:
            break
        else:
            # update the solution
            apply_move(copied.routes[x],
                       copied.routes[y],
                       mv, prob.weights)
            # steep improv
            steepest_improvement(copied.routes[x], prob.distance_matrix)
            steepest_improvement(copied.routes[y], prob.distance_matrix)
            # update the moves matrix
            mat = initial_calculate_savings(copied, prob)
    return copied


#cpdef solve(CVRPProblem prob, max_iterations):
    #"""solve the cvrp problem using tabu search"""
    
    ## initialize a solution using the random savings
    #initial_solution = generate_initial_solution(prob, 2)
    
    ## create the Tabu objects and parameters (tabu list and others)
    #TabuList tabu_list = TabuList(prob.num_clients, len(initial_solution.routes))
    #cdef int tabu_duration = ?
    #cdef int tabu_expiration

    ## remember the best solution
    #cdef Solution best_solution = initial_solution
    #cdef double best_score = initial_solution.score

    ## set the current solution and score
    #cdef Solution current_solution = initial_solution
    #cdef double current_score = initial_solution.score

    ## misc objects
    #p = progress_bar.ProgressBar("Tabu Search")
    #cdef tuple selected_move

    ## loop until termination is required
    #for current_iteration in range(max_iterations):

        ## update the progress bar
        #p.update(float(current_iteration) / max_iterations)

        ## choose feasable and admissible move in neighborhood
        #selected_move, delta = best_admissible(neighborhood(current_solution), tabu_list)

        ## if there are no admissible solutions in the neighborhood
        ##  a new solution is chosen from random savings
        #if(delta == np.inf):
            #current_solution, current_score, tabu_list, moves_matrix = restart(current_solution, current_score, tabu_list, moves_matrix)
        ## if the proposed move is admissible
        #else:
            ## assign the tabu status
            #tabu_expiration = current_iteration + tabu_duration
            #tabu_list.set_tabu(selected_move[0], selected_move[1], tabu_expiration)

            ## update the current solution
            #current_solution = update_solution(current_solution, selected_move)

    ## clean the progress bar
    #p.clean()

    #return best_solution, best_score
