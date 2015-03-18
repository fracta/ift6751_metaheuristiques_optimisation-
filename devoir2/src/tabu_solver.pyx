"""tabu solver for the CVRP, inspired by Taillard 1993"""


import numpy as np
cimport numpy as np

import clark_wright as cw

import copy

cimport cvrp
from cvrp cimport CVRPProblem

cimport routes
from routes cimport Route, steepest_improvement

cimport solution
from solution cimport Solution, get_solution_information

import progress_bar


cdef class Move:
    cdef public float value
    cdef public int c1_index
    cdef public int c2_index
    cdef public int r1_index
    cdef public int r2_index
    def __init__(self, float value, int c1, int c2, int r1, int r2):
        self.value=value
        self.c1_index=c1
        self.c2_index=c2
        self.r1_index=r1
        self.r2_index=r2

    def __repr__(self):
        return str(self)
    def __str__(self):
        return "{0}; ({1}->{2}) ({3}->{4})".format(self.value, self.c1_index, self.r1_index, self.c2_index, self.r2_index)


cpdef Solution generate_initial_solution(CVRPProblem prob, int k):
    """generate initial solution for the search
       (CVRPProblem problem, int k)"""
    cdef list routes = [Route([0, i, 0], prob.weights[i])
                                         for i in range(1, prob.num_clients+1)]
    cdef Solution result = Solution(cw.cw_parallel_random(routes,
                                                          prob.distance_matrix,
                                                          prob.vehicle_capacity,
                                                          k))
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
            # figure out which better interchange betwen (0, 1), (1, 0) and (1, 1)
            # best c2 -> route1, (0, 1) solution
            insertion = shift_0_1(route1, route2, distance_matrix, weights, vehicle_capacity)
            best = insertion
            # best c1 -> route2, (1, 0) solution
            deletion  = shift_1_0(route1, route2, distance_matrix, weights, vehicle_capacity)
            if deletion.value < best.value:
                best = deletion
            # best (1, 0) solution
            #swap = interchange_1_1(route1, route2, distance_matrix, weights, vehicle_capacity)
           # if swap.value < best.value:
            #    best = swap
            best_moves[index1, index2] = best
            best_moves[index2, index1] = best

    return best_moves



cpdef tuple least_insertion_cost(int client,
                                 Route route,
                                 np.ndarray distance_matrix,
                                 np.ndarray weights,
                                 double vehicle_capacity):
    cdef double least_distance = np.inf
    cdef int best_insertion_point = 0

    # if not feasable (weight is too much), return infinite distance
    if (weights[client] + route.weight) > vehicle_capacity:
        return (least_distance, best_insertion_point)

    cdef int left, right
    cdef double dist
    for index in range(len(route.nodes)-1):
        left = route.nodes[index]
        right = route.nodes[index+1]
        dist = distance_matrix[left, client] + distance_matrix[client, right]
        if dist < least_distance:
            least_distance = dist
            best_insertion_point = index
    return (least_distance, best_insertion_point)


cpdef Move shift_0_1(Route route1, Route route2, distance_matrix,  weights, vehicle_capacity):
    """c1 in route1 transfers to route2. reverse the routes to get (1, 0)"""
    #cdef double best_dist = np.inf
    #cdef int best_client = 0
    #cdef int best_insert = 0
    cdef double current_dist, dist_reduction
    cdef int current_client, previous_client, next_client, index
    cdef int current_insert

    cdef double best_dist = np.inf
    cdef int best_client = 0
    cdef int best_insert = 0

    for index in range(1, len(route1.nodes)-1):
        # assign the clients
        previous_client = route1.nodes[index-1]
        current_client = route1.nodes[index]
        next_client = route1.nodes[index+1]

        dist_reduction = distance_matrix[previous_client, current_client] + distance_matrix[current_client, next_client]
        current_dist, current_insert = least_insertion_cost(current_client, route2, distance_matrix, weights, vehicle_capacity)
        current_dist -= dist_reduction
        # assign if new optimal found
        if current_dist < best_dist:
            best_dist = current_dist
            best_client = current_client
            best_insert = current_insert

    return Move(best_dist, best_client, 0, best_insert, 0)


cpdef Move shift_1_0(Route route1, Route route2, distance_matrix,  weights, vehicle_capacity):
    """reuse shift_0_1 by inverting the routes and inverting the indices at the end"""
    cdef Move move = shift_0_1(route2, route1, distance_matrix, weights, vehicle_capacity)
    move.r1_index, move.r2_index = move.r2_index, move.r1_index  # invert the insertion point indices
    move.c1_index, move.c2_index = move.c2_index, move.c1_index  # invert the client indices
    return move


#cpdef interchange_1_1(Route route1, Route route2, distance_matrix, weights, vehicle_capacity):
    #"""lambda interchange (1, 1), more costly than the simple shifts"""
    #cdef list rt1 = copy.copy(route1.nodes)
    #cdef list rt2 = copy.copy(route2.nodes)

    #cdef int ind1, ind2, client1, client2
    #cdef double weight1 = route1.weight
    #cdef double weight2 = route2.weight

    #cdef double removal_saving, insertion_saving

    #for ind1 in range(1, len(c1[1:-1])):
        ## modify the route
        #removal_saving = distance_matrix
        #client1 = rt1.pop(ind1)
        #for ind2 in range(1, len(c2[1:-1])):
            ## 
            



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
    