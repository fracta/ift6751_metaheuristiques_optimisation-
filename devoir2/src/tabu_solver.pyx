"""tabu solver for the CVRP, inspired by Taillard 1993"""

import numpy as np
cimport numpy as np

import clark_wright as cw

cimport cvrp
from cvrp cimport CVRPProblem

cimport routes
from routes cimport Route

cimport interchange
from interchange cimport Move, least_insertion_cost, apply_move, removal_cost

import two_opt

cimport solution
from solution cimport Solution

import progress_bar

import timeit


cdef class TabuList:
    cdef np.ndarray matrix

    def __init__(self, int num_clients, int num_routes):
        """client 0 will be the null customer"""
        self.matrix = np.empty((num_clients + 1, num_routes), dtype=int)
        self.matrix.fill(np.iinfo(int).min)

    cpdef is_tabu(self, int client1, int route_index2, int client2, int route_index1, int current_iteration):
        """verify if the move is tabu"""
        # move is tabu if TABL(i, p) > iter and TABL(j, q) > iter
        if (self.matrix[client1, route_index2] > current_iteration and
            self.matrix[client2, route_index1] > current_iteration):
            return True
        else:
            return False


    cpdef set_tabu(self, int client1, int route2, int client2, int route1, int until_iteration):
        """set the move tabu until specified iteration"""
        self.matrix[client1, route2] = until_iteration
        self.matrix[client2, route1] = until_iteration


cpdef Move transfer_to_tabu(Solution sol,
                            int route_index1, int route_index2,
                            np.ndarray distance_matrix,
                            np.ndarray weights,
                            double vehicle_capacity,
                            int iteration,
                            TabuList tabulist):
    """c1 in route1 transfers to route2"""
    cdef Route route1 = sol.routes[route_index1]
    cdef Route route2 = sol.routes[route_index2]
    cdef int client1
    cdef int client2 = 0  # no insert back client
    cdef int insertion_point
    cdef double distance_difference
    cdef int best_client = 0
    cdef int best_insertion_point = 0
    cdef double best_distance_difference = np.inf

    for index in range(1, len(route1.nodes)-1):
        # assign the client
        client1 = route1.nodes[index]
        if not tabulist.is_tabu(client1, route_index2, client2, route_index1, iteration):
            # the distance decrease when the current client is taken out of the route1
            removal_difference = removal_cost(route1, index, distance_matrix)
            # the distance increase when the current client is inserted in route2
            insertion_difference, insertion_point = least_insertion_cost(client1,
                                                                         route2,
                                                                         distance_matrix,
                                                                         weights,
                                                                         vehicle_capacity)
            # assign
            distance_difference = removal_difference + insertion_difference
            # assign if new optimal found
            if distance_difference < best_distance_difference:
                best_client = client1
                best_insertion_point = insertion_point
                best_distance_difference = distance_difference
    return Move(best_distance_difference, best_client, best_insertion_point, 0, 0)


cpdef Move transfer_from_tabu(Solution sol,
                              int route_index1,
                              int route_index2,
                              np.ndarray distance_matrix,
                              np.ndarray weights,
                              double vehicle_capacity,
                              int iteration,
                              TabuList tabulist):
    """inverse of transfer to"""
    cdef Move move = transfer_to_tabu(sol,
                                      route_index1, route_index2,
                                      distance_matrix,
                                      weights,
                                      vehicle_capacity,
                                      iteration,
                                      tabulist)
    # invert the insertion point indices
    move.r1_index, move.r2_index = move.r2_index, move.r1_index
    # invert the client indices
    move.client1, move.client2 = move.client2, move.client1
    return move


cpdef Move best_client_interchange_tabu(Solution sol,
                                        int route_index1,
                                        int route_index2,
                                        np.ndarray distance_matrix,
                                        np.ndarray weights,
                                        double vehicle_capacity,
                                        int iteration,
                                        TabuList tabulist):
    """lambda interchange (1, 1), more costly than the simple shifts"""
    cdef Route route1 = sol.routes[route_index1]
    cdef Route route2 = sol.routes[route_index2]
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

            if not tabulist.is_tabu(client1, route_index2, client2, route_index1, iteration):
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


cpdef Move find_best_move_tabu(Solution sol,
                               int route_index1,
                               int route_index2,
                               np.ndarray distance_matrix,
                               np.ndarray weights,
                               double vehicle_capacity,
                               int iteration,
                               TabuList tabulist):
    """return the best move between the two routes"""
    cdef Move removal, insertion, interchange, best
    # route 1 -> route 2
    removal = transfer_to_tabu(sol,
                               route_index1,
                               route_index2,
                               distance_matrix,
                               weights,
                               vehicle_capacity,
                               iteration,
                               tabulist)
    best = removal
    # route 2 -> route 1
    insertion  = transfer_from_tabu(sol,
                                    route_index1,
                                    route_index2,
                                    distance_matrix,
                                    weights,
                                    vehicle_capacity,
                                    iteration,
                                    tabulist)
    if insertion.value < best.value:
        best = insertion
    # route 1 <-> route 2
    swap = best_client_interchange_tabu(sol,
                                        route_index1,
                                        route_index2,
                                        distance_matrix,
                                        weights,
                                        vehicle_capacity,
                                        iteration,
                                        tabulist)
    if swap.value < best.value:
        best = swap
    return best


cdef class MovesMatrixTabu:
    """modified moves matrix taking tabu into account"""
    cdef np.ndarray matrix
    def __init__(self, Solution sol,
                 np.ndarray distance_matrix,
                 np.ndarray weights,
                 double vehicle_capacity,
                 int iteration,
                 TabuList tabulist):
        """fill the best moves upper triangular matrix"""
        cdef int num_routes = len(sol.routes)
        self.matrix = np.empty((num_routes, num_routes), dtype=Move)
        cdef Move default = Move(np.inf, 0, 0, 0, 0)
        self.matrix.fill(default)

        for index1 in range(0, num_routes-1):
            for index2 in range((index1+1), num_routes):
                self.matrix[index1, index2] = find_best_move_tabu(sol,
                                                                  index1,
                                                                  index2,
                                                                  distance_matrix,
                                                                  weights,
                                                                  vehicle_capacity,
                                                                  iteration,
                                                                  tabulist)

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

    #cpdef update(self, Solution sol,
                 #int index1, int index2,
                 #np.ndarray distance_matrix,
                 #np.ndarray weights,
                 #double vehicle_capacity,
                 #int iteration,
                 #TabuList tabulist):
        #"""update all moves implying route 1 and route 2"""
        #cdef Move move
        #cdef Route route1, route2
        #cdef int num_routes = len(sol.routes)
        #cdef int i1, i2
        #for i1 in range(0, num_routes-1):
            #route1 = sol.routes[i1]
            #for i2 in range((i1+1), num_routes):
                #if i1==index1 or i1==index2 or i2==index1 or i2==index2:
                    #route2 = sol.routes[i2]
                    #move = find_best_move_tabu(sol,
                                               #i1,
                                               #i2,
                                               #distance_matrix,
                                               #weights,
                                               #vehicle_capacity,
                                               #iteration,
                                               #tabulist)
                    #self.matrix[i1, i2] = move
                    #self.matrix[i2, i1] = move

    #cpdef update_tabu(MovesMatrixTabu self,
                      #Solution sol,
                      #int client,
                      #np.ndarray distance_matrix,
                      #np.ndarray weights,
                      #double vehicle_capacity,
                      #int iteration,
                      #TabuList tabulist):
        #"""update all the moves implying the client"""
        #assert(client != 0), "invalid client"
        #cdef int route_index = -1
        #cdef int index

        ## find the route in which the client belongs
        #for index, route in enumerate(sol.routes):
            #if client in route.nodes:
                #route_index = index
                #break
        ## update all routes implicating the 
        #cdef int num_routes = len(sol.routes)
        #cdef int i1, i2
        #for i1 in range(0, num_routes-1):
            #for i2 in range((i1+1), num_routes):
                ## if the specified route is amongst them
                #if i1==route_index or i2==route_index:
                    #move = find_best_move_tabu(sol,
                                               #i1,
                                               #i2,
                                               #distance_matrix,
                                               #weights,
                                               #vehicle_capacity,
                                               #iteration,
                                               #tabulist)
                    #self.matrix[i1, i2] = move
                    #self.matrix[i2, i1] = move
        #return


# GENERATE INITIAL SOLUTION

cpdef Solution generate_initial_solution(CVRPProblem prob, int search_size=100):
    """generate initial solution with Clark & Wright savings"""
    cdef list routes = [Route([0, i, 0], prob.weights[i])
                                         for i in range(1, prob.num_clients+1)]
    cdef Solution sol = Solution(cw.cw_parallel(routes,
                                                prob.distance_matrix,
                                                prob.vehicle_capacity))
    # apply 2-opt steepest improvement
    cdef Route route
    for route in sol.routes:
        two_opt.steepest_improvement(route, prob.distance_matrix)

    return sol


cpdef Solution generate_new_solution(CVRPProblem prob, int k):
    """generate new solution with Clark & Wright random savings"""
    cdef list routes = [Route([0, i, 0], prob.weights[i])
                                         for i in range(1, prob.num_clients+1)]
    cdef Solution sol = Solution(cw.cw_parallel_random(routes,
                                                       prob.distance_matrix,
                                                       prob.vehicle_capacity,
                                                       k))

    # apply 2-opt steepest improvement
    cdef Route route
    for route in sol.routes:
        two_opt.steepest_improvement(route, prob.distance_matrix)

    return sol


cpdef solve(CVRPProblem prob, int time_limit, int patience, int k):
    """solve the cvrp problem using tabu search"""
    assert(time_limit > 0)
    assert(patience > 0)
    assert(k > 0)

    # initialize a solution using Clark & Wright savings
    cdef Solution sol = generate_initial_solution(prob)

    cdef Solution best_solution = None
    cdef double best_score = np.inf


    # create the Tabu list and possible moves matrix
    cdef TabuList tabulist = TabuList(prob.num_clients, len(sol.routes))
    cdef MovesMatrixTabu possible_moves

    # set the tabu duration with Osman's regression equation
    cdef int tabu_duration = max(7, (-40 + 9.6 * np.log(prob.num_clients * len(sol.routes))))
    cdef int tabu_expiration

    # misc objects
    p = progress_bar.ProgressBar("Tabu Search")
    cdef Move move
    cdef int x, y
    cdef int tabu_remov1, tabu_remov2
    cdef int current_patience = patience
    start = timeit.default_timer()
    cdef int current_iteration = 0
    while True:
        current_iteration += 1
        current = timeit.default_timer()
        if current - start > time_limit:
            break
        # update the status bar
        p.update(float(current - start) / time_limit)

        # choose feasable and admissible move in neighborhood
        possible_moves = MovesMatrixTabu(sol,
                                         prob.distance_matrix,
                                         prob.weights,
                                         prob.vehicle_capacity,
                                         current_iteration,
                                         tabulist)
        x, y = possible_moves.min()
        move = possible_moves.get(x, y)

        # apply the move
        apply_move(sol.routes[x], sol.routes[y], move, prob.weights)

        # optimize the routes touched by the move
        two_opt.steepest_improvement(sol.routes[x], prob.distance_matrix)
        two_opt.steepest_improvement(sol.routes[y], prob.distance_matrix)

        # add the inverse move to the tabulist
        tabu_expiration = current_iteration + tabu_duration
        tabulist.set_tabu(move.client1, x, move.client2, y, tabu_expiration)


        # check if the score is better and update
        sol.score = sol.get_distance(prob.distance_matrix)
        if sol.score < best_score:
            best_score = sol.score
            best_solution = sol.copy()
            last_optimal_found = current_iteration
            current_patience = patience

        # if there is no patience left, reinitialize a solution using random savings
        if current_patience <= 0:
            # initialize a new solution and Tabu List
            sol = generate_new_solution(prob, k)
            tabulist = TabuList(prob.num_clients, len(sol.routes))
            current_patience = patience
        # decrement the patience
        current_patience -= 1

    # clean the progress bar
    p.clean()
    return best_solution
