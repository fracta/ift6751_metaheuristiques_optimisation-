# cython: profile=True
"""genetic algorithm solver for the constrained vehicle routing problem (CVRP)"""

import numpy as np
cimport numpy as np

import clark_wright

cimport cvrp
from cvrp cimport CVRPProblem

import interchange

cimport routes
from routes cimport Route

from two_opt import steepest_improvement

cimport solution
from solution cimport Solution

import progress_bar


cpdef tuple find_inherited_routes(Solution solution1,
                                  Solution solution2,
                                  np.ndarray positions):
    """return the biggest union of non intersecting routes between solution 1 and 2"""

    cdef int to_select = np.floor(len(solution1.routes)/4.)
    # rearrange the routes by their angle relative to the depot and get values
    cdef list angles1 = solution.sort_routes_by_angle(solution1, positions)
    cdef list angles2 = solution.sort_routes_by_angle(solution2, positions)

    # select one route and its clockwise neighbors
    cdef list selected_routes1 = []
    cdef int index1 = np.random.randint(0, len(angles1))
    cdef set clients1 = set()
    cdef int ind
    for i in range(to_select):
        ind = (index1 + i) % len(angles1)
        # add the selected route
        selected_routes1.append(ind)
        # add the clients of these routes
        clients1 = clients1.union(set(solution1.routes[ind].nodes[1:-1]))

    # select routes from second parent who share no clients with
    # selected routes from first parent
    cdef list selected_routes2 = []
    cdef set set2
    for (index2, route2) in enumerate(solution2.routes):
        set2 = set(route2.nodes[1:-1])
        if set2.isdisjoint(clients1):
            selected_routes2.append(index2)

    # if there are more routes needed, select a subset
    if len(selected_routes2) > to_select:
        index2 = np.random.randint(0, len(selected_routes2))
        tmp = []
        for i in range(to_select):
            tmp.append(selected_routes2[(i+to_select)%len(selected_routes2)])
        selected_routes2 = tmp
    return (selected_routes1, selected_routes2)


cpdef set find_unserved_clients(list routes, int num_clients):
    """using sets, figure out which clients are still not in a route"""
    assert(num_clients > 0)
    cdef set clients = set(np.arange(1, num_clients+1))
    cdef set served_clients = set()
    for route in routes:
        served_clients = served_clients.union(set(route.nodes))
    return clients.difference(served_clients)


cpdef Solution crossover(CVRPProblem cvrp_problem,
                         Solution parent1, Solution parent2,
                         np.ndarray route_info):
    """"Optimised crossover genetic algoritm for capacited vehicle routing problem"
     by Nazif and Lee, 2012, modified with added Clark & Wright"""

    # select m / 2 routes with least discrepancy with the capacity limit from parent 1

    cdef tuple indices = find_inherited_routes(parent1, parent2, cvrp_problem.positions)
    cdef list inherited_routes = []
    cdef int index
    for index in indices[0]:
        inherited_routes.append(parent1.routes[index])
    for index in indices[1]:
        inherited_routes.append(parent2.routes[index])


    # let's reassemble the rest of the routes with savings :)
    cdef set unserved_clients = find_unserved_clients(inherited_routes, cvrp_problem.num_clients)
    cdef list new_routes = []
    for client in np.arange(1, cvrp_problem.num_clients+1):
        if client in unserved_clients:
            new_routes.append(Route([0, client, 0], cvrp_problem.weights[client]))
        else:
            new_routes.append(None)
    cdef list remaining_routes = clark_wright.cw_parallel(new_routes,
                                                          cvrp_problem.distance_matrix,
                                                          cvrp_problem.vehicle_capacity)
    for route in remaining_routes:
        steepest_improvement(route, cvrp_problem.distance_matrix)
    # concatenate the routes
    inherited_routes.extend(remaining_routes)
    return Solution(inherited_routes)


cpdef mutate(Solution sol, CVRPProblem prob):
    """mutate the solution using lambda interchange"""
    interchange.steepest_descent(sol,
                                 prob.distance_matrix,
                                 prob.weights,
                                 prob.vehicle_capacity,
                                 5)
    return


cpdef tuple select_2(int low, int high):
    """select 2 different random integers in the interval"""
    assert(high - 1 > low), "interval is nonsensical"
    cdef int first = np.random.randint(low, high)
    cdef int second = np.random.randint(low, high)
    while (first == second):
        second = np.random.randint(low, high)
    return (first, second)


cpdef list binary_tournament_selection(list population, int num_to_select):
    """binary tournament selection"""
    assert(num_to_select > 0)
    cdef list selected = []
    cdef int index1 = 0
    cdef int index2 = 0
    cdef int pop_size = len(population)
    for index in range(num_to_select):
        index1, index2 = select_2(0, pop_size)
        if population[index1].score < population[index2].score:
            selected.append(population[index1].copy())
        else:
            selected.append(population[index2].copy())
    return selected


cpdef list initialize_population(CVRPProblem cvrp_problem, int pop_size, int k):
    """use Clark & Wright with random choice (up to k worst) to initialize"""
    assert(pop_size > 0)
    assert(k > 0)

    # let's extract a few variables from the problem settings
    cdef np.ndarray distance_matrix = cvrp_problem.distance_matrix
    cdef np.ndarray weights = cvrp_problem.weights

    cdef list routes = [Route([0, i, 0], cvrp_problem.weights[i]) for i in range(1, cvrp_problem.num_clients+1)]
    cdef list solutions = []

    # add the "best" clark wright solution (the one that selects only the best saving)
    # the calculation should quite fast (about a quarter as small, maybe more)
    solutions.append(Solution(clark_wright.cw_parallel(routes, cvrp_problem.distance_matrix, cvrp_problem.vehicle_capacity)))

    # add now, until the population is filled
    bar = progress_bar.ProgressBar("Random Savings Initialization")
    for iteration in range(pop_size - 1):
        bar.update(float(iteration) / (pop_size - 2))
        solutions.append(Solution(clark_wright.cw_parallel_random(routes, distance_matrix, cvrp_problem.vehicle_capacity, k)))
    bar.clean()
    return solutions


cpdef double calculate_score(Solution sol, CVRPProblem cvrp_problem, double penalty=1000.):
    """calculate the fitness based on Graglia et al."""
    cdef np.ndarray information = sol.get_information(cvrp_problem.distance_matrix, cvrp_problem.weights)
    cdef double overcap = 0.
    cdef double total_distance = 0.
    cdef double score
    for (distance, capacity_used) in information:
        if capacity_used > cvrp_problem.vehicle_capacity:
            overcap += capacity_used - cvrp_problem.vehicle_capacity
        total_distance += distance
    score = (overcap * penalty) + total_distance
    return score


cpdef solve(CVRPProblem cvrp_problem,
            int population_size,
            int num_generations,
            double elitism = 0.1,
            double recombination_prob=0.65,
            double mutation_prob=0.05,
            int k=4):
    """solve the cvrp problem using a simple genetic algorithm"""
    # initialize the population using the k-savings
    assert(population_size > 0)
    assert(num_generations >= 0)
    assert(0 <= elitism <= 1)
    assert(0 <= recombination_prob <= 1)
    assert(0 <= mutation_prob <= 1)
    assert(k > 0)
    cdef list population = initialize_population(cvrp_problem,
                                                 population_size,
                                                 k)
    population.sort() # need to fix this
    # optimize the routes of all individuals using steepest improvement
    for solution in population:
        for route in solution.routes:
            steepest_improvement(route, cvrp_problem.distance_matrix)

    # some structures to remember the best individuals
    cdef list best_solutions = []
    cdef Solution parent1, parent2, child
    cdef int current_best_index
    cdef double current_best_score

    # start the progress meter
    bar = progress_bar.ProgressBar("Main loop")
    cdef double iteration = 0.

    cdef int elite_size = np.round(elitism * population_size)
    cdef int i
    cdef list children, parents, elite
    cdef np.ndarray p1_info

    for generation_index in range(num_generations):
        # output the progress bar
        iteration += 1.
        bar.update(iteration / num_generations)

        # score the solutions
        for index, solution in enumerate(population):
            # optimize the routes and assign the new scores
            solution.score = calculate_score(solution, cvrp_problem)

        # elitist selection
        population = sorted(population)
        elite = [sol.copy() for sol in population[: elite_size]]
        best_solutions.append(elite[0].copy())

        # selection process
        parents = binary_tournament_selection(population, (population_size-elite_size)*2)
        children = []
        for i in range(population_size-elite_size):
            parent1 = parents[i*2]
            parent2 = parents[(i*2)+1]

            # crossover
            if np.random.rand() < recombination_prob:
                p1_info = parent1.get_information(cvrp_problem.distance_matrix, cvrp_problem.weights)
                child = crossover(cvrp_problem, parent1, parent2, p1_info)
            else:
                child = parent1.copy()

            # mutation
            if np.random.rand() < mutation_prob:
                mutate(child, cvrp_problem)
            children.append(child)

        # replace the population by its children and the previous elite
        children.extend(elite)
        population = children

    # clean the progress bar
    bar.clean()

    # improve the solutions
    for sol in best_solutions:
        interchange.steepest_descent(sol,
                                     cvrp_problem.distance_matrix,
                                     cvrp_problem.weights,
                                     cvrp_problem.vehicle_capacity)
    return best_solutions
