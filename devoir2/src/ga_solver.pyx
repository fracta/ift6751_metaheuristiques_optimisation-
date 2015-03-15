"""genetic algorithm solver for the constrained vehicle routing problem (CVRP)"""

import numpy as np
cimport numpy as np

import clark_wright

cimport cvrp
from cvrp cimport CVRPProblem

cimport routes
from routes cimport Route, steepest_improvement

cimport solution
from solution cimport Solution, get_solution_information

import progress_bar


cpdef solution_union(Solution solution1, Solution solution2):
    """return the biggest union of non intersecting routes between solution 1 and 2"""
    cdef list intersection = []
    cdef int index1, index2
    for index1, route1 in enumerate(solution1.routes):
        tmp = (index1, set())
        set1 = set()
        for index2, route2 in enumerate(solution2.routes):
            if 
    pass


cpdef set find_unserved_clients(list routes, int num_clients):
    """using sets, figure out which clients are still not in a route"""
    assert(num_clients > 0)
    cdef set clients = set(np.arange(1, num_clients+1))
    cdef set served_clients = set()
    for route in routes:
        served_clients = served_clients.union(set(route.nodes))
    return clients.difference(served_clients)


cpdef Solution BRBAX(CVRPProblem cvrp_problem,
                     Solution parent1, Solution parent2,
                     np.ndarray route_info):
    """"Optimised crossover genetic algoritm for capacited vehicle routing problem"
     by Nazif and Lee, 2012, modified with added Clark & Wright"""

    # select m / 2 routes with least discrepancy with the capacity limit from parent 1
    cdef int to_select = np.round(len(parent1.routes)/2.)

    cdef np.ndarray capacity_difference = np.abs(np.subtract(route_info["weight"],
                                                             cvrp_problem.vehicle_capacity))
    cdef np.ndarray indices = np.argpartition(capacity_difference, to_select)[: to_select]
    #cdef np.ndarray indices = np.random.choice(np.arange(len(parent1.routes)), to_select)
    cdef list inherited_routes = []
    for index in indices:
        inherited_routes.append(parent1.routes[index])


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


cpdef mutate(Solution sol, int num_clients, np.ndarray weights):
    """insertion mutation operator (Graglia et al.)"""
    cdef int route1, route2, client, insert_position
    route1, route2 = select_2(0, len(sol.routes))
    if len(sol.routes[route1].nodes) <= 3:
        return # to avoid emptying routes

    # remove the client and adjust the weight of the route
    client = sol.routes[route1].nodes.pop(np.random.randint(1, len(sol.routes[route1].nodes)-1))

    # insert the client and adjust the weight of the route
    insert_position = np.random.randint(1, len(sol.routes[route2].nodes)-1)
    sol.routes[route2].nodes.insert(insert_position, client)
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
    bar = progress_bar.ProgressBar("Clark & Wright population initialization")
    for iteration in range(pop_size - 1):
        bar.update(float(iteration) / (pop_size - 2))
        solutions.append(Solution(clark_wright.cw_parallel(routes, distance_matrix, cvrp_problem.vehicle_capacity)))
    bar.clean()
    return solutions


cpdef double calculate_score(Solution sol, CVRPProblem cvrp_problem, double penalty=1000.):
    """calculate the fitness based on Graglia et al."""
    cdef np.ndarray information = get_solution_information(sol, cvrp_problem.distance_matrix, cvrp_problem.weights)
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
            double recombination_prob=0.65,
            double mutation_prob=0.1,
            int k=4):
    """solve the cvrp problem using a simple genetic algorithm"""
    # initialize the population using the k-savings
    cdef list population = initialize_population(cvrp_problem,
                                                 population_size,
                                                 k)

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
    bar = progress_bar.ProgressBar("genetic algorithm")
    cdef double iteration = 0.

    cdef int i
    cdef list parents
    cdef list children
    cdef np.ndarray p1_info

    for generation_index in range(num_generations):
        # output the progress bar
        iteration += 1.
        bar.update(iteration / num_generations)

        # find out the index of the best solution at this iteration
        current_best_index = 0
        current_best_score = np.inf

        # score the solutions
        for index, solution in enumerate(population):
            # optimize the routes and assign the new scores
            solution.score = calculate_score(solution, cvrp_problem)
            if solution.score < current_best_score:
                current_best_score = solution.score
                current_best_index = index
        best_solutions.append(population[current_best_index].copy())

        # selection process
        parents = binary_tournament_selection(population, population_size*2)
        children = []
        for i in range(population_size):
            parent1 = parents[i*2]
            parent2 = parents[(i*2)+1]

            # crossover
            if np.random.rand() < recombination_prob:
                p1_info = get_solution_information(parent1, cvrp_problem.distance_matrix, cvrp_problem.weights)
                child = BRBAX(cvrp_problem, parent1, parent2, p1_info)
            else:
                child = parent1.copy()

            # mutation
            if np.random.rand() < mutation_prob:
                mutate(child, cvrp_problem.num_clients, cvrp_problem.weights)
            children.append(child)
        # replace the population by its children
        population = children

    # clean the progress bar
    bar.clean()
    return population, best_solutions
