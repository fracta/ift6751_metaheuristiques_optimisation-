"""genetic algorithm solver for the constrained vehicle routing problem (CVRP)"""

import numpy as np
cimport numpy as np

import copy
import cvrp


cpdef int approx_num_vehicles(np.ndarray weights, double vehicle_capacity):
    """calculate approximate number of vehicles needed based on naive formula
    from Graglia et al."""
    tmp = sum(weights) / vehicle_capacity * 1.3
    return int (np.ceil (tmp))


###############################################################################
# ROUTE OBJECT

cdef class Route:
    """represents a route, sequence of integers.
       routes are separated by a chosen symbol
       in this case, the depot"""
    cpdef readonly np.ndarray nodes

    def __init__(self, np.ndarray nodes):
        assert(nodes[0] == 0)
        assert(nodes[-1]== 0)
        assert(len(nodes) > 1), "depot to depot routes are allowed"
        for i in range(1, len(nodes)-1):
            assert(i != 0)
        self.nodes = nodes

    def __getitem__(self, index):
        return self.nodes[index]

    def __len__(self):
        return len(self.nodes)
    def __str__(self):
        return str(self.nodes)
    def __repr__(self):
        return self.__str__()


cpdef get_route_information(Route route,
                            np.ndarray distance_matrix,
                            np.ndarray weights):
    """calculate the distance and the capacity used by the route"""
    cdef double distance = 0.
    cdef double capacity_used = 0.
    for (index, node) in enumerate(route.nodes[:-1]):
        # calculate the distance from this node to the next
        distance += distance_matrix[node][route.nodes[index+1]]
        capacity_used += weights[node]
    return (distance, capacity_used)



###############################################################################
# LOCAL OPTIMIZATION METHOD FOR ROUTE

cpdef two_opt(route, int ind1, int ind3):
    """2-opt procedure for local optimization"""
    assert(ind1 != ind3 and ind1 + 1 != ind3)
    assert(ind1 < ind3)
    cdef np.ndarray rev = route.nodes[ind1+1:ind3+1]
    rev = rev[::-1]
    route.nodes[ind1+1:ind3+1] = rev
    return


cpdef steepest_improvement(route, np.ndarray distance_matrix):
    """route reorganization optimization, greedy local search
       as described in: Solving the Vehicle Routing Problem with Genetic Algorithms,
       Áslaug Sóley Bjarnadóttir"""
    if len(route) < 5:
        # 2 nodes routes are empty, 3 and 4 are automatically optimal
        return
    cdef int ind1, ind3, n1, n2, n3, n4
    cdef int best_ind1 = 0
    cdef int best_ind3 = 0
    cdef double savings = 0.
    cdef double proposed_savings = 0.
    while True:  # iterate until there isn't any better local choice (2-opt)
        savings = 0.
        for ind1 in range(0, len(route)-2):
            for ind3 in range(ind1+2, len(route)-1):
                n1 = route[ind1]
                n2 = route[ind1 + 1]
                n3 = route[ind3]
                n4 = route[ind3+1]
                actual = distance_matrix[n1][n2] + distance_matrix[n3][n4]
                proposed = distance_matrix[n1][n3] + distance_matrix[n2][n4]
                proposed_savings = actual - proposed
                if proposed_savings > savings:
                    best_ind1 = ind1
                    best_ind3 = ind3
                    savings = proposed_savings
        if savings > 0.:
            two_opt(route, best_ind1, best_ind3)
        else:
            return
    return



###############################################################################
# ROUTES -> GENES, GENES -> ROUTES CONVERTERS

cpdef np.ndarray genes_to_routes(np.ndarray genes):
    """GENES -> ROUTES
       0 is used as separator between routes"""
    assert(genes[0] == 0)
    assert(genes[-1] == 0)
    cdef current_route = [0]
    cdef all_routes = []
    for client in genes[1:]:
        if client == 0:
            # end of the route
            current_route.append(0)
            all_routes.append(Route(np.array(current_route)))
            current_route = [0]
        else:
            current_route.append(client)
    return np.array(all_routes)


cpdef np.ndarray routes_to_genes(routes):
    """ROUTES -> GENES"""
    concatenated = np.array([0])
    for route in routes:
        concatenated = np.concatenate((concatenated, route[1:]))
    return concatenated



###############################################################################
# INDIVIDUALS USED IN THE GENETIC ALGORITHM


cdef inline bint richcmp_helper(int compare, int op):
    """Returns True/False for each compare operation given an op code.
    Compare should act similarly to Java's comparable interface"""
    if op == 2: # ==
        return compare == 0
    elif op == 3: # !=
        return compare != 0
    elif op == 0: # <
        return compare < 0
    elif op == 1: # <=
        return compare <= 0
    elif op == 4: # >
        return compare > 0
    elif op == 5: # >=
        return compare >= 0


cdef class Individual:
    """individuals upon which the evolution acts"""

    cdef readonly np.ndarray genes
    cdef readonly np.ndarray routes
    cdef public double score

    def __init__(self, np.ndarray genes, double score=-1):
        self.genes = genes
        self.routes = genes_to_routes(genes)
        self.score = score

    def __str__(self):
        return str(self.genes)
    def __repr__(self):
        return self.__str__()

    def __richcmp__(Individual self, Individual other not None, int op):
      """used to compare individuals for the ranking in the hall of fame"""
      cdef int compare
      cdef double v1 = self.score
      cdef double v2 = other.score
      if v1 > v2:
          compare = 1
      elif v1 < v2:
          compare = -1
      else:
          compare = 0
      return richcmp_helper(compare, op)

    def __copy__(self):
        return Individual(self.genes, self.score)


cpdef np.ndarray get_individual_information(Individual individual,
                                            np.ndarray distance_matrix,
                                            np.ndarray weights):
    """get both the capacity and the distance used by the route"""
    cdef np.ndarray information = np.zeros(len(individual.routes),
         dtype= ([("distance", np.float), ("capacity", np.float)]))

    for (index, route) in enumerate(individual.routes):
        information[index] = get_route_information(route, distance_matrix, weights)
    return information


cpdef optimize_routes(Individual individual,
                      np.ndarray distance_matrix):
    """optimize the routes using steepest improvement"""
    for route in individual.routes:
        steepest_improvement(route, distance_matrix)

    individual.genes = routes_to_genes(individual.routes)
    return


###############################################################################
# POPULATION USED IN THE GENETIC ALGORITHM

cdef class Population:
    """population of individuals"""
    cdef public np.ndarray individuals

    def __init__(self, np.ndarray individuals):
        self.individuals = individuals

    def __str__(self):
        ret = ""
        for ind in self.individuals:
            ret += str(ind)+"\n"
        return ret

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, index):
        return self.individuals[index]



###############################################################################
# CROSSOVER OPERATOR

cpdef Individual BRBAX(Individual parent1, Individual parent2,
                       np.ndarray route_info,
                       double target_capacity):
    """"Optimised crossover genetic algoritm for capacited vehicle routing problem"
     by Nazif and Lee, 2012"""

    # select m / 2 best routes from parent 1, by best we mean the ones having
    # that have the least discrepancy with the capacity limit
    cdef int to_select = np.round(len(parent1.routes)/2.)

    abs_capacity_difference = np.abs(np.subtract(route_info["capacity"], target_capacity))
    indices = np.argpartition(abs_capacity_difference, to_select)[: to_select]
    selected_routes = parent1.routes[indices]
    selected_genes = routes_to_genes(selected_routes)
    to_add = []
    # don't add separators for nothing
    separator_counter =  to_select + 1
    for client in parent2.genes:
        if client == 0 and separator_counter > 0:
            separator_counter -= 1
        elif client == 0 and separator_counter == 0:
            to_add.append(0)
        elif not client in selected_genes:
            to_add.append(client)
    child_genes = np.concatenate((selected_genes, to_add))
    return Individual(child_genes)



###############################################################################
# GENETIC OPERATORS

cpdef mutate(Individual ind):
    """insertion mutation operator (Graglia et al.)"""
    cdef int swap1 = np.random.randint(1, len(ind.genes)-1)
    cdef int swap2 = np.random.randint(1, len(ind.genes)-1)
    tmp = ind.genes[swap2]
    ind.genes[swap2] = ind.genes[swap1]
    ind.genes[swap1] = tmp
    # update routes
    ind.routes = genes_to_routes(ind.genes)
    return


cdef tuple select_2(int low, int high):
    """select 2 different random integers in the interval"""
    # high will be the LENGTH of the list
    assert(high - 1 > low), "interval is nonsensical"
    cdef int l = np.random.randint(low, high)
    cdef int h = np.random.randint(low, high)
    while (l == h):
        h = np.random.randint(low, high)
    return (l, h)


cpdef list binary_tournament_selection(Population population, int num_to_select):
    """binary tournament selection"""
    assert(num_to_select > 0)
    cdef list selected = []
    cdef int index1 = 0
    cdef int index2 = 0
    cdef int pop_size = len(population.individuals)
    for index in range(num_to_select):
        index1, index2 = select_2(0, pop_size)
        if population[index1].score < population[index2].score:
            selected.append(copy.copy(population[index1]))
        else:
            selected.append(copy.copy(population[index2]))
    return selected



###############################################################################
# GENETIC ALGORITHM HELPERS

cpdef Population initialize_population(int pop_size, int num_clients, int num_vehicles):
    """use numpy random permutation for the sequence of client visit"""
    pop = []
    for _ in range(pop_size):
        clients_and_splitters = np.concatenate((np.zeros(num_vehicles-1, dtype=int), np.arange(1, num_clients+1)))
        clients_and_splitters = np.random.permutation(clients_and_splitters)
        # the genes start at depot and end at depot
        genes = np.zeros(len(clients_and_splitters)+2, dtype=int)
        genes[1:-1] = clients_and_splitters
        pop.append(Individual(genes))
    return Population(np.array(pop))


cpdef double calculate_score(Individual ind,
                             double vehicle_capacity,
                             np.ndarray distance_matrix,
                             np.ndarray weights,
                             double penalty=1000.):
    """calculate the fitness based on Graglia et al."""
    route_info = get_individual_information(ind, distance_matrix, weights)
    cdef double overcap = 0.
    cdef double total_distance = 0.
    cdef double score
    for (distance, capacity_used) in route_info:
        #print "distance: {0} capacity: {1}".format(distance, capacity_used)
        if capacity_used > vehicle_capacity:
            overcap += capacity_used - vehicle_capacity
        total_distance += distance
    score = (overcap * penalty) + total_distance
    return score



###############################################################################
# THE GENETIC ALGORITHM LOOP

cpdef solve(problem,
            int population_size,
            int num_generations,
            int opt_step=75,
            double recombination_prob=0.65,
            double mutation_prob=0.1):
    """solve the cvrp problem using a simple genetic algorithm"""

    cdef int num_vehicles = approx_num_vehicles(problem.weights,
                                                problem.vehicle_capacity)
    cdef Population population = initialize_population(population_size,
                                                       problem.num_clients,
                                                       num_vehicles)
    cdef list best_individuals = []
    cdef Individual parent1, parent2, child

    # start the loop
    for generation_index in range(num_generations):
        if generation_index%25 == 0:
            print generation_index

        current_best_index = 0
        current_best_score = np.inf
        for index, individual in enumerate(population.individuals):
            # optimize the routes and assign the new scores
            if generation_index%opt_step==0:
                optimize_routes(individual, problem.distance_matrix)
            individual.score = calculate_score(individual,
                                               problem.vehicle_capacity,
                                               problem.distance_matrix,
                                               problem.weights)
            if individual.score < current_best_score:
                current_best_score = individual.score
                current_best_index = index
        best_individuals.append(copy.copy(population.individuals[current_best_index]))
        # selection process
        parents = binary_tournament_selection(population, population_size*2)
        children = []
        for i in range(population_size):
            parent1 = parents[i*2]
            parent2 = parents[(i*2)+1]

            # crossover, probability = 0.65
            if np.random.rand() < recombination_prob:
                p1_info = get_individual_information(parent1, problem.distance_matrix, problem.weights)
                child = BRBAX(parent1, parent2, p1_info, problem.vehicle_capacity)

            else:
                child = copy.copy(parent1)

            if np.random.rand() < mutation_prob:
                mutate(child)
            children.append(child)

        population = Population(np.array(children))

    tmp_scores = []
    for individual in population.individuals:
        # optimize the routes and assign the new scores
        optimize_routes(individual, problem.distance_matrix)
        individual.score = calculate_score(individual,
                                           problem.vehicle_capacity,
                                           problem.distance_matrix,
                                           problem.weights)
        tmp_scores.append(individual.score)
    # add the last generation's best to the best_individuals
    best_individuals.append(population.individuals[np.argmin(tmp_scores)])
    return population, best_individuals


# pseudo code view of the whole process from the article

# # INITIALIZATION STEP
# pop = initialize_population()
# elite = []
# evaluate(pop)
#
#   while not condition:
#     # SELECTION
#     selected = prob_bin_tournament_selection(pop)
#
#     # OPTIMIZED CROSSOVER OR NODE SWAP
#     crossed = []
#     for sol1, sol2 in selected:
#       if rand() < threshold: # either crossover or node swap
#         crossed.append(optimized_crossover(sol1, sol2))
#       else:
#         crossed.append(node_swap(sol1, sol2))
#       crossed.append(sol)
#     # 
#     crossed = optimised_crossover(selected) else node_swap(selected)
#
#     # MUTATION
#     pop = mutate(crossed)
#
#     # EVALUATE
#     evaluate(pop)
#
#     # ELITISM AND FILTRATION
#     pop = filtrate_pop(elite, pop)
#     elite = update_elite(elite, pop)
#
# return pop, elite


