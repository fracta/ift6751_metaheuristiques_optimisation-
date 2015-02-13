"""genetic algorithm solver for the constrained vehicule routing problem (CVRP)"""

import numpy as np
cimport numpy as np


# ROUTE RELATED STUFF----------------------------------------------------------
cdef class Route:
    """data representation of the route of a solution"""
    # fields
    cdef np.ndarray nodes
    # constructor
    def __init__(self, np.ndarray nodes):
        assert(nodes[0] == 0)
        assert(nodes[-1]== 0)
        for i in range(1, len(nodes)-1):
            assert(i != 0)
        assert(len(nodes) > 1)  # depot to depot routes are allowed
        self.nodes = nodes
    # getters
    def get_nodes(self):
        return self.nodes
    def __getitem__(self, index):
        return self.nodes[index]
    # representation
    def __len__(self):
        return len(self.nodes)
    def __str__(self):
        return str(self.nodes)
    def __repr__(self):
        return self.__str__()
    # optimization methods
    cpdef two_opt(self, int ind1, int ind3):
        """2-opt procedure for vertice exchange"""
        assert(ind1 != ind3 and ind1 + 1 != ind3)
        assert(ind1 < ind3)
        cdef np.ndarray rev = self.nodes[ind1+1:ind3+1]
        rev = rev[::-1]
        self.nodes[ind1+1:ind3+1] = rev
        return
    # Solving the Vehicle Routing Problem with Genetic Algorithms by
    # Áslaug Sóley Bjarnadóttir
    cpdef steepest_improvement(self, np.ndarray distance_matrix):
        """route reorganization optimization, greedy local search"""
        if len(self) < 5:
            # 2 nodes routes are empty, 3 and 4 are automatically optimal
            return
        cdef int ind1, ind3, n1, n2, n3, n4
        cdef int best_ind1 = 0
        cdef int best_ind3 = 0
        cdef double savings, proposed_savings
        while True:  # iterate until there isn't any better local choice (2-opt)
            savings = 0.
            for ind1 in range(0, len(self)-2):
                for ind3 in range(ind1+2, len(self)-1):
                    n1 = self[ind1]
                    n2 = self[ind1 + 1]
                    n3 = self[ind3]
                    n4 = self[ind3+1]
                    actual = distance_matrix[n1][n2] + distance_matrix[n3][n4]
                    proposed = distance_matrix[n1][n3] + distance_matrix[n2][n4]
                    proposed_savings = actual - proposed
                    if proposed_savings > savings:
                        best_ind1 = ind1
                        best_ind3 = ind3
                        savings = proposed_savings
            if savings > 0.:
                self.two_opt(best_ind1, best_ind3)
            else:
                return
        return

    # used for fitness evaluation, returns distance and capacity used
    cpdef get_information(self,
                          np.ndarray distance_matrix,
                          np.ndarray weights):
        """calculate the distance and the capacity used by the route"""
        cdef double distance = 0.
        cdef double capacity_used = 0.
        for (index, node) in enumerate(self.nodes[:-1]):
            # calculate the distance from this node to the next
            distance += distance_matrix[node][self.nodes[index+1]]
            capacity_used += weights[node]
        return (distance, capacity_used)



# DECODERS---------------------------------------------------------------------
cpdef genes_to_routes(np.ndarray genes):
    """translate the genes to a route given the vrp problem
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

cpdef np.ndarray routes_to_genes(np.ndarray routes):
    concatenated = np.array([0])
    for route in routes:
        concatenated = np.concatenate((concatenated, route[1:]))
    return concatenated





# MAIN OBJECTS FOR THE GENETIC ALGORITHM
cdef class Individual:
    """individuals upon which the evolution acts"""
    # fields
    cdef readonly int[:] genes
    cdef readonly Route [:] routes
    # constructor
    def __init__(self, np.ndarray genes):
        self.genes = genes
        self.routes = genes_to_routes(genes)
    # getters
    def get_routes(self):
        return self.routes
    def get_genes(self):
        return self.genes
    # string representation
    def __str__(self):
        return str(self.genes)
    def __repr__(self):
        return self.__str__()
    # route related methods
    cpdef optimize_routes(self, np.ndarray distance_matrix):
        for route in self.routes:
            route.steepest_improvement(distance_matrix)
        return
    cpdef get_information(self, np.ndarray distance_matrix, np.ndarray weights):
        """ get both the capacity and the distance used by the route"""
        cdef np.ndarray information = np.zeros(len(self.get_routes()), dtype= ([("distance", np.float), ("capacity", np.float)]))
        for (index, route) in enumerate(self.routes):
            information[index] = route.get_information(distance_matrix, weights)
        return information


cdef class Population:
    """population of individuals"""
    cdef np.ndarray individuals
    def __init__(self, np.ndarray individuals):
        self.individuals = individuals
    def get_individuals(self):
        return self.individuals
    def __str__(self):
        ret = ""
        for ind in self.individuals:
            ret += str(ind)+"\n"
        return ret
    def __repr__(self):
        return self.__str__()
    def __getitem__(self, index):
        return self.individuals[index]





# GENETIC OPERATORS------------------------------------------------------------
cpdef BRBAX(Individual parent1, Individual parent2, int num_separators):
    pass

"""solving the CVRP using the strategy outlined in "Optimised crossover genetic
algoritm for capacited vehicule routing problem" by Nazif and Lee, 2012 """

cpdef Population initialize_population(int pop_size, int num_clients, weights, max_capacity, seed=42):
    """use numpy random permutation for the sequence of client visit"""
    np.random.seed(seed)
    pop = []
    for i in range(pop_size):
        genes = np.random.permutation(np.arange(1, num_clients))
        route = genes_to_routes(genes)
        pop.append(Individual(genes, route))
    return Population(np.array(pop))



cpdef np.ndarray evaluate_population(Population pop, int pop_size, fitness_function):
    """apply the fitness function over all the individuals"""
    cdef np.ndarray scores = np.zeros(len(pop_size))
    for (i, ind) in enumerate(pop.get_individuals()):
        scores[i] = fitness_function(ind)
    return scores

cpdef np.ndarray selection(Population pop, int num_to_select, tournament_function):
    selected = []
    for _ in range(num_to_select):
        selected.append(tournament_function(pop))
    return np.array(selected)



  
#cpdef ga_optimize(int population_size,
                  #int num_generations,
                  #initialize_population,
                  #evaluate_solution,
                  #selection,
                  #crossover,
                  #mutate,
                  #apply_elitism):
    #cdef Population population = initialize_population()
    #cdef Population hall_of_fame
    
    #for generation_index in range(num_generations):
        
    
    #return


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