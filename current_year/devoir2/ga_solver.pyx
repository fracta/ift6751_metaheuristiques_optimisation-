"""genetic algorithm solver for the constrained vehicule routing problem (CVRP)"""

import numpy as np
cimport numpy as np


# MAIN OBJECTS FOR THE GENETIC ALGORITHM

cdef class Individual:
    """individuals upon which the evolution acts"""
    cdef readonly np.ndarray genes
    def __init__(self, np.ndarray genes):
        self.genes = genes
    def get_genes(self):
        return self.genes
    def __str__(self):
        return str(self.genes)
    def __repr__(self):
        return self.__str__()


cdef class Population:
    """population of individuals"""
    cdef np.ndarray individuals
    cdef int population_size
    def __init__(self, np.ndarray individuals):
        self.individuals = individuals
    
    def get_individuals(self):
        return self.individuals
    
    def get_pop_size(self):
        return self.population_size


cpdef double get_fitness(np.ndarray route, np.ndarray distance_matrix):
    cdef double distance = 0
    cdef int i
    for i in range(len(route)):
        distance += distance_matrix[route[i]][route[i+1]]
    return distance


# MAIN FUNCTIONS

cpdef np.ndarray genes_to_route(ind, weights, max_capacity):
    """translate the genes to a route given the vrp problem based on 
    the push forward heuristic (solomon)"""

    solution = [0]
    cdef double current_weight = 0
    cdef double client_weight
    cdef double total_weight
    for elem in ind.get_genes():
        print elem
        client_weight = weights[elem]
        total_weight = client_weight + current_weight
        # case 1: too big to fit on current solution
        if total_weight > max_capacity:
            current_weight = client_weight
            solution.append(0)
            solution.append(elem)
        # case 2: fits exactly right
        elif total_weight == max_capacity:
            current_weight = 0
            solution.append(elem)
            solution.append(0)
        # case 3: still space left
        else:
            current_weight = total_weight
            solution.append(elem)
    solution.append(0)
    return np.array(solution, dtype='i')


"""solving the CVRP using the strategy outlined in "Optimised crossover genetic
algoritm for capacited vehicule routing problem" by Nazif and Lee, 2012 """

cpdef Population initialize_population(int pop_size, int num_clients, seed=42):
    """use numpy random permutation for the sequence of client visit"""
    np.random.seed(seed)
    pop = []
    for i in range(pop_size):
        pop.append(Individual(np.random.permutation(num_clients)))
    return Population(np.array(pop))


cpdef swap_nodes

cpdef ga_optimize(int population_size,
                  int num_generations,
                  initialize_population,
                  evaluate_solution,
                  selection,
                  crossover,
                  mutate,
                  apply_elitism):
    cdef Population population = initialize_population()
    cdef Population hall_of_fame
    
    return


    
    

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