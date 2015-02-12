"""genetic algorithm solver for the constrained vehicule routing problem (CVRP)"""

import numpy as np
cimport numpy as np



# REPRESENTATION OF A SINGLE ROUTE (solution is multiple routes)
cdef class Route:
    """data representation of the routes of a solution"""
    cdef np.ndarray nodes
    def __init__(self, np.ndarray nodes):
        assert(nodes[0] == 0)
        assert(nodes[-1]== 0)
        for i in range(1, len(nodes)-1):
            assert(i != 0)
        assert(len(nodes) > 2)
        self.nodes = nodes
    
    cpdef two_opt(self, int ind1, int ind3):
        """2-opt procedure for vertice exchange"""
        assert(ind1 != ind3 and ind1 + 1 != ind3)
        assert(ind1 < ind3)
        cdef np.ndarray rev = self.nodes[ind1+1:ind3+1]
        rev = rev[::-1]
        self.nodes[ind1+1:ind3+1] = rev
        return
    
    cpdef respects_capacity(self, weights, max_capacity):
        """used to assert that the route created respects capacity limit"""
        cdef double total = 0.
        for client in self.nodes:
            total+=weights[client]
        if total <= max_capacity:
            return True
        else:
            return False
    
    cpdef get_distance(self, distance_matrix):
        """calculate the distance used in the route"""
        cdef double distance = 0
        cdef int i
        for i in range(len(self.nodes)-1):
            vertex = distance_matrix[self.nodes[i]][self.nodes[i+1]]
            distance += vertex
        return distance
    
    def get_nodes(self):
        return self.nodes
    def __getitem__(self, index):
        return self.nodes[index]
    def __len__(self):
        return len(self.nodes)
    def __str__(self):
        return str(self.nodes)
    def __repr__(self):
        return self.__str__()



# LOCAL SEARCH METHOD
cpdef steepest_improvement(Route route, np.ndarray distance_matrix):
    """route reorganization optimization, greedy local search"""
    if len(route) < 5:
        # 2 nodes are impossible, 3 and 4 are automatically optimal
        return
    cdef int ind1, ind3, n1, n2, n3, n4
    cdef double savings, proposed_savings
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
            route.two_opt(best_ind1, best_ind3)
        else:
            return
    return



# MAIN OBJECTS FOR THE GENETIC ALGORITHM
cdef class Individual:
    """individuals upon which the evolution acts"""
    cdef readonly np.ndarray genes
    cdef readonly np.ndarray routes
    def __init__(self, np.ndarray genes, np.ndarray routes):
        self.genes = genes
        self.routes = routes
    cpdef optimize_routes(self, np.ndarray distance_matrix):
        for route in self.routes:
            steepest_improvement(route, distance_matrix)
    cpdef double get_distance(self, np.ndarray distance_matrix):
        cdef double distance = 0
        for route in self.routes:
            distance += route.get_distance(distance_matrix)
        return distance
    def get_genes(self):
        return self.genes
    def __str__(self):
        return str(self.genes)
    def __repr__(self):
        return self.__str__()


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






# MAIN FUNCTIONS
cpdef np.ndarray genes_to_route(route, weights, max_capacity):
    """translate the genes to a route given the vrp problem based on 
    the push forward heuristic (solomon)"""

    solution = [0]
    all_routes = []
    cdef double current_weight = 0
    cdef double client_weight
    cdef double total_weight
    for elem in route:
        client_weight = weights[elem]
        total_weight = client_weight + current_weight
        # case 1: too big to fit on current solution
        if total_weight > max_capacity:
            current_weight = client_weight
            solution.append(0)
            all_routes.append(Route(np.array(solution)))
            solution = [0, elem]
        # case 2: fits exactly right
        elif total_weight == max_capacity:
            current_weight = 0
            solution.append(elem)
            solution.append(0)
            all_routes.append(Route(np.array(solution)))
            solution = [0]
        # case 3: still space left
        else:
            current_weight = total_weight
            solution.append(elem)
    if solution != [0]: # still a route in construction
        solution.append(0)
        all_routes.append(Route(np.array(solution)))
    return np.array(all_routes)


"""solving the CVRP using the strategy outlined in "Optimised crossover genetic
algoritm for capacited vehicule routing problem" by Nazif and Lee, 2012 """

cpdef Population initialize_population(int pop_size, int num_clients, weights, max_capacity, seed=42):
    """use numpy random permutation for the sequence of client visit"""
    np.random.seed(seed)
    pop = []
    for i in range(pop_size):
        genes = np.random.permutation(np.arange(1, num_clients))
        route = genes_to_route(genes, weights, max_capacity)
        pop.append(Individual(genes, route))
    return Population(np.array(pop))



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