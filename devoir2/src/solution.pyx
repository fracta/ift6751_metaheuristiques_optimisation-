"""representation of a solution to the cvrp problem"""

# used in both tabu and genetic algorithm solver


cimport numpy as np
import numpy as np


cdef class Solution:
    """a solution to the cvrp problem"""

    def __init__(self, list routes, double score=-1):
        self.routes = routes
        self.score = score

    def __str__(self):
        return str(self.routes)

    def __repr__(self):
        return self.__str__()

    def __richcmp__(Solution self, Solution other not None, int op):
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

    cpdef Solution copy(Solution self):
        return Solution([route.copy() for route in self.routes], self.score)


cpdef np.ndarray get_solution_information(Solution sol, np.ndarray distance_matrix, np.ndarray weights):
    """returns [(distance, weight) for route in sol.routes]"""
    cdef np.ndarray information = np.empty(len(sol.routes), dtype=[('dist', float),('weight', float)])
    cdef int index
    for index in range(len(sol.routes)):
        information[index] = get_information(sol.routes[index], distance_matrix, weights)
    return information
