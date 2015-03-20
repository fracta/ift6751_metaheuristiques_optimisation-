"""representation of a solution to the cvrp problem"""

# used in both tabu and genetic algorithm solver

import copy

cimport numpy as np
import numpy as np


cdef class Solution:
    """a solution to the cvrp problem"""

    def __init__(self, list routes, double score=-1):
        self.routes = copy.copy(routes)
        self.score = copy.copy(score)

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

    cpdef np.ndarray get_information(Solution self,
                                     np.ndarray distance_matrix,
                                     np.ndarray weights):
        """returns [(distance, weight) for route in sol.routes]"""
        cdef np.ndarray info = np.empty(len(self.routes), dtype=[('dist', float),('weight', float)])
        cdef Route route
        cdef int index
        cdef double dist, weight
        for index, route in enumerate(self.routes):
            dist = route.get_distance(distance_matrix)
            weight = route.get_weight(weights)
            info[index] = (dist, weight)
        return info


cpdef list get_centroids(Solution sol, np.ndarray positions):
    """get the centroid for each route of the solution"""
    cdef list centroids = []
    cdef double x_avg, y_avg
    cdef int num_clients
    for route in sol.routes:
        x_avg = 0.
        y_avg = 0.
        num_clients = len(route.nodes[1:-1])
        if num_clients == 0:
            centroids.append((0, 0))
        else:
            for client in route.nodes[1:-1]:
                position = positions[client]
                x_avg += position[0]
                y_avg += position[1]
            x_avg /= num_clients
            y_avg /= num_clients
            centroids.append((x_avg, y_avg))
    return centroids


cpdef double get_angle(double x1, double y1, double x2, double y2):
    """figure out the angle between two points (atan2 style)"""
    cdef double x_diff = x2 - x1
    cdef double y_diff = y2 - y1
    cdef double degree = np.degrees(np.arctan2(y_diff, x_diff))
    if degree < 0:
        degree = 360. + degree
    return degree


cpdef list sort_routes_by_angle(Solution sol, np.ndarray positions):
    """sort the routes of a solution by their centroid angle to the depot"""
    cdef list centroids = get_centroids(sol, positions)
    cdef list angles = []
    cdef double x1, y1
    for x1, y1 in centroids:
        angles.append(get_angle(positions[0][0], positions[0][1], x1, y1))
    sol.routes = [route for (angle, route) in sorted(zip(angles, sol.routes))]
    return sorted(angles)
