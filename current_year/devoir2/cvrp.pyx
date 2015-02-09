"""vehicule routing problem main description (general data holder)"""

import csv
import cython
import numpy as np
cimport numpy as np


cpdef list read_vrp(str file_problem_name):
    """file reader for the text file format """
    cdef list data = []
    with open(file_problem_name) as F:
        reader = csv.reader(F, delimiter=' ')
        for row in reader:
            data.append(row[1:])
    return data


cpdef np.ndarray get_distance_matrix(client_positions):
    """compute the euclidean distance matrix between all points"""
    cdef int num_clients = len(client_positions)
    cdef np.ndarray matrix = np.zeros((num_clients, num_clients))
    for i in range(num_clients):
        for j in range(i, num_clients):
            diff_x = client_positions[i][0] - client_positions[j][0]
            diff_x *= diff_x
            diff_y = client_positions[i][1] - client_positions[j][1]
            diff_y *= diff_y
            diff_y += diff_x
            diff_y = np.sqrt(diff_y)
            matrix[i][j] = diff_y
            matrix[j][i] = diff_y
    return matrix


cdef class CVRPProblem:
    """data for the constrained vrp problem"""

    cdef double vehicule_capacity
    cdef np.ndarray positions
    cdef np.ndarray weights
    cdef str problem_name
    cdef np.ndarray distance_matrix
    
    def __init__(self, data_table, problem_name=""):
        self.vehicule_capacity = float(data_table[0][1])
        num_clients = int(data_table[0][0])
        self.positions = np.zeros(num_clients+1, dtype=[("x", float), ("y", float)])
        self.weights = np.zeros(num_clients+1)
        self.positions[0] = (float(data_table[1][0]), float(data_table[1][1]))
        for (i, (x_coord, y_coord, capacity)) in enumerate(data_table[2:]):
            self.positions[i+1] = (float(x_coord), float(y_coord))
            self.weights[i+1] = float(capacity)
        self.problem_name = problem_name
        self.distance_matrix = get_distance_matrix(self.positions)
    
    def get_vehicule_capacity(self):
        return self.vehicule_capacity
    
    def get_positions(self):
        return self.positions
    
    def get_weights(self):
        return self.weights
    
    def get_distance_matrix(self):
        return self.distance_matrix
    
    def __str__(self):
        return "{0}\nConstrained Vehicule Routing Problem (CVRP)\n{1} clients, Q = {2}\n".format(self.problem_name, self.distance_matrix.shape[0]-1, self.vehicule_capacity)
    
    def __repr__(self):
        return self.__str__()

