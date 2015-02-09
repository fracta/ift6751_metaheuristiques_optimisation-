"""vehicule routing problem main description (general data holder)"""

import csv
import cython
import numpy as np
cimport numpy as np


cpdef list read_vrp(str file_name):
    """file reader for the text file format """
    cdef list data = []
    with open(file_name) as F:
        reader = csv.reader(F, delimiter=' ')
        for row in reader:
            data.append(row[1:])
    return data


cpdef np.ndarray get_distance_matrix(client_positions):
    """compute the euclidean distance matrix between all points"""
    cdef int num_clients = len(client_positions)
    print client_positions
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
    cdef int num_clients
    cdef double vehicule_capacity
    cdef tuple depot_position
    cdef np.ndarray client_positions
    cdef np.ndarray weights
    cdef str name
    cdef np.ndarray distances
    
    def __init__(self, data_table, name=""):
        self.num_clients = int(data_table[0][0])
        self.vehicule_capacity = float(data_table[0][1])
        self.depot_position = (float(data_table[1][0]), float(data_table[1][1]))
        self.client_positions = np.zeros(self.num_clients, dtype=[("x", float), ("y", float)])
        self.weights = np.zeros(self.num_clients)
        for (i, (x_coord, y_coord, capacity)) in enumerate(data_table[2:]):
            self.client_positions[i] = (float(x_coord), float(y_coord))
            self.weights[i] = float(capacity)
        self.name = name
        print self.client_positions
        self.distances = get_distance_matrix(self.client_positions)
    
    def get_num_clients(self):
        return self.num_clients
    
    def get_vehicule_capacity(self):
        return self.vehicule_capacity
    
    def get_depot_position(self):
        return self.depot_position
    
    def get_client_positions(self):
        return self.client_positions
    
    def get_weights(self):
        return self.weights
    
    def get_distances(self):
        return self.distances
    
    def __str__(self):
        return "{0}\nConstrained Vehicule Routing Problem (CVRP)\n{1} clients, Q = {2}\n".format(self.name, self.num_clients, self.vehicule_capacity)
    
    def __repr__(self):
        return self.__str__()

