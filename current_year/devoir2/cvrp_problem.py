"""vehicule routing problem main description (general data holder)"""

import csv
import cython
import numpy as np


def read_vrp(file_name):
    """file reader for the text file format """
    data = []
    with open(file_name) as F:
        reader = csv.reader(F, delimiter=' ')
        for row in reader:
            data.append(row[1:])
    return data


@cython.ccall
@cython.returns(np.ndarray)
@cython.locals(client_locations=np.ndarray, num_clients=cython.int, matrix=np.ndarray)
def get_distance_matrix(client_locations):
    """compute the euclidean distance matrix between all points"""
    num_clients = len(client_locations)
    matrix = np.zeros((num_clients, num_clients))
    for i in range(num_clients):
        for j in range(i, num_clients):
            diff_x = client_locations[i][0] - client_locations[j][0]
            diff_x *= diff_x
            diff_y = client_locations[i][1] - client_locations[j][1]
            diff_y *= diff_y
            diff_y += diff_x
            diff_y = np.sqrt(diff_y)
            matrix[i][j] = diff_y
            matrix[j][i] = diff_y
    return matrix


@cython.cclass
@cython.locals(num_clients=cython.int,
               vehicule_capacity=cython.double,
               depot_position=tuple,
               client_positions=np.ndarray,
               weights=np.ndarray,
               name=str,
               distances=np.ndarray)
class CVRPProblem(object):
    """data for the constrained vrp problem"""
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


