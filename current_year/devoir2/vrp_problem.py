"""vehicule routing problem parameters (general data holder)"""



import csv


def read_vrp(file_name):
    """ """
    data = []
    with open(file_name) as F:
        reader = csv.reader(F, delimiter=' ')
        for row in reader:
            data.append(row[1:])
    return data



class VrpProblem(object):
    
    def __init__(self, data_table):
        self.num_clients = data_table[0][0]
        self.vehicule_capacity = data_table[0][1]
        self.depot_position = (data_table[1][0], data_table_[1][1])
        self.client_positions = data_table[2:]
    
    def get_num_clients(self):
        return self.num_clients
    
    def get_vehicule_capacity(self):
        return self.vehicule_capacity
    
    def get_depot_position(self):
        return self.depot_position
    
    def get_client_positions(self):
        return self.client_positions


