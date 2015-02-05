"""reader for vehicule routing with capacity constraint text file format"""

import csv

def read_vrp(file_name):
    """ """
    data = []
    with open(file_name) as F:
        reader = csv.reader(F, delimiter=' ')
        for row in reader:
            data.append(row[1:])
    return data
