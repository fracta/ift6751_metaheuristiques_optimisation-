
# test for RNA finder
#%pylab
import pyximport
import numpy as np
pyximport.install(reload_support=True)

def read(f):
    data = []
    with open(f) as F:
        name = ""
        seq = ""
        subopts = []
        lines = F.readlines()
        for line in lines:
            if line.startswith(">"):  # new name
                if name != "":
                    data.append((name, seq, subopts))
                name = line.strip()[1:]  # remove ">" char
                sequence = ""
                subopts = []
            s = line.strip()[0]
            if s in "AUGC":  # new sequence
                seq = line.strip()
            elif s in ".()":  # new subopt
                splitted = line.strip().split()
                subopt, energy = splitted[0], float(splitted[1])
                subopts.append((subopt, energy))
            else:
                continue
    return data


data = read("DM1-33.txt")

db = '((..))((((((((((((..))))).)))))))'

