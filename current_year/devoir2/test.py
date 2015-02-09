import pyximport
import numpy as np
pyximport.install(reload_support=True)

import simple_ga, cvrp

data=cvrp.read_vrp("vrp_test.txt")
prob = cvrp.CVRPProblem(data)

route1 = simple_ga.Route (np.array([0,1,2,3,4,5,0]))
route2 = simple_ga.Route(np.array([0,1,3,2,5,4,0]))