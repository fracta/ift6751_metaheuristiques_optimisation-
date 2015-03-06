
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


#def show_problem(points):
    #"""returns a handle on a plt figure"""
    #depot = points[0]
    #clients = points[1:]
    #f = plt.figure()
    #ax = f.add_subplot(111)
    #ax.plot(depot['x'], depot['y'], 'o', color='r')
    #ax.plot(clients['x'], clients['y'], 'o', color='b')
    #return f


def show_routes(coordinates, routes):
    # plot the clients
    plt.scatter(coordinates[1:]['x'], coordinates[1:]['y']) # 
    plt.scatter(coordinates[0]['x'], coordinates[0]['y'], marker='s', color='r') # depot
    cmap = matplotlib.cm.Set1

    # plot the routes
    num_routes = len(routes)
    for (route_index, route) in enumerate(routes):
        path = np.zeros(len(route.nodes), dtype=[("x", float), ("y", float)])
        for (i, e) in enumerate(route.nodes):
            path[i] = coordinates[e]
        plt.plot(path['x'], path['y'], color=cmap(route_index / float(num_routes)))

    # show the plot
    plt.show()
    return


#def show_hull(coordinates, hull_coordinates):
    #f = plt.figure()
    #ax=f.add_subplot(111)
    
    #ax.scatter(coordinates['x'], coordinates['y'])
    #for i in hull_coordinates:
        #ax.scatter(i['x'], i['y'], color='green')
    #plt.show()


#def highlight_points(ax, x, y, color):
    #ax.scatter(x, y, 'o', color=color)
    #return
  
## from stackoverflow
#import numpy as np
#import matplotlib.pyplot as plt

#plt.axis([0, 1000, 0, 1])
#plt.ion()
#plt.show()

#for i in range(1000):
    #y = np.random.random()
    #plt.scatter(i, y)
    #plt.pause(0.01)