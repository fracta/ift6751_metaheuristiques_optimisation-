# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="davide"
__date__ ="$3-ago-2010 16.32.33$"

from copy import *
from customer import *

def __sorting(x):
    #print(x);
    #print(len(x));
    return len(x);

class Clark_and_Wright():

    def __init__(self,dima,elma,customers,allocate_truck,trucks_number,depot_index):
        self.__dima=dima;
        self.__elma=elma;
        self.__customers=customers;
        self.__allocate_truck=allocate_truck;
        self.__truck_number=trucks_number;
        self.__depot_index=depot_index;
        self.__savings={};
        self.__routes={};
        self.__max_truck=allocate_truck();
        self.__max_capacity=self.__max_truck.max_load;
        self.__start={};
        self.__end={};
        #print(customers);
        #customers=deepcopy(customers);
        #customers=customers-customers[depot_index];
        customers=list(customers);
        #print(customers);
        scan=range(len(customers));
        for i in scan:
            self.__routes[i+1]=[customers[i]];
            self.__end[customers[i]]=self.__start[customers[i]]=i+1;
            for j in scan:
                #print(i,j);
                saving=dima[i+1][depot_index]+dima[depot_index][j+1]-dima[i+1][j+1];
                couple=(i+1,j+1);
                #print(couple);
                if saving in self.__savings:
                    self.__savings[saving].append(couple);
                else:
                    self.__savings[saving]=[couple];
        #print(self.__savings);
        assert(len(self.__start)==len(self.__end)==len(self.__routes)==len(customers));
        
        
    def is_feasible(self,list):
        #print(list);
        depot_object=customers[self.__depot_index];
        total_load=0;
        total_time=depot_object.opening;
        prev_index=self.__depot_index;
        for index in list:
            customer=customers[index];
            total_load+=customer.demand;
            total_time+=self.__elma[prev_index][index];
            total_time=max(total_time,customer.opening)+download_time(index);
        total_time+=self.__elma[len(list)-1][self.__depot_index];
        
        return ((total_load<self.__max_capacity) and (total_time<depot_object.time_frame));
    
    def find_starting_solution(self):
        savings=self.__savings.keys();
        savings=sorted(savings,reverse=True);
        for saving in savings:
            #print('saving={0}:{1}'.format(saving,self.__savings[saving]));
            for couple in self.__savings[saving]:
                #print(couple);
                
                if((couple[0] in self.__start) and (couple[1] in self.__end)):
                    #print('{0}:{1} in start, {2}:{3} in end!'.format(couple[0],self.__start[couple[0]],couple[1],self.__end[couple[1]]));
                    route_2=self.__start[couple[0]];
                    route=self.__end[couple[1]];
                        
                elif((couple[0] in self.__end) and (couple[1] in self.__start)):
                    #print('{2}:{3} in start, {0}:{1} in end!'.format(couple[0],self.__end[couple[0]],couple[1],self.__start[couple[1]]));
                    route_2=self.__start[couple[1]];
                    route=self.__end[couple[0]];
                else:
                    #print('{0} and {1} not present!'.format(couple[0],couple[1]));
                    continue;
                    
                if((route not in self.__routes) or (route_2 not in self.__routes)):
                    continue;
                #print('routes exist');
                if(route==route_2):
                    #print('same route');
                    #print(self.__routes[route]);
                    continue;
                #print('routes are different');
                if self.is_feasible(self.__routes[route]+self.__routes[route_2]):
                    #print("combination_feasible!");
                    #print('{0}+{1}='.format(self.__routes[route],self.__routes[route_2]));
                    self.__routes[route]+=self.__routes[route_2];
                    #print('{0}:{1}'.format(route,self.__routes[route]));
                    del self.__routes[route_2];
                    #print('route {0} deleted'.format(route_2));
                
                    #print(self.__start);
                    elements=[];
                    for element,r in self.__start.items():
                        if r==route_2:
                            elements.append(element);
                    #print('elements={0}'.format(elements));
                    #print('start={0}'.format(self.__start));
                    for element in elements:
                        del self.__start[element];
                    #print('start={0}'.format(self.__start));
                    
                    elements=[];
                    #print('end={0}'.format(self.__end));        
                    for element,r in self.__end.items():
                        if r==route_2:
                            self.__end[element]=route;
                        elif r==route:
                            elements.append(element);
                    #print('elements={0}'.format(elements));
                    for element in elements:
                        del self.__end[element];
                    #print('end={0}'.format(self.__end));
        
        #print(len(self.__routes.keys()));
        routes=[];
        for key,route in self.__routes.items():
            routes.append(route);
        self.__routes=routes;
        
        if(len(self.__routes)>self.__truck_number):
            
            routes=sorted(self.__routes,key=globals()['__sorting'],reverse=True);
            k=0;
            t=self.__truck_number;
            #print(routes);
            for route in range(t,len(routes)):
                routes[k]+=routes[route];
                k=(k+1)%t;
            #print(routes);
            routes=routes[0:t];
            #print(routes);
            self.__routes=routes;
                 
        #print(len(self.__routes));
        solution=[];
        for route in self.__routes:
            tour={};
            tour['truck']=self.__max_truck;
            self.__max_truck=self.__allocate_truck();
            tour['route']=route;
            tour['new_tabu']={};
            tour['inserted']={};
            tour['deleted']={};
            tour['old_tabu']={};
            solution.append(tour);
        return solution;