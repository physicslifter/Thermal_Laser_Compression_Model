import numpy as np
import main

class diffEVAgent:
    def __init__(self, ID:str, parameters):
        self.ID=ID
        self.parameters=parameters

def differential_evolution(group:main.OptimizationGroup, 
                           maxiter:int,
                           population_size:int, 
                           optimization_data:main.OptimizationData,
                           bounds,
                           crossover_probability=0.9,
                           differential_weight=0.8):
    
    i=0
    my_func=group.optimizable_function
    #Check to ensure bounds dimensions agree with the parameter dimensions
    #While we have not yet reched the maximum number of iterations...
    
    #Initialize a population of agents
    for x in range(population_size):
        for param_bounds in bounds:
            param_value=np.random.uniform(param_bounds)
    
    while i <= maxiter:
        #initialize agents with random positions in the search space
        