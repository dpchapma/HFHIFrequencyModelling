#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 12:43:31 2021

@author: danielchapman
"""
##############################################################################      
############################################################################## 
##############################################################################      
##############################################################################  

## THIS SCRIPT EVALUATES GENETIC ALGORITHM                                  ##
## code modified from bluebrain EFEL github                                 ##

##############################################################################      
############################################################################## 
##############################################################################      
##############################################################################

#%% 
##############################################################################      
##############################################################################   
     
        ### Import Modules
        
##############################################################################      
############################################################################## 
import numpy as np
import matplotlib.pyplot as plt 
from NeuronClass import Cell
import random
import numpy
import deap
import deap.gp
import deap.benchmarks
from deap import base
from deap import creator
from deap import tools
from deap import algorithms
from VisualizeFinalPopData import Visualize_data
from scoop import futures

#%%
##############################################################################      
##############################################################################   
     
        ### SETUP ALGORITHM PARAMETERS
        
##############################################################################      
############################################################################## 
random.seed(1)

POP_SIZE = 500
# Number of offspring in every generation
OFFSPRING_SIZE = 500
# Number of generations
NGEN = 100
# select number of offspring for next gen
MU = 100
LAMBDA = OFFSPRING_SIZE
# Crossover probability
CXPB = 0.3
# Mutation probability, should sum to one together with CXPB
MUTPB = 0.7
# Eta parameter of cx and mut operators
ETA = 10.0

IND_SIZE = 13

##############################################################################      
##############################################################################   
     
        ### SET BOUNDS FOR EACH CHANNEL PARAMETER
        
##############################################################################      
############################################################################## 
# LOWER = [75, 50,20,110,130,0.5,20000,10000,-80,0.001,0.001,0.001,0.001,0.001,
#          0.00001,0.00001,0.00001,0.00001,0.00001,0.0000001,0.00001,0.000001,
#          250,40,0.00009,0.00001,300,30,0.01,0.01]
# UPPER = [125, 70,40,150,170,1.5,30000,14000,-60,0.05,0.05,0.05,0.05,0.05,0.0001,
#          0.0001,0.0001,0.0001,0.005,0.00001,0.0001,0.00005,300,60,0.0009,0.0005,
#          350,70,0.07,0.1]

# LOWER = [0.00001,0.00001,0.00001,0.00001,0.00001,0.0000001,0.00001,0.000001,
#           250,40,0.00009,0.00001,300,30,0.01,0.01]
# UPPER = [0.0001,0.0001,0.0001,0.0001,0.005,0.00001,0.0001,0.00005,300,60,0.0009,0.0005,
#           350,70,0.07,0.1]

# LOWER = [0.001,0.001,0.001,0.001,0.001,
#           0.00001,0.00001,0.00001,0.00001,0.00001,
#           250,40,0.00009,0.00001,300,30,0.01,0.01]
# UPPER = [0.05,0.05,0.05,0.05,0.05,0.0001,
#           0.0001,0.0001,0.0001,0.005,300,60,0.0009,0.0005,
#           350,70,0.07,0.1]
# LOWER = [0.00001,0.00001,0.00001,0.00001,0.00001,0.0000001,0.00001,0.000001,
#           250,40,0.00009,0.00001,300,30,0.01,0.01]
# UPPER = [0.0001,0.0001,0.0001,0.0001,0.005,0.00001,0.0001,0.00005,300,60,0.0009,0.0005,
#           350,70,0.07,0.1]

    
LOWER = [0.001,0.001,0.001,0.001,0.001,
          0.00005,     # soma_caL = 0.000075,
          0.0001,      # soma_caR = 0.00015,
          0.000009,    # soma_caN = 0.00001,
          0.0001,      # soma_caT = 0.00015,
          0.000075,     # soma_hbar = 0.0001,
          0.000003,    # soma_km = 0.000005,
          0.000075,     # soma_kdBG = 0.0001,
          0.0000075]    # soma_kca = 0.00001,

UPPER = [0.05,0.05,0.05,0.05,0.05,
         0.0001,     # soma_caL = 0.000075,
         0.0002,     # soma_caR = 0.00015,
         0.00002,    # soma_caN = 0.00001,
         0.0002,     # soma_caT = 0.00015,
         0.0002,     # soma_hbar = 0.0001,
         0.000007,    # soma_km = 0.000005,
         0.00025,     # soma_kdBG = 0.0001,
         0.000025]    # soma_kca = 0.00001,

# LOWER = [0.001,0.001,0.001,0.001,0.001]
# UPPER = [0.05,0.05,0.05,0.05,0.05]

##############################################################################      
##############################################################################   
     
        ### SET TOOLBOXES
        
##############################################################################      
############################################################################## 
FeatureNumber = 28

SELECTOR = "NSGA2"
creator.create("Fitness", base.Fitness, weights=[-1.0] * FeatureNumber)
creator.create("Individual", list, fitness=creator.Fitness)

def uniform(lower_list, upper_list, dimensions):
    """Fill array """

    if hasattr(lower_list, '__iter__'):
        return [random.uniform(lower, upper) for lower, upper in
                zip(lower_list, upper_list)]
    else:
        return [random.uniform(lower_list, upper_list)
                for _ in range(dimensions)]
    
    
    
# def pop_start(PathToFile):
#     pop_params = np.genfromtxt(PathToFile, delimiter=',')
#     pop_start
#     return pop_params[:,4]


toolbox = base.Toolbox()    
toolbox.register("uniformparams", uniform, LOWER, UPPER, IND_SIZE)
toolbox.register("map", futures.map)
# toolbox.register("pop_start",pop_start,'/Users/danielchapman/PythonDev/Final/Models/01.28.21Params100GenTwoStimCaParam.csv')


toolbox.register(
    "Individual",
    tools.initIterate,
    creator.Individual,
    toolbox.uniformparams)


# toolbox.register(
#     "Individual",
#     tools.initIterate,
#     creator.Individual,
#     toolbox.pop_start)

toolbox.register("population", tools.initRepeat, list, toolbox.Individual)

from EvaluateFunction import evaluate
toolbox.register("evaluate", evaluate)

toolbox.register(
    "mate",
    deap.tools.cxSimulatedBinaryBounded,
    eta=ETA,
    low=LOWER,
    up=UPPER)
toolbox.register("mutate", deap.tools.mutPolynomialBounded, eta=ETA,
                 low=LOWER, up=UPPER, indpb=0.1)
toolbox.register(
    "select",
    tools.selNSGA2)



first_stats = tools.Statistics(key=lambda ind: ind.fitness.values[0])
second_stats = tools.Statistics(key=lambda ind: ind.fitness.values[1])
third_stats = tools.Statistics(key=lambda ind: ind.fitness.values[2])
fourth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[3])
fifth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[4])
sixth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[5])
seventh_stats = tools.Statistics(key=lambda ind: ind.fitness.values[6])
eighth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[7])
ninth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[8])
tenth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[9])
eleventh_stats = tools.Statistics(key=lambda ind: ind.fitness.values[10])
twelfth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[11])
thirteenth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[12])
fourteenth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[13])
fifteenth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[14])
sixteenth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[15])
seventeenth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[16])
eighteenth = tools.Statistics(key=lambda ind: ind.fitness.values[17])
nineteenth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[18])
twentieth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[19])
twentyfirst_stats = tools.Statistics(key=lambda ind: ind.fitness.values[20])
twentysecond_stats = tools.Statistics(key=lambda ind: ind.fitness.values[21])
twentythird_stats = tools.Statistics(key=lambda ind: ind.fitness.values[22])
twentyfourth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[23])
twentyfifth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[24])
twentysixth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[25])
twentyseventh_stats = tools.Statistics(key=lambda ind: ind.fitness.values[26])
twentyeigth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[27])
# twentyninth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[28])
# thirtieth_stats = tools.Statistics(key=lambda ind: ind.fitness.values[29])

# stats = tools.MultiStatistics(SC_StimInt = first_stats, SC = second_stats, AP_amp = third_stats,
#                                 APHW = fourth_stats, 
#                                 APW = fifth_stats,
#                                 AHP_abs = sixth_stats,
#                                 AHP_slow = seventh_stats, AHP_slow_time = eighth_stats,
#                                 AHP_depth = ninth_stats, AHP_timeFP =  tenth_stats, 
#                                 voltage_base = eleventh_stats, SS_stimend = twelfth_stats,
#                                 SS = thirteenth_stats,
#                                 obj15 = fifteenth_stats,
#                                obj16 = sixteenth_stats, obj17 = seventeenth_stats, 
#                                obj18 = eighteenth, obj19 = nineteenth_stats, 
#                                obj20 = twentieth_stats, obj21 = twentyfirst_stats, 
#                                obj22 = twentysecond_stats, obj23 = twentythird_stats, 
#                                obj24 = twentyfourth_stats)
                               # , obj25 = twentyfifth_stats, 
                              # obj26 = twentysixth_stats, obj27 = twentyseventh_stats, 
                              # obj28 = twentyeigth_stats, obj29 = twentyninth_stats) 


                    
stats = tools.MultiStatistics(obj1 = first_stats, obj2 = second_stats,
                                obj3 = third_stats, 
                                obj4 = fourth_stats,
                                obj5 = fifth_stats,
                                obj6 = sixth_stats,
                                obj7 = seventh_stats,
                                obj8 = eighth_stats, obj9 =  ninth_stats, 
                                obj10 = tenth_stats, obj11 = eleventh_stats, 
                                obj12 = twelfth_stats,obj13 = thirteenth_stats,
                                obj14 = fourteenth_stats, obj15 = fifteenth_stats,
                                obj16 = sixteenth_stats, obj17 = seventeenth_stats, 
                                obj18 = eighteenth, obj19 = nineteenth_stats, 
                                obj20 = twentieth_stats, obj21 = twentyfirst_stats, 
                                obj22 = twentysecond_stats, obj23 = twentythird_stats, 
                                obj24 = twentyfourth_stats,
                                obj25 = twentyfifth_stats, 
                                obj26 = twentysixth_stats,
                                obj27 = twentyseventh_stats, 
                               obj28 = twentyeigth_stats)
                              # obj29 = twentyninth_stats) 

logbook = tools.Logbook()
pop = toolbox.population(n=MU)

# def best_sum(d):
#     "code adapted from https://www.nature.com/articles/s41467-017-02718-3#Sec10"
#     """ Return the lowest sum across columns
#     Parameters
#     ----------
#     d : (n, p) array
#         Array with `n` samples each with `p` feature errors
#     Returns
#     -------
#     float
#         Minimum feature error sum across samples in `d`
#     """
#     return np.sum(d, axis=1).min()

stats.register("min", numpy.min)
stats.register("avg",numpy.mean)

# stats.register("best", best_sum)

#%% have to run everything from a main to allow for parallel processing 
if __name__ == '__main__':
    # run the algorithm
    pop, logbook = algorithms.eaMuPlusLambda(
        pop,
        toolbox,
        MU,
        LAMBDA,
        CXPB,
        MUTPB,
        NGEN,
        stats,
        halloffame=None)
    
    # save the parameters
    params = np.empty([IND_SIZE,MU])
    for i in range(0,MU):
        params[0:IND_SIZE,i] = pop[i]
    np.savetxt("CloudCluster500GenTestParams.csv",params,delimiter=",")
    Gen = logbook.chapters['obj1'].select('gen')
    # save error data
    ErrorMins = np.empty([NGEN+1,FeatureNumber])
    ErrorAvg = np.empty([NGEN+1,FeatureNumber])
    for x in range(1,FeatureNumber):
        ErrorMins[:,x] = logbook.chapters['obj{}'.format(x)].select('min')
        ErrorAvg[:,x] = logbook.chapters['obj{}'.format(x)].select('avg')
    np.savetxt("CloudCluster500GenTestMins.csv",ErrorMins,delimiter=",")
    np.savetxt("CloudCluster500GenTestAvg.csv",ErrorAvg,delimiter=",")

    

    
    
    
    
    
    
    
    
    



    





