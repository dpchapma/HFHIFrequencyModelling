#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 12:39:26 2021

@author: danielchapman
"""
##############################################################################      
############################################################################## 
##############################################################################      
##############################################################################  

## THIS SCRIPT CONTAINS EVALUAT FUNCTION FOR GENETIC ALGORITHM              ##
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
from neuron import h
h.load_file('stdrun.hoc')
import efel
from NeuronClass import Cell
import numpy as np

#%%
##############################################################################      
##############################################################################   
     
        ### DEFINE FUNCTION
        
##############################################################################      
############################################################################## 



# def evaluate(individual, PathToData = '/Users/danielchapman/Downloads/July 17-18 rmTBI 24Hr/ShamCC120mANew.csv', 
              # PathToData2 = '/Users/danielchapman/Downloads/July 17-18 rmTBI 24Hr/ShamCC180mANew.csv'):
def evaluate(individual, PathToData = '/home/dpc53/StephData/ShamCC120mANew.csv', 
              PathToData2 = '/home/dpc53/StephData/ShamCC180mANew.csv'):
    """
    Evaluate a neuron model with parameters e_pas and g_pas, extracts
    eFeatures from resulting traces and returns a tuple with
    abs(voltage_base-target_voltage1) and
    abs(steady_state_voltage-target_voltage2)
    """
    from GetData import get_data 
    data_new = get_data(PathToData) # from 0.1 stim amp
    data_new2 = get_data(PathToData2) # from 0.2 stim amp
    print(data_new)
    h.v_init = data_new[12]
    # Create cell class with parameters set as deap "individual" classes

    # m = Cell(soma_ra = 100,
    # global_ra = 50,
    # Ra_end = 36,
    # Ra_dhalf = 142,
    # Ra_steep = 143,
    # cm = 1.2,
    # Rm = 36182,
    # Rm_end = 11730,
    # Vleak = -67,
    # gnaSoma = individual[0],
    # gnaSr = individual[1],
    # gkdr = individual[2],
    # gkap = individual[3],
    # gkad = individual[4],
    # soma_caL = individual[5],
    # soma_caR = individual[6],
    # soma_caN = individual[7],
    # soma_caT = individual[8],
    # soma_hbar = individual[9],
    # soma_km = individual[10],
    # soma_kdBG = individual[11],
    # soma_kca = individual[12],
    # dhalf = individual[13],
    # steep = individual[14],
    # gh_end = individual[15],
    # gh_soma = individual[16],
    # dlimit = individual[17],
    # dprox = individual[18],
    # dslope = individual[19],
    # okslope = individual[20])
    
    # m = Cell(soma_ra = 100,
    # global_ra = 50,
    # Ra_end = 36,
    # Ra_dhalf = 142,
    # Ra_steep = 143,
    # cm = 1.2,
    # Rm = 36182,
    # Rm_end = 11730,
    # Vleak = -67,
    # gnaSoma = individual[0],
    # gnaSr = individual[1],
    # gkdr = individual[2],
    # gkap = individual[3],
    # gkad = individual[4],
    # soma_caL = individual[5],
    # soma_caR = individual[6],
    # soma_caN = individual[7],
    # soma_caT = individual[8],
    # soma_hbar = individual[9],
    # soma_km = individual[10],
    # soma_kdBG = individual[11],
    # soma_kca = individual[12],
    # dhalf = individual[13],
    # steep = individual[14],
    # gh_end = individual[15],
    # gh_soma = individual[16],
    # dlimit = individual[17],
    # dprox = individual[18],
    # dslope = individual[19],
    # okslope = individual[20])
    
    # m = Cell(soma_ra = 100,
    # global_ra = 60,
    # Ra_end = 30,
    # Ra_dhalf = 130,
    # Ra_steep = 150,
    # cm = 1,
    # Rm = 30000,
    # Rm_end = 12000,
    # Vleak = -65,
    # gnaSoma = individual[0],
    # gnaSr = individual[1],
    # gkdr = individual[2],
    # gkap = individual[3],
    # gkad = individual[4],
    # soma_caL = 0.000075,
    # soma_caR = 0.00015,
    # soma_caN = 0.00001,
    # soma_caT = 0.00015,
    # soma_hbar = 0.0001,
    # soma_km = 0.000005,
    # soma_kdBG = 0.0001,
    # soma_kca = 0.00001,
    # dhalf = 280,
    # steep = 50,
    # gh_end = 0.0009,
    # gh_soma = 0.0001,
    # dlimit = 330,
    # dprox = 50,
    # dslope = 0.065,
    # okslope = 0.03)
    
    
    m = Cell(soma_ra = 100,
    global_ra = 50,
    Ra_end = 36,
    Ra_dhalf = 142,
    Ra_steep = 143,
    cm = 1.2,
    Rm = 36182,
    Rm_end = 11730,
    Vleak = -67,
    gnaSoma = individual[0],
    gnaSr = individual[1],
    gkdr = individual[2],
    gkap = individual[3],
    gkad = individual[4],
    soma_caL = individual[5],
    soma_caR = individual[6],
    soma_caN = individual[7],
    soma_caT = individual[8],
    soma_hbar = individual[9],
    soma_km = individual[10],
    soma_kdBG = individual[11],
    soma_kca = individual[12],
    dhalf = 280,
    steep = 50,
    gh_end = 0.0009,
    gh_soma = 0.0001,
    dlimit = 330,
    dprox = 50,
    dslope = 0.065,
    okslope = 0.03)
    

    from Run_sim import Run_sim
    SimData = Run_sim(StimAmp = 0.1, m = m)
    SimData2 = Run_sim(StimAmp = 0.18, m = m)
    
    
    
    # setup fitneess function 
    # Fitness = np.empty([2,12])
    # for i in range (0,12):
    #     Fitness[0,i] = abs((data_new[i]-SimData[i])/data_new[i])
    #     Fitness[1,i] = abs((data_new2[i]-SimData2[i])/data_new2[i])
    
    # FitnessValues = np.reshape(Fitness,(24,1))
    # FitnessValues = FitnessValues.tolist()
    
    # print(FitnessValues)
    # return FitnessValues 
    
    
    return abs((data_new[0]-SimData[0])/data_new[0]), abs((data_new2[0]-SimData2[0])/data_new2[0]),\
        abs((data_new[1]-SimData[1])/data_new[1]), abs((data_new2[1]-SimData2[1])/data_new2[1]),\
        abs((data_new[2]-SimData[2])/data_new[2]), abs((data_new2[2]-SimData2[2])/data_new2[2]),\
        abs((data_new[3]-SimData[3])/data_new[3]), abs((data_new2[3]-SimData2[3])/data_new2[3]), \
        abs((data_new[4]-SimData[4])/data_new[4]), abs((data_new2[4]-SimData2[4])/data_new2[4]), \
        abs((data_new[5]-SimData[5])/data_new[5]), abs((data_new2[5]-SimData2[5])/data_new2[5]), \
        abs((data_new[6]-SimData[6])/data_new[6]), abs((data_new2[6]-SimData2[6])/data_new2[6]), \
        abs((data_new[7]-SimData[7])/data_new[7]), abs((data_new2[7]-SimData2[7])/data_new2[7]), \
        abs((data_new[8]-SimData[8])/data_new[8]), abs((data_new2[8]-SimData2[8])/data_new2[8]), \
        abs((data_new[9]-SimData[9])/data_new[9]), abs((data_new2[9]-SimData2[9])/data_new2[9]), \
        abs((data_new[10]-SimData[10])/data_new[10]), abs((data_new2[10]-SimData2[10])/data_new2[10]), \
        abs((data_new[11]-SimData[11])/data_new[11]), abs((data_new2[11]-SimData2[11])/data_new2[11]), \
        abs((data_new[12]-SimData[12])/data_new[12]), abs((data_new2[12]-SimData2[12])/data_new2[12]), \
        abs((data_new[13]-SimData[13])/data_new[13]), abs((data_new2[13]-SimData2[13])/data_new2[13])




        
        









        
        
    