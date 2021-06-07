#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 15:32:47 2021

@author: danielchapman
"""
##############################################################################      
############################################################################## 
##############################################################################      
##############################################################################  

## THIS SCRIPT VISUALIZES OUTPUT FROM THE GENETIC ALGORITHM IN deap_efel.py ##

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
from NeuronClassCN import Cell
from neuron import h
import matplotlib.pyplot as plt 
import efel
h.load_file('stdrun.hoc')

#%%
#%%
##############################################################################      
##############################################################################   
     
        ### Define the function
        
##############################################################################      
############################################################################## 
def Visualize_data(PathToParams):
    parameters = np.genfromtxt(PathToParams, delimiter=',')



    for i in range (0,len(parameters[1,:])):
        # m = Cell(soma_ra = 100,
        #           global_ra = 50,
        #           Ra_end = 36,
        #           Ra_dhalf = 142,
        #           Ra_steep = 143,
        #           cm = 1.2,
        #           Rm = 36182,
        #           Rm_end = 11730,
        #           Vleak = -67,
        #           gnaSoma = parameters[0,i],
        #           gnaSr = parameters[1,i],
        #           gkdr = parameters[2,i],
        #           gkap = parameters[3,i],
        #           gkad = parameters[4,i],
        #           soma_caL = parameters[5,i],
        #           soma_caR = parameters[6,i],
        #           soma_caN = parameters[7,i],
        #           soma_caT = parameters[8,i],
        #           soma_hbar = parameters[9,i],
        #           soma_km = parameters[10,i],
        #           soma_kdBG = parameters[11,i],
        #           soma_kca = parameters[12,i],
        #           dhalf = parameters[13,i],
        #           steep = parameters[14,i],
        #           gh_end = parameters[15,i],
        #           gh_soma = parameters[16,i],
        #           dlimit = parameters[17,i],
        #           dprox = parameters[18,i],
        #           dslope = parameters[19,i],
        #           okslope = parameters[20,i], CN = i) 
        
        # m = Cell(soma_ra = 100,
        #           global_ra = 50,
        #           Ra_end = 36,
        #           Ra_dhalf = 142,
        #           Ra_steep = 143,
        #           cm = 1.2,
        #           Rm = 36182,
        #           Rm_end = 11730,
        #           Vleak = -67,
        #           gnaSoma = 0.05,
        #           gnaSr = 0.0475,
        #           gkdr = 0.05,
        #           gkap = 0.0488,
        #           gkad = 0.008973,
        #           soma_caL = parameters[0,i],
        #           soma_caR = parameters[1,i],
        #           soma_caN = parameters[2,i],
        #           soma_caT = parameters[3,i],
        #           soma_hbar = parameters[4,i],
        #           soma_km = parameters[5,i],
        #           soma_kdBG = parameters[6,i],
        #           soma_kca = parameters[7,i],
        #           dhalf = parameters[8,i],
        #           steep = parameters[9,i],
        #           gh_end = parameters[10,i],
        #           gh_soma = parameters[11,i],
        #           dlimit = parameters[12,i],
        #           dprox = parameters[13,i],
        #           dslope = parameters[14,i],
        #           okslope = parameters[15,i], CN = i)

        # m = Cell(soma_ra = 100,
        #          global_ra = 60,
        #          Ra_end = 30,
        #          Ra_dhalf = 130,
        #          Ra_steep = 150,
        #          cm = 1,
        #          Rm = 30000,
        #          Rm_end = 12000,
        #          Vleak = -65,
        #          gnaSoma = parameters[0,i],
        #          gnaSr = parameters[1,i],
        #          gkdr = parameters[2,i],
        #          gkap = parameters[3,i],
        #          gkad = parameters[4,i],
        #          soma_caL = 0.000075,
        #          soma_caR = 0.00015,
        #          soma_caN = 0.00001,
        #          soma_caT = 0.00015,
        #          soma_hbar = 0.0001,
        #          soma_km = 0.000005,
        #          soma_kdBG = 0.0001,
        #          soma_kca = 0.00001,
        #          dhalf = 280,
        #          steep = 50,
        #          gh_end = 0.0009,
        #          gh_soma = 0.0001,
        #          dlimit = 330,
        #          dprox = 50,
        #          dslope = 0.065,
        #          okslope = 0.03, CN = i)       
        
        
        m = Cell(soma_ra = 100,
                  global_ra = 50,
                  Ra_end = 36,
                  Ra_dhalf = 142,
                  Ra_steep = 143,
                  cm = 1.2,
                  Rm = 36182,
                  Rm_end = 11730,
                  Vleak = -67,
                  gnaSoma = parameters[0,i],
                  gnaSr = parameters[1,i],
                  gkdr = parameters[2,i],
                  gkap = parameters[3,i],
                  gkad = parameters[4,i],
                  soma_caL = parameters[5,i],
                  soma_caR = parameters[6,i],
                  soma_caN = parameters[7,i],
                  soma_caT = parameters[8,i],
                  soma_hbar = parameters[9,i],
                  soma_km = parameters[10,i],
                  soma_kdBG = parameters[11,i],
                  soma_kca = parameters[12,i],
                  dhalf = 280,
                  steep = 50,
                  gh_end = 0.0009,
                  gh_soma = 0.0001,
                  dlimit = 330,
                  dprox = 50,
                  dslope = 0.065,
                  okslope = 0.03, CN = i) 
        clamp = h.IClamp(m.soma(0.5)) 
        stim_start = 100
        stim_end = 260
        clamp.amp = 0.1
        clamp.delay = stim_start
        clamp.dur = 160
        
        # stim = h.IClamp(m.soma(0.5)) # neuron 2
        # Stim = 0.1
        # stim.amp = Stim
        t = h.Vector().record(h._ref_t) # record time steps  
        # stim.delay = 100 # delay of stimulation (ms)
        # stim.dur = 400 # duration of stimulation (ms)  
        # create vectors to record
        # voltage vectors
        soma_v = h.Vector().record(m.soma(0.5)._ref_v)
        # dend_v = h.Vector().record(m.apic[10](0.5)._ref_v)
        # action potential vectors
        soma_AP = h.APCount(m.soma(0.5))

        h.celsius = 20
        h.finitialize(-70) # initial potential 
        h.continuerun(600)
        f = plt.figure()
        plt.plot(t,soma_v,linewidth=0.5)
        print(soma_AP)
        trace = {}
        trace['T'] = t
        trace['V'] = soma_v
        trace['stim_start'] = [stim_start]
        trace['stim_end'] = [stim_end]
        traces = [trace]
    
    # setup feature extraction 
        features = efel.getFeatureValues(traces, ["Spikecount",
                                              "Spikecount_stimint",
                                              "AP_amplitude","AP1_amp",
                                              "AP2_amp","APlast_amp",
                                              "AP_duration_half_width",
                                              "AP_width","AHP_depth_abs",
                                              "AHP_depth_abs_slow",
                                              "AHP_slow_time","AHP_depth",
                                              "AHP_time_from_peak",
                                              "voltage_base",
                                              "steady_state_voltage_stimend",
                                              "steady_state_voltage",
                                              "decay_time_constant_after_stim"
                                              ])
    



        Spikecount = features[0]["Spikecount"][0]
        decay_time = features[0]["decay_time_constant_after_stim"][0]
        print(Spikecount)
        print(decay_time)
    plt.show()
    
    # Stim = 0.18
    # # plot at different injection levels
    # f = plt.figure()
    # i = 0 
    # for x in np.linspace(0,Stim,9):
    #     clamp.amp = x 
    #     h.celsius = 20
    #     h.finitialize(-70) # initial potential 
    #     h.continuerun(600)
    #     plt.plot(t,soma_v)
    #     somaAPC[i] = soma_AP.n
    #     i = i + 1
    #     plt.show()

    #     # plot frequency curve
    #     HzFactor = (1000/clamp.dur)
    #     f2 = plt.figure()
    #     plt.plot(APCX, somaAPC*HzFactor,'bo')


#%%

# # visualize the parameters
PathToParams = '/Users/danielchapman/PythonDev/Final/Models/Results/CloudCluster500GenTestParams.csv'
Visualize_data(PathToParams)




