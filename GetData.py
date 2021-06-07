#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 14:53:25 2021

@author: danielchapman
"""

##############################################################################      
############################################################################## 
##############################################################################      
##############################################################################  

## THIS SCRIPT GET'S FEATURES FROM STEPH'S DATA TO IMPORT INTO EVALUATE     ##
## FUNCTION                                                                 ##
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

import efel
# import csv
import numpy as np
# import matplotlib.pyplot as plt 

    

#%%
#%%
##############################################################################      
##############################################################################   
     
        ### Define the function
        
##############################################################################      
############################################################################## 
def get_data(PathToFile):
    # Get data file 
    data_new = np.genfromtxt(PathToFile,delimiter=',')
    
    # setup stim parameters
    stim_start = 58.75
    stim_end = 218.75
    
    #Put into new variables for time and voltage
    time = data_new[:,0]
    voltage = data_new[:,1]
    
    # setup trace for efel input 
    trace = {}
    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]
    traces = [trace]
    
    # setup feature extraction 
    features = efel.getFeatureValues(traces, ["inv_time_to_first_spike",
                                              "inv_last_ISI","Spikecount",
                                              "Spikecount_stimint",
                                              "AP_amplitude","AP1_amp",
                                              "AP2_amp","APlast_amp",
                                              "voltage_base","AHP_depth_abs",
                                              "AHP_depth_abs_slow",
                                              "AHP_slow_time","AHP_depth",
                                              "AHP_time_from_peak",
                                              "AP_duration_half_width",
                                              "AP_width",
                                              "steady_state_voltage_stimend",
                                              "steady_state_voltage",
                                              "decay_time_constant_after_stim",
                                             ])
    
    # set variables for extracted features
    inv_time_to_first_spike = features[0]["inv_time_to_first_spike"][0]
    inv_last_ISI = features[0]["inv_last_ISI"][0]
    Spikecount = features[0]["Spikecount"][0]
    Spikecount_stimint = features[0]["Spikecount_stimint"][0]
    AP_amplitude = features[0]["AP_amplitude"][0]
    AP1_amp = features[0]["AP1_amp"][0]
    AP2_amp = features[0]["AP2_amp"][0]
    APlast_amp = features[0]["APlast_amp"][0]
    voltage_base = features[0]["voltage_base"][0]
    AHP_depth_abs = features[0]["AHP_depth_abs"][0]
    AHP_depth_abs_slow = features[0]["AHP_depth_abs_slow"][0]
    AHP_slow_time = features[0]["AHP_slow_time"][0]
    AHP_depth = features[0]["AHP_depth"][0]
    AHP_time_from_peak = features[0]["AHP_time_from_peak"][0]
    AP_duration_half_width = features[0]["AP_duration_half_width"][0]
    AP_width = features[0]["AP_width"][0]
    steady_state_voltage_stimend = features[0]["steady_state_voltage_stimend"][0]
    steady_state_voltage = features[0]["steady_state_voltage"][0]
    decay_time_constant_after_stim = features[0]["decay_time_constant_after_stim"][0]
    # ohmic_input_resistance= features[0]["ohmic_input_resistance"][0]
    
    # returen the features 
    # return inv_time_to_first_spike,inv_last_ISI,Spikecount,Spikecount_stimint,\
    # AP_amplitude,AP1_amp,AP2_amp,APlast_amp,voltage_base,AHP_depth_abs,\
    # AHP_depth_abs_slow,AHP_slow_time,AHP_depth,AHP_time_from_peak,\
    # AP_duration_half_width, AP_width,steady_state_voltage_stimend,\
    # steady_state_voltage,decay_time_constant_after_stim
    
    print(decay_time_constant_after_stim)
    
    return Spikecount,Spikecount_stimint,\
    AP_amplitude, voltage_base,AHP_depth_abs,\
    AHP_depth_abs_slow,AHP_slow_time,AHP_depth,AHP_time_from_peak,\
    AP_duration_half_width, AP_width,steady_state_voltage_stimend,\
    steady_state_voltage, decay_time_constant_after_stim


# example 
# data = get_data('/Users/danielchapman/Downloads/Dan Engram 12.15/ShamCC.csv')
# plt.plot(data[:,0],data[:,1])

PathToFile ='/Users/danielchapman/Desktop/SampleEphysFolder/CCIV/2021_02_15_0021.abf'








