#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 12:34:22 2021

@author: danielchapman
"""

##############################################################################      
############################################################################## 
##############################################################################      
##############################################################################  

## THIS SCRIPT RUNS SIM ON MY CELL AND GETS SAME FEATURES WE GET FROM REAL  ##
## CELL. TAKES STIMULUS AMPLITUDE AS ARGUMENT                               ##                           ##

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

#%%
##############################################################################      
##############################################################################   
     
        ### Define the function
        
##############################################################################      
############################################################################## 

def Run_sim(StimAmp, m):

    # setup stimulation paraemeters 
    clamp = h.IClamp(m.soma(0.5)) 
    stim_start = 110
    stim_end = 270
    clamp.amp = StimAmp
    clamp.delay = stim_start
    clamp.dur = 160
    
    # recording paremeters 
    voltage = h.Vector()
    voltage.record(m.soma(0.5)._ref_v)
    h.celsius = 20
    time = h.Vector()
    time.record(h._ref_t)
    h.tstop = 300
    h.run()
    
    # setup trace for efel input for feature extraction
    trace = {}
    trace['T'] = time
    trace['V'] = voltage
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

    
    # inv_time_to_first_spike = features[0]["inv_time_to_first_spike"][0]
    # inv_last_ISI = features[0]["inv_last_ISI"][0]
    Spikecount = features[0]["Spikecount"][0]
    Spikecount_stimint = features[0]["Spikecount_stimint"][0]
    print(Spikecount)
    
    AP_amplitude = features[0]["AP_amplitude"]
    if AP_amplitude is None:
        AP_amplitude = 0
    else:
        AP_amplitude = features[0]["AP_amplitude"][0]
        
    # AP1_amp = features[0]["AP1_amp"][0]
    # AP2_amp = features[0]["AP2_amp"][0]
    
    # APlast_amp = features[0]["APlast_amp"]
    # if APlast_amp is None:
    #     APlast_amp = 0
    # else:
    #     APlast_amp = features[0]["APlast_amp"] 
    
    AHP_depth_abs_slow = features[0]["AHP_depth_abs_slow"]
    if AHP_depth_abs_slow is None:
        AHP_depth_abs_slow = 0
    else: 
        AHP_depth_abs_slow = features[0]["AHP_depth_abs_slow"][0]     
    
    AHP_depth_abs = features[0]["AHP_depth_abs"]
    if AHP_depth_abs is None:
        AHP_depth_abs = 0
    else:
        AHP_depth_abs = features[0]["AHP_depth_abs"][0]
    
    AHP_slow_time = features[0]["AHP_slow_time"]
    if AHP_slow_time is None:
        AHP_slow_time = 0
    else:
        AHP_slow_time = features[0]["AHP_slow_time"][0]
    
    AHP_depth = features[0]["AHP_depth"]
    if AHP_depth is None:
        AHP_depth = 0
    else:
        AHP_depth = features[0]["AHP_depth"][0]
    
    AHP_time_from_peak = features[0]["AHP_time_from_peak"]
    if AHP_time_from_peak is None:
        AHP_time_from_peak = 0
    else:
        AHP_time_from_peak = features[0]["AHP_time_from_peak"][0]
    
    
    AP_duration_half_width = features[0]["AP_duration_half_width"]
    if AP_duration_half_width is None:
        AP_duration_half_width = 0
    else:
        AP_duration_half_width = features[0]["AP_duration_half_width"][0]
    
    AP_width = features[0]["AP_width"]
    if AP_width is None:
        AP_width = 0
    else:
        AP_width = features[0]["AP_width"][0]
        
        
    voltage_base = features[0]["voltage_base"][0]
    steady_state_voltage_stimend = features[0]["steady_state_voltage_stimend"][0]
    steady_state_voltage = features[0]["steady_state_voltage"][0]
    decay_time_constant_after_stim = features[0]["decay_time_constant_after_stim"][0]
    # ohmic_input_resistance= features[0]["ohmic_input_resistance"][0]
    
    # return the features
    return Spikecount,Spikecount_stimint,\
    AP_amplitude, voltage_base,AHP_depth_abs,\
    AHP_depth_abs_slow,AHP_slow_time,AHP_depth,AHP_time_from_peak,\
    AP_duration_half_width, AP_width,steady_state_voltage_stimend,\
    steady_state_voltage, decay_time_constant_after_stim
    
    
     
    
    
    
    
    
    ##
    
    
    
    
    
    