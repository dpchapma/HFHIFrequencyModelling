#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 17:00:21 2021

@author: danielchapman
"""

import numpy as np
from NeuronClassCN import Cell
from neuron import h
import matplotlib.pyplot as plt 
h.load_file('stdrun.hoc')


def Visualize_data(PathToParams):
    parameters = np.genfromtxt(PathToParams, delimiter=',')



    for i in range (0,len(parameters[1,:])):
        m = Cell(soma_ra = parameters[0,i],
                 global_ra =parameters[1,i],
                 Ra_end = parameters[2,i],
                 Ra_dhalf = parameters[3,i],
                 Ra_steep = parameters[4,i],
                 cm = parameters[5,i],
                 Rm = parameters[6,i],
                 Rm_end =parameters[7,i],
                 Vleak = parameters[8,i],
                 gnaSoma = 0,
                 gnaSr = 0,
                 gkdr = 0,
                 gkap = 0,
                 gkad = 0,
                 soma_caL = 0,
                 soma_caR = 0,
                 soma_caN = 0,
                 soma_caT = 0,
                 soma_hbar = 0.0000000001,
                 soma_km = 0,
                 soma_kdBG = 0,
                 soma_kca = 0,
                 dhalf = 0.00000000001,
                 steep = 0.00000000001,
                 gh_end = 0.00000000001,
                 gh_soma = 0.0000000001,
                 dlimit = 0,
                 dprox = 0,
                 dslope = 0,
                 okslope = 0, CN = i)  
        stim = h.IClamp(m.soma(0.5)) # neuron 2
        Stim = 0.04
        stim.amp = Stim
        t = h.Vector().record(h._ref_t) # record time steps  
        stim.delay = 100 # delay of stimulation (ms)
        stim.dur = 400 # duration of stimulation (ms)  
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
        plt.plot(t,soma_v)
        print(soma_AP)

    plt.show()

    # ap python vectors
    # somaAPC = np.empty(9)
    # APCX = np.linspace(0,stim*1000,9)
    



    # # plot at different injection leve
    # f = plt.figure()
    # i = 0 
    # for x in np.linspace(0,Stim,9):
    #     stim.amp = x 
    #     h.celsius = 20
    #     h.finitialize(-70) # initial potential 
    #     h.continuerun(600)
    #     plt.plot(t,soma_v)
    #     somaAPC[i] = soma_AP.n
    #     i = i + 1
    #     plt.show()

    #     # plot frequency curve
    #     HzFactor = (1000/stim.dur)
    #     f2 = plt.figure()
    #     plt.plot(APCX, somaAPC*HzFactor,'bo')


