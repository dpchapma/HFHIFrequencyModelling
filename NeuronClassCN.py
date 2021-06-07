#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 08:24:03 2021

@author: danielchapman
"""

##############################################################################      
############################################################################## 
##############################################################################      
##############################################################################  

## THIS SCRIPT CONTAINS THE CELL CLASS USING MORPHOLOGY, INSERTING CHANNELS,##
## AND TAKES PARAMETERS LISTED IN THE __init__                              ##

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
from neuron import h,load_mechanisms,hclass
import numpy as np
h.load_file('stdrun.hoc')


#%% import the morphology 
##############################################################################      
##############################################################################   
     
        ### Import morohology and mechanisms 
        
##############################################################################      
############################################################################## 
######### UNCOMMENT OUT THIS LINE IN THE FIRST TIME RUNNING THE CODE #########
######### TO LOAD THE MECHANISM, THEN YOU CAN COMMENT IT OUT SO      #########
######### THE KERNEL DOENSN' NEED TO BE RESTARTED EACH TIME IT'S RUN #########
######### OR JUST RESTART THE KERNEL OR ELSE AN ERROR WILL BE THROWN #########
load_mechanisms('/Users/danielchapman/PythonDev/Final/Mechanism/')
h.load_file('/Users/danielchapman/PythonDev/Final/Morphology/Neuron3Oblique.hoc')


#%% 
##############################################################################      
##############################################################################   
     
        ### Create cell class from NEURON template and setup biophysical 
        ### framework 
        
##############################################################################      
############################################################################## 
class Cell(hclass(h.Cell)):
    conductance_based = True
    parameter_names = []
    # def __init__(self,**parameters): # intialize the cell
    def __init__(self,soma_ra,global_ra,Ra_end,Ra_dhalf,Ra_steep,cm,
                 Rm,Rm_end,Vleak,gnaSoma,gnaSr,gkdr,gkap,gkad,soma_caL,
                 soma_caR,soma_caN,soma_caT,soma_hbar,soma_km,soma_kdBG,
                 soma_kca,dhalf,steep,gh_end,gh_soma,
                 dlimit,dprox,dslope,okslope,CN):
        self.source_section = self.soma 
        self.source = self.source_section(0.5)._ref_v
        self.recording_time = False # for pyNN
        self.spike_times = h.Vector(0) # for pyNN
        
        
        ### setup self parameters to be passed in during optimization 
        self.soma_ra = soma_ra    
        self.global_ra = global_ra
        self.Ra_end = Ra_end
        self.Ra_dhalf = Ra_dhalf
        self.Ra_steep = Ra_steep
        self.cm = cm
        self.Rm = Rm
        self.Rm_end = Rm_end
        self.Vleak = Vleak
        self.gnaSoma = gnaSoma
        self.gnaSr = gnaSr
        self.gkdr = gkdr
        self.gkap = gkap
        self.gkad = gkad
        self.soma_caL = soma_caL
        self.soma_caR = soma_caR 
        self.soma_caN = soma_caN
        self.soma_caT = soma_caT
        self.soma_hbar = soma_hbar
        self.soma_km = soma_km
        self.soma_kdBG = soma_kdBG
        self.soma_kca = soma_kca 
        self.dhalf = dhalf
        self.steep = steep
        self.gh_end = gh_end
        self.gh_soma = gh_soma
        self.dlimit = dlimit
        self.dprox = dprox
        self.dslope = dslope
        self.okslope = okslope 
        self.cellnumber = CN
        
        self._setup_biophysics() # call biophysics  
        # self.excitatory = h.ExpSyn(0.5, sec=self.soma[0])
        # self.inhibitory = h.ExpSyn(0.5, sec=self.soma[0])
        # for syn in self.excitatory, self.inhibitory:
        #     syn.tau = 2.0
        # self.excitatory.e = 0.0
        # self.inhibitory.e = -70.0
        self.v_init = -70.0 # for pyNN
        # self.parameter_names = ('soma_ra') # for pyNN
        self.traces = {} # for pyNN
        self.recording_time = False # for pyNN
        self.rec = h.NetCon(self.source, None, sec=self.source_section) # for pyNN
        
    def record_v(self, active): # for pyNN
        if active:
            self.vtrace = h.Vector()
            self.vtrace.record(self.source_section(0.5)._ref_v)
        if not self.recording_time:
            self.record_times = h.Vector()
            self.record_times.record(h._ref_t)
            self.recording_time += 1
        else:
            self.vtrace = None
            self.recording_time -= 1
        if self.recording_time == 0:
            self.record_times = None
    
    def record(self, active): # for pyNN 
        if active:
            self.rec = h.NetCon(self.source, None, sec=self.source_section)            
            self.rec.record(self.spike_times)
        
    def memb_init(self, v_init=None): # for pyNN 
        if v_init:
            self.v_init = v_init
            assert self.v_init is not None, "cell is a %s" 
            self.__class__.__name__
        for seg in self.soma:
            seg.v = self.v_init
      
    def _setup_biophysics(self): # setup biophysics
##############################################################################      
##############################################################################   
     
        ### Let's start with intrinsic parameters 
        
##############################################################################      
##############################################################################  
            
        soma_ra = self.soma_ra
        # soma_ra = 100 # axial resistance for soma
        global_ra = self.global_ra
        # global_ra = 60 # axial resistance for dendrites
        # Ra_end = 30 # end Ra for dends increading increases frequency slope vice versa
        Ra_end = self.Ra_end
        # Ra_dhalf = 130 # 
        Ra_dhalf = self.Ra_dhalf
        # Ra_steep = 150 # 
        Ra_steep = self.Ra_steep
        # dend_ra = 43

        # cm = 1 # membrane capacitance
        cm = self.cm
        # Rm = 30000 # membrane resistance (used for passive conductance)
        Rm = self.Rm
        # Rm_end = 12000 # membrane resitance at end of dendrites
        Rm_end = self.Rm_end
        # Vleak = -65 # reversal potential for passive
        Vleak = self.Vleak
       
##############################################################################      
##############################################################################   
     
        ### Major voltage gated Na and K currents, all channels from Spruston
        
##############################################################################      
##############################################################################  
        # gnaSoma = 0.015 # NaV conductance in soma 
        gnaSoma = self.gnaSoma
        # gnaSr = gnaSoma*0.8 # NaV conductance in dendrites 
        gnaSr = self.gnaSr
        
        # gkdr = 0.3*gnaSoma # delayed rectifier K current expressed as perecentage of gNa
        gkdr = self.gkdr 
        # setgk = 0.15*gnaSoma

        # gkap = setgk # proximal A type current conductance
        gkap = self.gkap
        # gkad = setgk # distal A type current conductance 
        gkad = self.gkad 
        
        
##############################################################################      
##############################################################################   
     
        ### Accessory currents all channels from aged 
        
##############################################################################      
##############################################################################  
        # soma_caL = 0.000075 # L type calcium 
        soma_caL = self.soma_caL
        # soma_caR = 0.00015 # R type calcium 
        soma_caR = self.soma_caR
        # soma_caN = 0.00001 # N type calcium 
        soma_caN = self.soma_caN 
        # soma_caT = 0.00015 # T type calcium 
        soma_caT = self.soma_caT
        # soma_hbar = 0.0001 # funny current
        soma_hbar = self.soma_hbar
        # soma_km = 0.000005 # m type potassium current
        soma_km = self.soma_km
        # mykca_init = 0.0001
        # soma_kdBG = 0.0001 # somatic d type potassium current
        soma_kdBG = self.soma_kdBG
        # kdBG_init = 0.001  # altering d type currrent in dendrites
        # kdBG_init = self.kdBG_init
        # soma_kca = 0.00001 # calcium dependent k current 
        soma_kca = self.soma_kca 
        # kir_gkbar = 0.0000001 # inward rectifying potassium from DG paper
        
##############################################################################      
##############################################################################   
     
        ### Setup parametrs by which to vary currents in dendrites
        
##############################################################################      
##############################################################################  
        
        # make funny current in distal dendrites 9x somatic funny (Aged)
        # dhalf = 280 # half distance for funny current cutoff
        h_dhalf = self.dhalf
        # steep = 50 # slope for funny current
        h_steep = self.steep
        # gh_end = soma_hbar*9 # make distal funny current 9x somatic
        gh_end = self.gh_end
        # gh_soma = soma_hbar # just for conistency from aged paper 
        gh_soma = self.gh_soma
        
        # Switch from proximal to distal A type current and vary by distance 
        # dlimit=330	    	# cut-off for increase of A-type density 
        dlimit = self.dlimit
        # dprox=50         	# distance to switch from proximal to distal type 
        dprox = self.dprox 
        # dslope=0.065       # slope of A-type density 
        dslope = self.dslope 

        # okslope = setokslope # oblique potassium channel gradient 
        okslope = self.okslope 
        
        CN = self.cellnumber
##############################################################################      
##############################################################################   
     
        ### Insert channels into soma and set conductances
        
##############################################################################      
##############################################################################  
        self.soma.insert('pas') # standard neuron mechanism
        self.soma.Ra = soma_ra
        self.soma.cm = cm
        self.soma.insert('nax') # from spruston et al.
        self.soma.insert('kdr') # from spruston et al.
        self.soma.insert('kap') # from spruston et al.
        self.soma.insert('kad') # from spruston et al.
        # self.soma.insert('kir') # from DG paper
        self.soma.insert('cal') # from aged 
        self.soma.insert('can') # from aged 
        self.soma.insert('car') # from aged
        self.soma.insert('cat') # from aged
        self.soma.insert('cad') # from aged
        self.soma.insert('cadL') # from aged
        self.soma.insert('cadN') # from aged
        self.soma.insert('h') # from aged
        self.soma.insert('kdBG') # from aged
        self.soma.insert('km') # from aged
        self.soma.insert('kca') # from aged
        # self.soma.insert('mykca') # from aged

        for seg in self.soma:     
            seg.pas.g = 1/Rm
            seg.pas.e = Vleak
            seg.nax.gbar = gnaSoma
            # seg.nax.ar2 = 0.8
            seg.kdr.gkdrbar = gkdr
            seg.kap.gkabar = gkap
            seg.kad.gkabar = 0
            seg.cal.gcalbar = soma_caL
            seg.can.gcalbar = soma_caN
            seg.car.gcabar = soma_caR/10
            seg.cat.gcatbar = soma_caT/10
            # seg.kir.vhalfl = -90
            # seg.kir.gkbar = kir_gkbar
            seg.h.gbar = soma_hbar
            seg.h.vhalf = -82
            seg.km.gbar = soma_km
            seg.kdBG.gbar = soma_kdBG
            seg.kca.gbar = soma_kca
            # seg.mykca.gkbar = mykca_init
            seg.eca = 137
            # seg.ena = 55
            # seg.ek = -80

##############################################################################      
##############################################################################   
     
        ### Insert channels into basal dendrites
        
##############################################################################      
##############################################################################           
        DendLength = len(h.Cell[CN].dend)
        for x in range(0,DendLength):
            self.dend[x].insert('pas') # insert passive mechanism into dends        
            self.dend[x].insert('nax') # from Spruston et al. 
            self.dend[x].insert('kdr') # delayed rectifier potassium current
            self.dend[x].insert('kap') # proximal potassium
            self.dend[x].insert('kad') # distal potassium
            # self.dend[x].insert('kir') # inward rectifying potassium
            self.dend[x].insert('kdBG') # from aged 
            self.dend[x].insert('cal') # from aged 
            self.dend[x].insert('can') # from aged 
            self.dend[x].insert('car') # from aged
            self.dend[x].insert('cat') # from aged
            self.dend[x].insert('cad') # from aged
            self.dend[x].insert('cadL') # from aged
            self.dend[x].insert('cadN') # from aged
            self.dend[x].insert('h') # from aged
            self.dend[x].insert('kca') # from aged
            self.dend[x].insert('km') # from aged
            
            # find distance from soma of dend
            xdist = h.distance(self.soma(0.5), self.dend[x](0.5)) # neuron 2
            
            # Vary Ra sigmoidally along dendrites
            # Ra_end = 35
            dhalf = Ra_dhalf
            steep = Ra_steep
            self.dend[x].Ra = global_ra + (Ra_end - global_ra)/(1 + np.exp((dhalf-xdist)/steep)) # axial resistance
            # self.dend[x].Ra = dend_ra
            print(self.dend[x].Ra)
            
            # Vary Rm sigmoidally with distance
            # Rm_end =12000 # most distal Rm value
            dhalf = 200 # distance for inflection
            steep = 50 # slope
            Rm_sigmoid = Rm + (Rm_end - Rm)/(1+np.exp((dhalf-xdist)/steep))
            # print(Rm_sigmoid)
            seg.pas.g = 1/Rm_sigmoid
            # seg.pas.g = 1/Rm
            self.dend[x].cm = cm
            
            # set capicitance and passive conductance to make distal membrane
            # less excitable using Rm sigmoid and spinefactor for Cm
            # if xdist < spinelimit:
            #     self.dend[x].cm = cm
            #     for seg in self.dend[x]:
            #         seg.pas.g = 1/Rm_sigmoid
            # else: 
            #     self.dend[x].cm = spinefactor*cm
            #     for seg in self.dend[x]:
            #         seg.pas.g = spinefactor/Rm_sigmoid
            
            
            # use proximal potassium for membrane closer than dprox and distal for 
            # membrane segments further than dprox 
            if xdist > dprox:
                for seg in self.dend[x]:
                    seg.kad.gkabar = gkad*(1+xdist*dslope)
                    seg.kap.gkabar = 0
                    seg.kdBG.gbar = soma_kdBG
            else:
                for seg in self.dend[x]:
                    seg.kap.gkabar = gkap*(1+xdist*dslope)
                    seg.kad.gkabar = 0
                    # seg.kdBG.gbar = kdBG_init
        
            #### from aged ####                
#/* Inserting LVA Ca++ T-type channels along the apical trunk in
# a linearly increasing manner, for xdist > 100 um 
# */
            caT_distal_distance = 300 # distance for maximum conductance
            caT_proximal_distance = 100 # distance in dendrites for maxiumu cond.
            fr = xdist/caT_distal_distance
            if xdist < caT_proximal_distance:
                for seg in self.dend[x]:
                    seg.cat.gcatbar = soma_caT/10
            elif xdist < caT_distal_distance:
                for seg in self.dend[x]:
                    seg.cat.gcatbar = soma_caT*(1+2.4*fr)
            else:
                for seg in self.dend[x]:
                    seg.cat.gcatbar = soma_caT*3.4

                ### from aged ###
# /* Inserting HVAm Ca++ R-type and N-type, and HVA L-type channesls along
# the apical trunk. The R-type current is distributed as T-type current, in 
# a linearly increasing manner, for xdist < 300 um. 
# The L-type current is distributed in a 
# linearly decreasing conductance for distances xdist  < 100 um
# while the N-type current is distributed in a
# fixed conductance, as does L-type for distances xdist  > 100 um
# */

            caR_distal_distance = 300
            fr = xdist/caR_distal_distance
            if xdist < 50:
                for seg in self.dend[x]:
                    seg.car.gcabar = soma_caR/2
            elif xdist < caR_distal_distance:
                for seg in self.dend[x]:
                    seg.car.gcabar = soma_caR*(1+5*fr)
            else:
                for seg in self.dend[x]:
                    seg.car.gcabar = soma_caR*6


            caL_distal_distance = 100
            if xdist < caL_distal_distance:
                for seg in self.dend[x]:
                    seg.cal.gcalbar = soma_caL*(1-2.*fr/3.)
                    seg.can.gcalbar = soma_caN*(1-2.*fr/3.)
            else:
                for seg in self.dend[x]:
                    seg.cal.gcalbar = soma_caL/6.
                    seg.can.gcalbar = soma_caN/4.
   
            ### set conductance values for various channels in basal dendrites
            for seg in self.dend[x]:
                seg.pas.e = Vleak
                seg.kdr.gkdr = gkdr
                seg.nax.gbar = gnaSr
                seg.cal.gcalbar = soma_caL
                seg.can.gcalbar = soma_caN
                seg.car.gcabar = soma_caR/10
                seg.cat.gcatbar = soma_caT/10
                # seg.kir.vhalfl = -90
                # seg.kir.gkbar = kir_gkbar
                seg.eca = 137
                seg.kca.gbar = soma_kca
                seg.ena = 55
                seg.ek = -80
      
       
                              
#### from aged #### 
# /* To make the distal trunk h-current conductance, g_h, about 7
# times higher (at 300 um) than the somatic value vis-a-vis Magee
# J., J. of Neuroscience 18(19) 7613-7624, 1998, we vary I_h
# conductance sigmoidally along the apical trunk.
# */
                seg.h.gbar = gh_soma + (gh_end - gh_soma)/(1.0 + np.exp((dhalf-xdist)/steep))
        
##############################################################################      
##############################################################################   
     
        ### Insert channels into apical trunk
        
##############################################################################      
############################################################################## 
        self.Trunk = h.Cell[CN].Trunk
        for sec in self.Trunk:
            sec.insert('pas')
            sec.insert('nax')
            sec.insert('kdr')
            sec.insert('kap')
            sec.insert('kad')
            sec.insert('h')
            # sec.insert('kir') # inward rectifying potassium 
            sec.insert('cal') # from aged 
            sec.insert('can') # from aged 
            sec.insert('car') # from aged
            sec.insert('cat') # from aged
            sec.insert('cad') # from aged
            sec.insert('cadL') # from aged
            sec.insert('cadN') # from aged
            sec.insert('kca')
            sec.insert('km')
            sec.insert('kdBG')
            # seg.ek = -80
            # seg.ena = 55
            
            # find distance from soma 
            xdist = h.distance(self.soma(0.5), sec(0.5)) # neuron 2
            
            # set conductance value of parameters
            for seg in sec:
                seg.pas.e = Vleak
                seg.kdr.gkdr = gkdr
                seg.nax.gbar = gnaSr
                seg.h.gbar = gh_soma + (gh_end - gh_soma)/(1.0 + np.exp((dhalf-xdist)/steep))
                seg.eca = 137
                seg.kca.gbar = soma_kca
                # seg.kir.gkbar = kir_gkbar
                # seg.kir.vhalfl = -90
                # seg.kdBG.gbar = soma_kdBG
                seg.km.gbar = soma_km
            
            # vary Ra sigmoidally along trunk
            # Ra_end = 35
            dhalf = Ra_dhalf
            steep = Ra_steep
            sec.Ra = global_ra + (Ra_end - global_ra)/(1 + np.exp((dhalf-xdist)/steep)) # axial resistance
            # sec.Ra = dend_ra
            
            # vary Rm sigmoidally along trunk
            # Rm_end =12000
            dhalf = 200
            steep = 50
            Rm_sigmoid = Rm + (Rm_end - Rm)/(1+np.exp((dhalf-xdist)/steep))
            
            for seg in sec:   
                seg.pas.g = 1/Rm_sigmoid
                # seg.pas.g = 1/Rm
            sec.cm = cm
            # set capicitance and passive conductance to make distal membrane
            # less excitable using Rm sigmoid and spinefactor for Cm
            # if xdist < spinelimit:
            #     sec.cm = cm
            #     for seg in sec:
            #         seg.pas.g = 1/Rm_sigmoid
            # else: 
            #     sec.cm = spinefactor*cm
            #     for seg in sec:
            #         seg.pas.g = spinefactor/Rm_sigmoid
 
            # use proximal potassium for membrane closer than dprox and distal for 
            # membrane segments further than dprox 
            if xdist > dprox:
                for seg in sec:
                    seg.kad.gkabar = gkad*(1+xdist*dslope)
                    seg.kap.gkabar = 0
                    seg.kdBG.gbar = soma_kdBG
            else:
                for seg in sec:
                    seg.kap.gkabar = gkap*(1+xdist*dslope)
                    seg.kad.gkabar = 0
                    # seg.kdBG.gbar = kdBG_init
    
    
            #### from aged ####                
# /* Inserting LVA Ca++ T-type channels along the apical trunk in
# a linearly increasing manner, for xdist > 100 um 
# */
            caT_distal_distance = 300 # distance for maximum conductance
            caT_proximal_distance = 100 # distance in dendrites for maxiumu cond.
            fr = xdist/caT_distal_distance
            if xdist < caT_proximal_distance: 
                for seg in sec:
                    seg.cat.gcatbar = soma_caT/10
            elif xdist < caT_distal_distance:
                for seg in sec:
                    seg.cat.gcatbar = soma_caT*(1+2.4*fr)
            else:
                for seg in sec:
                    seg.cat.gcatbar = soma_caT*3.4
                
            
                            ### from aged ###
#/* Inserting HVAm Ca++ R-type and N-type, and HVA L-type channesls along
# the apical trunk. The R-type current is distributed as T-type current, in 
# a linearly increasing manner, for xdist < 300 um. 
# The L-type current is distributed in a 
# linearly decreasing conductance for distances xdist  < 100 um
# while the N-type current is distributed in a
# fixed conductance, as does L-type for distances xdist  > 100 um
# */
            caR_distal_distance = 300
            fr = xdist/caR_distal_distance
            if xdist < 50:      
                for seg in sec:
                    seg.car.gcabar = soma_caR/2
            elif xdist < caR_distal_distance:
                for seg in sec:
                    seg.car.gcabar = soma_caR*(1+5*fr)
            else:
                for seg in sec:
                    seg.car.gcabar = soma_caR*6

            caL_distal_distance = 100
            if xdist < caL_distal_distance:
                for seg in sec:
                    seg.cal.gcalbar = soma_caL*(1-2.*fr/3.)
                    seg.can.gcalbar = soma_caN*(1-2.*fr/3.)
            else:
                for seg in sec:
                    seg.cal.gcalbar = soma_caL/6.
                    seg.can.gcalbar = soma_caN/4.
            
##############################################################################      
##############################################################################   
     
        ### Insert channels into apical obliques (non-trunk)
        
##############################################################################      
##############################################################################  
        self.NonTrunk = h.Cell[CN].NonTrunk
        for sec in self.NonTrunk:
            sec.insert('pas')
            sec.insert('nax')
            sec.insert('kdr')
            sec.insert('kap')
            sec.insert('kad')
            sec.insert('h')
            # sec.insert('kir') # inward rectifying potassium 
            sec.insert('cal') # from aged 
            sec.insert('can') # from aged 
            sec.insert('car') # from aged
            sec.insert('cat') # from aged
            sec.insert('cad') # from aged
            sec.insert('cadL') # from aged
            sec.insert('cadN') # from aged
            sec.insert('kca')
            sec.insert('km')
            
            dhalf = h_dhalf
            steep = h_steep
            
            # sett conductances for various channels 
            for seg in sec:
                seg.pas.e = Vleak
                seg.kdr.gkdr = gkdr
                seg.nax.gbar = gnaSr
                seg.h.gbar = gh_soma + (gh_end - gh_soma)/(1.0 + np.exp((dhalf-xdist)/steep))
                seg.kca.gbar = soma_kca
                # seg.kdBG.gbar = soma_kdBG
                seg.km.gbar = soma_km
                          
            
            # find distance from soma 
            xdist = h.distance(self.soma(0.5), sec(0.5)) 
            
            # Vary Ra sigmoidally with disttance 
            # Ra_end = 35
            dhalf = Ra_dhalf
            steep = Ra_steep
            sec.Ra = global_ra + (Ra_end - global_ra)/(1 + np.exp((dhalf-xdist)/steep)) # axial resistance
            # sec.Ra = dend_ra
            
            # Vary Rm sigmoidally with distance 
            # Rm_end =12000
            dhalf = 200
            steep = 50
            Rm_sigmoid = Rm + (Rm_end - Rm)/(1+np.exp((dhalf-xdist)/steep))
            for seg in sec:   
                seg.pas.g = 1/Rm_sigmoid
                # seg.pas.g = 1/Rm
            sec.cm = cm
            # set capicitance and passive conductance to make distal membrane
            # less excitable using Rm sigmoid and spinefactor for Cm
            # if xdist < spinelimit:
            #     sec.cm = cm
            #     for seg in sec:
            #         seg.pas.g = 1/Rm_sigmoid
            # else: 
            #     sec.cm = spinefactor*cm
            #     for seg in sec:
            #         seg.pas.g = spinefactor/Rm_sigmoid
 
                    ### from spruston###
#  forsec obliqueList {
# 				for (x) {
# 					odist = distance(x) 	// odist is the distance of each segment along the oblique branch
# 					pdist_k = dlimit
# 					gkabar_kap(x) = gkap*(1+pdist_k*dslope+odist*okslope)
# 					gkabar_kad(x) = gkad*(1+pdist_k*dslope+odist*okslope)
# 				}
# 		}
            # use proximal potassium for membrane closer than dprox and distal for 
            # membrane segments further than dprox 
            odist = xdist ## need to actually find this???
            pdist_k = dlimit
            for seg in sec:
                seg.kad.gkabar = gkad*(1+pdist_k*dslope+odist*okslope)
                seg.kap.gkabar = gkad*(1+pdist_k*dslope+odist*okslope)         
            
                        #### from aged ####                
# /* Inserting LVA Ca++ T-type channels along the apical trunk in
# a linearly increasing manner, for xdist > 100 um 
# */
            caT_distal_distance = 300 # distance for maximum conductance
            caT_proximal_distance = 100 # distance in dendrites for maxiumu cond.
            fr = xdist/caT_distal_distance
            if xdist < caT_proximal_distance: 
                for seg in sec:
                    seg.cat.gcatbar = soma_caT/10
            elif xdist < caT_distal_distance:
                for seg in sec:
                    seg.cat.gcatbar = soma_caT*(1+2.4*fr)
            else:
                for seg in sec:
                    seg.cat.gcatbar = soma_caT*3.4
                
            
                            ### from aged ###
# /* Inserting HVAm Ca++ R-type and N-type, and HVA L-type channesls along
# the apical trunk. The R-type current is distributed as T-type current, in 
# a linearly increasing manner, for xdist < 300 um. 
# The L-type current is distributed in a 
# linearly decreasing conductance for distances xdist  < 100 um
# while the N-type current is distributed in a
# fixed conductance, as does L-type for distances xdist  > 100 um
# */
            caR_distal_distance = 300
            fr = xdist/caR_distal_distance
            if xdist < 50:      
                for seg in sec:
                    seg.car.gcabar = soma_caR/2
            elif xdist < caR_distal_distance:
                for seg in sec:
                    seg.car.gcabar = soma_caR*(1+5*fr)
            else:
                for seg in sec:
                    seg.car.gcabar = soma_caR*6


            caL_distal_distance = 100
            if xdist < caL_distal_distance:
                for seg in sec:
                    seg.cal.gcalbar = soma_caL*(1-2.*fr/3.)
                    seg.can.gcalbar = soma_caN*(1-2.*fr/3.)
            else:
                for seg in sec:
                    seg.cal.gcalbar = soma_caL/6.
                    seg.can.gcalbar = soma_caN/4.