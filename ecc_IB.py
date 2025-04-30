#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:16:50 2024

@author: mattiastefano
"""

from neuron import h
from matplotlib import pyplot 

class PYRS: 
    def __init__ (self): #, label):   
        #self._gid=gid  , gid
        self.soma = soma = h.Section(name='soma',cell=self)

        self.initsoma()

        #self.label = label 
    def initsoma (self):
        soma = self.soma 
        soma.nseg = 1
        soma.diam = 96
        soma.L = 96
        soma.cm = 1
        soma.Ra = 100
    
        soma.insert ('hh2')
        soma.ek_hh2 = -100 # K+ current reversal potential (mV)
        soma.ena_hh2 = 50 # Na+ current reversal potential (mV)
        soma.gnabar_hh2 = 0.05  # Sodium conductance in S/cm2
        soma.gkbar_hh2 = 0.005  # Potassium conductance in S/cm2
        soma.vtraub_hh2 = -55
    
        soma.insert ('pas')
        soma.e_pas = -85 #-70
        soma.g_pas = 1e-5
    
        soma.insert ('im')		# M current 
        soma.gkbar_im = 3e-5	
        soma.tauM_im = 608
        
        soma.insert ('ical')		# IL current 
        soma.gcabar_ical= 2e-4 #1e-4	## N.B. nel paper Fig.5 è 2e-4 	    
