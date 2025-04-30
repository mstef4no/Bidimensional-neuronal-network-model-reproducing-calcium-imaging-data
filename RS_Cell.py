#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 16:55:45 2024

@author: mattiastefano
"""

# Regular Spiking pyramidal neuron Ref. Fig.2a Popischil 2008 
from neuron import h
from matplotlib import pyplot

class RS_CELL:
    def __init__ (self): #, label):   
        #self._gid=gid  , gid
        self.soma = soma = h.Section(name='soma',cell=self)
        
        self.initsoma()
    
        #self.label = label 
    def initsoma (self): 
        soma = self.soma 
        soma.nseg = 1
        soma.diam = 61.4
        soma.L = 61.4 
        soma.cm = 1
        soma.Ra = 100 
        """        
        soma.insert ('hh2')
        soma.ek_hh2 = -90 # K+ current reversal potential (mV)
        soma.ena_hh2 = 50 # Na+ current reversal potential (mV)
        soma.gnabar_hh2 = 0.05  # Sodium conductance in S/cm2
        soma.gkbar_hh2 = 0.006  # Potassium conductance in S/cm2
        soma.vtraub_hh2 = -56.2
        """ 
        soma.insert ('hhers')
        soma.ek_hhers = -90 #-100 # K+ current reversal potential (mV)
        soma.ena_hhers = 50 # Na+ current reversal potential (mV)
        soma.gnabar_hhers = 0.056 #0.05  # Sodium conductance in S/cm2
        soma.gkdrbar_hhers =  0.006 #0.005  # Potassium conductance in S/cm2
        soma.vt_hhers = -56.2 #-55
        soma.taumax_hhers = 608
        soma.gmbar_hhers = 7.5*1e-5 



        soma.insert('pas')
        soma.e_pas = -70.3 
        soma.g_pas = 2.05*1e-5 


        """
        soma.insert('im') 
        soma.gkbar_im = 7.5*1e-5 
        """

        #def stampa_label(self): 
        #    return f" Label: {self.label}\n soma.L: {self.soma.L} \n soma.nseg: {self.soma.nseg}"  



