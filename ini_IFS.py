#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:02:18 2024

@author: mattiastefano
"""

from neuron import h
from matplotlib import pyplot


class INFS:
    def __init__ (self): #, label):   
        #self._gid=gid  , gid
        self.soma = soma = h.Section(name='soma',cell=self)

        self.initsoma()
    
    #self.label = label 
    def initsoma (self):
        soma = self.soma 
        soma.nseg = 1
        soma.diam = 56.9 
        soma.L =56.9
        soma.cm = 1
        soma.Ra = 100 
        
        soma.insert ('hhers') 
        soma.ek_hhers = -90
        soma.gkdrbar_hhers = 0.0039 # 5*1e-3 
        soma.ena_hhers = 50 
        soma.gnabar_hhers = 0.058 #50*1e-3 
        soma.vt_hhers = -57.9
        soma.gmbar_hhers = 7.87 * 1e-5 
        soma.taumax_hhers = 502 
        
        soma.insert ('pas') 
        soma.e_pas = -70.4 
        soma.g_pas = 3.8 * 1e-5 