#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 14:29:35 2018

@author: chris
"""
import numpy as np
import math as m

# Defining error propagation function for Pendulum
def errorprop_pendulum(L, Lerr, T, Terr):
    
        """
    Computes the error on gravitational acceleration g for a
    'Pendulum Experiment'

    INPUT parameters
    -----------
    L : Length
    Lerr= Uncertainty on L
    T : time array of shape [n_samples, 1]
    Terr: Uncertainty on time mmeasurements
    
    OUTPUT returns
    ---------
    sigma_g : error on 'g' as a numpy array
    """
    
    sigma_g= np.sqrt( ((2*np.pi / T)**4) * Lerr**2 + ((-2*L * ((2*np.pi)**2) / T**3)**2) * Terr**2 )

    return: sigma_g
    
