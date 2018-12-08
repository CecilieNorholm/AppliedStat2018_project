# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 16:39:08 2018

@author: Cecilie
"""

import numpy as np
import math as m

# Defining error propagation function for Pendulum
def gcalc_pendulum(L, T):
    
        """
    Computes the gravitational acceleration g for a
    'Pendulum Experiment'

    INPUT parameters
    -----------
    L : Length
    T : time array of shape [n_samples, 1]
    
    OUTPUT returns
    ---------
    g_pendulum : gravitational acceleration 'g' as a numpy array
    """
    
    g_pendulum = L * (2*np.pi / T)**2

    return: g_pendulum