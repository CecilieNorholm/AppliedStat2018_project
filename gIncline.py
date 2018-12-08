# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 16:41:48 2018

@author: Cecilie
"""

import numpy as np
import math as m

def g_Incline(a,t,d,d_rail,D_ball,theta,Delta_theta,sigma_t,sigma_d,sigma_rail,sigma_ball,sigma_theta):
    """
    Computes the error on gravitational acceleration g for a
    'Ball Incline Experiment'

    INPUT parameters
    -----------
    a : acceleration array of shape [n_samples, 1]
    t : time array of shape [n_samples, 1]
    d : distance between gates array of shape [n_samples, 1]
    d_rail : diameter of rail as float
    D_ball : diameter of ball as float
    theta : angle in radians as float
    Delta_theta : difference angle in radians as float

    OUTPUT returns
    ---------
    g_incline : gravitational acceleration 'g' as a numpy array
    """

    # make sure that we have numpy arrays; also
    # reshape the array X to ensure that we have
    # a multidimensional numpy array (ndarray) np.array(a).reshape((a.shape[0], -1))
    a = np.array(a).reshape((len(a),1))
    d = np.array(d).reshape((len(d),1))
    t = np.array(t).reshape((len(t),1))

    #Rewriting parameters to ease computation
    delta = D_ball**2 / (D_ball**2 - d_rail**2)
    phi = Delta_theta

    #Calculating the gravitational acceleration, g
    g_Incline = ( a / np.sin(theta + phi) ) * ( 1 + (2/5)*delta ) 
 
    return g_Incline