# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 16:41:48 2018

@author: Cecilie
"""

import numpy as np
import math as m

def g_Incline(a,d_rail,D_ball,theta):
    """
    Computes the error on gravitational acceleration g for a
    'Ball Incline Experiment'

    INPUT parameters
    -----------
    a : acceleration float
    d_rail : diameter of rail as float
    D_ball : diameter of ball as float
    theta : angle in radians as float

    OUTPUT returns
    ---------
    g_incline : gravitational acceleration 'g' as a float
    """

    theta = np.radians(theta)


    #Rewriting parameters to ease computation
    delta = D_ball**2 / (D_ball**2 - d_rail**2)

    #Calculating the gravitational acceleration, g
    g_Incline = ( a / np.sin(theta) ) * ( 1 + (2/5)*delta )

    return g_Incline
