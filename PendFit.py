#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 20:52:55 2018

@author: chris
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 10:59:44 2018

@author: chris
"""

import numpy as np                                     # Matlab like syntax for linear algebra and functions
import matplotlib.pyplot as plt                        # Plots and figures like you know them from Matlab
from iminuit import Minuit                             # The actual fitting tool, better than scipy's
from probfit import BinnedLH, Chi2Regression, Extended, UnbinnedLH # Helper tool for fitting
import sys
import seaborn as sns
from scipy import stats
from scipy.special import erfc

# Read in data from files and put them into arrays
sys.path.append('/External_Functions')
from ExternalFunctions import nice_string_output, add_text_to_ax # useful functions to print fit results on figure
def linear(x,a,b):
    y=a*x+b
    return y
def PendFit(timer_dat):
    """
    Fitting function for pendulum data, 
    
    timer_dat: Time series of pendulum period measurements
    
    """
    
    # Set data and plot
    y=timer_dat
    n=np.linspace(1,len(y),len(y))
    x=n
    
    # Perform and plot Fit
    chi2_object = Chi2Regression(linear,x, y)
    minuit = Minuit(chi2_object, pedantic=False, a=0,b=0) #   
    minuit.migrad()  # perform the actual fit
    
    chi2 = minuit.fval
    T=minuit.values['a']
    sigmaT=minuit.errors['a']
    
    return T, sigmaT, chi2