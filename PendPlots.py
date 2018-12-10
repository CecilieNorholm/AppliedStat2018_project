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
        


def PendPlots(timer_dat):
    """
    Plotting and Fitting function for pendulum data, 
    
    timer_dat: Time series of pendulum period measurements
    
    """
    
    # Set data and plot
    fig_fit, ax_fit = plt.subplots(figsize=(12, 8))
    y=timer_dat
    n=np.linspace(1,len(y),len(y))
    x=n
    ax_fit.plot(x,y, 'o')
    ax_fit.set(xlabel="Measurement number",ylabel="Time elapsed (s)")
    
    # Perform and plot Fit
    chi2_object = Chi2Regression(linear,x, y)
    minuit = Minuit(chi2_object, pedantic=False, a=0,b=0) #   
    minuit.migrad()  # perform the actual fit
    chi2 = minuit.fval
    
    N_var = 2                     # Number of variables (p0:p3)
    N_dof = len(x) - N_var   # Number of degrees of freedom
    chi2_prob = stats.chi2.sf(chi2, N_dof) # The chi2 probability given N_DOF degrees of freedom
    
    xaxis = n
    yaxis = linear(xaxis, *minuit.args)
    ax_fit.plot(xaxis,yaxis)
    ax_fit.set(xlabel="Measurement number",ylabel="Time elapsed (s)")
    
    d = {'Chi2':     chi2,
         'ndf':      N_dof,
         'Prob':     chi2_prob,
         'Haeldning': [minuit.values['a'], minuit.errors['a']],
         'Skaering': [minuit.values['b'], minuit.errors['b']]
        }
    
    string = nice_string_output(d, extra_spacing=2, decimals=3)
    add_text_to_ax(0.02, 0.97, string, ax_fit, fontsize=14)
        
    # Residuals    
    y=y-minuit.values["b"]
    res=[]    
    for i in range(len(y)):
        res.append(y[i]-y[i-1])
    res[0]=y[0]
    res=[x-minuit.values["a"] for x in res]
    
    res2=[x**2 for x in res]
    eT=np.sqrt(np.sum(res2)/len(x))

     
    # Plot residuals around 0
    fig_res, ax_res = plt.subplots(figsize=(8, 6))
    ax_res.errorbar(x, res, eT, fmt='k_', ecolor='k', elinewidth=1, capsize=2, capthick=1)
        
#    xaxis = np.linspace(0, 25, 25)
#    yaxis = np.linspace (-1,1,20)
#    ax_res.plot(xaxis,yaxis)
    
    T=minuit.values['a']
    sigmaT=minuit.errors['a']
    
    return T,sigmaT, res
    