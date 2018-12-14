#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 10:59:44 2018

@author: chris
"""

import numpy as np                                     # Matlab like syntax for linear algebra and functions
import matplotlib.pyplot as plt  
from matplotlib import ticker                      # Plots and figures like you know them from Matlab
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

def get_bincenter_and_counts_in_range(hist, xmin=None, xmax=None):
    
    if xmin is None:
        xmin = np.min(hist)
    if xmax is None:
        xmax = np.max(hist)
    
    counts, bin_edges, _ = hist
    bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
    mask1 = (xmin < bin_centers) & (bin_centers <= xmax) 
    mask2 = counts > 0
    mask_final = mask1 & mask2
    return bin_centers[mask_final], counts[mask_final], np.sqrt(counts[mask_final])

def gauss_pdf(x, mu, sigma):
    """Normalized Gaussian"""
    return 1 / np.sqrt(2 * np.pi) / sigma * np.exp(-(x - mu) ** 2 / 2. / sigma ** 2)

#def PendPlots(timer_dat):
    """
    Plotting and Fitting function for pendulum data, 
    
    timer_dat: Time series of pendulum period measurements
    
    """
    
    # Set data and plot
timer_R1=[]
tmp_timer_R1=[]
infiles = ["data/Timer_Dat/timer_R1.dat"]
for infile in infiles:
        n, tmp_timer_R1=np.loadtxt(infile, skiprows=0, unpack=True)
        timer_R1=np.append(timer_R1, tmp_timer_R1)        
  
timer_dat=timer_R1
y=timer_dat
y=[x-y[0] for x in y]
y=np.array(y[1:len(y)])
n=np.linspace(1,len(y),len(y))
x=n

 # Perform liner fit of measurements
chi2_object= Chi2Regression(linear,x, y)
minuit = Minuit(chi2_object, pedantic=False, a=0,b=0) #   
minuit.migrad()  # perform the actual fit
    
# Compute residuals
res=[]    
for i in range(len(y)):
    res.append(y[i]-y[i-1])
res[0]=y[0]
res=[x-minuit.values["a"] for x in res]


#Fit residuals to gauss
minL = min(res)-0.1
maxL = max(res)+0.1
bins= 10


fig_gauss, ax_gauss = plt.subplots(figsize=(14,8))

hist_gauss = ax_gauss.hist(res, bins=bins, range=(minL, maxL), histtype='step', label='Binned Residuals')
ax_gauss.set(xlabel='Residuals [s]', ylabel='Frequency', title='Gaussian fit of residuals')

res_x, res_y, res_sy = get_bincenter_and_counts_in_range(hist_gauss, minL, maxL)

res=np.array(res)
chi2_res = Chi2Regression(gauss_pdf, res_x, res_y, res_sy) 
minuit_res = Minuit(chi2_res, pedantic=False, mu=res.mean(), sigma=res.std(ddof=1)) 
minuit_res.migrad();
chi2_res = minuit_res.fval

xaxisres = np.linspace(minL, maxL, 100) 
yaxisres = gauss_pdf(xaxisres, *minuit_res.args) 
ax_gauss.plot(xaxisres, yaxisres, '-', label='Fit')


sigma=minuit_res.values["sigma"]


#Fit measurements with error from residuals

sigma=np.ones(len(x))*sigma
chi2_object_T = Chi2Regression(linear, x, y, error=sigma)
minuit_T = Minuit(chi2_object_T, pedantic=False, a=3.31,b=0.01) #   
minuit_T.migrad()  # perform the actual fit
chi2_T = minuit_T.fval

N_var = 2                     # Number of variables (p0:p3)
N_dof = len(x) - N_var   # Number of degrees of freedom
chi2_prob = stats.chi2.sf(chi2_T, N_dof) # The chi2 probability given N_DOF degrees of freedom

eT=minuit_T.values["b"]
 
N_var = 2                     # Number of variables (p0:p3)
N_dof = len(x) - N_var   # Number of degrees of freedom
chi2_prob = stats.chi2.sf(chi2_T, N_dof) # The chi2 probability given N_DOF degrees of freedom

fig_main, ax_main = plt.subplots(figsize=(10,7))
ax_main.plot(x,y, "o")
ax_main.set(xlabel="Measurement number",ylabel="Time elapsed (s)")

xaxis = n
yaxis = linear(xaxis, *minuit.args)
ax_main.plot(xaxis,yaxis)

d = {'Chi2':     chi2_T,
     'ndf':      N_dof,
     'Prob':     chi2_prob,
     'Haeldning': [minuit.values['a'], minuit.errors['a']],
     'Skaering': [minuit.values['b'], minuit.errors['b']]
    }

string = nice_string_output(d, extra_spacing=2, decimals=3)
add_text_to_ax(0.02, 0.97, string, ax_main, fontsize=14)
    
res=res.tolist()

# Plot residuals around 0
fig_res, ax_res = plt.subplots(figsize=(8, 6))
ax_res.errorbar(x, res, yerr=eT, xerr=None, fmt='k_', ecolor='k', elinewidth=1, capsize=2, capthick=1)
    
#    xaxis = np.linspace(0, 25, 25)
#    yaxis = np.linspace (-1,1,20)
#    ax_res.plot(xaxis,yaxis)
minorloc = ticker.AutoMinorLocator()
fig_test, ax_test = plt.subplots(figsize=(12, 10))
ax_test.plot(x, y, 'ko', markersize=7)
ax_test.plot(xaxis, yaxis, 'r-', linewidth=1.5)
ax_test.text(2,77,'Fit Results', fontsize=18, color='r')#, horizontalalignment='center')
ax_test.yaxis.set_minor_locator(minorloc)
#ax_test.xaxis.set_minor_locator(minorloc)
ax_test.tick_params(which='major', direction='in', length=7, width=2, labelsize='large')
ax_test.tick_params(which='minor', direction='in', length=4, width=2, labelsize='medium')
ax_test.set_xlabel('Number of measurements', fontsize=14)
ax_test.set_ylabel('Elapsed time [s]', fontsize=14)
ax_test.set_ylim(bottom=-32, top=85)

L, B, W, H = [0.125, 0.14, 0.775, 0.2]
ax_test2 = fig_test.add_axes([L, B, W, H])
ax_test2.plot(x, res, 'o', markersize=6)
ax_test2.errorbar(x, res, yerr=eT, xerr=None, fmt='o', ecolor='k', color='xkcd:royal blue')
ax_test2.plot((1,24),(0,0), 'k-')
ax_test2.plot((1,24), (eT,eT), 'b--')
ax_test2.plot((1,24), (-eT,-eT), 'b--')
ax_test2.axes.get_xaxis().set_visible(False)
ax_test2.yaxis.tick_right()
ax_test2.yaxis.set_label_position("right")
ax_test2.yaxis.set_minor_locator(minorloc)
ax_test2.tick_params(which='major', direction='in', length=7, width=2, labelsize='large')
ax_test2.tick_params(which='minor', direction='in', length=4, width=2, labelsize='medium')
ax_test2.set_ylabel('Fit residuals [s]', fontsize=14)
 
   # return 
    