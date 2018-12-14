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

# Define your PDF / model 
def gauss_pdf(x, mu, sigma):
    """Normalized Gaussian"""
    return 1 / np.sqrt(2 * np.pi) / sigma * np.exp(-(x - mu) ** 2 / 2. / sigma ** 2)


#def PendFit(timer_dat):
    """
Fitting function for pendulum data, 

timer_dat: Time series of pendulum period measurements

"""
timer_C2=np.array([])
#
timer_dat=[]
tmp_timer_dat=[]
infiles = ["data/Timer_Dat/timer_R1.dat"]
for infile in infiles:
        n, tmp_timer_dat=np.loadtxt(infile, skiprows=0, unpack=True)
        timer_C1=np.append(timer_dat, tmp_timer_dat)
timer_dat=timer_C1

# Set data and plot
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


fig_gauss, ax_gauss = plt.subplots(figsize=(10,8))

hist_gauss = ax_gauss.hist(res, bins=bins, range=(minL, maxL), histtype='step', label='Binned Residuals', linewidth=5, color='xkcd:royal blue')
#ax_gauss.set(xlabel='Residuals [s]', ylabel='Frequency')#, title='Gaussian fit of residuals')
ax_gauss.set_ylabel('Frequency', fontsize=20, fontweight='bold')
ax_gauss.set_xlabel('Residuals [s]', fontsize=20, fontweight='bold')
ax_gauss.tick_params(direction='in', labelsize=18, length=7, width=2)

res_x, res_y, res_sy = get_bincenter_and_counts_in_range(hist_gauss, minL, maxL)

res=np.array([res])
chi2_res = Chi2Regression(gauss_pdf, res_x, res_y, res_sy) 
minuit_res = Minuit(chi2_res, pedantic=False, mu=res.mean(), sigma=res.std(ddof=1)) 
minuit_res.migrad();
chi2_res = minuit_res.fval

xaxisres = np.linspace(minL, maxL, 100) 
yaxisres = gauss_pdf(xaxisres, *minuit_res.args) 
ax_gauss.plot(xaxisres, yaxisres, '-', label='Fit', color='r', linewidth=5)
plt.savefig('Residuals_Gaussfit_Pendulum.png')

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

T=minuit_T.values["a"]
eT=minuit_T.values["b"]
    
    #return T, eT, chi2_T, chi2_prob