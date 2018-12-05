### APPSTAT PROJEKT 


import numpy as np                                     # Matlab like syntax for linear algebra and functions
import matplotlib.pyplot as plt                        # Plots and figures like you know them from Matlab
from iminuit import Minuit                             # The actual fitting tool, better than scipy's
from probfit import BinnedLH, Chi2Regression, Extended, UnbinnedLH # Helper tool for fitting
import sys
from scipy import stats
from scipy.special import erfc

# Read in data from files and put them into arrays

infiles = ["data/data_DBall.txt",
           "data/data_DRail.txt",
           "data/data_Gates.txt",
           "data/data_AngADoor.txt",
           "data/data_AngAWall.txt",
           "data/data_AngBDoor.txt",
           "data/data_AngBWall.txt",
           "data/data_AngBWall.txt", 
                   
           
           
           "data/data_DBalltxt"]


