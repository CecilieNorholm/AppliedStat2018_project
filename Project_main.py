### APPSTAT PROJEKT 


import numpy as np                                     # Matlab like syntax for linear algebra and functions
import matplotlib.pyplot as plt                        # Plots and figures like you know them from Matlab
from iminuit import Minuit                             # The actual fitting tool, better than scipy's
from probfit import BinnedLH, Chi2Regression, Extended, UnbinnedLH # Helper tool for fitting
import sys
from scipy import stats
from scipy.special import erfc

# Read in data from files and put them into arrays


blinded = False
if blinded:
    blinding = r.normal(0, 0.1)      # I add a constant (Gaussian with +-10cm) to remain "blind"
else:
    blinding = 0

DBall=np.array([])
eDBall=np.array([])
DRail=np.array([])
eDRail=np.array([])
Gates=np.array([])
eGates=np.array([])
AngADoor=np.array([])
eAngADoor=np.array([])
AngAWall=np.array([])
eAngAWall=np.array([])
AngBDoor=np.array([])
eAngBDoor=np.array([])
AngBWall=np.array([])
eAngBWall=np.array([])
Hook=np.array([])
eHook=np.array([])
Pendul=np.array([])
ePendul=np.array([])
StringL=np.array([])
eStringL=np.array([])
StringR=np.array([])
eStringR=np.array([])

infiles = ["data/data_DBall.txt"]
for infile in infiles:
    tmp_DBall, tmp_eDBall = np.loadtxt(infile, skiprows=2, unpack=True)
    DBall=np.append(DBall, tmp_DBall)
    eDBall=np.append(eDBall, tmp_eDBall)
    
infiles = ["data/data_DRail.txt"] 
for infile in infiles:
    tmp_DRail, tmp_eDRail = np.loadtxt(infile, skiprows=2, unpack=True)
    DRail=np.append(DRail, tmp_DRail)
    eDRail=np.append(eDRail, tmp_eDRail)
   
infiles = ["data/data_Gates.txt"] 
for infile in infiles:
    tmp_Gates, tmp_eGates = np.loadtxt(infile, skiprows=2, unpack=True)
    Gates=np.append(Gates, tmp_Gates)
    eGates=np.append(eGates, tmp_eGates)
   
infiles = ["data/data_AngADoor.txt"] 
for infile in infiles:
    tmp_AngADoor, tmp_eAngADoor = np.loadtxt(infile, skiprows=2, unpack=True)
    AngADoor=np.append(AngADoor, tmp_AngADoor)
    eAngADoor=np.append(eAngADoor, tmp_eAngADoor)
    
infiles = ["data/data_AngAWall.txt"] 
for infile in infiles:
    tmp_AngAWall, tmp_eAngAWall = np.loadtxt(infile, skiprows=2, unpack=True)
    AngAWall=np.append(AngAWall, tmp_AngAWall)
    eAngAWall=np.append(eAngAWall, tmp_eAngAWall)    

infiles = ["data/data_AngBDoor.txt"] 
for infile in infiles:
    tmp_AngBDoor,  tmp_eAngBDoor = np.loadtxt(infile, skiprows=2, unpack=True)
    AngBDoor =np.append(AngBDoor, tmp_AngBDoor)
    eAngBDoor=np.append(eAngBDoor, tmp_eAngBDoor) 

infiles = ["data/data_AngBWall.txt"] 
for infile in infiles:
    tmp_AngBWall,  tmp_eAngBWall = np.loadtxt(infile, skiprows=2, unpack=True)
    AngBWall =np.append(AngBWall, tmp_AngBWall)
    eAngBWall=np.append(eAngBWall, tmp_eAngBWall) 
    
infiles = ["data/data_Hook.txt"] 
for infile in infiles:
    tmp_Hook,  tmp_eHook = np.loadtxt(infile, skiprows=2, unpack=True)
    Hook =np.append(Hook, tmp_Hook)
    eHook=np.append(eHook, tmp_eHook) 
    
infiles = ["data/data_Pendul.txt"] 
for infile in infiles:
    tmp_Pendul,  tmp_ePendul = np.loadtxt(infile, skiprows=2, unpack=True)
    Pendul=np.append(Pendul, tmp_Pendul)
    ePendul=np.append(ePendul, tmp_ePendul)         

infiles = ["data/data_StringL.txt"] 
for infile in infiles:
    tmp_StringL,  tmp_eStringL = np.loadtxt(infile, skiprows=2, unpack=True)
    StringL =np.append(StringL, tmp_StringL)
    eStringL=np.append(eStringL, tmp_eStringL) 

infiles = ["data/data_StringR.txt"] 
for infile in infiles:
    tmp_StringR,  tmp_eStringR = np.loadtxt(infile, skiprows=2, unpack=True)
    StringR =np.append(StringR, tmp_StringR)
    eStringR=np.append(eStringR, tmp_eStringR)     

  
    