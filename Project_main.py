### APPSTAT PROJEKT 
#Test chr 222

import numpy as np                                     # Matlab like syntax for linear algebra and functions
import matplotlib.pyplot as plt                        # Plots and figures like you know them from Matlab
from iminuit import Minuit                             # The actual fitting tool, better than scipy's
from probfit import BinnedLH, Chi2Regression, Extended, UnbinnedLH # Helper tool for fitting
import sys
import seaborn as sns
from scipy import stats
from scipy.special import erfc

from PendFit import PendFit
from PendPlots import PendPlots

# Read in data from files and put them into arrays
sys.path.append('/External_Functions')
from ExternalFunctions import nice_string_output, add_text_to_ax # useful functions to print fit results on figure

# Blinding
blinded = False
if blinded:
    blinding = r.normal(0, 0.1)      # I add a constant (Gaussian with +-10cm) to remain "blind"
else:
    blinding = 0
   
##### Lengths and angles read in 
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

infiles = ["data/Lengths/data_DBall.txt"]
for infile in infiles:
    tmp_DBall, tmp_eDBall = np.loadtxt(infile, skiprows=2, unpack=True)
    DBall=np.append(DBall, tmp_DBall)
    eDBall=np.append(eDBall, tmp_eDBall)
infiles = ["data/Lengths/data_DRail.txt"] 
for infile in infiles:
    tmp_DRail, tmp_eDRail = np.loadtxt(infile, skiprows=2, unpack=True)
    DRail=np.append(DRail, tmp_DRail)
    eDRail=np.append(eDRail, tmp_eDRail)
infiles = ["data/Lengths/data_Gates.txt"] 
for infile in infiles:
    tmp_Gates, tmp_eGates = np.loadtxt(infile, skiprows=2, unpack=True)
    Gates=np.append(Gates, tmp_Gates)
    eGates=np.append(eGates, tmp_eGates)
infiles = ["data/Lengths/data_AngADoor.txt"] 
for infile in infiles:
    tmp_AngADoor, tmp_eAngADoor = np.loadtxt(infile, skiprows=2, unpack=True)
    AngADoor=np.append(AngADoor, tmp_AngADoor)
    eAngADoor=np.append(eAngADoor, tmp_eAngADoor)
infiles = ["data/Lengths/data_AngAWall.txt"] 
for infile in infiles:
    tmp_AngAWall, tmp_eAngAWall = np.loadtxt(infile, skiprows=2, unpack=True)
    AngAWall=np.append(AngAWall, tmp_AngAWall)
    eAngAWall=np.append(eAngAWall, tmp_eAngAWall)    
infiles = ["data/Lengths/data_AngBDoor.txt"] 
for infile in infiles:
    tmp_AngBDoor,  tmp_eAngBDoor = np.loadtxt(infile, skiprows=2, unpack=True)
    AngBDoor =np.append(AngBDoor, tmp_AngBDoor)
    eAngBDoor=np.append(eAngBDoor, tmp_eAngBDoor) 
infiles = ["data/Lengths/data_AngBWall.txt"] 
for infile in infiles:
    tmp_AngBWall,  tmp_eAngBWall = np.loadtxt(infile, skiprows=2, unpack=True)
    AngBWall =np.append(AngBWall, tmp_AngBWall)
    eAngBWall=np.append(eAngBWall, tmp_eAngBWall) 
infiles = ["data/Lengths/data_Hook.txt"] 
for infile in infiles:
    tmp_Hook,  tmp_eHook = np.loadtxt(infile, skiprows=2, unpack=True)
    Hook =np.append(Hook, tmp_Hook)
    eHook=np.append(eHook, tmp_eHook) 
infiles = ["data/Lengths/data_Pendul.txt"] 
for infile in infiles:
    tmp_Pendul,  tmp_ePendul = np.loadtxt(infile, skiprows=2, unpack=True)
    Pendul=np.append(Pendul, tmp_Pendul)
    ePendul=np.append(ePendul, tmp_ePendul)         
infiles = ["data/Lengths/data_StringL.txt"] 
for infile in infiles:
    tmp_StringL,  tmp_eStringL = np.loadtxt(infile, skiprows=2, unpack=True)
    StringL =np.append(StringL, tmp_StringL)
    eStringL=np.append(eStringL, tmp_eStringL) 
infiles = ["data/Lengths/data_StringR.txt"] 
for infile in infiles:
    tmp_StringR,  tmp_eStringR = np.loadtxt(infile, skiprows=2, unpack=True)
    StringR =np.append(StringR, tmp_StringR)
    eStringR=np.append(eStringR, tmp_eStringR)     

### Calculating length of the pendulum
# String - hook + pendulum/2
# Combined string length is computed first
String_L_array = np.hstack((StringL,StringR))
String_Lerr_array = np.hstack((eStringL,eStringR))
String_L_combined =  (sum(String_L_array)) / 8
String_Lerr_combined = sum(np.sqrt((String_Lerr_array**2) / 8**2))
Pendul = Pendul*0.01 # Converting pendulum length from cm to m

# Total length of pendulum is calculated
Pendulum_L = String_L_combined - Hook + (Pendul / 2)
Pendulum_Lerr = np.sqrt(String_Lerr_combined**2 + eHook**2 + (ePendul**2 / 2**2))

# Combined length of 4 pendulum-length and corresponding errors is computed
Pendulum_L_combined = sum(Pendulum_L) / 4
Pendulum_Lerr_combined = sum(np.sqrt(Pendulum_Lerr**2 / 4**2))



######### Pendulum time series read in 
timer_C1=np.array([])
timer_C2=np.array([])
timer_Cr1=np.array([])
timer_Cr2=np.array([])
timer_R1=np.array([])
timer_R2=np.array([])
timer_Z1=np.array([])
timer_Z2=np.array([])

infiles = ["data/Timer_Dat/timer_C1.dat"]
for infile in infiles:
        n, tmp_timer_C1=np.loadtxt(infile, skiprows=0, unpack=True)
        timer_C1=np.append(timer_C1, tmp_timer_C1)
infiles = ["data/Timer_Dat/timer_C2.dat"]
for infile in infiles:
        n, tmp_timer_C2=np.loadtxt(infile, skiprows=0, unpack=True)
        timer_C2=np.append(timer_C2, tmp_timer_C2)
infiles = ["data/Timer_Dat/timer_Cr1.dat"]
for infile in infiles:
        n, tmp_timer_Cr1=np.loadtxt(infile, skiprows=0, unpack=True)
        timer_Cr1=np.append(timer_Cr1, tmp_timer_Cr1)        
infiles = ["data/Timer_Dat/timer_Cr2.dat"]
for infile in infiles:
        n, tmp_timer_Cr2=np.loadtxt(infile, skiprows=0, unpack=True)
        timer_Cr2=np.append(timer_Cr2, tmp_timer_Cr2)        
infiles = ["data/Timer_Dat/timer_R1.dat"]
for infile in infiles:
        n, tmp_timer_R1=np.loadtxt(infile, skiprows=0, unpack=True)
        timer_R1=np.append(timer_R1, tmp_timer_R1)        
infiles = ["data/Timer_Dat/timer_R2.dat"]
for infile in infiles:
        n, tmp_timer_R2=np.loadtxt(infile, skiprows=0, unpack=True)
        timer_R2=np.append(timer_R2, tmp_timer_R2)        
infiles = ["data/Timer_Dat/timer_Z1.dat"]
for infile in infiles:
        n, tmp_timer_Z1=np.loadtxt(infile, skiprows=0, unpack=True)
        timer_Z1=np.append(timer_Z1, tmp_timer_Z1)        
infiles = ["data/Timer_Dat/timer_Z2.dat"]
for infile in infiles:
        n, tmp_timer_Z2=np.loadtxt(infile, skiprows=0, unpack=True)
        timer_Z2=np.append(timer_Z2, tmp_timer_Z2)
 
PendPlots(timer_R2)   
    
# Perioden med usikkerheder regnes    
T_comb=[]
sigmaT_comb=[]
T, sigmaT, Res=PendFit(timer_C1);
T_comb.append(T), sigmaT_comb.append(sigmaT)
T, sigmaT, Res=PendFit(timer_C2);
T_comb.append(T), sigmaT_comb.append(sigmaT)
T, sigmaT, Res=PendFit(timer_Cr1);
T_comb.append(T), sigmaT_comb.append(sigmaT)
T, sigmaT, Res=PendFit(timer_Cr2);
T_comb.append(T), sigmaT_comb.append(sigmaT)
T, sigmaT, Res=PendFit(timer_R1);
T_comb.append(T), sigmaT_comb.append(sigmaT)
T, sigmaT, Res=PendFit(timer_R2);
T_comb.append(T), sigmaT_comb.append(sigmaT)
T, sigmaT, Res=PendFit(timer_Z1);
T_comb.append(T), sigmaT_comb.append(sigmaT)
T, sigmaT, Res=PendFit(timer_Z2);
T_comb.append(T), sigmaT_comb.append(sigmaT)

T_comb=np.array([T_comb])
T_mean=T_comb.mean()
T=T_mean

sigmaT_comb=np.array([sigmaT_comb])
#sigmaT= ??? RMS of residuals ?

# g regnes
g=4*np.pi**2*Pendulum_L_combined/T**2
print(g)


