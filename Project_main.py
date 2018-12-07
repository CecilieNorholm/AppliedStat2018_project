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

# Read in data from files and put them into arrays
sys.path.append('/External_Functions')
from ExternalFunctions import nice_string_output, add_text_to_ax # useful functions to print fit results on figure


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
 
                 
def linear(x,a,b):
    y=a*x+b
    return y
        

N_var = 2                     # Number of variables (p0:p3)
N_dof = len(x) - N_var   # Number of degrees of freedom
from scipy import stats
chi2_prob = stats.chi2.sf(chi2, N_dof) # The chi2 probability given N_DOF degrees of freedom

fig, ax = plt.subplots(figsize=(12, 8))
y=timer_Cr1
x=n
ax.plot(x,y, 'o')
ax.set(xlabel="Measurement number",ylabel="Time elapsed (s)")


chi2_object = Chi2Regression(linear,x, y)
minuit = Minuit(chi2_object, pedantic=False, a=0,b=0) #   
minuit.migrad()  # perform the actual fit
chi2 = minuit.fval


xaxis = n
yaxis = linear(xaxis, *minuit.args)
ax.plot(xaxis,yaxis)
ax.set(xlabel="Measurement number",ylabel="Time elapsed (s)")

d = {'Chi2':     chi2,
     'ndf':      N_dof,
     'Prob':     chi2_prob,
     'Haeldning': [minuit.values['a'], minuit.errors['a']],
     'Skaering': [minuit.values['b'], minuit.errors['b']]
    }

string = nice_string_output(d, extra_spacing=2, decimals=3)
add_text_to_ax(0.02, 0.97, string, ax, fontsize=14)


fig
    