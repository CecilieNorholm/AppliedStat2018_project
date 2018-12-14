import numpy as np
import math as m
import sys
import csv
from iminuit import Minuit
from matplotlib import pyplot as plt
from probfit import Chi2Regression
sys.path.append('../External_Functions')
from ExternalFunctions import nice_string_output, add_text_to_ax # useful functions to print fit results on figure
from scipy import stats
from gIncline import*
from SigmaIncline import*


#LOADING DATA

#DIAMETER OF BALL AND ERROR
Dball, eball = np.loadtxt('Lengths/data_DBall.txt', unpack=True, skiprows=2)
Rball = Dball[0]
eRball = eball[0]
Zball = Dball[1]
eZball = eball[1]
Cball = Dball[2]
eCball = eball[2]
Crball = Dball[3]
eCrball = eball[3]

#LENGHT OF RAIL AND ERROR
Drail, erail = np.loadtxt('Lengths/data_DRail.txt', unpack=True, skiprows=2)
Rrail = Drail[0]
eRrail = erail[0]
Zrail = Drail[1]
eZrail = erail[1]
Crail = Drail[2]
eCrail = erail[2]
Crrail = Drail[3]
eCrrail = erail[3]

#load acceleration
accA_c, acc_errorA_c = np.loadtxt('cecilie_acc_a.txt', unpack=True)
accB_c, acc_errorB_c = np.loadtxt('cecilie_acc_b.txt', unpack=True)
accA_cr, acc_errorA_cr = np.loadtxt('christian_acc_a.txt', unpack=True)
accB_cr, acc_errorB_cr = np.loadtxt('christian_acc_b.txt', unpack=True)
accA_r, acc_errorA_r = np.loadtxt('rune_acc_a.txt', unpack=True)
accB_r, acc_errorB_r = np.loadtxt('rune_acc_b.txt', unpack=True)
accA_z, acc_errorA_z = np.loadtxt('zlatko_acc_a.txt', unpack=True)
accB_z, acc_errorB_z = np.loadtxt('zlatko_acc_b.txt', unpack=True)

#load angles
angA_c, angB_c  = np.loadtxt('Angles/Ang_C.txt', unpack=True)
angA_error_c, angB_error_c = np.loadtxt('Angles/Ang_eC.txt', unpack=True)
angA_cr, angB_cr  = np.loadtxt('Angles/Ang_Cr.txt', unpack=True)
angA_error_cr, angB_error_cr = np.loadtxt('Angles/Ang_eCr.txt', unpack=True)
angA_r, angB_r  = np.loadtxt('Angles/Ang_R.txt', unpack=True)
angA_error_r, angB_error_r = np.loadtxt('Angles/Ang_eR.txt', unpack=True)
angA_z, angB_z  = np.loadtxt('Angles/Ang_Z.txt', unpack=True)
angA_error_z, angB_error_z = np.loadtxt('Angles/Ang_eZ.txt', unpack=True)

phi_c = angA_c - angB_c
phi_cr = angA_cr - angB_cr
phi_r = angA_r - angB_r
phi_z = angA_z - angB_z



#compute g
acc_c = (accA_c + accB_c) / 2.
ang_c = (angA_c + angB_c) / 2.

g_c = g_Incline(acc_c,Crail,Cball,ang_c)

# # TEST
# print(f"mean acc {g_c}")
# gA_c = g_Incline(accA_c,Rrail,Rball,angA_c)
# print(f"acc on A side {gA_c}")
# gB_c = g_Incline(accB_c,Rrail,Rball,angB_c)
# print(f"acc on B side {gB_c}")
# print(f"average acc og A and B side {(gA_c+gB_c)/2}")
# # END OF TEST



acc_cr = (accA_cr + accB_cr) / 2.
ang_cr = (angA_cr + angB_cr) / 2.
g_cr = g_Incline(acc_cr,Crrail,Crball,ang_cr)


acc_r = (accA_r + accB_r) / 2.
ang_r = (angA_r + angB_r) / 2.
g_r = g_Incline(acc_r,Rrail,Rball,ang_r)

# TEST
print(f"mean acc {g_r}")
gA_r = g_Incline(accA_r,Rrail,Rball,angA_r)
print(f"acc on A side {gA_r}")
gB_r = g_Incline(accB_r,Rrail,Rball,angB_r)
print(f"acc on B side {gB_r}")
print(f"average acc og A and B side {(gA_r+gB_r)/2}")
# END OF TEST


acc_z = (accA_z + accB_z) / 2.
ang_z = (angA_z + angB_z) / 2.
g_z = g_Incline(acc_z,Zrail,Zball,ang_z)

# print(g_c)
# print(g_cr)
# print(g_r)
# print(g_z)


#COMPUTE MEAN ERROR ON ACC. AND ANGLE
sigma_a_c = np.sqrt( (0.5 * acc_errorA_c)**2 + (0.5 * acc_errorB_c)**2 )
sigma_theta_c = np.sqrt( (0.5 * angA_error_c)**2 + (0.5 * angB_error_c)**2 )

sigma_a_cr = np.sqrt( (0.5 * acc_errorA_cr)**2 + (0.5 * acc_errorB_cr)**2 )
sigma_theta_cr = np.sqrt( (0.5 * angA_error_cr)**2 + (0.5 * angB_error_cr)**2 )

sigma_a_r = np.sqrt( (0.5 * acc_errorA_r)**2 + (0.5 * acc_errorB_r)**2 )
sigma_theta_r = np.sqrt( (0.5 * angA_error_r)**2 + (0.5 * angB_error_r)**2 )

sigma_a_z = np.sqrt( (0.5 * acc_errorA_z)**2 + (0.5 * acc_errorB_z)**2 )
sigma_theta_z = np.sqrt( (0.5 * angA_error_z)**2 + (0.5 * angB_error_z)**2 )



#compute error g
#CECILE
error_g_c = SigmaIncline(acc_c,Crail,Cball,ang_c,sigma_a_c,eCrail,eCball,sigma_theta_c)
print(error_g_c)

#CHRISTIAN
error_g_cr = SigmaIncline(acc_cr,Crrail,Crball,ang_cr,sigma_a_cr,eCrrail,eCrball,sigma_theta_cr)
print(error_g_cr)

#RUNE
error_g_r = SigmaIncline(acc_r,Rrail,Rball,ang_r,sigma_a_r,eRrail,eRball,sigma_theta_r)
print(error_g_r)

#ZLATKO
error_g_z = SigmaIncline(acc_z,Zrail,Zball,ang_z,sigma_a_z,eZrail,eZball,sigma_theta_z)
print(error_g_z)

#COLLECT IN LIST
g_full = [g_c, g_cr, g_r, g_z]
print(g_full)
error_full = [error_g_c, error_g_cr, error_g_r, error_g_z]

#CONVERT TO NUMPY ARRAY
g_full = np.asarray(g_full)
error_full = np.asarray(error_full)


g_full_mean = np.average(g_full, weights = error_full)
#error_full_mean = np.mean(error_full) / np.sqrt(4)

error_full_mean = np.sqrt( 1. / np.sum(1./error_full**2) )

print(f"Our measuremnt : {g_full_mean} m/s^2")
print(f"Our error : {error_full_mean} m/s^2")

g_dtu = 9.81546301

error_dtu = 10**(-8)


person = [1,2,3,4]
name = ['Cecilie','Christian','Rune','Zlatko']

fig1, ax = plt.subplots(figsize=(10,6))

ax.errorbar(person,g_full,error_full,fmt='ro',ecolor ='k',elinewidth=1,capsize=2,capthick=1)
ax.axhline(y=g_full_mean, color='r', linestyle='-',label = "Measurement in First Lab for Ball on Incline Experiment $g = (9.60 \pm 0.06)$ ${\\rm m/s^2}$")
ax.axhline(y=g_full_mean+error_full_mean, color='r', linestyle=':')
ax.axhline(y=g_full_mean-error_full_mean, color='r', linestyle=':')
ax.axhline(y=g_dtu, color='b', linestyle='--',label = "Measurement by DTU at CPH University $g = (9.81546301 \pm 0.00000001)$  ${\\rm m/s^2}$")
ax.set_xlabel('Individual')
ax.set_xticks(person)
ax.set_xticklabels(name)
ax.set_ylabel("Gravitational acceleration ${\\rm [m/s^2]}$")
ax.set_title("Resulting Gravity Measurements in First Lab")
ax.legend()
fig1.tight_layout()
fig1.savefig("Ball_Incline_g.pdf",dpi=600)
fig1
