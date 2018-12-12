import numpy as np
import math as m
import sys
import csv
from iminuit import Minuit
from matplotlib import pyplot as plt
from probfit import Chi2Regression
sys.path.append('../../External_Functions')
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

print(accA_c-accB_c)
print(accA_cr-accB_cr)
print(accA_r-accB_r)
print(accA_z-accB_z)


#load angles
angA_c, angB_c  = np.loadtxt('Angles/Ang_C.txt', unpack=True)
angA_error_c, angB_error_c = np.loadtxt('Angles/Ang_eC.txt', unpack=True)
angA_cr, angB_cr  = np.loadtxt('Angles/Ang_Cr.txt', unpack=True)
angA_error_cr, angB_error_cr = np.loadtxt('Angles/Ang_eCr.txt', unpack=True)
angA_r, angB_r  = np.loadtxt('Angles/Ang_R.txt', unpack=True)
angA_error_r, angB_error_r = np.loadtxt('Angles/Ang_eR.txt', unpack=True)
angA_z, angB_z  = np.loadtxt('Angles/Ang_Z.txt', unpack=True)
angA_error_z, angB_error_z = np.loadtxt('Angles/Ang_eZ.txt', unpack=True)

#compute g
acc_c = (accA_c + accB_c) / 2.
ang_c = (angA_c + angB_c) / 2.
g_c = g_Incline(acc_c,Crail,Cball,ang_c)

acc_cr = (accA_cr + accB_cr) / 2.
ang_cr = (angA_cr + angB_cr) / 2.
g_cr = g_Incline(acc_cr,Crrail,Crball,ang_cr)

acc_r = (accA_r + accB_r) / 2.
ang_r = (angA_r + angB_r) / 2.
g_r = g_Incline(acc_r,Rrail,Rball,ang_r)

acc_z = (accA_z + accB_z) / 2.
ang_z = (angA_z + angB_z) / 2.
g_z = g_Incline(acc_z,Zrail,Zball,ang_z)

#
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

# #compute error g
# #CECILE
# error_g_c = SigmaIncline(acc_c,Crail,Cball,ang_c,sigma_a_c,eCrail,eCball,sigma_theta_c)
# print(error_g_c)
#
# #CHRISTIAN
# error_g_cr = SigmaIncline(acc_cr,Crrail,Crball,ang_cr,sigma_a_cr,eCrrail,eCrball,sigma_theta_cr)
# print(error_g_cr)
#
# #RUNE
# error_g_r = SigmaIncline(acc_r,Rrail,Rball,ang_r,sigma_a_r,eRrail,eRball,sigma_theta_r)
# print(error_g_r)
#
# #ZLATKO
# error_g_z = SigmaIncline(acc_z,Zrail,Zball,ang_z,sigma_a_z,eZrail,eZball,sigma_theta_z)
# print(error_g_z)
