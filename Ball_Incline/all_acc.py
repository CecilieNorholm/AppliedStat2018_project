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

##################### RUNE ########################
#Loading data
d = np.loadtxt('data_Gates.txt',skiprows=2)
tA = np.loadtxt('Incline_sorted_data/Gate_timings/gates_R_a.txt',delimiter=',')
tB = np.loadtxt('Incline_sorted_data/Gate_timings/gates_R_b.txt',delimiter=',')
sigma_d = d[0:5,1]
sigma_t = np.array([0.0001,0.0001,0.0001,0.0001,0.0001])
D_ball = 0.01270 #diameter of ball

#Declaring constants
Npoints = 5
Nvar = 1
Ndof_fit = Npoints - Nvar

#Initialize list full acc and error for each run
acc_array_a = []
acc_array_b = []
error_array_a = []
error_array_b = []

#Declaring Fitting function (const. acceleration)
def fit_func(t,a,v0,s0):
    return 0.5 * a * t**2 + v0*t+s0

#Computes acceleration on 'A' side of every run
for i in range(5):

    d_a = d[0:5,0] #CHANGE THIS FOR EVERY PERSON!!!!!!!!!!
    print(d_a)
    t_a = tA[i,0:5]

    chi2_object = Chi2Regression(fit_func, t_a, d_a,error = sigma_d)
    minuit = Minuit(chi2_object, pedantic=False, a = 9.8, v0 = 16, s0 = 2)
    minuit.migrad();
    chi2_trans = minuit.fval

    acc_a = minuit.values['a']
    acc_array_a.append(acc_a)

    error_acc_a = minuit.errors['a']
    error_array_a.append(error_acc_a) #error on every full run

#PLOTTING A
#OBS - REMEMBER TITLE AND AXIS AND LIMITS
x = np.linspace(0,2,100)
fig1, ax = plt.subplots(figsize=(10,6))
ax.errorbar(t_a,d_a,sigma_d,fmt='ro',ecolor ='k',elinewidth=1,capsize=2,capthick=1)
ax.plot(x,fit_func(x,*minuit.args),'--b',label='Chi2Fit')
ax.set_title('Eqation of Motion with Constant Acceleration for the Ball on Incline Experiment')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Distance [m]')
ax.set_xlim(0,1.5)
a = minuit.values['a']
sigma_a = minuit.errors['a']
Chi2_fit = minuit.fval
Prob_fit = stats.chi2.sf(Chi2_fit, Ndof_fit)
info = {'a':   [a, sigma_a],
        'Chi2':     Chi2_fit,
        'ndf':      Ndof_fit,
        'Prob':     Prob_fit,
        }
text = nice_string_output(info, extra_spacing=2, decimals=3)
add_text_to_ax(0.02, 0.95, text,ax, fontsize=14)
fig1.tight_layout()
fig1.savefig("EoM_Ball_Incline.pdf",dpi=600)
fig1

#Computes acceleration on 'B' side of every run
for i in range(5):

    d_b = d[0:5,0]
    t_b = tB[i,0:5]

    chi2_object = Chi2Regression(fit_func, t_b, d_b,error = sigma_d)
    minuit = Minuit(chi2_object, pedantic=False, a = 9.8, v0 = 16, s0 = 2)
    minuit.migrad();
    chi2_trans = minuit.fval

    acc_b = minuit.values['a']
    acc_array_b.append(acc_b)

    error_acc_b = minuit.errors['a']
    error_array_b.append(error_acc_b) #error on every full run

# #PLOTTING B
# #OBS - REMEMBER TITLE AND AXIS AND LIMITS
# x = np.linspace(0,2,100)
# fig2, ax = plt.subplots(figsize=(10,6))
# ax.errorbar(t_b,d_b,sigma_d,fmt='ro',ecolor ='k',elinewidth=1,capsize=2,capthick=1)
# ax.plot(x,fit_func(x,*minuit.args),'-r',label='Chi2Fit')
# ax.set_title('Constant Acceleration - Ball Incline Experiment B side')
# ax.set_xlabel('Time [s]')
# ax.set_ylabel('Distance [m]')
# ax.set_xlim(0,1.5)
# a = minuit.values['a']
# sigma_a = minuit.errors['a']
# Chi2_fit = minuit.fval
# Prob_fit = stats.chi2.sf(Chi2_fit, Ndof_fit)
# info = {'a':   [a, sigma_a],
#         'Chi2':     Chi2_fit,
#         'ndf':      Ndof_fit,
#         'Prob':     Prob_fit,
#         }
# text = nice_string_output(info, extra_spacing=2, decimals=3)
# add_text_to_ax(0.02, 0.95, text,ax, fontsize=14)
# fig2.tight_layout()
# fig2

#Convert to numpy array
acc_array_a = np.asarray(acc_array_a)
acc_array_b = np.asarray(acc_array_b)
error_array_a = np.asarray(error_array_a)
error_array_b = np.asarray(error_array_b)

#Computes the mean acceleration and error on mean
mean_acc_a = np.mean(acc_array_a)
mean_acc_b = np.mean(acc_array_b)
error_mean_acc_a = np.mean(error_array_a) / np.sqrt(len(error_array_a))
error_mean_acc_b = np.mean(error_array_b) / np.sqrt(len(error_array_b))

#change to list and combine mean and error in one list
acc_full_a = [mean_acc_a,error_mean_acc_a]
acc_full_b = [mean_acc_b,error_mean_acc_b]


# csvfile = "/Users/zsaldic/Documents/nbi/AppStat2018/Project/Ball_Incline/rune_acc_a.txt"
# with open(csvfile,"w") as output:
#     writer = csv.writer(output,lineterminator='\n')
#     for val in acc_full_a:
#         writer.writerow([val])
#
# csvfile = "/Users/zsaldic/Documents/nbi/AppStat2018/Project/Ball_Incline/rune_acc_b.txt"
# with open(csvfile,"w") as output:
#     writer = csv.writer(output,lineterminator='\n')
#     for val in acc_full_b:
#         writer.writerow([val])
#==============================================================================
#################### ZLATKO ########################
#Loading data
d = np.loadtxt('data_Gates.txt',skiprows=2)
tA = np.loadtxt('Incline_sorted_data/Gate_timings/gates_Z_a.txt',delimiter=',')
tB = np.loadtxt('Incline_sorted_data/Gate_timings/gates_Z_b.txt',delimiter=',')
sigma_d = d[5:10,1]
sigma_t = np.array([0.0001,0.0001,0.0001,0.0001,0.0001])
D_ball = 0.01270 #diameter of ball

#Declaring constants
Npoints = 5
Nvar = 1
Ndof_fit = Npoints - Nvar

#Initialize list full acc and error for each run
acc_array_a = []
acc_array_b = []
error_array_a = []
error_array_b = []

#Declaring Fitting function (const. acceleration)
def fit_func(t,a,v0,s0):
    return 0.5 * a * t**2 + v0*t+s0

#Computes acceleration on 'A' side of every run
for i in range(5):

    d_a = d[5:10,0] #CHANGE FOR ZLATKO
    t_a = tA[i,0:5]

    chi2_object = Chi2Regression(fit_func, t_a, d_a,error = sigma_d)
    minuit = Minuit(chi2_object, pedantic=False, a = 9.8, v0 = 16, s0 = 2)
    minuit.migrad();
    chi2_trans = minuit.fval

    acc_a = minuit.values['a']
    acc_array_a.append(acc_a)

    error_acc_a = minuit.errors['a']
    error_array_a.append(error_acc_a) #error on every full run

# #PLOTTING A
# #OBS - REMEMBER TITLE AND AXIS AND LIMITS
# x = np.linspace(0,2,100)
# fig1, ax = plt.subplots(figsize=(10,6))
# ax.errorbar(t_a,d_a,sigma_d,fmt='ro',ecolor ='k',elinewidth=1,capsize=2,capthick=1)
# ax.plot(x,fit_func(x,*minuit.args),'-r',label='Chi2Fit')
# ax.set_title('Constant Acceleration - Ball Incline Experiment A side')
# ax.set_xlabel('Time [s]')
# ax.set_ylabel('Distance [m]')
# ax.set_xlim(0,1.5)
# a = minuit.values['a']
# sigma_a = minuit.errors['a']
# Chi2_fit = minuit.fval
# Prob_fit = stats.chi2.sf(Chi2_fit, Ndof_fit)
# info = {'a':   [a, sigma_a],
#         'Chi2':     Chi2_fit,
#         'ndf':      Ndof_fit,
#         'Prob':     Prob_fit,
#         }
# text = nice_string_output(info, extra_spacing=2, decimals=3)
# add_text_to_ax(0.02, 0.95, text,ax, fontsize=14)
# fig1.tight_layout()
# fig1

#Computes acceleration on 'B' side of every run
for i in range(5):

    d_b = d[5:10,0]
    t_b = tB[i,0:5]

    chi2_object = Chi2Regression(fit_func, t_b, d_b,error = sigma_d)
    minuit = Minuit(chi2_object, pedantic=False, a = 9.8, v0 = 16, s0 = 2)
    minuit.migrad();
    chi2_trans = minuit.fval

    acc_b = minuit.values['a']
    acc_array_b.append(acc_b)

    error_acc_b = minuit.errors['a']
    error_array_b.append(error_acc_b) #error on every full run

# #PLOTTING B
# #OBS - REMEMBER TITLE AND AXIS AND LIMITS
# x = np.linspace(0,2,100)
# fig2, ax = plt.subplots(figsize=(10,6))
# ax.errorbar(t_b,d_b,sigma_d,fmt='ro',ecolor ='k',elinewidth=1,capsize=2,capthick=1)
# ax.plot(x,fit_func(x,*minuit.args),'-r',label='Chi2Fit')
# ax.set_title('Constant Acceleration - Ball Incline Experiment B side')
# ax.set_xlabel('Time [s]')
# ax.set_ylabel('Distance [m]')
# ax.set_xlim(0,1.5)
# a = minuit.values['a']
# sigma_a = minuit.errors['a']
# Chi2_fit = minuit.fval
# Prob_fit = stats.chi2.sf(Chi2_fit, Ndof_fit)
# info = {'a':   [a, sigma_a],
#         'Chi2':     Chi2_fit,
#         'ndf':      Ndof_fit,
#         'Prob':     Prob_fit,
#         }
# text = nice_string_output(info, extra_spacing=2, decimals=3)
# add_text_to_ax(0.02, 0.95, text,ax, fontsize=14)
# fig2.tight_layout()
# fig2

#Convert to numpy array
acc_array_a = np.asarray(acc_array_a)
acc_array_b = np.asarray(acc_array_b)
error_array_a = np.asarray(error_array_a)
error_array_b = np.asarray(error_array_b)

#Computes the mean acceleration and error on mean
mean_acc_a = np.mean(acc_array_a)
mean_acc_b = np.mean(acc_array_b)
error_mean_acc_a = np.mean(error_array_a) / np.sqrt(len(error_array_a))
error_mean_acc_b = np.mean(error_array_b) / np.sqrt(len(error_array_b))

#change to list and combine mean and error in one list
acc_full_a = [mean_acc_a,error_mean_acc_a]
acc_full_b = [mean_acc_b,error_mean_acc_b]

# csvfile = "/Users/zsaldic/Documents/nbi/AppStat2018/Project/Ball_Incline/zlatko_acc_a.txt"
# with open(csvfile,"w") as output:
#     writer = csv.writer(output,lineterminator='\n')
#     for val in acc_full_a:
#         writer.writerow([val])
#
# csvfile = "/Users/zsaldic/Documents/nbi/AppStat2018/Project/Ball_Incline/zlatko_acc_b.txt"
# with open(csvfile,"w") as output:
#     writer = csv.writer(output,lineterminator='\n')
#     for val in acc_full_b:
#         writer.writerow([val])
#==============================================================================

#################### CHRISTIAN ########################
#Loading data
d = np.loadtxt('data_Gates.txt',skiprows=2)
tA = np.loadtxt('Incline_sorted_data/Gate_timings/gates_Cr_a.txt',delimiter=',')
tB = np.loadtxt('Incline_sorted_data/Gate_timings/gates_Cr_b.txt',delimiter=',')
sigma_d = d[10:15,1]
sigma_t = np.array([0.0001,0.0001,0.0001,0.0001,0.0001])
D_ball = 0.01270 #diameter of ball

#Declaring constants
Npoints = 5
Nvar = 1
Ndof_fit = Npoints - Nvar

#Initialize list full acc and error for each run
acc_array_a = []
acc_array_b = []
error_array_a = []
error_array_b = []

#Declaring Fitting function (const. acceleration)
def fit_func(t,a,v0,s0):
    return 0.5 * a * t**2 + v0*t+s0

#Computes acceleration on 'A' side of every run
for i in range(5):

    d_a = d[10:15,0]
    t_a = tA[i,0:5]

    chi2_object = Chi2Regression(fit_func, t_a, d_a,error = sigma_d)
    minuit = Minuit(chi2_object, pedantic=False, a = 9.8, v0 = 16, s0 = 2)
    minuit.migrad();
    chi2_trans = minuit.fval

    acc_a = minuit.values['a']
    acc_array_a.append(acc_a)

    error_acc_a = minuit.errors['a']
    error_array_a.append(error_acc_a) #error on every full run

# #PLOTTING A
# #OBS - REMEMBER TITLE AND AXIS AND LIMITS
# x = np.linspace(0,2,100)
# fig1, ax = plt.subplots(figsize=(10,6))
# ax.errorbar(t_a,d_a,sigma_d,fmt='ro',ecolor ='k',elinewidth=1,capsize=2,capthick=1)
# ax.plot(x,fit_func(x,*minuit.args),'-r',label='Chi2Fit')
# ax.set_title('Constant Acceleration - Ball Incline Experiment A side')
# ax.set_xlabel('Time [s]')
# ax.set_ylabel('Distance [m]')
# ax.set_xlim(0,1.5)
# a = minuit.values['a']
# sigma_a = minuit.errors['a']
# Chi2_fit = minuit.fval
# Prob_fit = stats.chi2.sf(Chi2_fit, Ndof_fit)
# info = {'a':   [a, sigma_a],
#         'Chi2':     Chi2_fit,
#         'ndf':      Ndof_fit,
#         'Prob':     Prob_fit,
#         }
# text = nice_string_output(info, extra_spacing=2, decimals=3)
# add_text_to_ax(0.02, 0.95, text,ax, fontsize=14)
# fig1.tight_layout()
# fig1

#Computes acceleration on 'B' side of every run
for i in range(5):

    d_b = d[10:15,0] #CHANGE FOR CHRISTIAN
    t_b = tB[i,0:5]

    chi2_object = Chi2Regression(fit_func, t_b, d_b,error = sigma_d)
    minuit = Minuit(chi2_object, pedantic=False, a = 9.8, v0 = 16, s0 = 2)
    minuit.migrad();
    chi2_trans = minuit.fval

    acc_b = minuit.values['a']
    acc_array_b.append(acc_b)

    error_acc_b = minuit.errors['a']
    error_array_b.append(error_acc_b) #error on every full run

# #PLOTTING B
# #OBS - REMEMBER TITLE AND AXIS AND LIMITS
# x = np.linspace(0,2,100)
# fig2, ax = plt.subplots(figsize=(10,6))
# ax.errorbar(t_b,d_b,sigma_d,fmt='ro',ecolor ='k',elinewidth=1,capsize=2,capthick=1)
# ax.plot(x,fit_func(x,*minuit.args),'-r',label='Chi2Fit')
# ax.set_title('Constant Acceleration - Ball Incline Experiment B side')
# ax.set_xlabel('Time [s]')
# ax.set_ylabel('Distance [m]')
# ax.set_xlim(0,1.5)
# a = minuit.values['a']
# sigma_a = minuit.errors['a']
# Chi2_fit = minuit.fval
# Prob_fit = stats.chi2.sf(Chi2_fit, Ndof_fit)
# info = {'a':   [a, sigma_a],
#         'Chi2':     Chi2_fit,
#         'ndf':      Ndof_fit,
#         'Prob':     Prob_fit,
#         }
# text = nice_string_output(info, extra_spacing=2, decimals=3)
# add_text_to_ax(0.02, 0.95, text,ax, fontsize=14)
# fig2.tight_layout()
# fig2

#Convert to numpy array
acc_array_a = np.asarray(acc_array_a)
acc_array_b = np.asarray(acc_array_b)
error_array_a = np.asarray(error_array_a)
error_array_b = np.asarray(error_array_b)

#Computes the mean acceleration and error on mean
mean_acc_a = np.mean(acc_array_a)
mean_acc_b = np.mean(acc_array_b)
error_mean_acc_a = np.mean(error_array_a) / np.sqrt(len(error_array_a))
error_mean_acc_b = np.mean(error_array_b) / np.sqrt(len(error_array_b))

#change to list and combine mean and error in one list
acc_full_a = [mean_acc_a,error_mean_acc_a]
acc_full_b = [mean_acc_b,error_mean_acc_b]

# csvfile = "/Users/zsaldic/Documents/nbi/AppStat2018/Project/Ball_Incline/christian_acc_a.txt"
# with open(csvfile,"w") as output:
#     writer = csv.writer(output,lineterminator='\n')
#     for val in acc_full_a:
#         writer.writerow([val])
#
# csvfile = "/Users/zsaldic/Documents/nbi/AppStat2018/Project/Ball_Incline/christian_acc_acc_b.txt"
# with open(csvfile,"w") as output:
#     writer = csv.writer(output,lineterminator='\n')
#     for val in acc_full_b:
#         writer.writerow([val])
##==============================================================================
#
##################### CECILIE ########################
#Loading data
d = np.loadtxt('data_Gates.txt',skiprows=2)
tA = np.loadtxt('Incline_sorted_data/Gate_timings/gates_C_a.txt',delimiter=',')
tB = np.loadtxt('Incline_sorted_data/Gate_timings/gates_C_b.txt',delimiter=',')
sigma_d = d[15:21,1]
sigma_t = np.array([0.0001,0.0001,0.0001,0.0001,0.0001])
D_ball = 0.01270 #diameter of ball

#Declaring constants
Npoints = 5
Nvar = 1
Ndof_fit = Npoints - Nvar

#Initialize list full acc and error for each run
acc_array_a = []
acc_array_b = []
error_array_a = []
error_array_b = []

#Declaring Fitting function (const. acceleration)
def fit_func(t,a,v0,s0):
    return 0.5 * a * t**2 + v0*t+s0

#Computes acceleration on 'A' side of every run
for i in range(5):

    d_a = d[15:21,0]
    t_a = tA[i,0:5]

    chi2_object = Chi2Regression(fit_func, t_a, d_a,error = sigma_d)
    minuit = Minuit(chi2_object, pedantic=False, a = 9.8, v0 = 16, s0 = 2)
    minuit.migrad();
    chi2_trans = minuit.fval

    acc_a = minuit.values['a']
    acc_array_a.append(acc_a)

    error_acc_a = minuit.errors['a']
    error_array_a.append(error_acc_a) #error on every full run

# #PLOTTING A
# #OBS - REMEMBER TITLE AND AXIS AND LIMITS
# x = np.linspace(0,2,100)
# fig1, ax = plt.subplots(figsize=(10,6))
# ax.errorbar(t_a,d_a,sigma_d,fmt='ro',ecolor ='k',elinewidth=1,capsize=2,capthick=1)
# ax.plot(x,fit_func(x,*minuit.args),'-r',label='Chi2Fit')
# ax.set_title('Constant Acceleration - Ball Incline Experiment A side')
# ax.set_xlabel('Time [s]')
# ax.set_ylabel('Distance [m]')
# ax.set_xlim(0,1.5)
# a = minuit.values['a']
# sigma_a = minuit.errors['a']
# Chi2_fit = minuit.fval
# Prob_fit = stats.chi2.sf(Chi2_fit, Ndof_fit)
# info = {'a':   [a, sigma_a],
#         'Chi2':     Chi2_fit,
#         'ndf':      Ndof_fit,
#         'Prob':     Prob_fit,
#         }
# text = nice_string_output(info, extra_spacing=2, decimals=3)
# add_text_to_ax(0.02, 0.95, text,ax, fontsize=14)
# fig1.tight_layout()
# fig1

#Computes acceleration on 'B' side of every run
for i in range(5):

    d_b = d[15:21,0]
    t_b = tB[i,0:5]

    chi2_object = Chi2Regression(fit_func, t_b, d_b,error = sigma_d)
    minuit = Minuit(chi2_object, pedantic=False, a = 9.8, v0 = 0.2, s0 = 0.1)
    minuit.migrad();
    chi2_trans = minuit.fval

    acc_b = minuit.values['a']
    acc_array_b.append(acc_b)

    error_acc_b = minuit.errors['a']
    error_array_b.append(error_acc_b) #error on every full run

# ##PLOTTING B
# x = np.linspace(0,2,100)
# fig2, ax = plt.subplots(figsize=(10,6))
# ax.errorbar(t_b,d_b,sigma_d,fmt='ro',ecolor ='k',elinewidth=1,capsize=2,capthick=1)
# ax.plot(x,fit_func(x,*minuit.args),'-r',label='Chi2Fit')
# ax.set_title('Constant Acceleration - Ball Incline Experiment B side')
# ax.set_xlabel('Time [s]')
# ax.set_ylabel('Distance [m]')
# ax.set_xlim(0,1.5)
# a = minuit.values['a']
# sigma_a = minuit.errors['a']
# Chi2_fit = minuit.fval
# Prob_fit = stats.chi2.sf(Chi2_fit, Ndof_fit)
# info = {'a':   [a, sigma_a],
#         'Chi2':     Chi2_fit,
#         'ndf':      Ndof_fit,
#         'Prob':     Prob_fit,
#         }
# text = nice_string_output(info, extra_spacing=2, decimals=3)
# add_text_to_ax(0.02, 0.95, text,ax, fontsize=14)
# fig2.tight_layout()
# fig2

#Convert to numpy array
acc_array_a = np.asarray(acc_array_a)
acc_array_b = np.asarray(acc_array_b)
error_array_a = np.asarray(error_array_a)
error_array_b = np.asarray(error_array_b)

#Computes the mean acceleration and error on mean
mean_acc_a = np.mean(acc_array_a)
mean_acc_b = np.mean(acc_array_b)
error_mean_acc_a = np.mean(error_array_a) / np.sqrt(len(error_array_a))
error_mean_acc_b = np.mean(error_array_b) / np.sqrt(len(error_array_b))

#change to list and combine mean and error in one list
acc_full_a = [mean_acc_a,error_mean_acc_a]
acc_full_b = [mean_acc_b,error_mean_acc_b]

# csvfile = "/Users/zsaldic/Documents/nbi/AppStat2018/Project/Ball_Incline/cecilie_acc_a.txt"
# with open(csvfile,"w") as output:
#     writer = csv.writer(output,lineterminator='\n')
#     for val in acc_full_a:
#         writer.writerow([val])
#
# csvfile = "/Users/zsaldic/Documents/nbi/AppStat2018/Project/Ball_Incline/cecilie_acc_b.txt"
# with open(csvfile,"w") as output:
#     writer = csv.writer(output,lineterminator='\n')
#     for val in acc_full_b:
#         writer.writerow([val])
