import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.io import FortranFile
import glob
from scipy.ndimage import zoom

size = 12
med_size = 13
big_size = 14

plt.rc('font', size=size)
plt.rc('axes', titlesize=size)
plt.rc('axes', labelsize=med_size)
plt.rc('xtick', labelsize=size)
plt.rc('ytick', labelsize=size)
plt.rc('legend', fontsize=size)
plt.rc('figure', titlesize=big_size)
plt.rcParams['figure.figsize'] = (4.5, 3)
plt.rcParams['figure.autolayout'] = True

cm_subsection = np.linspace(0.0, 1.0, 4) 
colors = [ cm.plasma(x) for x in cm_subsection ]

path = 'output/2d_'

reses = ['1e7', '5e6', '3e6', '2e6', '1e6', '5e5', '3e5', '2e5', '1e5', '5e4', 
         '3e4', '2e4', '1e4', '5e3', '3e3', '2e3', '1e3', '5e2']

nr = len(reses)

Id_f = np.zeros(nr)
Vd_f = np.zeros(nr)
Id_f2 = np.zeros(nr)
Vd_f2 = np.zeros(nr)

#fig1 = plt.figure(1)
#plt.xlabel('Time')
#plt.ylabel('Discharge Current')

#fig2 = plt.figure(2)
#plt.xlabel('Time')
#plt.ylabel('Discharge Voltage')

for i in range(nr):
    res = reses[i]
    path = 'output/res_' + res

    t = np.fromfile(path + '/time.dat', dtype=float)
    Id = np.fromfile(path + '/id.dat', dtype=float)
    Vd = np.fromfile(path + '/vd.dat', dtype=float)/1e5
    
    Id_f[i] = Id[-1]
    Vd_f[i] = Vd[-1]
    
    path = 'output/2d_res_' + res

    t = np.fromfile(path + '/time.dat', dtype=float)
    Id = np.fromfile(path + '/id.dat', dtype=float)
    Vd = np.fromfile(path + '/vd.dat', dtype=float)/1e5
    
    Id_f2[i] = Id[-1]
    Vd_f2[i] = Vd[-1]
    
    #plt.figure(1)
    #plt.plot(t[:len(Id)], -Id, label = res)
    
    #plt.figure(2)
    #plt.plot(t[:len(Vd)], Vd, label = res)
    

fig3 = plt.figure()
plt.plot(-Id_f, Vd_f, label='1d')
plt.plot(-Id_f2, Vd_f2, label='2d')
plt.ylabel('Discharge Voltage')
plt.xlabel('Discharge Current')
plt.xscale('log')
plt.legend()


plt.show()
    
    
    

