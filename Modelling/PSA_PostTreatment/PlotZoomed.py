'''
This script get the Modelica simulation files using the buildingspy library.
This allows to interpolate over the results to get the boundary conditions,i.e., DNI, T_su, T_ex as if they came from the experimental data
for the SF EuroThrough collectors and plot it for the whole day.
It plots the zoomed simulation results from some days

'''

from __future__ import division
import sys,os
import buildingspy.simulate.Simulator as si
from buildingspy.io.outputfile import Reader
from buildingspy.io.postprocess import Plotter
from scipy import optimize
import numpy as np
from modelicares import SimRes
import matplotlib.pyplot as plt
from scipy import optimize
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)



# --------------    SAVE GRAPHS   -----------------------------
SAVE_FIGURES = 'False'
# INSERT Here the directory path to save the figures
DIRECTORY_FIGURES = '/Users/adriano/Dropbox/GitHub/PSA_SFERAII/Modelling/ResultsFigure/' #'C:/Users/Admin/Documents/GitHub/PSA_SFERAII/Modelling/Results/'


# --------------    Select the day of experiments  -----------------
Days = ['20160629', '20160630',  '20160701',  '20160704',  '20160705', '20160706']
StartModTime = [46000,       40000,        37000,      43000,          43000,     41000]
StopModTime =  [58000,       54000,        53000,      57000,          56000,     50000]

# Define simulation file
FileSimulation = 'PTTL_SF_'


# Define Path where Modelica Results are stored:
ModResults = '/Users/adriano/Dropbox/GitHub/PSA_SFERAII/Modelling/ModelicaResults/' #'C:/Users/Admin/Desktop/SFERA_II_Sim/'
TimeSim = []
TimeExp = []
NN = []
T_su = []
T_ex = []
T_ex_exp = []
Delta_T = []
DNI = []
m_wf = []
T_amb = []
v_wind = []

DNIExp = []
T_suExp = []
T_exExp = []
m_wfExp = []

vector = [ 0, 1]
for KKK in range(len(StartModTime)):
    resultFile=(ModResults+Days[KKK]+'/'+FileSimulation+".mat")
    r=Reader(resultFile, "dymola")
    DNI.append(r.values('DNI.y'))
    NN.append(r.values('EuroTrough.N'))
    T_su.append(r.values('SensTsu.fluidState.T'))
    T_ex.append(r.values('SensTex.fluidState.T'))
    TimeSim.append(T_ex[KKK][0])
    T_ex_exp.append(r.values('t_htf_ex.y'))
    Delta_T.append(r.values('DeltaT.Delta'))
    m_wf.append(r.values('m_dot_htf.y'))

    # Create a time vector which values every 5 seconds

    TimeExp.append(np.linspace(StartModTime[KKK],StopModTime[KKK],num=int((StopModTime[KKK]-StartModTime[KKK])/5)))

    # Interpolate the Modelica results over the Time vector
    DNIExp.append(Plotter.interpolate(TimeExp[KKK], DNI[KKK][0], DNI[KKK][1]))
    T_suExp.append(Plotter.interpolate(TimeExp[KKK],T_su[KKK][0],T_su[KKK][1]-273.15))
    T_exExp.append(Plotter.interpolate(TimeExp[KKK],T_ex_exp[KKK][0],T_ex_exp[KKK][1]-273.15))
    m_wfExp.append(Plotter.interpolate(TimeExp[KKK],m_wf[KKK][0],m_wf[KKK][1]))


T_lim_min = [170, 150, 140, 280,  200, 150]
T_lim_max = [300, 400, 300, 400,  320, 350]

NUMBER = [r'I',r'II',r'III',r'IV',r'V']


print max(DNIExp[3])
print max(m_wfExp[3])
print '--'
print max(DNIExp[4])
print max(m_wfExp[4])
print '--'
print max(DNIExp[0])
print max(m_wfExp[0])

# Over the day different nodes
# Simulation time

NUMBER = [r'$I$',r'$II$',r'$III$',r'$IV$',r'$V$']
ZoomTimeShift = [9000,2,3,4000,8500]
Title = [r'a) MFE: Oil mass flow change',r'b) TE: Oil inlet temperature change','c) SBE: solar beam change']
##  ---------   PLOTS ZOOM OF THE PREVIOUS DAYS TO SHOW THE DIFFERENT EXPERIMENTAL RESULTS ------------
fig = plt.figure()
fig.set_size_inches(8.27,11.69)
fig.subplots_adjust(hspace=.21,right=.52,left=.06,top = .97)
ax1 = fig.add_subplot(3,1,1)
lns1 = ax1.plot(TimeExp[3]-StartModTime[3]-ZoomTimeShift[3],T_suExp[3],'o',markeredgewidth=0.05,markeredgecolor='dodgerblue',markerfacecolor='None',label=r'T$_\mathrm{su}$ exp')  #alpha=.2,color='dodgerblue',markeredgecolor='none',label=r'T$_\mathrm{su}$ exp')
lns2 = ax1.plot(TimeExp[3]-StartModTime[3]-ZoomTimeShift[3],T_exExp[3],'o',markeredgewidth=0.05,markeredgecolor='red',markerfacecolor='None',label=r'T$_\mathrm{ex}$ exp')
lns3 = ax1.plot(T_ex[3][0]-StartModTime[3]-ZoomTimeShift[3],T_ex[3][1]-273.15,color='black',label=r'T$_\mathrm{ex}$ sim')
lns4 = ax1.plot(TimeExp[3]-StartModTime[3]-ZoomTimeShift[3],T_exExp[3]-1,'--',color='black',label=r'$\pm$~1~K')
ax1.plot(TimeExp[3]-StartModTime[3]-ZoomTimeShift[3],T_exExp[3]+1,'--',color='black')
plt.ylabel(r'Temperature [$^{\circ}C$]')
# Plot DNI
ax2 = ax1.twinx()# for an array of equal element
lns5 = ax2.plot(TimeExp[3]-StartModTime[3]-ZoomTimeShift[3],DNIExp[3]/max(DNIExp[3]),'o', markeredgewidth=0.05,markerfacecolor='None',markeredgecolor='darkgoldenrod',label = r'DNI')
lns6 = ax2.plot(TimeExp[3]-StartModTime[3]-ZoomTimeShift[3],m_wfExp[3]/max(m_wfExp[3]),'o', markeredgewidth=0.05,markerfacecolor='None',markeredgecolor='black',label = r'$\dot{m}_\mathrm{wf}$')
plt.ylabel(r'DNI - $\dot{m}_\mathrm{wf}$ [-]')
lns = lns1+lns2+lns3+lns4+lns5+lns6
labs = [l.get_label() for l in lns]
ax1.legend(lns,labs,loc = 'upper right',fancybox='True',shadow='True',labelspacing=0.1,fontsize=10,numpoints=1)
ax1.set_ylim(T_lim_min[3],T_lim_max[3])
bbox_props = dict(boxstyle = "round",fc = 'white',ec ="0",alpha = 0.5)
ax1.text(200,T_lim_max[3]-20, r'Day:'+NUMBER[3], weight = 'bold',ha ="center", va ="center",style='italic',
bbox=bbox_props)
ax2.set_ylim(0,1.6)
ax1.set_xlim(0,2000)
ax1.set_title(Title[0])
#bbox_props = dict(boxstyle = "circle",fc = 'white',ec ="0",alpha = 0.5)
#ax1.text(100,T_lim_max[3]-5, NUMBER[3], weight = 'bold',ha ="center", va ="center",style='italic',
#bbox=bbox_props)
plt.tight_layout()
plt.grid()
#


ax1 = fig.add_subplot(3,1,2)
lns1 = ax1.plot(TimeExp[4]-StartModTime[4]-ZoomTimeShift[4],T_suExp[4],'o',markeredgewidth=0.05,markeredgecolor='dodgerblue',markerfacecolor='None',label=r'T$_\mathrm{su}$ exp')  #alpha=.2,color='dodgerblue',markeredgecolor='none',label=r'T$_\mathrm{su}$ exp')
lns2 = ax1.plot(TimeExp[4]-StartModTime[4]-ZoomTimeShift[4],T_exExp[4],'o',markeredgewidth=0.05,markeredgecolor='red',markerfacecolor='None',label=r'T$_\mathrm{ex}$ exp')
lns3 = ax1.plot(T_ex[4][0]-StartModTime[4]-ZoomTimeShift[4],T_ex[4][1]-273.15,color='black',label=r'T$_\mathrm{ex}$ sim')
lns4 = ax1.plot(TimeExp[4]-StartModTime[4]-ZoomTimeShift[4],T_exExp[4]-1,'--',color='black',label=r'$\pm$~1~K')
ax1.plot(TimeExp[4]-StartModTime[4]-ZoomTimeShift[4],T_exExp[4]+1,'--',color='black')
plt.ylabel(r'Temperature [$^{\circ}C$]')
# Plot DNI
ax2 = ax1.twinx()# for an array of equal element
lns5 = ax2.plot(TimeExp[4]-StartModTime[4]-ZoomTimeShift[4],DNIExp[4]/max(DNIExp[4]),'o', markeredgewidth=0.05,markerfacecolor='None',markeredgecolor='darkgoldenrod',label = r'DNI')
lns6 = ax2.plot(TimeExp[4]-StartModTime[4]-ZoomTimeShift[4],m_wfExp[4]/max(m_wfExp[4]),'o', markeredgewidth=0.05,markerfacecolor='None',markeredgecolor='black',label = r'$\dot{m}_\mathrm{wf}$')
plt.ylabel(r'DNI - $\dot{m}_\mathrm{wf}$ [-]')
lns = lns1+lns2+lns3+lns4+lns5+lns6
labs = [l.get_label() for l in lns]
ax1.legend(lns,labs,loc = 'upper right',fancybox='True',shadow='True',labelspacing=0.1,fontsize=10,numpoints=1)
ax1.set_ylim(T_lim_min[4],T_lim_max[4])
ax2.set_ylim(0,1.6)
ax1.set_xlim(0,2000)
ax1.set_title(Title[1])
bbox_props = dict(boxstyle = "round",fc = 'white',ec ="0",alpha = 0.5)
ax1.text(200,T_lim_max[4]-20, r'Day:'+NUMBER[4], weight = 'bold',ha ="center", va ="center",style='italic',
bbox=bbox_props)
plt.tight_layout()
plt.grid()


ax1 = fig.add_subplot(3,1,3)
lns1 = ax1.plot(TimeExp[0]-StartModTime[0]-ZoomTimeShift[0],T_suExp[0],'o',markeredgewidth=0.05,markeredgecolor='dodgerblue',markerfacecolor='None',label=r'T$_\mathrm{su}$ exp')  #alpha=.2,color='dodgerblue',markeredgecolor='none',label=r'T$_\mathrm{su}$ exp')
lns2 = ax1.plot(TimeExp[0]-StartModTime[0]-ZoomTimeShift[0],T_exExp[0],'o',markeredgewidth=0.05,markeredgecolor='red',markerfacecolor='None',label=r'T$_\mathrm{ex}$ exp')
lns3 = ax1.plot(T_ex[0][0]-StartModTime[0]-ZoomTimeShift[0],T_ex[0][1]-273.15,color='black',label=r'T$_\mathrm{ex}$ sim')
lns4 = ax1.plot(TimeExp[0]-StartModTime[0]-ZoomTimeShift[0],T_exExp[0]-1,'--',color='black',label=r'$\pm$~1~K')
ax1.plot(TimeExp[0]-StartModTime[0]-ZoomTimeShift[0],T_exExp[0]+1,'--',color='black')
plt.ylabel(r'Temperature [$^{\circ}C$]')
# Plot DNI
ax2 = ax1.twinx()# for an array of equal element
lns5 = ax2.plot(TimeExp[0]-StartModTime[0]-ZoomTimeShift[0],DNIExp[0]/max(DNIExp[0]),'o', markeredgewidth=0.05,markerfacecolor='None',markeredgecolor='darkgoldenrod',label = r'DNI')
lns6 = ax2.plot(TimeExp[0]-StartModTime[0]-ZoomTimeShift[0],m_wfExp[0]/max(m_wfExp[0]),'o', markeredgewidth=0.05,markerfacecolor='None',markeredgecolor='black',label = r'$\dot{m}_\mathrm{wf}$')
plt.ylabel(r'DNI - $\dot{m}_\mathrm{wf}$ [-]')
lns = lns1+lns2+lns3+lns4+lns5+lns6
labs = [l.get_label() for l in lns]
ax1.legend(lns,labs,loc = 'upper right',fancybox='True',shadow='True',labelspacing=0.1,fontsize=10,numpoints=1)
ax1.set_ylim(T_lim_min[0],T_lim_max[0])
ax2.set_ylim(0,1.6)
ax1.set_xlim(0,2000)
ax1.set_title(Title[2])
bbox_props = dict(boxstyle = "round",fc = 'white',ec ="0",alpha = 0.5)
ax1.text(200,T_lim_max[0]-20, r'Day:'+NUMBER[0], weight = 'bold',ha ="center", va ="center",style='italic',
bbox=bbox_props)
ax1.set_xlabel(r'Time[sec]')
plt.tight_layout()
plt.grid()

if SAVE_FIGURES == 'True': fig.savefig(DIRECTORY_FIGURES+'/'+'_MassFlowChange_Unc.pdf',dpi=1000)
plt.show()