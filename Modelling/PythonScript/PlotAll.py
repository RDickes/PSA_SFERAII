'''
This script get the Modelica simulation files using the buildingspy library.
This allows to interpolate over the results to get the boundary conditions,i.e., DNI, T_su, T_ex as if they came from the experimental data
for the SF EuroThrough collectors and plot it for the whole day.
For now I am plotting in the same graph only the first 5 days of simulation

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


# FUNCTIONS
# Function to plot uncertainty as a shadow
def pdense(ax,x, y, sigma, M=1000):
    """ Plot probability density of y with known stddev sigma
    """
    assert len(x) == len(y) and len(x) == len(sigma)
    N = len(x)
    # TODO: better y ranging
    ymin, ymax = min(y - 2 * sigma), max(y + 2 * sigma)
    yy = np.linspace(ymin, ymax, M)
    a = [np.exp(-((Y - yy) / s) ** 2) / s for Y, s in zip(y, sigma)]
    A = np.array(a)
    A = A.reshape(N, M)
    plt.imshow(-A.T, cmap='gray', aspect='auto',
               origin='lower', extent=(min(x), max(x), ymin, ymax))

# --------------    SAVE GRAPHS   -----------------------------
SAVE_FIGURES = 'True'
DIRECTORY_FIGURES = 'C:\Users\susanna\Documents\GitHub\PSA_SFERAII\Modelling\ResultsFigure/'

# --------------    Select the day of experiments  -----------------
Days = ['20160629', '20160630',  '20160701',  '20160704',  '20160705', '20160706']
StartModTime = [46000,       40000,        37000,      43000,          43000,     41000]
StopModTime =  [58000,       54000,        53000,      57000,          56000,     50000]

FileSimulation = 'PTTL_SF_'


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
    resultFile=('C:\Users\susanna\Documents\GitHub\PSA_SFERAII\Modelling\ModelicaResults/'+Days[KKK]+'/'+FileSimulation+".mat")
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

    TimeExp.append(np.linspace(StartModTime[KKK],StopModTime[KKK],num=int((StopModTime[KKK]-StartModTime[KKK])/10)))

    # Interpolate the Modelica results over the Time vector
    DNIExp.append(Plotter.interpolate(TimeExp[KKK], DNI[KKK][0], DNI[KKK][1]))
    T_suExp.append(Plotter.interpolate(TimeExp[KKK],T_su[KKK][0],T_su[KKK][1]-273.15))
    T_exExp.append(Plotter.interpolate(TimeExp[KKK],T_ex_exp[KKK][0],T_ex_exp[KKK][1]-273.15))
    m_wfExp.append(Plotter.interpolate(TimeExp[KKK],m_wf[KKK][0],m_wf[KKK][1]))




T_lim_min = [170, 150, 140, 280,  200, 150]
T_lim_max = [300, 400, 300, 400,  320, 350]

NUMBER = [r'$I$',r'$II$',r'$III$',r'$IV$',r'$V$']
fig = plt.figure()
fig.set_size_inches(8.27,11.69)
fig.subplots_adjust(hspace=.21,right=.52,left=.06,top = .97)
for kk in range(len(TimeSim)-1):
    ax1 = fig.add_subplot(5,1,kk+1,)
    lns1 = ax1.plot(TimeExp[kk]-StartModTime[kk]-500,T_suExp[kk],'o',markeredgewidth=0.05,markeredgecolor='dodgerblue',markerfacecolor='None',label=r'T$_\mathrm{su}$ exp')  #alpha=.2,color='dodgerblue',markeredgecolor='none',label=r'T$_\mathrm{su}$ exp')
    lns2 = ax1.plot(TimeExp[kk]-StartModTime[kk]-500,T_exExp[kk],'o',markeredgewidth=0.05,markeredgecolor='red',markerfacecolor='None',label=r'T$_\mathrm{ex}$ exp')
    lns3 = ax1.plot(T_ex[kk][0]-StartModTime[kk]-500,T_ex[kk][1]-273.15,color='black',label=r'T$_\mathrm{ex}$ sim')
    if T_lim_min[kk]==170:ax1.plot([600,600],[T_lim_min[0],T_lim_max[0]],linestyle='--',color='black')
    if T_lim_min[kk]==170:ax1.plot([780,780],[T_lim_min[0],T_lim_max[0]],linestyle='--',color='black')
    plt.ylabel(r'Temperature [$^{\circ}C$]')
    # Plot DNI
    ax2 = ax1.twinx()# for an array of equal element
    lns4 = ax2.plot(TimeExp[kk]-StartModTime[kk]-500,DNIExp[kk]/max(DNIExp[kk]),'o', markeredgewidth=0.05,markerfacecolor='None',markeredgecolor='darkgoldenrod',label = r'DNI')
    lns5 = ax2.plot(TimeExp[kk]-StartModTime[kk]-500,m_wfExp[kk]/max(m_wfExp[kk]),'o', markeredgewidth=0.05,markerfacecolor='None',markeredgecolor='black',label = r'$\dot{m}_\mathrm{wf}$')
    plt.ylabel(r'DNI - $\dot{m}_\mathrm{wf}$ [-]')
    ax2.set_ylim(0,1.6)
    lns = lns1+lns2+lns3+lns4+lns5
    labs = [l.get_label() for l in lns]
    ax1.legend(lns,labs,loc = 'upper center',fancybox='True',shadow='True',labelspacing=0.1,fontsize=8,numpoints=1)
    ax1.set_ylim(T_lim_min[kk],T_lim_max[kk])
    bbox_props = dict(boxstyle = "round",fc = 'white',ec ="0",alpha = 0.5)
    ax1.text(1000,T_lim_max[kk]-18, r'Day:'+NUMBER[kk], weight = 'bold',ha ="center", va ="center",style='italic',
    bbox=bbox_props)
    plt.xlim(0,StopModTime[kk]-StartModTime[kk]-500)
    plt.grid()
    if T_lim_max[kk]==320: ax1.set_xlabel(r'Time[sec]')
    plt.tight_layout()
plt.tight_layout()
if SAVE_FIGURES == 'True':fig.savefig(DIRECTORY_FIGURES+'/'+'_FirstFiveDays.pdf',dpi=1000)
plt.show()







'''
# In this part of the script I load the Modelica results with the ModelicaRes library and I calculate the maximum and minimum
of the different values. I used this to build the table in Chapter 5 where I show the working conditions range during the experimental campaign

for KKK in range(len(Days)):
    # ----------- DEFINE SIMULATION TIME ------------
                 # 2016.06.29   #2016.06.30   #2016.07.01  #2016.07.04    #2016.07.05    #2016.07.06
    Delta_Time = StopModTime[KKK] - StartModTime[KKK]
    StoreModResult = 'C:/Users/Admin/Desktop/SFERA_II_Sim/'+Days[KKK]
    # Define simulation results name
    FileSimulation = 'PTTL_SF_'
    # Set Time for plot
    TimeNoInit = 500                            # Set a value to phase out the initialization line
    StartPlotTime = [StartModTime[KKK]+TimeNoInit]
    StopPlotTime = [StopModTime[KKK]]
    Data_Sim = StoreModResult +'/'+ FileSimulation+'.mat'
    sim = SimRes(Data_Sim)
    initialTime = (sim['Time'].IV())
    FinalTime=(sim['Time'].FV())
    TimeSim.append(((sim['Time'].values(t=(StartPlotTime[0],StopPlotTime[0]))-(StartPlotTime[0]))))
    NN.append(sim['EuroTrough.N'].values(t=(StartPlotTime[0])))
    Tpt_su.append((sim['SensTsu.fluidState.T'].values(t=(StartPlotTime[0],StopPlotTime[0])))-273.15)
    Tpt_ex.append((sim['SensTex.fluidState.T'].values(t=(StartPlotTime[0],StopPlotTime[0])))-273.15)
    Tpt_ex_exp.append((sim['t_htf_ex.y'].values(t=(StartPlotTime[0],StopPlotTime[0])))-273.15)
    Delta_T.append((sim['DeltaT.Delta'].values(t=(StartModTime[KKK]+600,StartModTime[KKK]+800))))
    DNI.append((sim['DNI.y'].values(t=(StartPlotTime[0],StopPlotTime[0]))))
    m_wf.append((sim['m_dot_htf.y'].values(t=(StartPlotTime[0],StopPlotTime[0]))))
    T_amb.append((sim['T_amb.y'].values(t=(StartPlotTime[0],StopPlotTime[0])))-273.15)
    v_wind.append((sim['v_wind.y'].values(t=(StartPlotTime[0],StopPlotTime[0]))))
T_lim_min = [170, 150, 140, 280,  200, 150]
T_lim_max = [300, 400, 300, 400,  320, 350]

M_max = []
M_min = []
T_su_MAX = []
T_su_MIN = []
T_ex_MAX = []
T_ex_MIN = []
p_MAX = []
p_MIN = []
DNI_MAX = []
DNI_MIN = []
T_amb_MAX = []
T_amb_MIN = []
v_MAX = []
v_MIN = []
for kk in range(len(T_lim_min)):
    M_max.append(max(m_wf[kk]))
    M_min.append(min(m_wf[kk]))
    T_su_MAX.append(max(Tpt_su[kk]))
    T_su_MIN.append(min(Tpt_su[kk]))
    T_ex_MAX.append(max(Tpt_ex_exp[kk]))
    T_ex_MIN.append(min(Tpt_ex_exp[kk]))
    DNI_MAX.append(max(DNI[kk]))
    DNI_MIN.append(min(DNI[kk]))
    T_amb_MAX.append(max(T_amb[kk]))
    T_amb_MIN.append(min(T_amb[kk]))
    v_MAX.append(max(v_wind[kk]))
    v_MIN.append(min(v_wind[kk]))

print '       --------'
print 'Mass flow max', max(M_max)
print 'Mass flow min', min(M_min)
print '       --------'
print 'T inlet max:', max(T_su_MAX)
print 'T inlet min:', min(T_su_MIN)
print  '         --------'
print 'T outlet max:', max(T_ex_MAX)
print 'T outlet min:', min(T_ex_MIN)
print  '         --------'
print 'DNI Max:', max(DNI_MAX)
print 'DNI Min:', min(DNI_MIN)
print  '         --------'
print 'T amb Max', max(T_amb_MAX)
print 'T amb Min', min(T_amb_MIN)
print '          --------'
print 'wind max', max(v_MAX)
print 'wind min', min(v_MIN)

#Sigma_T = []
#for kk in len(T_lim_min):
Sigma_T = []
for aa in range(len(T_lim_min)):
    Sigma_T.append(np.empty(len(Tpt_ex_exp[aa])))
    Sigma_T[aa].fill(0.9)

print len(Sigma_T)


# -------------------- PLOT RESULTS --------------------------------------------------------
fig = plt.figure()
fig.set_size_inches(8.27,11.69)
fig.subplots_adjust(hspace=.21,right=.52,left=.06,top = .97)

for kk in range(len(TimeSim)-1):
    ax1 = fig.add_subplot(5,1,kk+1)
    ax1.plot(TimeSim[kk],Tpt_su[kk],'o',color='dodgerblue',markeredgecolor='grey',label=r'T$_\mathrm{su}$ exp')
    ax1.plot(TimeSim[kk],Tpt_ex[kk],color='red',label=r'T$_\mathrm{ex}$ sim')
    ax1.plot(TimeSim[kk],Tpt_ex_exp[kk],'o',alpha=.5,color='black',markeredgecolor='grey',label=r'T$_\mathrm{ex}$ exp')
    ax1.legend(loc = 'best',fancybox='True',shadow='True',labelspacing=0.1,fontsize=12)
    plt.xlabel(r'Time[sec]')
    plt.ylabel(r'Temperature [$^{\circ}C$]')
    # Plot DNI
    plt.tight_layout()
    ax2 = ax1.twinx()# for an array of equal element
    ax2.plot(TimeSim[kk],DNI[kk]/max(DNI[kk]),'o', alpha=.5,markeredgecolor='grey',label = r'DNI', color = 'darkgoldenrod')
    ax2.plot(TimeSim[kk],m_wf[kk]/max(m_wf[kk]),'o', alpha=.5,markeredgecolor='grey',label = r'$\dot{m}_\mathrm{wf}$', color = 'black')
    plt.ylabel(r'DNI - $\dot{m}_\mathrm{wf}$ [-]')
    #for xx in ax2.get_yticklabels():
    #    xx.set_color('darkgoldenrod')
    plt.tight_layout()
    #ax2.set_ylim(600,1000)
    ax1.set_ylim(T_lim_min[kk],T_lim_max[kk])
    plt.xlim(0,TimeSim[kk][-1])
    plt.grid()
plt.tight_layout()
if SAVE_FIGURES == 'True': fig.savefig(DIRECTORY_FIGURES+'/'+'_FirstFiveDays.pdf')
plt.show()

'''