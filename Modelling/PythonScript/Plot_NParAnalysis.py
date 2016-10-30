'''
This script get the Modelica simulation files using the buildingspy library.
And plot it for the Node analysis performed on the forth day of experimental data
20160704_NodeAnalysis
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
SAVE_FIGURES = 'True'
DIRECTORY_FIGURES = 'C:\Users\susanna\Documents\GitHub\PSA_SFERAII\Modelling\ResultsFigure/'

# --------------    Select the day of experiments  -----------------

StartModTime = [43000]
StopModTime =  [57000]

FileSimulation = 'PTTL_SF_'
Nodes = ["N2",'N5','N10','N20','N50']

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

E_htf = []
vector = [ 0, 1]
for KKK in range(len(Nodes)):
    resultFile=('C:\Users\susanna\Documents\GitHub\PSA_SFERAII\Modelling\ModelicaResults\20160704_NodeAnalysis/'+FileSimulation+Nodes[KKK]+".mat")
    r=Reader(resultFile, "dymola")
    DNI.append(r.values('DNI.y'))
    NN.append(r.values('EuroTrough.N'))
    T_su.append(r.values('SensTsu.fluidState.T'))
    T_ex.append(r.values('SensTex.fluidState.T'))
    TimeSim.append(T_ex[KKK][0])
    T_ex_exp.append(r.values('t_htf_ex.y'))
    Delta_T.append(r.values('DeltaT.Delta'))
    m_wf.append(r.values('m_dot_htf.y'))
    E_htf.append(r.integral('EuroTrough.Summary.Q_htf_tot'))

    TimeExp.append(np.linspace(StartModTime[0],StopModTime[0],num=int((StopModTime[0]-StartModTime[0])/10)))

    # Interpolate the Modelica results over the Time vector
    DNIExp.append(Plotter.interpolate(TimeExp[KKK], DNI[KKK][0], DNI[KKK][1]))
    T_suExp.append(Plotter.interpolate(TimeExp[KKK],T_su[KKK][0],T_su[KKK][1]-273.15))
    T_exExp.append(Plotter.interpolate(TimeExp[KKK],T_ex_exp[KKK][0],T_ex_exp[KKK][1]-273.15))
    m_wfExp.append(Plotter.interpolate(TimeExp[KKK],m_wf[KKK][0],m_wf[KKK][1]))

    #r.integral('preHea.port.Q_flow')


# Compute percentage difference between total energy absorbed with the different CVs taking the 50 CVs as the reference

PDTE = []
for aaa in range(len(Nodes)):
    PDTE.append(abs(E_htf[aaa]-E_htf[4])/E_htf[4]*100)

print PDTE


# Set up computational time vector for bar plots
TimeVector = [10.3,13.3,24.4,32.4,60.4]
RealTime = (StopModTime[0]-StartModTime[0])

PerTVec = []
for xx in range(len(TimeVector)):
    PerTVec.append(TimeVector[xx]/RealTime*100)

NumOfBars = len(TimeVector)   # Number of bars

position = np.arange(NumOfBars)+0.5  # The bar center on the y axis
WidthBar = 0.5






NodeNum = [' CVs:2',' CVs:5',' CVs:10',' CVs:20',' CVs:50']
colors = ['dodgerblue','magenta','green','black','yellow']


fig = plt.figure()
fig.set_size_inches(11.69,8.27)
fig.subplots_adjust(hspace=.21,right=.52,left=.06,top = .97)
ax1 = fig.add_subplot(2,1,1)
lns1 = ax1.plot(TimeExp[0]-StartModTime[0]-500,T_suExp[0],'o',markeredgewidth=0.05,markeredgecolor='dodgerblue',markerfacecolor='None',label=r'T$_\mathrm{su}$ exp')
lns2 = ax1.plot(TimeExp[0]-StartModTime[0]-500,T_exExp[0],'o',markeredgewidth=0.05,markeredgecolor='red',markerfacecolor='None',label=r'T$_\mathrm{ex}$ exp')
# Plot DNI and m dot
ax2 = plt.twinx()# for an array of equal element
lns3 = ax2.plot(TimeExp[0]-StartModTime[0]-500,DNIExp[0]/max(DNIExp[0]),'o', markeredgewidth=0.05,markerfacecolor='None',markeredgecolor='darkgoldenrod',label = r'DNI')
lns4 = ax2.plot(TimeExp[0]-StartModTime[0]-500,m_wfExp[0]/max(m_wfExp[0]),'o', markeredgewidth=0.05,markerfacecolor='None',markeredgecolor='black',label = r'$\dot{m}_\mathrm{wf}$')
#Plot simulation results
lns5 = ax1.plot(T_ex[0][0]-StartModTime[0]-500,T_ex[0][1]-273.15,color=colors[0],label=r'T$_\mathrm{ex}$ sim'+NodeNum[0])
lns6 = ax1.plot(T_ex[1][0]-StartModTime[0]-500,T_ex[1][1]-273.15,color=colors[1],label=r'T$_\mathrm{ex}$ sim'+NodeNum[1])
lns7 = ax1.plot(T_ex[2][0]-StartModTime[0]-500,T_ex[2][1]-273.15,color=colors[2],label=r'T$_\mathrm{ex}$ sim'+NodeNum[2])
lns8 = ax1.plot(T_ex[3][0]-StartModTime[0]-500,T_ex[3][1]-273.15,color=colors[3],label=r'T$_\mathrm{ex}$ sim'+NodeNum[3])
lns9 = ax1.plot(T_ex[4][0]-StartModTime[0]-500,T_ex[4][1]-273.15,color=colors[4],label=r'T$_\mathrm{ex}$ sim'+NodeNum[4])
# Set labels
ax1.set_ylabel(r'Temperature [$^{\circ}C$]')
ax2.set_ylabel(r'DNI - $\dot{m}_\mathrm{wf}$ [-]')
ax1.set_xlabel(r'Time [sec]')
#Legend
lns = lns1+lns2+lns3+lns4+lns5+lns6+lns7+lns8+lns9
labs = [l.get_label() for l in lns]
plt.legend(lns,labs,loc = 'upper right',fancybox='True',shadow='True',labelspacing=0.1,fontsize=10,numpoints=1)
# Text
bbox_props = dict(boxstyle = "round",fc = 'white',ec ="0",alpha = 0.5)
ax1.text(1500,400-15, r'Day:$IV$', weight = 'bold',ha ="center", va ="center",style='italic',
bbox=bbox_props)
# Limits
ax1.set_ylim(280,400)
ax2.set_ylim(0,1.6)
ax1.set_xlim(0,StopModTime[0]-StartModTime[0]-500)
# Set title
ax1.set_title(r'a) SF CVs parametric analysis',fontsize=12)
plt.grid()
plt.tight_layout()
#
#
#
#
ax1 = fig.add_subplot(2,1,2)
ax1.barh(position,PerTVec,align ='center',color='b')
ax1.set_xlabel(r'PCE [$\%$]')
ax1.set_yticks(position)
ax1.set_yticklabels((r'2',r'5',r'10',r'20',r'50'))
ax1.set_ylabel(r'$\#$ of CVs in SF FV models')
plt.grid(True)
plt.tight_layout()
ax1.set_title(r'b) Computational effort by discretized volumes',fontsize=12)
fig.savefig(DIRECTORY_FIGURES+'/'+'_SF_NodesParaALL.pdf',dpi=1000)
plt.show()