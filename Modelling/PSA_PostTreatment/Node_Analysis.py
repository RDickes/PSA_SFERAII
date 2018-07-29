'''
This script get the Modelica simulation files using the buildingspy library.
And plot the results for the different nodes...

'''
from __future__ import division

import sys,os
import buildingspy.simulate.Simulator as si
from buildingspy.io.outputfile import Reader
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
FolderFig = '2016.06.29/'

# ---------    LOAD EXPERIMENTAL RESULTS    ---------
DirectoryExpData = 'C:\Users\susanna\Documents\GitHub\PSA_SFERAII\ExperimentalData/'
file = '2016_06_29_DATA.csv'
# data as dataframe
df = pd.read_csv(DirectoryExpData+file)



# ---------------------  MODELICA MODEL ---------------
# Define directory to store results
ResultDirectory = 'C:\Users\susanna\Documents\GitHub\PSA_SFERAII\Modelling\ModelicaResults/'
# Define simulation results name
FileSimulation = 'PTTL_SF_N'#'PTTL_SF_'

# SET TIME FOR SIMULATION
StartModTime = [46000]
StopModTime = [58000]
Delta_Time = StopModTime[0] - StartModTime[0]

SIMULATE = 'False'
# SET TIME FOR PLOT
TimeNoInit = 500                        # Set a value to phase out the initialization line
StartPlotTime = [StartModTime[0]+TimeNoInit]
StopPlotTime = [StopModTime[0]]

Nodes = [2,5,10,20,50]

if SIMULATE =='True':
    # LOAD MODELICA MODEL
    model ='package_PSA_SFERAII.Simulations.PTTL_SF_basic'
    s = si.Simulator(model,'dymola',outputDirectory= ResultDirectory,packagePath='C:\Users\susanna\Documents\GitHub\PSA_SFERAII\Modelling/package_PSA_SFERAII')
    s.setNumberOfIntervals(500)
    s.setSolver('Dassl')
    for kk in range(len(Nodes)):
        s.addParameters({'EuroTrough.N':Nodes[kk]})
        s.addParameters({'EuroTrough.eps6': 0.943771684737})
        s.setResultFile(FileSimulation+str(Nodes[kk]))
        s.setStartTime(StartModTime[0])
        s.setStopTime(StopModTime[0])
        s.printModelAndTime()
        s.simulate()
else:
    s = 'None'

# IMPORT MODELICA RESULT USING MODELICA RES
File_Sim = []
for kk in range(len(Nodes)):
    File_Sim.append(FileSimulation+str(Nodes[kk]))
TimeSim = []
NN =[]
Tpt_su =[]
Tpt_ex = []
Delta_T = []
DNI = []
m_wf = []
T_amb = []
for kk in range(len(File_Sim)):
    Data_Sim = ResultDirectory + File_Sim[kk]+'.mat'
    sim = SimRes(Data_Sim)
    initialTime = (sim['Time'].IV())
    FinalTime=(sim['Time'].FV())
    TimeSim.append(((sim['Time'].values(t=(StartPlotTime[0],StopPlotTime[0]))-(StartPlotTime[0]))))
    NN.append(sim['EuroTrough.N'].values(t=(StartPlotTime[0])))
    Tpt_su.append((sim['SensTsu.fluidState.T'].values(t=(StartPlotTime[0],StopPlotTime[0])))-273.15)
    Tpt_ex.append((sim['SensTex.fluidState.T'].values(t=(StartPlotTime[0],StopPlotTime[0])))-273.15)
    Delta_T.append((sim['DeltaT.Delta'].values(t=(StartModTime[0]+600,StartModTime[0]+800))))
    DNI.append((sim['DNI.y'].values(t=(StartPlotTime[0],StopPlotTime[0]))))
    m_wf.append((sim['m_dot_htf.y'].values(t=(StartPlotTime[0],StopPlotTime[0]))))
    T_amb.append((sim['T_amb.y'].values(t=(StartPlotTime[0],StopPlotTime[0])))-273.15)



# ORGANIZE VARIABLES IN LIST
#VariableExp_I = [Time, T_eva_su_X, p_exp_su_X, p_exp_ex_X, P_exp_X,  m_orc_X,T_exp_su_X]
VariableSim = [TimeSim, Tpt_su,  Tpt_ex,   Delta_T,  DNI,  m_wf, T_amb]
# FOR LOOP TO CHANGE NUMBER OF NODE 2-5-10-20-30

# PLOT THE DATA
fig = plt.figure()
fig.set_size_inches(8.27,4.15)
fig.subplots_adjust(hspace=.21,right=.52,left=.06,top = .97)
ax1 = fig.add_subplot(1,1,1)
ax1.plot(df['Seconds']-StartPlotTime[0],df['TA066'],label=r'T$_\mathrm{pt,ex}$ exp',color='red',linestyle='-.')
ax1.plot(df['Seconds']-StartPlotTime[0],df['TA060'],label=r'T$_\mathrm{pt,su}$ exp',color='dodgerblue',linestyle='-.')
for kk in range(len(TimeSim)):
    ax1.plot(TimeSim[kk],Tpt_ex[kk],label=r'N '+str(int(NN[kk])))
#ax1.set_title(r'a) MassFlow',fontsize=fontsizes[0])
ax1.legend(loc = 'best',fancybox='True',shadow='True',labelspacing=0.1,fontsize=13)
plt.xlim(0,TimeSim[0][-1])
plt.ylim(190,300)
plt.ylabel(r'Temperature [$^{\circ}$C]')
plt.xlabel(r'Time [$^{sec}$]')
plt.tight_layout()
if SAVE_FIGURES == 'True': fig.savefig(DIRECTORY_FIGURES+FolderFig+'_Node.pdf')
plt.show()

