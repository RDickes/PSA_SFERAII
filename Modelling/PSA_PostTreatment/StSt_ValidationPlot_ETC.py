

'''
This script get the Modelica simulation files using the buildingspy library.
This allows to interpolate over the results to get the boundary conditions,i.e., DNI, T_su, T_ex as if they came from the experimental data
for the SF EuroThrough collectors and plot the steady state points
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
DIRECTORY_FIGURES = '/Users/adriano/Dropbox/GitHub/PSA_SFERAII/Modelling/ResultsFigure/'


# --------------    Define vectors for the different day of experiments  -----------------
Days = ['20160629', '20160630',  '20160701',  '20160704',  '20160705', '20160706']
StartModTime = [46000,       40000,        37000,      43000,          43000,     41000]
StopModTime =  [58000,       54000,        53000,      57000,          56000,     50000]
# Define simulation file
FileSimulation = 'PTTL_SF_'


# Define name sheet of excel with steady-state value
sheetnames= ['Day I 20160629','Day II 20160630','Day III 20160701','Day IV 20160704','Day V 20160705']

# Directory of excel with start and end time for steady-state data for each day
Directory = '/Users/adriano/Dropbox/GitHub/PSA_SFERAII/ExperimentalData/'
file = 'StSt_ETC.xls'
Path = Directory+file
DATA = {sheet: pd.read_excel(Path,sheetname=sheet) for sheet in sheetnames}



#print DATA['Day IV 20160704']


# Vector for simulation results
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


for KKK in range(len(StartModTime)):
    resultFile=('/Users/adriano/Dropbox/GitHub/PSA_SFERAII/Modelling/ModelicaResults/'+Days[KKK]+'/'+FileSimulation+".mat")
    r=Reader(resultFile, "dymola")
    DNI.append(r.values('DNI.y'))
    NN.append(r.values('EuroTrough.N'))
    T_su.append(r.values('SensTsu.fluidState.T'))
    T_ex.append(r.values('SensTex.fluidState.T'))
    TimeSim.append(T_ex[KKK][0])
    T_ex_exp.append(r.values('t_htf_ex.y'))
    Delta_T.append(r.values('DeltaT.Delta'))
    m_wf.append(r.values('m_dot_htf.y'))


    # Create a time vector which values every 10 seconds
    TimeExp.append(np.linspace(StartModTime[KKK],StopModTime[KKK],num=int((StopModTime[KKK]-StartModTime[KKK])/10)))

    # Interpolate the Modelica results over the Time vector
    DNIExp.append(Plotter.interpolate(TimeExp[KKK], DNI[KKK][0], DNI[KKK][1]))
    T_suExp.append(Plotter.interpolate(TimeExp[KKK],T_su[KKK][0],T_su[KKK][1]-273.15))
    T_exExp.append(Plotter.interpolate(TimeExp[KKK],T_ex_exp[KKK][0],T_ex_exp[KKK][1]-273.15))
    m_wfExp.append(Plotter.interpolate(TimeExp[KKK],m_wf[KKK][0],m_wf[KKK][1]))

##########################################################################
#      Steady-state values for model and experiments
##########################################################################
# Get index for Modelica model

#Define function to get to closer value of index
def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx


indexStart = []
indexStop = []
# Get starting and ending index
for kkk in range(len(sheetnames)):
    indexStart.append(find_closest(T_ex[kkk][0],StartModTime[kkk]+DATA[sheetnames[kkk]]['time start']))
    indexStop.append(find_closest(T_ex[kkk][0],StartModTime[kkk]+DATA[sheetnames[kkk]]['time stop']))

# Get the averaged value of Temperature
T_stst_model = []
for aa in range(len(sheetnames)):
    for bb in range(len(indexStart[aa])):
        T_stst_model.append(sum(T_ex[aa][1][indexStart[aa][bb]:indexStop[aa][bb]]) / len(T_ex[aa][1][indexStart[aa][bb]:indexStop[aa][bb]]) - 273.15)



# Get the index for Experimental data
T_stst_exp = []
for aa in range(len(sheetnames)):
    for bb in range(len(indexStart[aa])):
        T_stst_exp.append(np.mean(T_exExp[aa][int((DATA[sheetnames[aa]]['time start'][bb])/10):int((DATA[sheetnames[aa]]['time stop'][bb])/10)]))





########################################################################################################################
###########################                      PLOT   ----      Steady state validation                        ###########################
########################################################################################################################
FontSizes = [30]
markerSizes = [16]

fig = plt.figure()
fig.set_size_inches(4.5,4.5)
fig.subplots_adjust(hspace=.24,right=.99,left=.07,top = .96)
ax1 = fig.add_subplot(1,1,1)
ax1.plot(T_stst_exp,T_stst_model,linestyle="None",color = 'firebrick', marker = '^',label = r'$T_\mathrm{ex,ETC}$ [$^{\circ}$C]')
alpha_min = 150
alpha_max = 360
ax1.plot([alpha_min,alpha_max],[alpha_min,alpha_max],label = '$X =Y$', color='red')
ax1.plot([0,alpha_max],[0,alpha_max+3],label = '$\pm 3$K', color='black')
ax1.plot([0,alpha_max],[0,alpha_max-3], color='black')
ax1.set_ylabel(r'$T_\mathrm{ex,ETC,model} [^{\circ}C]$')
ax1.set_xlabel(r'$T_\mathrm{ex,ETC,exp} [^{\circ}C]$')
plt.grid()
ax1.legend(loc = 'lower right',numpoints = 1,fancybox='True',shadow='True',labelspacing=0.1)
plt.xlim(alpha_min,alpha_max)
plt.ylim(alpha_min,alpha_max)
plt.tight_layout()
if SAVE_FIGURES == 'True':fig.savefig('/Users/adriano/Dropbox/GitHub/PSA_SFERAII/Modelling/ResultsFigure/StSt_Validation.pdf')

plt.show()
