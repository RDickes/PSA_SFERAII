from __future__ import division
__author__ = 'Admin'
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
    plt.imshow(-A.T, cmap='inferno', aspect='auto',
               origin='lower', extent=(min(x), max(x), ymin, ymax))


# --------------    SAVE GRAPHS   -----------------------------
SAVE_FIGURES = 'False'
DIRECTORY_FIGURES = 'C:\Users\susanna\Documents\GitHub\PSA_SFERAII\Modelling\ResultsFigure/'

# --------------    Select the day of experiments  -----------------
Days = ['20160629', '20160630',  '20160701',  '20160704',  '20160705', '20160706']
KKK = 5
# --------------------  SET FLAG  --------------
OPTIMIZE = 'False'   # To true to optimize eps6
SIMULATE = 'True'    # To true to simulate the modelica data
                     # Both to false to plot the results !!!! MAKE SURE THE SIMULATION FILE EXIST IN THE DEFINED DIRECTORY!!!
# --------------     LOAD EXPERIMENTAL RESULTS  -----------------
DirectoryExpData = 'C:\Users\susanna\Documents\GitHub\PSA_SFERAII\ExperimentalData/'
file = ['2016_06_29_DATA.csv','2016_06_30_DATA.csv','2016_07_01_DATA.csv','2016_07_04_DATA.csv','2016_07_05_DATA.csv','2016_07_06_DATA.csv']
# data as dataframe
df = pd.read_csv(DirectoryExpData+file[KKK])

# --------------------- DIRECTORY TO STORE  MODELICA MODEL ---------------
StoreModResult = 'C:\Users\susanna\Documents\GitHub\PSA_SFERAII\Modelling\ModelicaResults/'+Days[KKK]
# Define simulation results name
FileSimulation = 'PTTL_SF_'

# ----------- DEFINE SIMULATION TIME ------------
             # 2016.06.29   #2016.06.30   #2016.07.01  #2016.07.04    #2016.07.05    #2016.07.06
StartModTime = [46000,       40000,        37000,      43000,          43000,     41000]
StopModTime =  [58000,       54000,        53000,      57000,          56000,     50000]
Delta_Time = StopModTime[KKK] - StartModTime[KKK]

#----------- LOAD AND TRANSLATE MODELICA MODEL           ------------
model ='package_PSA_SFERAII.Simulations.PTTL_SF_basic_'+Days[KKK]


PACKAGE = os.path.join('C:\Users\susanna\Documents\GitHub\PSA_SFERAII\Modelling/', 'package_PSA_SFERAII')
s = si.Simulator(model,'dymola',outputDirectory= StoreModResult,packagePath=PACKAGE)

# Set Time for plot
TimeNoInit = 500                            # Set a value to phase out the initialization line
StartPlotTime = [StartModTime[KKK]+TimeNoInit]
StopPlotTime = [StopModTime[KKK]]


if OPTIMIZE == 'True':
    s.setNumberOfIntervals(500)
    s.setSolver('Dassl')
    s.translate()
    def optimizeCoeff(In):
        eps = In[0]
        print 'Current eps:', eps
        s.addParameters({'EuroTrough.eps6': eps})
        s.setResultFile(FileSimulation)
        s.setStartTime(StartModTime[0])
        s.setStopTime(StopModTime[0])
        s.printModelAndTime()
        #s.showGUI(show=True)
        #s.exitSimulator(exitAfterSimulation=False)
        s.simulate_translated()
        Data_Sim = StoreModResult +'/'+ FileSimulation+'.mat'
        sim = SimRes(Data_Sim)
        #initialTime = (sim['Time'].IV())
        #FinalTime=(sim['Time'].FV())
        #TimeSim=(((sim['Time'].values(t=(StartModTime[0],StopModTime[0]))-StartModTime[0])))
        #NN = (sim['EuroTrough.N'].values(t=(StartModTime[0])))
        #Tpt_su = ((sim['SensTsu.fluidState.T'].values(t=(StartModTime[0],StopModTime[0])))-273.15)
        #Tpt_ex = ((sim['SensTex.fluidState.T'].values(t=(StartModTime[0],StopModTime[0])))-273.15)
        Delta_T = ((sim['DeltaT.Delta'].values(t=(StartModTime[0]+600,StartModTime[0]+660))))
        print 'Residue:'+str(sum(Delta_T**2))
        return Delta_T

    opt,cov,infodict,mesg,ier  = optimize.leastsq(optimizeCoeff,(0.9),ftol=10,full_output=True)

    perr = np.sqrt(np.diag(cov))

    print perr

elif SIMULATE == 'True':
    # PLOT RESULTS
    s.addParameters({'EuroTrough.eps6': 0.943771684737})
    s.setResultFile(FileSimulation)
    s.setStartTime(StartModTime[KKK])
    s.setStopTime(StopModTime[KKK])
    s.setNumberOfIntervals(500)
    s.setSolver('Dassl')
    s.printModelAndTime()
    #s.showGUI(show=True)
    #s.exitSimulator(exitAfterSimulation=False)
    s.simulate()
    #resultFile=(ResultDirectory+'/'+FileSimulation+".mat")
    #r=Reader(resultFile, "dymola")
    #print r


# --------------------    IMPORT MODELICA RESULT USING MODELICA RES  -----------------------------
Data_Sim = StoreModResult +'/'+ FileSimulation+'.mat'
sim = SimRes(Data_Sim)
initialTime = (sim['Time'].IV())
FinalTime=(sim['Time'].FV())
TimeSim=(((sim['Time'].values(t=(StartPlotTime[0],StopPlotTime[0]))-(StartPlotTime[0]))))
NN = (sim['EuroTrough.N'].values(t=(StartPlotTime[0])))
Tpt_su = ((sim['SensTsu.fluidState.T'].values(t=(StartPlotTime[0],StopPlotTime[0])))-273.15)
Tpt_ex = ((sim['SensTex.fluidState.T'].values(t=(StartPlotTime[0],StopPlotTime[0])))-273.15)
Delta_T = ((sim['DeltaT.Delta'].values(t=(StartModTime[KKK]+600,StartModTime[KKK]+800))))
DNI = ((sim['DNI.y'].values(t=(StartPlotTime[0],StopPlotTime[0]))))
m_wf = ((sim['m_dot_htf.y'].values(t=(StartPlotTime[0],StopPlotTime[0]))))
T_amb = ((sim['T_amb.y'].values(t=(StartPlotTime[0],StopPlotTime[0])))-273.15)


Sigma_T=np.empty(len(df['TA066'].values)); Sigma_T.fill(0.9) # for an array of equal element



# ----------     PLOTS RESULTS  -------------------------
# TRITTICO DNI, m_dot, Temperatures
fig = plt.figure()
fig.set_size_inches(8.27,11.69)
fig.subplots_adjust(hspace=.21,right=.52,left=.06,top = .97)
ax1 = fig.add_subplot(3,1,1)
ax1.set_title(r'a) DNI')
ax1.plot(TimeSim,DNI,label='Exp Data')
plt.ylabel(r'DNI [W m$^{-2}$]')
ax1.legend(loc = 'best',fancybox='True',shadow='True',labelspacing=0.1)
plt.xlim(0,TimeSim[-1])
plt.tight_layout()
plt.grid()
#
ax1 = fig.add_subplot(3,1,2)
ax1.set_title(r'b) Htf mass flow')
ax1.plot(TimeSim,m_wf,label='Exp Data')
ax1.legend(loc = 'best',fancybox='True',shadow='True',labelspacing=0.1)
plt.ylabel(r'Mass flow rate [kg s$^{-1}$]')
plt.xlim(0,TimeSim[-1])
plt.tight_layout()
plt.grid()
#
ax1 = fig.add_subplot(3,1,3)
ax1.set_title(r'c) Htf Temperatures')
ax1.plot(TimeSim,Tpt_ex,label='T$_\mathrm{pt,ex}$ sim',color='orangered')
ax1.plot(df['Seconds']-StartPlotTime[0],df['TA066'],label=r'T$_\mathrm{pt,ex}$ exp',color='red',linestyle='-.')
ax1.plot(df['Seconds']-StartPlotTime[0],df['TA060'],label=r'T$_\mathrm{pt,su}$ exp',color='dodgerblue',linestyle='-.')
#pdense(ax1,df['Seconds'].values-StartPlotTime[0], df['TA066'].values, Sigma_T, M=1000)
ax1.legend(loc = 'best',fancybox='True',shadow='True',labelspacing=0.1)
plt.ylabel(r'Temperature [$^{\circ}$C]')
plt.xlim(0,TimeSim[-1])
ax1.set_ylim(180,350)
plt.tight_layout()
plt.grid()
if SAVE_FIGURES == 'True': fig.savefig(DIRECTORY_FIGURES+Days[KKK]+'/'+'_Trittico.pdf')

# PLOT Temperatures and DNI double Y
fig = plt.figure()
fig.set_size_inches(4.15,4.15)
fig.subplots_adjust(hspace=.21,right=.52,left=.06,top = .97)
ax1 = fig.add_subplot(1,1,1)
ax1.plot(TimeSim,Tpt_ex,label='T$_\mathrm{pt,ex}$ sim',color='orangered')
ax1.plot(df['Seconds']-StartPlotTime[0],df['TA066'],label=r'T$_\mathrm{pt,ex}$ exp',color='red',linestyle='-.')
ax1.plot(df['Seconds']-StartPlotTime[0],df['TA060'],label=r'T$_\mathrm{pt,su}$ exp',color='dodgerblue',linestyle='-.')
#ax1.set_title(r' EuroThrough Validation')
ax1.legend(loc = 'best',fancybox='True',shadow='True',labelspacing=0.1)
# Plot DNI
ax2 = ax1.twinx()
lns6 = ax2.plot(df['Seconds'].values-StartPlotTime[0],df['IA028'].values,'--', label = 'DNI', color = 'darkgoldenrod')
plt.ylabel('DNI',color = 'darkgoldenrod')
for xx in ax2.get_yticklabels():
    xx.set_color('darkgoldenrod')
plt.tight_layout()
ax2.set_ylim(600,1000)
ax1.set_ylim(180,350)
plt.xlim(0,TimeSim[-1])
plt.grid()
plt.tight_layout()
if SAVE_FIGURES == 'True': fig.savefig(DIRECTORY_FIGURES+Days[KKK]+'/'+'_T_DNI.pdf')
#
#fig = plt.figure()
#fig.set_size_inches(4.15,4.15)
#fig.subplots_adjust(hspace=.21,right=.52,left=.06,top = .97)
#ax1 = fig.add_subplot(1,1,1)
#ax1.plot(Delta_T,label='Exp Data')
#ax1.set_title(r'a) Delta T')
#ax1.legend(loc = 'best',fancybox='True',shadow='True',labelspacing=0.1)
#plt.tight_layout()
#plt.grid()
plt.show()