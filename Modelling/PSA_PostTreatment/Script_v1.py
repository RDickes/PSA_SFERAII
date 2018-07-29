from __future__ import division
__author__ = 'Admin'
from datetime import datetime
import re
import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI
import matplotlib.dates as md
import matplotlib.pyplot as plt
# Latex caption in the figures:
from modelicares import SimRes
from decimal import *

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


KKK = 0
DIRECTORY_FIGURES = 'C:/Users/Admin/Documents/GitHub/PSA_SFERAII/Modelling/Results/'
# --------------    Select the day of experiments  -----------------
Days = ['20160629', '20160630',  '20160701',  '20160704',  '20160705', '20160706']
Files = ['29-06-2016_II.txt', '30-06-2016_bis.txt', '01-07-2016_bis.txt', '04-07-2016_bis.txt', '05-07-2016_bis.txt', '06-07-2016_bis.txt']
StartTime = [46000,       40000,        37000,      43000,          43000,     41000]
StopTime =  [58000,       54000,        53000,      57000,          56000,     50000]
Directory ='C:/Users/Admin/Documents/GitHub/PSA_SFERAII/ExperimentalData/RAW_DATA/'


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
for kk in range(len(Days)):
    File = Files[kk]
    df = pd.read_table(Directory+'/'+File,low_memory=False)
    Time_h = df['Hora'].values
    # CONVERT hours to seconds
    Time_s_list = []
    for t in Time_h:
        times = map(int, re.split(r"[:,]",t))
        Time_s_list.append(times[0]*3600+times[1]*60+times[2])
    Time_s = np.asarray(Time_s_list)
    df["Seconds"] = Time_s

    df = df[df['Seconds']>StartTime[kk]]
    df = df[df['Seconds']<StopTime[kk]]

    M_max.append(max(df['FA032'].values))
    M_min.append(min(df['FA032'].values))
    T_su_MAX.append(max(df['TA060'].values))
    T_su_MIN.append(min(df['TA060'].values))
    #  df.loc[df['Value'].idxmax()]
    T_ex_MAX.append(max(df['TA066'].values))
    T_ex_MIN.append(min(df["TA066"].values))

    p_MAX.append(max(df['PA052'].values))
    p_MIN.append(min(df['PA052'].values))

    DNI_MAX.append(max(df['IA028'].values))
    DNI_MIN.append(min(df['IA028'].values))

    T_amb_MAX.append(max(df['TA029'].values))
    T_amb_MIN.append(min(df['TA029'].values))

    v_MAX.append(max(df['SA026'].values))
    v_MIN.append(min(df['SA026'].values))



print  '--------'
print 'Max mass flow rate:', max(M_max)
print 'Min mass flow rate:',  min(M_min)
print  '         --------'
print 'Max Oil T inlet:', max(T_su_MAX)
print 'Min Oil T inlet:', min(T_su_MIN)
print  '         --------'
print 'Max Oil T outlet:', max(T_ex_MAX)
print 'Min Oil T outlet:', min(T_ex_MIN)
print  '         --------'
print 'Max pressure:', max(p_MAX)
print 'Min pressure:', min(p_MIN)



File = Files[KKK]
df = pd.read_table(Directory+'/'+File,low_memory=False)

Day = Days[KKK]

# TIME IN hours 'h:m:s'
Time_h = df['Hora'].values
# CONVERT hours to seconds
Time_s_list = []
for t in Time_h:
    times = map(int, re.split(r"[:,]",t))
    Time_s_list.append(times[0]*3600+times[1]*60+times[2])
Time_s = np.asarray(Time_s_list)
df["Seconds"] = Time_s


# CONVERT String hours to datetime object
Time_hours = []
for kk in range(len(df['Hora'].values)):
    xx = datetime.strptime(df['Hora'].values[kk],"%H:%M:%S")
    Time_hours.append(xx)

# Temperature inlet and outlet
T_cooler_ex = df['TA012'].values
T_pt_su = df['TA060'].values
T_pt_ex = df['TA066'].values
T_avg = (T_pt_ex+ T_pt_su)/2

# Pressures
p_pt_su = df['PA052'].values
p_tank = df['PA021'].values

# Flows
V_wf_su = df['FA023'].values
m_wf_su = df['FA032'].values

# Ambient condition
DNI = df['IA028'].values
Tamb = df['TA029'].values
df['WD089'] = df['WD089'] /3.6  # km/h to m/s
df['WD088'] = df['WD088'] /3.6 # km/h to m/s


print df["Seconds"].values[-1]

# PLOT ALL DAY FIGURES
fig = plt.figure()
fig.set_size_inches(8.27,7.79)
fig.subplots_adjust(hspace=.24,right=.99,left=.07,top = .96)
ax2 = fig.add_subplot(2,2,1)
ax2.plot(Time_s,T_pt_su,label=r'T$_{pt,su}$')
ax2.plot(Time_s,T_pt_ex,label=r'T$_{pt,ex}$')
ax2.set_ylabel(r'$T [^{\circ}C]$')
ax2.legend(loc = 'best',numpoints = 1,fancybox='True',shadow='True',labelspacing=0.1)
#
ax2 = fig.add_subplot(2,2,2)
ax2.plot(Time_s,p_pt_su,'ro',label=r'p$_{pt,su}$')
ax2.plot(Time_s,p_tank,'o',label=r'p$_{tank}$')
ax2.set_ylabel(r'$p [bar]$')
ax2.legend(loc = 'best',numpoints = 1,fancybox='True',shadow='True',labelspacing=0.1)
#
ax2 = fig.add_subplot(2,2,3)
ax2.plot(Time_s/3600,V_wf_su,label=r'V$_{pt,su}$')
ax2.set_ylabel(r'$[m^3 s^{-1}]$')
ax2.legend(loc = 'best',numpoints = 1,fancybox='True',shadow='True',labelspacing=0.1)
#
ax2 = fig.add_subplot(2,2,4)
ax2.plot(Time_s,m_wf_su,label=r'$\dot{m}_{pt,su}$')
ax2.set_ylabel(r'$[kg s^{-1}]$')
ax2.legend(loc = 'best',numpoints = 1,fancybox='True',shadow='True',labelspacing=0.1)
#
# DNI and Temperatures and mass flow
fig = plt.figure()
fig.set_size_inches(8.27,7.79)
fig.subplots_adjust(hspace=.24,right=.99,left=.07,top = .96)
ax2 = fig.add_subplot(3,1,1)
ax2.plot(Time_hours,DNI,label=r'DNI')
ax2.set_ylabel(r'DNI [W m$^{-2}$]')
ax2.set_ylim([0,1000])
ax2.set_xlim([datetime(1900, 1, 1, 10, 30, 00),datetime(1900, 1, 1, 15, 30, 00)])
ax2.legend(loc = 'best',numpoints = 1,fancybox='True',shadow='True',labelspacing=0.1)
#
ax2 = fig.add_subplot(3,1,2)
ax2.plot(Time_hours,T_pt_ex,label=r'T$_{pt,ex}$',color='red')
ax2.plot(Time_hours,T_pt_su,label=r'T$_{pt,su}$',color='blue')
ax2.set_ylabel(r'$T [^{\circ}C]$')
ax2.set_xlim([datetime(1900, 1, 1, 10, 30, 00),datetime(1900, 1, 1, 15, 30, 00)])
ax2.legend(loc = 'best',numpoints = 1,fancybox='True',shadow='True',labelspacing=0.1)
#
ax2 = fig.add_subplot(3,1,3)
ax2.plot(Time_hours,m_wf_su,label=r'$\dot{m}_{pt,su}$')
ax2.set_ylabel(r'Mass flow $[kg s^{-1}]$')
ax2.set_xlabel(r'$Time [hours]$')
ax2.legend(loc = 'best',numpoints = 1,fancybox='True',shadow='True',labelspacing=0.1)
ax2.set_xlim([datetime(1900, 1, 1, 10, 30, 00),datetime(1900, 1, 1, 15, 30, 00)])
plt.tight_layout()
#fig.savefig(DIRECTORY_FIGURES+Days[KKK]+'/'+'_ExpData.pdf')
plt.show()


#In case you want to slice the data
#SLICE DATA FRAME
#Start_time = 40000
#Stop_time = 50000
#df = df.loc[:,'Seconds']>Start_time


# MAKE DATA FRAME SUITABLE FOR COPY PASTE INTO MODELICA SPLINE
df['TA029'] = df['TA029'] +273.15
df['TA060'] = df['TA060'] + 273.15
df['TA066'] = df['TA066'] + 273.15
df['WD089'] = df['WD089'] /3.6
df['Incidencia'] = df['Incidencia'] *np.pi/180




df['TIME_DNI'] = df["Seconds"].astype(str) + ',' + df['IA028'].astype(str) + ';'
df['TIME_Theta'] = df["Seconds"].astype(str) + ',' + df['Incidencia'].astype(str) + ';'
df['TIME_Vwind'] = df["Seconds"].astype(str) + ',' + df['WD089'].astype(str) + ';'
df['TIME_Tamb'] = df["Seconds"].astype(str) + ',' + df['TA029'].astype(str) + ';'
df['TIME_Thtfsu'] = df["Seconds"].astype(str) + ',' + df['TA060'].astype(str) + ';'
df['TIME_Thtfex'] = df["Seconds"].astype(str) + ',' + df['TA066'].astype(str) + ';'
df['TIME_M_dot_htf'] = df["Seconds"].astype(str) + ',' + df['FA032'].astype(str) + ';'
df['TIME_p_htf'] = df["Seconds"].astype(str) + ',' + df['PA052'].astype(str) + ';'

#df.to_csv('C:/Users/Admin/Documents/GitHub/PSA_SFERAII/ExperimentalData/'+Day+'_DATA.csv')




