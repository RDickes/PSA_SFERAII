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


Directory ='C:/Users/Admin/Documents/GitHub/PSA_SFERAII/ExperimentalData/RAW_DATA/'
File = '06-07-2016_bis.txt'
df = pd.read_table(Directory+'/'+File,low_memory=False)

Day = '2016_07_06'

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
ax2 = fig.add_subplot(4,1,1)
ax2.plot(Time_s,DNI,label=r'DNI')
ax2.set_ylabel(r'$ [W m^{-2}]')
ax2.legend(loc = 'best',numpoints = 1,fancybox='True',shadow='True',labelspacing=0.1)
#
ax2 = fig.add_subplot(4,1,2)
ax2.plot(Time_s,T_pt_su,label=r'T$_{pt,su}$')
ax2.set_ylabel(r'$T [^{\circ}C]$')
ax2.legend(loc = 'best',numpoints = 1,fancybox='True',shadow='True',labelspacing=0.1)
#
ax2 = fig.add_subplot(4,1,3)
ax2.plot(Time_s,T_pt_ex,label=r'T$_{pt,ex}$')
ax2.set_ylabel(r'$T [^{\circ}C]$')
ax2.legend(loc = 'best',numpoints = 1,fancybox='True',shadow='True',labelspacing=0.1)
#
ax2 = fig.add_subplot(4,1,4)
ax2.plot(Time_s,m_wf_su,label=r'$\dot{m}_{pt,su}$')
ax2.set_ylabel(r'$[kg s^{-1}]$')
ax2.legend(loc = 'best',numpoints = 1,fancybox='True',shadow='True',labelspacing=0.1)
plt.show()


#In case you want to slice the data
#SLICE DATA FRAME
#Start_time = 40000
#Stop_time = 50000
#df = df.loc[:,'Seconds']>Start_time


# MAKE DATA FRAME SUITABLE FOR COPY PASTE INTO MODELICA SPLINE
df['TA029'] = df['TA029'] +273.15
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

df.to_csv('C:/Users/Admin/Documents/GitHub/PSA_SFERAII/ExperimentalData/'+Day+'_DATA.csv')