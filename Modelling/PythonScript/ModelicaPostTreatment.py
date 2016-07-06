from __future__ import division
__author__ = 'Admin'

from CoolProp.CoolProp import PropsSI
from os import listdir
from os.path import isfile, join
import csv
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.cm as cm
from uncertainties import unumpy
import csv


# GET EXPERIMENTAL DATA
Directory ='C:/Users/Admin/Documents/GitHub/PSA_SFERAII/ExperimentalData/RAW_DATA/'
File = '29-06-2016_bis.txt'
df = pd.read_table(Directory+'/'+File,low_memory=False)

# GET MODELICA DATA
Directory_Sim = "C:/Users/Admin/Desktop/SFERA_II_Sim/"
File_Sim = ["StepDown4Hz"]
StartModTime = [0]
StopModTime = [1500]

TimeSim = []
T_eva_su =[]
T_exp_su =[]
m_orc = []
T_cond_sf_su = []
P_exp = []
p_exp_ex = []
T_cond_ex = []
T_tank_ex = []
p_exp_su =[]
for kk in range(len(File_Sim)):
    Data_Sim = Directory_Sim + File_Sim[kk]+'.mat'
    sim = SimRes(Data_Sim)
    initialTime = (sim['Time'].IV())
    FinalTime=(sim['Time'].FV())
    TimeSim.append(((sim['Time'].values(t=(StartModTime[0],StopModTime[0]))-StartModTime[0])))
    T_eva_su.append((sim['Sensor_T_eva_su.T'].values(t=(StartModTime[0],StopModTime[0])))-273.15)
    T_exp_su.append((sim['Sensor_T_exp_su.T'].values(t=(StartModTime[0],StopModTime[0])))-273.15)
    m_orc.append((sim['Pump.M_dot'].values(t=(StartModTime[0],StopModTime[0]))))
    T_cond_sf_su.append((sim['Condenser.inletSf.T'].values(t=(StartModTime[0],StopModTime[0]))-273.15))
    p_exp_su.append((sim['Expander.InFlow.p'].values(t=(StartModTime[0],StopModTime[0]))/1e5))
    P_exp.append((sim['Expander.W_dot'].values(t=(StartModTime[0],StopModTime[0]))))
    p_exp_ex.append((sim['Expander.OutFlow.p'].values(t=(StartModTime[0],StopModTime[0]))/1e5))
    T_cond_ex.append((sim['Sensor_T_cond_ex.T'].values(t=(StartModTime[0],StopModTime[0]))))
    T_tank_ex.append((sim['Pump.fluidState.T'].values(t=(StartModTime[0],StopModTime[0]))-273.15))