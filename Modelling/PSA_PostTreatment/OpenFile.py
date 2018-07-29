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


# --------------    Select the day of experiments  -----------------
Days = ['20160629', '20160630',  '20160701',  '20160704',  '20160705', '20160706']
StartModTime = [46000,       40000,        37000,      43000,          43000,     41000]
StopModTime =  [58000,       54000,        53000,      57000,          56000,     50000]

# Define simulation file
FileSimulation = 'PTTL_SF_'


# Define Path where Modelica Results are stored:
ModResults = '/Users/adriano/Google Drive/GitHub/PSA_SFERAII/Modelling/ModelicaResults/' #'C:/Users/Admin/Desktop/SFERA_II_Sim/'



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



#resultFile='/Users/adriano/Google Drive/GitHub/PSA_SFERAII/Modelling/ModelicaResults/20160629/PTTL_SF_.mat'
#r=Reader(resultFile, "dymola")
#DNI.append(r.values('DNI.y'))
#NN.append(r.values('EuroTrough.N'))
#T_su.append(r.values('SensTsu.fluidState.T'))
#T_ex.append(r.values('SensTex.fluidState.T'))
#TimeSim.append(T_ex[KKK][0])
##T_ex_exp.append(r.values('t_htf_ex.y'))
#Delta_T.append(r.values('DeltaT.Delta'))
#m_wf.append(r.values('m_dot_htf.y'))