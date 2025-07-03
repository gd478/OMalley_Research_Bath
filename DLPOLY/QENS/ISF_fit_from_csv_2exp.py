from polypy import read as rd
from polypy import msd as msd
from polypy import utils as ut
from polypy import write as wr
from tqdm import tqdm
import numpy as np
import scipy.special as spcl
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import seaborn as sns
import time
import csv


#-----------------------------------------#
#input values                   #
#-----------------------------------------#

print_step = 0.1 #ps


#-----------------------------------------#
#Read csv file                   #
#-----------------------------------------#

file = "ISF_trans_test_64.csv"
print('Reading csv file...')
F_QT_list_Qsplit=[]
with open(file,"r") as csv_file:
    csv_reader = list(csv.reader(csv_file,delimiter=","))
    Q = [float(x) for x in csv_reader[1][1:-1]]
    for i in range(len(Q)):
        temp=[]
        for lines in csv_reader[2:]:
            temp.append(float(lines[i+1]))
        F_QT_list_Qsplit.append(temp)
print(Q)        
#-----------------------------------------#
#EISF fitting                             #
#-----------------------------------------# 
print("Fitting ISF functions...")


#Define function for fitting such that sum of pre-expos + base == 1 and that they are positive 
def func(FQT, Base, ratio, expo, expo_2):
    pre_total = 1 - Base
    pre = pre_total * ratio
    pre_2 = pre_total * (1-ratio)
    return pre*(np.exp(-expo*FQT)) + pre_2*(np.exp(-expo_2*FQT)) + Base

params_list = []
guess = np.array([1,0,0,0])
for q_val in F_QT_list_Qsplit:
    xdata = np.arange(len(q_val))*print_step
    ydata = q_val
    params, params_covariance = optimize.curve_fit(func, xdata, ydata, p0=guess)
    params_list.append(params)
    guess = params
    
#-----------------------------------------#
#write output files                       #
#-----------------------------------------# 

output_fit = " timestep, 0.483618819, 0.607871243, 0.729146627, 0.847027383, 0.960857079, 1.070239178, 1.174515119, 1.273198728, 1.36594491, 1.452158832, 1.531578505, \n"
timestep = 0
while timestep < len(F_QT_list_Qsplit[0]):
    Q_value = 0
    output_fit += str(timestep) + ","
    while Q_value < len(F_QT_list_Qsplit): 
        output_fit += str(func_2((timestep*print_step),params_list[Q_value][0],params_list[Q_value][1],params_list[Q_value][2],params_list[Q_value][3])) + ","
        Q_value += 1
    output_fit += "\n"
    timestep += 1 

with open("ISF_trans_test_64_fit.csv","w") as OutputFile:
    OutputFile.write(output_fit)
output_2 = "params, \n"
output_2 += "Q value, pre-exponential, exponent, baseline, \n"

output_1 = "params, \n"
output_1 += "Q value, pre-exponential 1, exponent 1, pre-exponential 2, exponent 2, Baseline \n"
for i, q in enumerate(params_list):
    output_1 += str(Q[i]) + "," + str((1-q[0])*q[1]) + "," + str(q[2])+ "," + str((1-q[0])*(1-q[1])) + "," + str(q[3]) + "," + str(q[0]) + "\n" 


with open("ISF_params_fromcsv.csv","w") as OutputFile:
    OutputFile.write(output_1)    
    
    

print("Complete.")
