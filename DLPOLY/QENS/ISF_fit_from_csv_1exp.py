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

file = "ISF.csv"
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
def func(FQT, Base, expo):
    return (1-Base)*(np.exp(-expo*FQT)) + Base

params_list = []
guess = np.array([1,1])
for q_val in F_QT_list_Qsplit:
    xdata = np.arange(len(q_val))*print_step
    ydata = q_val
    params, params_covariance = optimize.curve_fit(func, xdata, ydata, p0=guess)
    params_list.append(params)
    guess = params

    
#-----------------------------------------#
#write output files                       #
#-----------------------------------------# 


output_1 = "params, \n"
output_1 += "Q value, pre-exponential 1, exponent 1, Baseline \n"
for i, q in enumerate(params_list):
    output_1 += str(Q[i]) + "," + str(q[1])+ "," + str(q[0]) + "\n" 

with open("ISF_params_fromcsv.csv","w") as OutputFile:
    OutputFile.write(output_1)    