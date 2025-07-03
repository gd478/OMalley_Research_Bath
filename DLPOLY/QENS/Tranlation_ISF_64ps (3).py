from socket import has_ipv6
from polypy import read as rd
from polypy import msd as msd
from polypy import utils as ut
from polypy import write as wr
from tqdm import tqdm
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import seaborn as sns
import time

#-----------------------------------------#
#Inputs based on your system#
#Change these numbers!!!!!#
#-----------------------------------------#
loading = 9
print_step = 0.1 #ps
origin_len_ps = 64 #ps
num_origin = 1
num_H = 3
box = [50.5280, 50.5280, 52.3720]
atom_list = ["H","H_N"]


x_box = box[0]
y_box = box[1]
z_box = box[2]
#-----------------------------------------#
#Read trajectories file#
#-----------------------------------------#

file = "HISTORY"
#File must be in DL_POLY HISTORY format
print('Reading HISTORY file...')
data = rd.read_history(file,atom_list)
coords = data["trajectories"]
print('Reading HISTORY Complete.')


H_1 = coords[0::3]
H_2 = coords[1::3]
H_3 = coords[2::3]

H_list = [H_1, H_2, H_3]

#-----------------------------------------#
#Calculate rotation#
#-----------------------------------------#

#calculate the self-part of the intermediate scattering function
Q = [0.483618819, 0.607871243, 0.729146627, 0.847027383, 0.960857079, 1.070239178, 1.174515119, 1.273198728, 1.36594491, 1.452158832, 1.531578505]
F_QT_origin_list = []
F_QT_list_Qsplit = []
dist = 0
mol_F_QT = 0

print('Calculating ISF...')
for qi in tqdm(Q):
    F_QT_int_list = []
    for origin in range(0, int(data["timesteps"]-(origin_len_ps/print_step)), int((data["timesteps"]-(origin_len_ps/print_step))/num_origin)):
        F_QT_origin_list = []
        F_QT_origin_list.append(1)
        if origin == 0:
            for t in range(1, data["timesteps"]):
                F_QT = 0
                F_QT_avg = 0
                for molecule in range(0, loading):
                    for H_i in H_list:
                        sqr_dist = 0
                        for xyz in range(0,3):
                            if xyz == 0:
                                if H_i[molecule + loading*(origin+t-1)][xyz] - H_i[molecule + loading*(origin+t)][xyz] > x_box/2:
                                    H_i[molecule + loading*(origin+t)][xyz] = H_i[molecule + loading*(origin+t)][xyz] + x_box
                                    x_SD = (H_i[molecule + loading*(origin)][xyz] - H_i[molecule + loading*(origin+t)][xyz])**2
                                elif H_i[molecule + loading*(origin+t-1)][xyz] - H_i[molecule + loading*(origin+t)][xyz] <= -x_box/2:
                                    H_i[molecule + loading*(origin+t)][xyz] = H_i[molecule + loading*(origin+t)][xyz] - x_box
                                    x_SD = (H_i[molecule + loading*(origin)][xyz] - H_i[molecule + loading*(origin+t)][xyz])**2
                                else:
                                    x_SD = (H_i[molecule + loading*(origin)][xyz] - H_i[molecule + loading*(origin+t)][xyz])**2
                            elif xyz == 1:
                                if H_i[molecule + loading*(origin+t-1)][xyz] - H_i[molecule + loading*(origin+t)][xyz] > y_box/2:
                                    H_i[molecule + loading*(origin+t)][xyz] = H_i[molecule + loading*(origin+t)][xyz] + y_box
                                    y_SD = (H_i[molecule + loading*(origin)][xyz] - H_i[molecule + loading*(origin+t)][xyz])**2
                                elif H_i[molecule + loading*(origin+t-1)][xyz] - H_i[molecule + loading*(origin+t)][xyz] <= -y_box/2:
                                    H_i[molecule + loading*(origin+t)][xyz] = H_i[molecule + loading*(origin+t)][xyz] - y_box
                                    y_SD = (H_i[molecule + loading*(origin)][xyz] - H_i[molecule + loading*(origin+t)][xyz])**2
                                else:
                                    y_SD = (H_i[molecule + loading*(origin)][xyz] - H_i[molecule + loading*(origin+t)][xyz])**2
                            elif xyz == 2:
                                if H_i[molecule + loading*(origin+t-1)][xyz] - H_i[molecule + loading*(origin+t)][xyz] > z_box/2:
                                    H_i[molecule + loading*(origin+t)][xyz] = H_i[molecule + loading*(origin+t)][xyz] + z_box
                                    z_SD = (H_i[molecule + loading*(origin)][xyz] - H_i[molecule + loading*(origin+t)][xyz])**2
                                elif H_i[molecule + loading*(origin+t-1)][xyz] - H_i[molecule + loading*(origin+t)][xyz] <= -z_box/2:
                                    H_i[molecule + loading*(origin+t)][xyz] = H_i[molecule + loading*(origin+t)][xyz] - z_box
                                    z_SD = (H_i[molecule + loading*(origin)][xyz] - H_i[molecule + loading*(origin+t)][xyz])**2
                                else:
                                    z_SD = (H_i[molecule + loading*(origin)][xyz] - H_i[molecule + loading*(origin+t)][xyz])**2
                        sqr_dist += x_SD + y_SD + z_SD
                        F_QT += ((np.sin(qi*(np.sqrt(sqr_dist))))/(qi*(np.sqrt(sqr_dist))))
                F_QT_avg = (F_QT/(loading*num_H))
                if t <= (origin_len_ps/print_step):
                    F_QT_origin_list.append(F_QT_avg)
        else:
            for t in range(1, int(origin_len_ps/print_step)):
                F_QT = 0
                F_QT_avg = 0
                if t > (origin_len_ps/print_step): break
                for molecule in range(0, loading):
                    for H_i in H_list:
                        sqr_dist = 0
                        x_SD = (H_i[molecule + loading*(origin)][0] - H_i[molecule + loading*(origin+t)][0])**2
                        y_SD = (H_i[molecule + loading*(origin)][1] - H_i[molecule + loading*(origin+t)][1])**2
                        z_SD = (H_i[molecule + loading*(origin)][2] - H_i[molecule + loading*(origin+t)][2])**2
                        sqr_dist += x_SD + y_SD + z_SD
                        F_QT += ((np.sin(qi*(np.sqrt(sqr_dist))))/(qi*(np.sqrt(sqr_dist))))
                F_QT_avg = (F_QT/(loading*num_H))
                F_QT_origin_list.append(F_QT_avg)    
        F_QT_int_list.append(F_QT_origin_list)

    F_QT_list = []
    for t in range(0,int(origin_len_ps/print_step)):
        F_QT_val = 0
        for origin in F_QT_int_list:
            print(t)
            F_QT_val += origin[t]
        F_QT_val = (F_QT_val/len(F_QT_int_list))
        F_QT_list.append(F_QT_val)      
    F_QT_list_Qsplit.append(F_QT_list)

######  WRITING OUTPUT FOR ISF
output = "self-part of the intermediate scattering function, first row is Q values, \n"
output += " timestep, 0.483618819, 0.607871243, 0.729146627, 0.847027383, 0.960857079, 1.070239178, 1.174515119, 1.273198728, 1.36594491, 1.452158832, 1.531578505, \n"
timestep = 0
while timestep < (origin_len_ps/print_step):
    Q_value = 0
    output += str(timestep) + ","
    while Q_value < 11: 
        output += str(F_QT_list_Qsplit[Q_value][timestep]) + ","
        Q_value += 1
    output += "\n"
    timestep += 1  
with open("ISF_trans_test_64.csv","w") as OutputFile:
    OutputFile.write(output)



#-----------------------------------------#
#EISF fitting#
#-----------------------------------------# 
print("Fitting ISF functions...")

def func_1(FQT,pre,expo,base,):
    return pre*(np.exp(-expo*FQT)) + base
    
def func_2(FQT, Base, ratio, expo, expo_2):
    pre_total = 1 - Base
    pre = pre_total * ratio
    pre_2 = pre_total * (1-ratio)
    return pre*(np.exp(-expo*FQT)) + pre_2*(np.exp(-expo_2*FQT)) + Base


params_list = []
guess = [0,1,0.05,1.5]
for q_val in F_QT_list_Qsplit:
    xdata = np.arange(len(q_val))*print_step
    ydata = q_val
    if ydata[-1] <= 0:
        B_line_lim = 0.0000000001
    else:
        B_line_lim = ydata[-1]
    params, params_covariance = optimize.curve_fit(func_2, xdata, ydata, bounds=([0,0,0,0],[B_line_lim,np.inf,np.inf,np.inf]),p0=guess)  
    guess = [0]
    guess.extend(params[1:4])
    params_list.append(params)


#-----------------------------------------#
#write output files#
#-----------------------------------------# 


output_fit = " timestep, 0.483618819, 0.607871243, 0.729146627, 0.847027383, 0.960857079, 1.070239178, 1.174515119, 1.273198728, 1.36594491, 1.452158832, 1.531578505, \n"
timestep = 0
while timestep < origin_len_ps/print_step:
    Q_value = 0
    output_fit += str(timestep) + ","
    while Q_value < 11: 
        output_fit += str(func_2((timestep/print_step),params_list[Q_value][0],params_list[Q_value][1],params_list[Q_value][2],params_list[Q_value][3])) + ","
        Q_value += 1
    output_fit += "\n"
    timestep += 1 

  

 
with open("ISF_trans_test_64_fit_2.csv","w") as OutputFile:
    OutputFile.write(output_fit)
output_2 = "params, \n"
output_2 += "Q value, pre-exponential, exponent, baseline, \n"



output_1 = "Q value, pre-exponential 1, exponent 1, pre-exponential 2, exponent 2, Baseline \n"
for i, q in enumerate(params_list):
    output_1 += str(Q[i]) + "," + str((1-q[0])*q[1]) + "," + str(q[2])+ "," + str((1-q[0])*(1-q[1])) + "," + str(q[3]) + "," + str(q[0]) + "\n" 


with open("ISF_trans_test_64_parms_2.csv","w") as OutputFile:
    OutputFile.write(output_1)    

print("Complete.")