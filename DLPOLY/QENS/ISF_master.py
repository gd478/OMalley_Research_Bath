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
#Read trajectories file#
#-----------------------------------------#

file = "HISTORY"
#File must be in DL_POLY HISTORY format
print('Reading HISTORY file...')
atom_list = ["H_Cm","H_Om"] #List of atom names you wish to use (for ISF fitting this should be all the protons in your diffusive molecule)
data = rd.read_history(file,atom_list)
real_coords = data["trajectories"]
print('Reading HISTORY Complete.')

#-----------------------------------------#
#Inputs based on your system#
#Change these numbers!!!!!#
#-----------------------------------------#
loading = 48   # number of molecules
print_step = 0.1 #ps This has to equal how often your history is printing in ps
origin_len_ps = 64 #ps Should equal the max time resolution of the spectrometer you want to simulate
num_origin = 30  # comp cost variable
num_H = 4 #Number of H per molecule
origin_steps = int(origin_len_ps/print_step)

##CONVERTING FRAC_COORDS TO CARTESIAN
a_vector=np.array([40.18,0.0000000000,0.0000000000])
b_vector=np.array([0,39.476,0.0000000000])
c_vector=np.array([0.0000000000,0.0000000000,52.568])

#convert traj into fraction traj
my_frac_coords=real_coords*np.array([1/np.linalg.norm(a_vector),1/np.linalg.norm(b_vector),1/np.linalg.norm(c_vector)])
#a_unit_vector=a_vector/np.linalg.norm(a_vector)
#b_unit_vector=b_vector/np.linalg.norm(b_vector)
#c_unit_vector=c_vector/np.linalg.norm(c_vector)
unit_cell=np.vstack([a_vector,b_vector,c_vector]).T
#unit_cell_inv=np.linalg.inv(unit_cell)
coords = np.matmul(unit_cell,my_frac_coords.T).T
#real_coords = np.matmul(unit_cell_inv,real_coords.T).T
max_coords = np.matmul(unit_cell,np.array([1,1,1]).T).T


#CALCING MAX CARTESIAN LENGTHS TO CHECK FOR PARTICLES CROSSING BOANDARIES
x_box = np.matmul(unit_cell,np.array([1,0,0]).T).T[0]
y_box = np.matmul(unit_cell,np.array([0,1,0]).T).T[1]
z_box = np.matmul(unit_cell,np.array([0,0,1]).T).T[2]


#Assigin each H to its own list (you will need to make a list for each H) - original script is made for methanol with 4 H atoms.
H_1 = coords[0::4]
H_2 = coords[1::4]
H_3 = coords[2::4]
H_4 = coords[3::4]

H_list = [H_1,H_2,H_3,H_4]
#-----------------------------------------#
#Calculate rotation#
#Nothing below here should need changing unless you want new Q values#
#If you want new Qs then change the Q variable list and you will need to change the headers for the outputs in the file output section#
#I will probably automate this one day if I get time.#
#-----------------------------------------#

#calculate the self-part of the intermediate scattering function
Q = [0.2435, 0.4055, 0.5644, 0.7190, 0.8682, 1.0108, 1.1457, 1.2719, 1.3884, 1.4944, 1.5891, 1.6717, 1.7417, 1.7982]
F_QT_origin_list = []
F_QT_list_Qsplit = []
dist = 0
mol_F_QT = 0

print('Calculating ISF...')
for qi in tqdm(Q):
    F_QT_int_list = []
    for origin in range(0, int(data["timesteps"]-origin_steps), int((data["timesteps"]-origin_steps)/num_origin)):
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
                if t <= origin_steps:
                    F_QT_origin_list.append(F_QT_avg)
        else:
            for t in range(1, origin_steps):
                F_QT = 0
                F_QT_avg = 0
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
    for t in range(0,origin_steps):
        F_QT_val = 0
        for origin in F_QT_int_list:
            F_QT_val += origin[t]
        F_QT_val = (F_QT_val/len(F_QT_int_list))
        F_QT_list.append(F_QT_val)      
    F_QT_list_Qsplit.append(F_QT_list)
    
#-----------------------------------------#
#Write ISF Output#
#-----------------------------------------# 

output = "self-part of the intermediate scattering function, first row is Q values, \n"
output += " timestep, 0.2435, 0.4055, 0.5644, 0.7190, 0.8682, 1.0108, 1.1457, 1.2719, 1.3884, 1.4944, 1.5891, 1.6717, 1.7417, 1.7982, \n"
timestep = 0
while timestep < origin_steps:
    Q_value = 0
    output += str(timestep) + ","
    while Q_value < len(Q): 
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
guess = [(F_QT_list_Qsplit[0][-1]-0.02),1,0.05,1.5]
for q_val in F_QT_list_Qsplit:
    xdata = np.arange(len(q_val))*print_step
    ydata = q_val
    guess[0]=ydata[-1]-0.02
    params, params_covariance = optimize.curve_fit(func_2, xdata, ydata, bounds=([0,0,0,0],[ydata[-1],1,20,20]),p0=guess)  
    guess = [0]
    guess.extend(params[1:4])
    params_list.append(params)


#-----------------------------------------#
#write fitting output files#
#-----------------------------------------# 


output_fit = " timestep, 0.2435, 0.4055, 0.5644, 0.7190, 0.8682, 1.0108, 1.1457, 1.2719, 1.3884, 1.4944, 1.5891, 1.6717, 1.7417, 1.7982, \n"
timestep = 0
while timestep < origin_steps:
    Q_value = 0
    output_fit += str(timestep) + ","
    while Q_value < len(Q): 
        output_fit += str(func_2((timestep*print_step),params_list[Q_value][0],params_list[Q_value][1],params_list[Q_value][2],params_list[Q_value][3])) + ","
        Q_value += 1
    output_fit += "\n"
    timestep += 1 

 
with open("ISF_trans_test_64_fit.csv","w") as OutputFile:
    OutputFile.write(output_fit)
output_2 = "params, \n"
output_2 += "Q value, pre-exponential, exponent, baseline, \n"


output_1 = "Q value, pre-exponential 1, exponent 1, pre-exponential 2, exponent 2, Baseline \n"
for i, q in enumerate(params_list):
    output_1 += str(Q[i]) + "," + str((1-q[0])*q[1]) + "," + str(q[2])+ "," + str((1-q[0])*(1-q[1])) + "," + str(q[3]) + "," + str(q[0]) + "\n" 


with open("ISF_trans_test_64_parms.csv","w") as OutputFile:
    OutputFile.write(output_1)    

print("Complete.")