from polypy import read as rd
from polypy import msd as msd
from polypy import utils as ut
from polypy import write as wr
import numpy as np


#-----------------------------------------#
#Read trajectories file#
#-----------------------------------------#

file = "HISTORY"
#File must be in DL_POLY HISTORY format

#Methanol = C_m, 3* H_Cm, O_m, H_Om  
atom_list = ["C_m","H_Cm","O_m","H_Om"]
data = rd.read_history(file,atom_list)

#-----------------------------------------#
#Calculate centre of mass for each step#
#-----------------------------------------#

#data["trajectories"] is first split into lists of atoms
#all atoms of same type are listed in the order they are printed in the history
#i.e. C(1) and then c(2) from the first step are printed and then c(1) and c(2) from the 2nd step

mass = {'H': 1.007825, 'C': 12.01, 'O': 15.9994, 'N': 14.0067, 'F': 18.998403}

#split up into atoms (based on ordering in CONFIG)
Cm = data["trajectories"][0::6]
H1 = data["trajectories"][1::6]
H2 = data["trajectories"][2::6]
H3 = data["trajectories"][3::6]
Om = data["trajectories"][4::6]
Ho = data["trajectories"][5::6]

#multiply each list by its mass
Cm_multi = [i * mass['C'] for i in Cm] 
H1_multi = [i * mass['H'] for i in H1] 
H2_multi = [i * mass['H'] for i in H2]
H3_multi = [i * mass['H'] for i in H3]
Om_multi = [i * mass['O'] for i in Om]
Ho_multi = [i * mass['H'] for i in Ho]

#calculate center of mass for each molecule
#model system has 128 molecules

timestep = 0
molecule_coords = []
while timestep < data["timesteps"]:
    molecule = 0
    while molecule <=127: 
        c_of_m = []
        xyz = 0
        while xyz <= 2:
            c_of_m.append(
            ((Cm_multi[molecule +(128*timestep)][xyz]) +
            (H1_multi[molecule +(128*timestep)][xyz]) +
            (H2_multi[molecule +(128*timestep)][xyz]) +
            (H3_multi[molecule +(128*timestep)][xyz]) +
            (Om_multi[molecule +(128*timestep)][xyz]) +
            (Ho_multi[molecule +(128*timestep)][xyz]))/32.0419
            )
            xyz += 1
        molecule_coords.append(c_of_m)
        molecule += 1
    timestep += 1


#-----------------------------------------#
#Calculate rotation#
#-----------------------------------------#

#calculate vector from centre of mass to H i.e. H - molecule coords
timestep = 0
rot_coords = []
while timestep < data["timesteps"]:
    molecule = 0
    while molecule <=127:
        xyz = 0
        rot_coord = []
        while xyz <=2:
            rot_coord.append(
            H1[molecule + (128*timestep)][xyz] - molecule_coords[molecule + (128*timestep)][xyz]
            )
            xyz += 1
        rot_coords.append(rot_coord)
        molecule += 1
    timestep += 1

#calculate angle change
v1 = []
v2 = []
angles = []
timestep = 1
while timestep < data["timesteps"]:
    molecule = 0
    while molecule <=127:
        xyz=0
        while xyz <= 2:
            v1.append(rot_coords[molecule][xyz])
            v2.append(rot_coords[molecule + (128*timestep)][xyz])
            xyz += 1
        v1_u = v1 / np.linalg.norm(v1)
        v2_u = v2 / np.linalg.norm(v2)
        angles.append(
#np.rad2deg if need degrees#
        (np.arccos((np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))))
        )        
        molecule += 1
    v1 = []
    v2 = []
    timestep += 1


#calculate mean squared angular displacement
ang_sqrd = [a*a for a in angles]
avg_ang_sqrd = []
timestep = 0
while timestep < (data["timesteps"]-1):
    molecule = 0
    mean_ang_sqrd_step = 0
    while molecule <=127:
        mean_ang_sqrd_step += ang_sqrd[molecule + (128*timestep)]
        molecule += 1
    mean_ang_sqrd_step = (mean_ang_sqrd_step/128)
    avg_ang_sqrd.append(mean_ang_sqrd_step)
    timestep += 1

#calculate the self-part of the intermediate scattering function
Q = [0.01, 0.1, 0.2435, 0.4055, 0.5644, 0.7190, 0.8682, 1.0108, 1.1457, 1.2719, 1.3884, 1.4944, 1.5891, 1.6717, 1.7417, 1.7982]
print_step = 0.1
origin_len_ps = 50
F_QT_origin_list = []
F_QT_list = []
dist = 0
mol_F_QT = 0
qi = 0
while qi < 16:
    origin = 0
    F_QT_int_list = []
    while origin < 50:
        t = 0
        F_QT_origin_list = []
        while t <= origin_len_ps/print_step:
            molecule = 0
            F_QT = 0
            while molecule <=127:
                xyz = 0
                sqr_dist = 0
                while xyz <=2:
                    sqr_dist += (rot_coords[molecule +(128*(t+origin))][xyz] - rot_coords[molecule + 128*origin][xyz])**2
                    xyz += 1
                dist = np.sqrt(sqr_dist)
                mol_F_QT = ((np.sin(Q[qi]*dist))/(Q[qi]*dist))  
                F_QT += mol_F_QT
                molecule += 1
            F_QT = (F_QT/128)
            F_QT_origin_list.append(F_QT)
            t += 1
        F_QT_int_list.append(F_QT_origin_list)
        origin += 1
    
    t = 0
    while t <= origin_len_ps/print_step:
        origin = 0
        F_QT_val = 0
        while origin < 50:
            F_QT_val += F_QT_int_list[origin][t]
            origin += 1
        F_QT_val = (F_QT_val/(origin+1))
        F_QT_list.append(F_QT_val)
        t += 1    
    qi += 1   
    
#-----------------------------------------#
#write output files#
#-----------------------------------------#

output1 = "x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8, \n"
timestep = 0
while timestep < data["timesteps"]:
    molecule = 0
    while molecule <=127:
        xyz=0
        while xyz <= 2:
            output1 += str(rot_coords[molecule + (128*timestep)][xyz]) + ',' 
            xyz += 1
        molecule += 1
    output1 += ' \n'
    timestep += 1

output2 = "mean squard angle displacement, \n"
timestep = 0
while timestep < (data["timesteps"]-1):
    output2 += str(avg_ang_sqrd[timestep]) + "\n"
    timestep += 1    

output3 = "self-part of the intermediate scattering function, first row is Q values, \n"
output3 += " timestep, 0.01, 0.1, 0.2435, 0.4055, 0.5644, 0.7190, 0.8682, 1.0108, 1.1457, 1.2719, 1.3884, 1.4944, 1.5891, 1.6717, 1.7417, 1.7982, \n"
timestep = 0
while timestep <= origin_len_ps/print_step:
    Q_value = 0
    output3 += str(timestep) + ","
    while Q_value < 16: 
        output3 += str(F_QT_list[timestep + Q_value*int(origin_len_ps/print_step+1)]) + ","
        Q_value += 1
    output3 += "\n"
    timestep += 1  


with open("rot_raw_output.csv","w") as OutputFile:
    OutputFile.write(output1)
    
with open("rot_msad.csv","w") as OutputFile:
    OutputFile.write(output2)
  
with open("ISF.csv","w") as OutputFile:
    OutputFile.write(output3)




