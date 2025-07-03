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

num_mols = 128
timestep = 0
molecule_coords = []
while timestep < data["timesteps"]:
    molecule = 0
    while molecule < num_mols: 
        c_of_m = []
        xyz = 0
        while xyz <= 2:
            c_of_m.append(
            ((Cm_multi[molecule +(num_mols*timestep)][xyz]) +
            (H1_multi[molecule +(num_mols*timestep)][xyz]) +
            (H2_multi[molecule +(num_mols*timestep)][xyz]) +
            (H3_multi[molecule +(num_mols*timestep)][xyz]) +
            (Om_multi[molecule +(num_mols*timestep)][xyz]) +
            (Ho_multi[molecule +(num_mols*timestep)][xyz]))/32.0419 #Total molecular mass
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
    while molecule <num_mols:
        xyz = 0
        rot_coord = []
        while xyz <=2:
            rot_coord.append(
            H1[molecule + (num_mols*timestep)][xyz] - molecule_coords[molecule + (num_mols*timestep)][xyz]
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
    while molecule < num_mols:
        xyz=0
        while xyz <= 2:
            v1.append(rot_coords[molecule][xyz])
            v2.append(rot_coords[molecule + (num_mols*timestep)][xyz])
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
    while molecule < num_mols:
        mean_ang_sqrd_step += ang_sqrd[molecule + (num_mols*timestep)]
        molecule += 1
    mean_ang_sqrd_step = (mean_ang_sqrd_step/num_mols)
    avg_ang_sqrd.append(mean_ang_sqrd_step)
    timestep += 1


#-----------------------------------------#
#write output files#
#-----------------------------------------#

output1 = "x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8, \n"
timestep = 0
while timestep < data["timesteps"]:
    molecule = 0
    while molecule < num_mols:
        xyz=0
        while xyz <= 2:
            output1 += str(rot_coords[molecule + (num_mols*timestep)][xyz]) + ',' 
            xyz += 1
        molecule += 1
    output1 += ' \n'
    timestep += 1

output2 = "mean squard angle displacement, \n"
timestep = 0
while timestep < (data["timesteps"]-1):
    output2 += str(avg_ang_sqrd[timestep]) + "\n"
    timestep += 1    

with open("rot_raw_output.csv","w") as OutputFile:
    OutputFile.write(output1)
    
with open("rot_msad.csv","w") as OutputFile:
    OutputFile.write(output2)
  




