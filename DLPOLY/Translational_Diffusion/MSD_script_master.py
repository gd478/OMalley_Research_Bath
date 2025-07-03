from polypy import read as rd
from polypy import msd as msd
from polypy import utils as ut
from polypy import write as wr
import numpy as np
from tqdm import tqdm


#-----------------------------------------#
#Read history file#
#-----------------------------------------#

file = "HISTORY"
atom_list = ["C_m"]

data = rd.read_history(file,atom_list)

frac_coords = data["frac_trajectories"]
real_coords = data["trajectories"]




#-----------------------------------------#
#Calculate and output MSD#
#-----------------------------------------#

num_origin = 1500 #Number of time origins
loading = 64 #Loading of molecules
msd_len = 500 # Length of the MSD you require.
multi =  1 #Number of steps between each origin - you should aim to space this out over the whole simulation if your number of origins is not the maximum.


##CONVERTING FRAC_COORDS TO CARTESIAN
#https://chemistry.stackexchange.com/questions/136836/converting-fractional-coordinates-into-cartesian-coordinates-for-crystallography
#This should look just like the top of your DLPOLY CONFIG file. Copy this in for your system
a_vector=np.array([54.7000000000,0.0000000000,0.0000000000])
b_vector=np.array([-27.3499998866,47.3715896525,0.0000000000])
c_vector=np.array([0.0000000000,0.0000000000,59.0680000000])

#convert traj into fraction traj
my_frac_coords=real_coords*np.array([1/np.linalg.norm(a_vector),1/np.linalg.norm(b_vector),1/np.linalg.norm(c_vector)])
unit_cell=np.vstack([a_vector,b_vector,c_vector]).T
coords = np.matmul(unit_cell,my_frac_coords.T).T
max_coords = np.matmul(unit_cell,np.array([1,1,1]).T).T


#CALCING LENGTHS TO CHECK FOR PARTICLES CROSSING BOANDARIES
x_box = np.matmul(unit_cell,np.array([1,0,0]).T).T[0]
y_box = np.matmul(unit_cell,np.array([0,1,0]).T).T[1]
z_box = np.matmul(unit_cell,np.array([0,0,1]).T).T[2]

msd_list = []

for origin in tqdm(range(0, num_origin*multi, multi)):
    msd_int_list = [0]
    if origin == 0:
        for t in range(1, data["timesteps"]):
            tsd = 0
            for molecule in range(0, loading):
                x_SD = 0
                y_SD = 0
                z_SD = 0
                for xyz in range(0,3):
                    if xyz == 0:
                        if coords[molecule + loading*(origin+t-1)][xyz] - coords[molecule + loading*(origin+t)][xyz] > x_box/2:
                            coords[molecule + loading*(origin+t)][xyz] = coords[molecule + loading*(origin+t)][xyz] + x_box
                            x_SD = (coords[molecule + loading*(origin)][xyz] - coords[molecule + loading*(origin+t)][xyz])**2
                        elif coords[molecule + loading*(origin+t-1)][xyz] - coords[molecule + loading*(origin+t)][xyz] <= -x_box/2:
                            coords[molecule + loading*(origin+t)][xyz] = coords[molecule + loading*(origin+t)][xyz] - x_box
                            x_SD = (coords[molecule + loading*(origin)][xyz] - coords[molecule + loading*(origin+t)][xyz])**2
                        else:
                            x_SD = (coords[molecule + loading*(origin)][xyz] - coords[molecule + loading*(origin+t)][xyz])**2
                    elif xyz == 1:
                        if coords[molecule + loading*(origin+t-1)][xyz] - coords[molecule + loading*(origin+t)][xyz] > y_box/2:
                            coords[molecule + loading*(origin+t)][xyz] = coords[molecule + loading*(origin+t)][xyz] + y_box
                            y_SD = (coords[molecule + loading*(origin)][xyz] - coords[molecule + loading*(origin+t)][xyz])**2
                        elif coords[molecule + loading*(origin+t-1)][xyz] - coords[molecule + loading*(origin+t)][xyz] <= -y_box/2:
                            coords[molecule + loading*(origin+t)][xyz] = coords[molecule + loading*(origin+t)][xyz] - y_box
                            y_SD = (coords[molecule + loading*(origin)][xyz] - coords[molecule + loading*(origin+t)][xyz])**2
                        else:
                            y_SD = (coords[molecule + loading*(origin)][xyz] - coords[molecule + loading*(origin+t)][xyz])**2
                    elif xyz == 2:
                        if coords[molecule + loading*(origin+t-1)][xyz] - coords[molecule + loading*(origin+t)][xyz] > z_box/2:
                            coords[molecule + loading*(origin+t)][xyz] = coords[molecule + loading*(origin+t)][xyz] + z_box
                            z_SD = (coords[molecule + loading*(origin)][xyz] - coords[molecule + loading*(origin+t)][xyz])**2
                        elif coords[molecule + loading*(origin+t-1)][xyz] - coords[molecule + loading*(origin+t)][xyz] <= -z_box/2:
                            coords[molecule + loading*(origin+t)][xyz] = coords[molecule + loading*(origin+t)][xyz] - z_box
                            z_SD = (coords[molecule + loading*(origin)][xyz] - coords[molecule + loading*(origin+t)][xyz])**2
                        else:
                            z_SD = (coords[molecule + loading*(origin)][xyz] - coords[molecule + loading*(origin+t)][xyz])**2
                tsd += x_SD + y_SD + z_SD
            msd_int_list.append(tsd/loading)
    else:
        for t in range(1, msd_len):
            tsd = 0
            for molecule in range(0, loading):
                x_SD = (coords[molecule + loading*(origin)][0] - coords[molecule + loading*(origin+t)][0])**2
                y_SD = (coords[molecule + loading*(origin)][1] - coords[molecule + loading*(origin+t)][1])**2
                z_SD = (coords[molecule + loading*(origin)][2] - coords[molecule + loading*(origin+t)][2])**2
                tsd += x_SD + y_SD + z_SD
            msd_int_list.append(tsd/loading)
    msd_list.append(msd_int_list)
 

#-----------------------------------------#
#For reference calculate and output MSD#
#-----------------------------------------#

timestep = 1
msd_data = msd.msd(data, timestep)
output = "0"

for step in msd_data['msd']:
    output += '\n' + str(step)

with open("MSD_raw_old_output.csv","w") as OutputFile:
    OutputFile.write(output)
    
#-----------------------------------------#
#Do normal msd and write output#
#-----------------------------------------#

MSD_output = 'MSD \n' 
for step in msd_list[0]:
    MSD_output += str(step) + '\n'  

with open("MSD_new_non-ortho.csv","w") as OutputFile:
    OutputFile.write(MSD_output)
     
#-----------------------------------------#
#Calculate origin msd and write output#
#-----------------------------------------#

origin = 0
O_MSD = []
t = 0
while t < msd_len:
    OMSD = 0
    for list in msd_list:
        OMSD += list[t]
    O_MSD.append(OMSD/len(msd_list))    
    t+=1

origin_output = "MSD \n"
for step in O_MSD:
    origin_output += str(step) + '\n'  

with open("Origin_new_non-ortho.csv","w") as OutputFile:
    OutputFile.write(origin_output)



