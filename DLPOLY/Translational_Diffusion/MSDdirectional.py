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
atom_list = ["O_w"]

data = rd.read_history(file, atom_list)

frac_coords = data["frac_trajectories"]
real_coords = data["trajectories"]

#-----------------------------------------#
#Calculate and output MSD#
#-----------------------------------------#

num_origin = 301 #Number of time origins
loading = 69 #Loading of molecules
msd_len = 100 # Length of the MSD you require.
multi = 1 #Number of steps between each origin - you should aim to space this out over the whole simulation if your number of origins is not the maximum.

## CONVERTING FRAC_COORDS TO CARTESIAN
a_vector = np.array([50.760869738868, 0.0000000000, 0.0000000000])
b_vector = np.array([0.00000, 50.760869738868, 0.0000000000])
c_vector = np.array([0.0000000000, 0.0000000000, 52.613368230782])

my_frac_coords = real_coords * np.array([1 / np.linalg.norm(a_vector), 1 / np.linalg.norm(b_vector), 1 / np.linalg.norm(c_vector)])
unit_cell = np.vstack([a_vector, b_vector, c_vector]).T
coords = np.matmul(unit_cell, my_frac_coords.T).T
max_coords = np.matmul(unit_cell, np.array([1, 1, 1]).T).T

# CALCULATING LENGTHS TO CHECK FOR PARTICLE CROSSING BOUNDARIES
x_box = np.matmul(unit_cell, np.array([1, 0, 0]).T).T[0]
y_box = np.matmul(unit_cell, np.array([0, 1, 0]).T).T[1]
z_box = np.matmul(unit_cell, np.array([0, 0, 1]).T).T[2]

msd_list = []

for origin in tqdm(range(0, num_origin * multi, multi)):
    msd_int_list = [[0, 0, 0, 0]]  # Stores X, Y, Z, and total MSD
    if origin == 0:
        for t in range(1, data["timesteps"]):
            x_SD = 0
            y_SD = 0
            z_SD = 0
            for molecule in range(0, loading):
                for xyz in range(3):
                    if xyz == 0:
                        if coords[molecule + loading * (origin + t - 1)][xyz] - coords[molecule + loading * (origin + t)][xyz] > x_box / 2:
                            coords[molecule + loading * (origin + t)][xyz] += x_box
                        elif coords[molecule + loading * (origin + t - 1)][xyz] - coords[molecule + loading * (origin + t)][xyz] <= -x_box / 2:
                            coords[molecule + loading * (origin + t)][xyz] -= x_box
                        x_SD += (coords[molecule + loading * (origin)][xyz] - coords[molecule + loading * (origin + t)][xyz]) ** 2
                    elif xyz == 1:
                        if coords[molecule + loading * (origin + t - 1)][xyz] - coords[molecule + loading * (origin + t)][xyz] > y_box / 2:
                            coords[molecule + loading * (origin + t)][xyz] += y_box
                        elif coords[molecule + loading * (origin + t - 1)][xyz] - coords[molecule + loading * (origin + t)][xyz] <= -y_box / 2:
                            coords[molecule + loading * (origin + t)][xyz] -= y_box
                        y_SD += (coords[molecule + loading * (origin)][xyz] - coords[molecule + loading * (origin + t)][xyz]) ** 2
                    elif xyz == 2:
                        if coords[molecule + loading * (origin + t - 1)][xyz] - coords[molecule + loading * (origin + t)][xyz] > z_box / 2:
                            coords[molecule + loading * (origin + t)][xyz] += z_box
                        elif coords[molecule + loading * (origin + t - 1)][xyz] - coords[molecule + loading * (origin + t)][xyz] <= -z_box / 2:
                            coords[molecule + loading * (origin + t)][xyz] -= z_box
                        z_SD += (coords[molecule + loading * (origin)][xyz] - coords[molecule + loading * (origin + t)][xyz]) ** 2
            msd_int_list.append([x_SD / loading, y_SD / loading, z_SD / loading, (x_SD + y_SD + z_SD) / loading])
    else:
        for t in range(1, msd_len):
            x_SD = 0
            y_SD = 0
            z_SD = 0
            for molecule in range(0, loading):
                x_SD += (coords[molecule + loading * (origin)][0] - coords[molecule + loading * (origin + t)][0]) ** 2
                y_SD += (coords[molecule + loading * (origin)][1] - coords[molecule + loading * (origin + t)][1]) ** 2
                z_SD += (coords[molecule + loading * (origin)][2] - coords[molecule + loading * (origin + t)][2]) ** 2
            msd_int_list.append([x_SD / loading, y_SD / loading, z_SD / loading, (x_SD + y_SD + z_SD) / loading])
    msd_list.append(msd_int_list)

#-----------------------------------------#
#For reference calculate and output MSD#
#-----------------------------------------#

# Create a single output file for raw MSD with 4 columns: X, Y, Z, and total MSD
raw_msd_output = "MSD_X, MSD_Y, MSD_Z, Total_MSD\n"  # Header with space alignment
for step in msd_list[0]:
    raw_msd_output += "{:<12.6f}, {:<12.6f}, {:<12.6f}, {:<12.6f}\n".format(step[0], step[1], step[2], step[3])

with open("MSD_raw_new.csv", "w") as OutputFile:
    OutputFile.write(raw_msd_output)

#-----------------------------------------#
#Calculate origin msd and write output#
#-----------------------------------------#

O_MSD_X = []
O_MSD_Y = []
O_MSD_Z = []
O_MSD_Total = []

t = 0
while t < msd_len:
    OMSD_X = 0
    OMSD_Y = 0
    OMSD_Z = 0
    OMSD_Total = 0
    for list in msd_list:
        OMSD_X += list[t][0]  # X MSD
        OMSD_Y += list[t][1]  # Y MSD
        OMSD_Z += list[t][2]  # Z MSD
        OMSD_Total += list[t][3]  # Total MSD
    O_MSD_X.append(OMSD_X / len(msd_list))
    O_MSD_Y.append(OMSD_Y / len(msd_list))
    O_MSD_Z.append(OMSD_Z / len(msd_list))
    O_MSD_Total.append(OMSD_Total / len(msd_list))
    t += 1

# Write origin MSD output to a single CSV file with 4 columns: X, Y, Z, and total MSD
origin_output = "MSD_X, MSD_Y, MSD_Z, Total_MSD\n"  # Header with space alignment
for step in range(len(O_MSD_X)):
    origin_output += "{:<12.6f}, {:<12.6f}, {:<12.6f}, {:<12.6f}\n".format(O_MSD_X[step], O_MSD_Y[step], O_MSD_Z[step], O_MSD_Total[step])

with open("Origin_MSD_new_non-ortho.csv", "w") as OutputFile:
    OutputFile.write(origin_output)
