from polypy import read as rd
from polypy import msd as msd
from polypy import utils as ut
from polypy import write as wr
import numpy as np
from scipy import optimize
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import seaborn as sns
import time

#-----------------------------------------#
#Read history file#
#-----------------------------------------#

file = "HISTORY"
atom_list = ["O_w", "H_b"]

t1 = time.time()
data = rd.read_history(file,atom_list)
t2 = time.time()
print('Time to read HISTORY: ' + str(int((t2-t1))) + 's')

num_water = 1692 #define the number of water/adorbate molecules
num_H = 256 #define the number of protons/adsorption sites

#-----------------------------------------#
#define Ni and Pi#
#-----------------------------------------#
dist = 3.0 # define distance for contact 
box = [48.5152,48.5152,48.5152] #define dimenstion of simulation cell 

def dist_calc(coords_1, coords_2):
#Function which calculates the distance between two points - accounting for periodic boundaries.
    sqr_dist = 0
    for xyz_1, xyz_2, box_xyz in zip(coords_1,coords_2,box):
        xyz_1 = xyz_1 + box_xyz/2
        xyz_2 = xyz_2 + box_xyz/2
        if abs(xyz_1 - xyz_2) < (box_xyz-dist):
            sqr_dist += (xyz_1 - xyz_2)**2
        elif abs(xyz_1 - xyz_2) >= (box_xyz-dist):
            sqr_dist += (abs(xyz_1 - xyz_2)- box_xyz)**2 
        else:
            sqr_dist += (xyz_1 - xyz_2)**2
    return np.sqrt(sqr_dist)


def Ni_calc(coords_1, coords_2):
#Function which tests if two points are within a set cutoff distance
    if dist_calc(coords_1, coords_2) <= dist:
        return 1
    else:
        return 0
        
def test_func(t,tau):
    return np.exp(-t/tau)   
 
#-----------------------------------------#
#Sort atom lists into timesteps#
#-----------------------------------------#   

O_list = []
H_list = []

for timestep in range(0, data["timesteps"]):
    O_list.append(data["trajectories"][(timestep*(num_water+num_H)):(timestep*(num_water+num_H) + num_water)])
    H_list.append(data["trajectories"][(timestep*(num_water+num_H)+ num_water):(timestep*(num_water+num_H)+(num_water+num_H))])

#-----------------------------------------#
#Calulate Ct#
#-----------------------------------------# 

num_origin = 200 # define the number of time origins desired
len_origin = 200 # define the length of each origin

new_contact_mol_list = []
origin_ct_list = []


for origin in tqdm(range(0,num_origin)):
    timestep = 0
    ct_list = []
    o = origin*int((len(O_list)-len_origin)/num_origin)
    for timestep_O, timestep_H in zip(O_list[o:(o+len_origin)],H_list[o:(o+len_origin)]):
        if timestep == 0:
            contact_mol_list = []
            for m, molecule in enumerate(timestep_O):
                for s, site in enumerate(timestep_H):
                    if Ni_calc(molecule, site) == 1:
                        contact_mol_list.append([m,s])
            ct_list.append(1)
            initial_len = len(contact_mol_list)
            timestep += 1
        else:
            for pair in contact_mol_list:
                if Ni_calc(timestep_O[pair[0]], timestep_H[pair[1]]) == 1:
                        new_contact_mol_list.append([pair[0],pair[1]])
            ct_list.append(len(new_contact_mol_list)/initial_len)
            contact_mol_list = new_contact_mol_list
            new_contact_mol_list = []
            timestep += 1
    origin_ct_list.append(ct_list)

origined_ct = []
for timestep in range(0,(len_origin)):
    average_ct = 0
    for origin in origin_ct_list:
        average_ct += origin[timestep]
    origined_ct.append(average_ct/(len(origin_ct_list)))


#-----------------------------------------#
#Fit curve to data and plot#
#-----------------------------------------#

x_axis = np.arange(len(origined_ct))

params, params_covariance = optimize.curve_fit(test_func, x_axis,origined_ct)

t3 = time.time()

print('Tau = ' + str(params[0]))
print('Total runtime ' + str(int(t3-t1)) + 's')
#graph plotting stuff try seaborn
#read up matplot lib rcparams
#Need to add axis labels and tune the graph up a bit.

sns.set_context('paper',rc={"font.size":20,"axes.titlesize":20,"axes.labelsize":5})
sns.despine()
sns.color_palette('tab10')
plt.tight_layout()
plt.rc('font',family='serif')
plt.scatter(x_axis, origined_ct, label='Simulation', marker='x', s=3)
plt.plot(x_axis, test_func(x_axis, params[0]),label='Exponential fit')
plt.xlim(0, 100)
plt.ylabel('C(t)',fontsize=18)
plt.xlabel('Time / ps',fontsize=18)
plt.text(40,0.5, r'$\tau = $' + str(f'{params[0]:.6}')+ "ps",fontsize=18 )
plt.legend(loc='best',fontsize=14)
plt.savefig('contact_correlation.svg')

#-----------------------------------------#
#write output#
#-----------------------------------------#

output = "step , Ct , fit , tau parameter = ," + str(params[0]) + "," + "\n"
for s,step in enumerate(origined_ct):
    output += str(s) + "," + str(step) + "," + str(test_func(s, params[0])) + "\n"  

with open("contact_correlation.csv","w") as OutputFile:
    OutputFile.write(output)



