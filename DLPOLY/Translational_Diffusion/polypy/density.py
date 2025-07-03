import numpy as np
from polypy import read as rd
from polypy import utils as ut


class Density():
    '''The Density class calculates the total number of atoms in one
    dimensional planes, two dimensional boxes and a combination of
    both. This class is ideal for observing how particle density
    changes within a structure.

    Parameters
    ----------
    data : obj
        Class object containing the trajectory
    atom_type : list (optional)
        atoms to be analysed
    '''
    def __init__(self, data, atom_type=None):
        self.data = data
        self.atom_type = atom_type
        if len(np.unique(self.data['label'])) > 1 and self.atom_type is None:
            print("Multiple atom types detected - Splitting Coordinates")
        elif len(np.unique(self.data['label'])) > 1:
            self.data = rd.get_atom(self.data, self.atom_type)
        self.rcplvs, self.lengths = ut.calculate_rcplvs(data['lv'][-1])
     #   self.frac_coords = ut.cart_2_frac(self.data['trajectories'], self.lengths, self.rcplvs)
     #   self.frac_coords = data['fractional']


    def one_dimensional_density(self, histogram_width=0.1, direction="x"):
        '''Calculate the particle density within one dimensional
        slices of a structure.

        Parameters
        ----------
        histogram_width : float (optional)
            Width of histograms, perpendicular to a given direction.
        direction : string (optional)
            Direction perpendicular to the histograms

        Returns
        -------
        x : array like
            values corresponding to the histograms coordinates.
        histograms : array like
            total number of species in each histograms.
        '''
        if direction == "x":
            val = 0
        elif direction == "y":
            val = 1
        elif direction == "z":
            val = 2
        c = self.data['frac_trajectories'][:, val]
        x = np.ceil(self.lengths[val] / histogram_width).astype(int)
        bin_array = np.zeros(x)
        c.tolist()
        for j in range(0, self.data['frac_trajectories'][:, val].size):
            plane = 0
            plane = (c[j] * x).astype(int)
            bin_array[plane] = bin_array[plane] + 1
        x = (np.arange(0, bin_array.size))
        x = x * histogram_width
        return x, bin_array

    def two_dimensional_density(self, box=0.1, direction="x"):
        '''Calculate the atomic number density within two dimensional
        boxes in a structure.

        Parameters
        ----------
        box : float (optional)
            Edge length of the boxes
        direction : string (optional)
            Direction pointing into the boxes

        Returns
        -------
        x : array like
            x axis coordinates
        y : array like
            y axis coordinates
        z : array like
            Total number of species in each box as a function of x and y.
        '''
        if direction == "x":
            val = [1, 2]
        elif direction == "y":
            val = [0, 2]
        elif direction == "z":
            val = [0, 1]
        xc = self.data['frac_trajectories'][:, val[0]]
        yc = self.data['frac_trajectories'][:, val[1]]
        x = np.ceil(self.lengths[val][0] / box).astype(int)
        y = np.ceil(self.lengths[val][1] / box).astype(int)
        if x < y:
            x, y = y, x
            xc, yc = yc, xc
        bin_array = np.zeros(((y), (x)))
        xc = xc.tolist()
        yc = yc.tolist()

        for j in range(0, self.data['frac_trajectories'][:, val[0]].size):

            xbox = 0
            ybox = 0
            xbox = (xc[j] * x).astype(int)
            ybox = (yc[j] * y).astype(int)
            bin_array[ybox, xbox] = bin_array[ybox, xbox] + 1

        x = (np.arange(0, x)) * box
        y = (np.arange(0, y)) * box
        z = bin_array + 0.001

        return x, y, z

    def one_dimensional_density_sb(self, ul, ll, direction="x"):
        '''Calculate the total number of a given species within a
        1D slice of a structure.

        Parameters
        ----------
        ul : float
            Upper bin limit
        ll : float
            Lower bin limit
        direction : string (optional)
            direction perpendicular to the slice

        Returns
        -------
        plane : int
            Total number of species within specified bin.
        '''
        if direction == "x":
            val = 0
        elif direction == "y":
            val = 1
        elif direction == "z":
            val = 2

        c = self.coords[:, val]
        c.tolist()
        plane = 0

        for j in range(0, self.coords[:, val].size):

            if c[j] > ll and c[j] < ul:
                plane = plane + 1

        return plane

    def one_and_two_dimension_overlay(self, box=0.1, direction="x"):
        '''Combination of the one and two dimensional density functions.
        Calculates both total number of atoms in one and two dimensions
        and overlays them in one plot.

        Parameters
        ----------
        box : float (optional)
            size of the boxes and bins
        direction : string (optional)
            direction perpendicular to the planes

        Returns
        -------
        x : array like
            X axis coordinates.
        y : array like
            Y axis coordinates.
        histograms : array like
            Total number of species within each histogram.
        z : array like
            Total number of species within each box.
        '''
        if direction == "x":
            val = [1, 2]
        elif direction == "y":
            val = [0, 2]
        elif direction == "z":
            val = [0, 1]

        xc = self.coords[:, val[0]]
        xc = xc + (np.average(self.lv[:, [val[0]]]) / 2)
        yc = self.coords[:, val[1]]
        yc = yc + (np.average(self.lv[:, [val[1]]]) / 2)
        x = ut.bin_choose((np.amax(self.lv[:, val[0]])), box) + 1
        y = ut.bin_choose((np.amax(self.lv[:, val[1]])), box) + 1
        if x < y:
            x, y = y, x
            xc, yc = yc, xc
        histograms = np.zeros((x))
        z = np.zeros(((y), (x)))
        xc = xc.tolist()
        yc = yc.tolist()

        for j in range(0, self.coords[:, val[0]].size):

            xbox = 0
            ybox = 0
            xbox = ut.bin_choose(xc[j], box)
            ybox = ut.bin_choose(yc[j], box)
            histograms[xbox] = histograms[xbox] + 1
            z[ybox, xbox] = z[ybox, xbox] + 1

        x = np.arange((x))
        y = np.arange((y))
        x = ((x * box)) - (np.average(self.lv[:, [val[0]]]) / 2)
        y = ((y * box))
        z = z + 0.001
        histograms = histograms + 0.001

        return x, y, z, histograms


def regional_residence_time(data, ul, ll, direction, timestep):
    if direction == "x":
        val = 0
    elif direction == "y":
        val = 1
    elif direction == "z":
        val = 2

    times = np.array([])
    for i in range(data['natoms']):
        trajectory = rd.get_trajectory(data, i)
        in_region = False
        time_in_box = 0
        previously_in = False
        for j in range(trajectory[:,0].size):
           # print(trajectory[i,val])

            if trajectory[j,val] > ll and trajectory[j,val] < ul:
                if previously_in == True:  
                    time_in_box = time_in_box + 1

                if in_region == False:
                    previously_in = True

                in_region = True            
            elif trajectory[j,val] < ll or trajectory[j,val] > ul:
                if in_region == True:
                    in_region = False
                    previously_in = False
                    if time_in_box > 10:
                        times = np.append(times, time_in_box)
                    time_in_box = 0
            
           # print(trajectory[j,val], in_region, previously_in, j, i, time_in_box)
    print(times)
    steps = np.average(times)
    times = steps * timestep
    return steps, times