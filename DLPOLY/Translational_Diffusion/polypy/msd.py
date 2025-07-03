import sys as sys
import numpy as np
from polypy import utils as ut
from scipy.constants import codata

kb = codata.value('Boltzmann constant')
ev = codata.value('electron volt')
ev = -ev


def two_dimension_square_distance(distance, n, direction):
    '''Calculate the MSD for a series of distances

    Parameters
    ----------
    distance : array like
        Distance between atomic coordinates
    n : integer
        1 = 2D array, 0 = 1D array
    direction : str
        which directions
    Returns
    -------
    msd : array like
        squared displacement
    '''
    if direction == 'x':
        val = [0, 1]
    elif direction == 'y':
        val = [0, 2]
    else:
        val = [1, 2]

    if n == 1:
        msd = (distance[:, val[0]] ** 2) + (
               distance[:, val[1]] ** 2)
    elif n == 0:
        msd = (distance[val[0]] ** 2) + (
               distance[val[1]] ** 2)
    return msd


def square_distance(distance, n):
    '''Calculate the MSD for a series of distances

    Parameters
    ----------
    distance : array like
        Distance between atomic coordinates
    n : integer
        1 = 2D array, 0 = 1D array

    Returns
    -------
    msd : array like
        squared displacement
    '''
    if n == 1:
        msd = (distance[:, 0] ** 2) + (
               distance[:, 1] ** 2) + (
               distance[:, 2] ** 2)
    elif n == 0:
        msd = (distance[0] ** 2) + (
               distance[1] ** 2) + (
               distance[2] ** 2)
    return msd


def run_msd(trajectories, lv, timesteps, natoms, start, timestep):
    '''MSD calculator

    Parameters
    ----------
    trajectories : array like
        atomic coordinates
    lv : array like
        Lattive Vectors
    timesteps : int
        Total Number of Timesteps
    natoms : int
        Total Number of Atoms
    start : int
        Total number of trajectory loops
    timestep : int
        Timestep of the simulation

    Returns
    -------
    msd_data : dictionary
        Dictionary containing 3D msd, 1D msd in the x, y, z directions
        and the time.
    '''
    trajectories = np.asarray(trajectories)

    msd = np.zeros(timesteps-1)
    xmsd = np.zeros(timesteps-1)
    ymsd = np.zeros(timesteps-1)
    zmsd = np.zeros(timesteps-1)
    xymsd = np.zeros(timesteps-1)
    xzmsd = np.zeros(timesteps-1)
    yzmsd = np.zeros(timesteps-1)
    time = np.zeros(timesteps-1)

    r0 = trajectories[start-1]
    rOd = trajectories[start-1]
    for j in range((start), timesteps):
        r1 = trajectories[j]
        distance_new = r1 - r0
        r1.tolist()
        rOd.tolist()
        if distance_new.size > 3:
            n = 1
            for k in range(0, distance_new[:, 0].size):
                for i in range(0, 3):
                    cross, r_new = ut.pbc(r1[k, i], rOd[k, i], i)
                    if cross is True:
                        r1[k, i] = r_new
                        distance_new[k, i] = r_new - r0[k, i]
        else:
            n = 0
            r1 = r1.flatten()
            rOd = rOd.flatten()
            r0 = r0.flatten()
            distance_new = distance_new.flatten()

            for i in range(0, 3):
                cross, r_new = ut.pbc(r1[i], rOd[i], i)
               # print(j, r0[i], rOd[i], r1[i], r_new, distance_new[i])

                if cross is True:
                    #print("Pay Attention")
                    r1[i] = r_new
                    distance_new[i] = r_new - r0[i]
        if n == 0:
            distance = distance_new.flatten()
        else:
            distance = distance_new

        r1 = np.asarray(r1)
        rOd = np.asarray(rOd)
        rOd = r1
        distance_n = np.array([])
        if distance.size > 3:
            for i in range(distance[:,0].size):
                d = np.matmul(lv[j], distance[i])
                distance_n = np.append(distance_n, d)
            distance_n = np.reshape(distance_n, (int((distance_n.size / 3)), 3) )

        else:
            d = np.matmul(lv[j], distance)
            distance_n = np.append(distance_n, d)
        msd_new = square_distance(distance_n, n)
        xy_msd = two_dimension_square_distance(distance_n, n, 'x')
        xz_msd = two_dimension_square_distance(distance_n, n, 'y')
        yz_msd = two_dimension_square_distance(distance_n, n, 'z')

        msd_new = np.average(msd_new)
        xy_msd = np.average(xy_msd)
        xz_msd = np.average(xz_msd)
        yz_msd = np.average(yz_msd)


        msd[j-start] = msd_new
        xymsd[j-start] = xy_msd
        xzmsd[j-start] = xz_msd
        yzmsd[j-start] = yz_msd

        time[j-start] = ((j- start) * timestep)
 
        if n == 1:
            xmsd[j-start] = (np.average((distance_n[:, 0] ** 2)))
            ymsd[j-start] = (np.average((distance_n[:, 1] ** 2)))
            zmsd[j-start] = (np.average((distance_n[:, 2] ** 2)))
        elif n == 0:
            xmsd[j-start] = (np.average((distance_n[0] ** 2)))
            ymsd[j-start] = (np.average((distance_n[1] ** 2)))
            zmsd[j-start] = (np.average((distance_n[2] ** 2)))

        msd_data = {'msd': msd,
                    'xymsd': xymsd,
                    'xzmsd': xzmsd,
                    'yzmsd': yzmsd,
                    'xmsd': xmsd,
                    'ymsd': ymsd,
                    'zmsd': zmsd,
                    'time': time}
    return msd_data


def check_trajectory(trajectory, xc, lv, timesteps, timestep, ul, ll, runs):
    '''From an assigned bin determine if any part of a trajectory crosses
    a given 1D bin.

    Parameters
    ----------
    trajectory : array like
        Trajectories
    xc : array like
        Coordinates for one dimension
    lv : array like
        Lattice vectors
    timesteps : int
        Total Number of Timesteps
    timestep : float
        Timestep of simulation
    ul : float
        Upper Bin limit
    ll : float
        Lower Bin Limit
    runs : float
        Number of trajectory sweeps

    Return
    ------
    dco : float
        Diffusion Coefficient for a given atom within a given bin
    '''
    ib = False
    count = 0
    trajectory_slice = np.array([])
    vecs = []
    msd = np.array([0])
    xmsd = np.array([0])
    ymsd = np.array([0])
    zmsd = np.array([0])
    time = np.array([0])
    do = np.array([])
    dx = np.array([])
    dy = np.array([])
    dz = np.array([])
    conductivity_count = 0
    for i in range(0, xc.size):
        if xc[i] > ll and xc[i] < ul:
            ib = True
            count = count + 1
            conductivity_count = conductivity_count + 1
            trajectory_slice = np.append(trajectory_slice, trajectory[i])
            vecs.append(lv[i])

        elif xc[i] < ll or xc[i] > ul:
            if count > 100 and ib is True:

                trajectory_slice = np.split(trajectory_slice,
                                            (trajectory_slice.size / 3))

                for i in range(0, runs):
                    start = i + 5
                    msd_data = run_msd(trajectory_slice,
                                       vecs,
                                       count,
                                       1,
                                       start,
                                       timestep)
                    msd = np.append(msd, msd_data['msd'])
                    xmsd = np.append(xmsd, msd_data['xmsd'])
                    ymsd = np.append(ymsd, msd_data['ymsd'])
                    zmsd = np.append(zmsd, msd_data['zmsd'])
                    time = np.append(time, msd_data['time'])
                    d = ut.linear_regression(msd_data['time'],
                                             msd_data['msd'])[0]
                    xd = ut.linear_regression(msd_data['time'],
                                              msd_data['xmsd'])[0]
                    yd = ut.linear_regression(msd_data['time'],
                                              msd_data['ymsd'])[0]            
                    zd = ut.linear_regression(msd_data['time'],
                                              msd_data['zmsd'])[0]

                    d = ut.three_d_diffusion_coefficient(d)
                    xd = ut.one_d_diffusion_coefficient(xd)
                    yd = ut.one_d_diffusion_coefficient(yd)
                    zd = ut.one_d_diffusion_coefficient(zd)
                    do = np.append(do, d)
                    dx = np.append(dx, xd)
                    dy = np.append(dy, yd)
                    dz = np.append(dz, zd)

                count = 0
                trajectory_slice = np.array([])
                vecs = []
            else:
                ib = False
                trajectory_slice = np.array([])
                vecs = []
                count = 0

    if count > 100 and ib is True:

        trajectory_slice = np.split(trajectory_slice,
                                    (trajectory_slice.size / 3))
        for i in range(0, runs):
            start = i + 5
            msd_data = run_msd(trajectory_slice,
                               vecs,
                               count,
                               1,
                               start,
                               timestep)
            msd = np.append(msd, msd_data['msd'])
            xmsd = np.append(xmsd, msd_data['xmsd'])
            ymsd = np.append(ymsd, msd_data['ymsd'])
            zmsd = np.append(zmsd, msd_data['zmsd'])
            time = np.append(time, msd_data['time'])
            d = ut.linear_regression(msd_data['time'],
                                     msd_data['msd'])[0]
            xd = ut.linear_regression(msd_data['time'],
                                     msd_data['xmsd'])[0]
            yd = ut.linear_regression(msd_data['time'],
                                     msd_data['ymsd'])[0]            
            zd = ut.linear_regression(msd_data['time'],
                                     msd_data['zmsd'])[0]

            d = ut.three_d_diffusion_coefficient(d)
            xd = ut.one_d_diffusion_coefficient(xd)
            yd = ut.one_d_diffusion_coefficient(yd)
            zd = ut.one_d_diffusion_coefficient(zd)
            do = np.append(do, d)
            dx = np.append(dx, xd)
            dy = np.append(dy, yd)
            dz = np.append(dz, zd)

        count = 0

    msd_data = {'time': time,
                'msd': msd,
                'xmsd': xmsd,
                'ymsd': ymsd,
                'zmsd': zmsd}

    return msd_data, conductivity_count, do, dx, dy, dz


def msd(data, timestep):
    '''Function that runs all of the parts of the MSD calcualtion.

    Parameters
    ----------
    data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    timestep : float
        Simulation timestep.
    conductivity : bool
        True/False True - calculate conductivity.
    temperature : int
        Temperature of the simulation - needed for conductivity.

    Returns
    -------
    msd_data : dictionary
        Dictionary containing 3D msd, 1D msd in the x, y, z directions
        and the time.
    '''
    if data['timesteps'] == 1:
        print("ERROR: - Only one timestep has been found")
    if data['timesteps'] < 100:
        print("WARNING: Small number of timesteps - Poor statistics likely")
    if len(np.unique(data['label'])) > 1:
        print("ERROR: MSD can only handle one atom type. Exiting")
        sys.exit(0)

    trajectories = np.split(data['frac_trajectories'], data['timesteps'])
    msd_data = run_msd(trajectories, data['lv'],
                       data['timesteps'],
                       data['natoms'],
                       1,
                       timestep)
    return msd_data


def smooth_msd(data, timestep, runs=None):
    '''MSD Launcher for a Smoothed MSD calc.

    Parameters
    ----------
    data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    timestep : float
        simulation timestep
    runs : int
        How many sweeps across the trajectory

    Returns
    -------
    Outputs diffusion info
    '''
    if runs is None:
        runs = 5

    smsd = np.array([])
    sxymsd = np.array([])
    sxzmsd = np.array([])
    syzmsd = np.array([])
    sxmsd = np.array([])
    symsd = np.array([])
    szmsd = np.array([])
    stime = np.array([])

    trajectories = np.split(data['frac_trajectories'], data['timesteps'])

    for i in range(1, runs):
        start = i * 10
        msd_data = run_msd(trajectories, data['lv'],
                           data['timesteps'],
                           data['natoms'],
                           start,
                           timestep)
        smsd = np.append(smsd, msd_data['msd'])
        sxymsd = np.append(sxymsd, msd_data['xymsd'])
        sxymsd = np.append(sxzmsd, msd_data['xzmsd'])
        sxymsd = np.append(syzmsd, msd_data['yzmsd'])
        sxmsd = np.append(sxmsd, msd_data['xmsd'])
        symsd = np.append(symsd, msd_data['ymsd'])
        szmsd = np.append(szmsd, msd_data['zmsd'])
        stime = np.append(stime, msd_data['time'])

    smsd_data = {'time': stime, 
                 'msd': smsd, 
                 'xymsd': sxymsd, 
                 'xzmsd': sxzmsd, 
                 'yzmsd': syzmsd, 
                 'xmsd': sxmsd,
                 'ymsd': symsd, 
                 'zmsd': szmsd}

    msd_data = {'time': ut.smooth_msd_data(smsd_data['time'],
                                           smsd_data['msd'])[0],
                'msd':  ut.smooth_msd_data(smsd_data['time'],
                                           smsd_data['msd'])[1],
                'xymsd': ut.smooth_msd_data(smsd_data['time'],
                                           smsd_data['xymsd'])[1],
                'xzmsd': ut.smooth_msd_data(smsd_data['time'],
                                           smsd_data['xzmsd'])[1],
                'yzmsd': ut.smooth_msd_data(smsd_data['time'],
                                           smsd_data['yzmsd'])[1],                
                'xmsd': ut.smooth_msd_data(smsd_data['time'],
                                           smsd_data['xmsd'])[1],
                'ymsd': ut.smooth_msd_data(smsd_data['time'],
                                           smsd_data['ymsd'])[1],
                'zmsd': ut.smooth_msd_data(smsd_data['time'],
                                           smsd_data['zmsd'])[1]}
    return msd_data


def plane_msd(data, timestep, runs=None, ul=None, ll=None,
              direction=None):
    '''Calculate an MSD value within a area of a structure.

    Parameters
    ----------
    data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    timestep : float
        Simulation timestep.
    runs : int
        Number of trajectory sweeps.
    ul : float
        Upper bin limit.
    ll : float
        Lower bin limit.
    direction : str
        Direction normal to slices.

    Returns
    -------
    File containing the result of the calc
    '''
    if runs is None:
        runs = 1
    if ul is None:
        sys.exit(0)
    if ll is None:
        sys.exit(0)
    if direction is None:
        direction = "x"
    if direction == "x":
        val = 0
    elif direction == "y":
        val = 1
    elif direction == "z":
        val = 2
    msd = np.array([])
    xmsd = np.array([])
    ymsd = np.array([])
    zmsd = np.array([])
    time = np.array([])

    xc = np.reshape(data['trajectories'][:, val], ((data['timesteps']),
                    data['natoms']))

    trajectories = np.split(data['frac_trajectories'], data['timesteps'])
    trajectories = np.asarray(trajectories)
    bin_atoms = np.array([])
    diffusion = np.array([])
    xdiffusion = np.array([])
    ydiffusion = np.array([])
    zdiffusion = np.array([])

    for i in range(0, (data['natoms'])):

        msd_dat, conductivity_count, do, dx, dy, dz = check_trajectory(trajectories[:, i], xc[:, i], data['lv'],
                              data['timesteps'], timestep, ul, ll, runs)
        diffusion = np.append(diffusion, do)
        xdiffusion = np.append(xdiffusion, dx)
        ydiffusion = np.append(ydiffusion, dy)
        zdiffusion = np.append(zdiffusion, dz)

        bin_atoms = np.append(bin_atoms, conductivity_count)
        msd = np.append(msd, msd_dat['msd'])
        xmsd = np.append(xmsd, msd_dat['xmsd'])
        ymsd = np.append(ymsd, msd_dat['ymsd'])
        zmsd = np.append(zmsd, msd_dat['zmsd'])
        time = np.append(time, msd_dat['time'])
    
    diffusion = np.average(diffusion)
    xdiffusion = np.average(xdiffusion)
    ydiffusion = np.average(ydiffusion)
    zdiffusion = np.average(zdiffusion)

    timesteps, counts = np.unique(time, return_counts=True)
 
    msd_data = {'time': time,
                  'msd':  msd,
                  'xmsd': xmsd,
                  'ymsd': ymsd,
                  'zmsd': zmsd,
                  'timesteps' : timesteps, 
                  'stats': counts}

    transport_data = {'Diffusion': diffusion,
                      'XDiffusion': xdiffusion,
                      'YDiffusion': ydiffusion,
                      'ZDiffusion': zdiffusion,
                      'TotalAtoms': np.sum(bin_atoms) / data['timesteps']}



    return msd_data, transport_data
