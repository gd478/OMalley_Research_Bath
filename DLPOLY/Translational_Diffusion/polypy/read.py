import os as os
import sys as sys
import numpy as np
from polypy import utils as ut

def read_history(file, atom_list):
    '''Read a DL_POLY HISTORY file.

    Parameters
    ----------
    file : string
        Filename to be read
    atom_list : list
        List of atoms to be read

    Returns
    -------
    data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    '''
    if os.path.isfile(file):
        trajectories = []
        frac_trajectories = []
        rcplv = []
        atname = []
        lv = []
        timesteps = 0
        count = 0
        c = 0
        name = False
        tstep = False
        history = open(file, 'r')
        current_lv = []

        for line in history:
            x = line.split()
            if c == 3:
                c = 0
                tstep = False
                current_lv = np.asarray(current_lv, dtype=float)
                rcplvs, lengths = ut.calculate_rcplvs(current_lv)
                lv.append(current_lv)
                rcplv.append(rcplvs)
                current_lv = []
            if c < 3 and tstep is True:
                if c == 4:
                  #  lv.insert(0, line.split())
                    current_lv.insert(0, line.split())
                else:
                   # lv.append(line.split())
                    current_lv.append(line.split())
                c = c + 1
            if name:
                name = False
                trajectories.append(line.split())
                frac = np.matmul( rcplvs, np.asarray(line.split(), dtype=float))
                frac = np.mod( frac, 1 )
                frac_trajectories.append(frac)

            if x[0] in atom_list:
                atname.append(x[0])
                name = True
                count = count + 1
            if x[0] == "timestep":
                timesteps = timesteps + 1
                tstep = True

        trajectories = np.asarray(trajectories, dtype=float)
        frac_trajectories = np.asarray(frac_trajectories, dtype=float)

        atname = np.asarray(atname, dtype=str)

        lv = np.asarray(lv, dtype=float)
        #rcplvs = np.asarray(rcplvs, dtype=float)

        natoms = count / timesteps
        natoms = int(natoms)

       # lv = np.split(lv, timesteps)

        data = {'label': atname,
                'trajectories': trajectories,
                'frac_trajectories': frac_trajectories,
                'rcplvs': rcplvs,
                'lv': lv,
                'timesteps': timesteps,
                'natoms': natoms}

    else:
        print("File cannot be found")
        sys.exit(0)

    if natoms == 0:
        print("No Atoms of specified type exist within the HISTORY file")
        sys.exit(0)

    history.close()
    return data


def read_archive(file, atom_list):
    '''Read a DL_MONTE ARCHIVE file
    Parameters
    ----------
    file : str
        Name of the dlmonte ARCHIVE file
    atom_list : list
        list of atoms types to be read

    Returns
    -------
    data : dict
        Dictionary containing the atom labels, trajectories,
        lattice vectors, timesteps and number of atoms
    '''
    if os.path.isfile(file):
        trajectories = []
        atname = []
        lv = []
        natoms_at_timestep = []
        count = 0
        c = 0
        atom_count = 0
        timesteps = 1
        skipline = 1
        name = False
        tstep = True
        archive = open(file, 'r')
        config_label = archive.readline()
        config_label = config_label.split()
        for line in archive:
            x = line.split()
            if c == 3:
                c = 0
                skipline = 0
                tstep = False
            if c < 3 and tstep is True and skipline == 2:
                lv.append(line.split())
                c = c + 1
            if name:
                name = False
                trajectories.append(line.split())
                atom_count = atom_count + 1
            if x[0] in atom_list:
                atname.append(x[0])
                name = True
                count = count + 1
            if x[0] == config_label[0]:
                timesteps = timesteps + 1
                tstep = True
                skipline = 1
                natoms_at_timestep.append(atom_count)
                atom_count = 0
            elif tstep is True:
                skipline = 2

        trajectories = np.asarray(trajectories, dtype=float)
        atname = np.asarray(atname, dtype=str)
        lv = np.asarray(lv, dtype=float)
        natoms = count / timesteps
        natoms = int(natoms)
        vec = np.array([])
        lv = np.split(lv, timesteps)

        for i in range(0, timesteps):

            vec = np.append(vec, (lv[i].sum(axis=0)))

        lv = np.reshape(vec, (timesteps, 3))
        data = {'label': atname,
                'trajectories': trajectories,
                'lv': lv,
                'timesteps': timesteps,
                'natoms': natoms,
                "atoms_per_timestep": natoms_at_timestep}

    else:
        print("File cannot be found")
        sys.exit(0)

    if natoms == 0:
        print("No Atoms of specified type exist within the ARCHIVE file")
        sys.exit(0)

    archive.close()
    return data


def read_config(file, atom_list):
    '''Read a DL_POLY CONFIG file

    Parameters
    ----------
    file : string
        Filename to be read
    atom : list
        List of atoms to be read

    Returns
    -------
    data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    '''
    if os.path.isfile(file):
        coords = []
        config = open(file, 'r')
        name = False
        count = 0
        lv = []
        atname = []
        frac_trajectories  = []
        title = config.readline()
        stuff = config.readline()

        for i in range(0, 3):
            l = config.readline()
            lv.append(l.split())
        lv = np.asarray(lv, dtype="float")
        print(lv)
        rcplvs, lengths = ut.calculate_rcplvs(lv)
        for line in config:
            x = line.split()
            if name:
                name = False
                coords.append(line.split())
                frac = np.matmul( rcplvs, np.asarray(line.split(), dtype=float))
                frac = np.mod( frac, 1 )
                frac_trajectories.append(frac)
            if x[0] in atom_list:
                atname.append(x[0])
                name = True
                count = count + 1

        coords = np.asarray(coords, dtype=float)
        atname = np.asarray(atname, dtype=str)
        natoms = int(count)
    #    vec = lv.sum(axis=0)

        data = {'label': atname,
                'trajectories': coords,
                'fractional': frac_trajectories,
                'lv': lv,
                'timesteps': 1,
                'natoms': natoms}
    else:
        print("File cannot be found")
        sys.exit(0)
    if natoms == 0:
        print("No Atoms of specified type exist within the CONFIG file")
        sys.exit(0)

    config.close()
    return data


def get_atom(data, atom):
    """Reads through the dictionary returned from the reading functions
    and returns a given atom.

    Parameters
    ----------
    data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    atom : str
        atom type to be found

    Returns
    -------
    atom_data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    """
    if len(np.unique(data['label'])) == 1:
        return data
    frac_coords = []
    coords = []
    count = 0
    for i in range(0, data['label'].size):
        if data['label'][i] == atom:
            coords.append(data['trajectories'][i])
            frac_coords.append(data['frac_trajectories'][i])
            count = count + 1

    coords = np.asarray(coords)
    frac_coords = np.asarray(frac_coords)
    natoms = int(count / data['timesteps'])
    atom_data = {'label': atom,
                 'trajectories': coords,
                 'frac_trajectories': frac_coords,
                 'lv': data['lv'],
                 'timesteps': data['timesteps'],
                 'natoms': natoms}

    return atom_data


def get_config(data, timestep):
    """Returns single config from entire trajectory

    Parameters
    ----------
    data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    timestep : int
        Timestep of desired config

    Returns
    -------
    array like
        Coordinates of desired config
    """
    configs = np.split(data['trajectories'], data['timesteps'])
    return configs[timestep]

def get_trajectory(data, atom):
    configs = np.split(data['trajectories'], data['timesteps'])
    configs = np.asarray(configs)
    return configs[:,atom]
