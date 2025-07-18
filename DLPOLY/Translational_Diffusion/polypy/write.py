import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



def line_plot(x, y, xlab, ylab, output, set_style="default", palette="tab10",
              figsize=None):
    '''Plots the system volume vs timestep.

    Parameters
    ----------
    x : array like
        X axis
    y : array like
        Y axis
    xlab : str (optional)
        X label
    ylab : str (optional)
        Y label
    output : str (optional)
        Output filename
    set_style : str (optional)
        Plot style
    palette : str (optional)
        Color palette
    figsize : bool (optional)
        Size of plot
    '''
    sns.palette = palette
    plt.style.use(set_style)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    ax.plot(x, y)
    ax.set_xlabel(xlab, fontsize=13)
    ax.set_ylabel(ylab, fontsize=13)
    ax.tick_params(labelsize=12)
    if output:
        plt.savefig(output, dpi=600)
    plt.tight_layout()
    plt.show()
    plt.close()


def msd_plot(msd_data, set_style="default", palette="tab10",
             figsize=None, output=None):
    '''
    MSDPlot - Plot MSD

    Parameters
    ----------
    msd_data  : Dictionary {'msd': msd, 'xmsd': xmsd,
                'ymsd': ymsd, 'zmsd': zmsd, 'time': time}
    set_style : str (optional)
        Plot style
    palette : str (optional)
        Color palette
    figsize : bool (optional)
        Size of plot
    output : str (optional)
        Output filename
    '''
    sns.palette = palette
    plt.style.use(set_style)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.set_ylim(ymin=0, ymax=np.amax(msd_data['msd']))
    ax.set_xlim(xmin=0, xmax=np.amax(msd_data['time']))
    ax.plot(msd_data['time'], msd_data['msd'], label="XYZMSD")
    ax.plot(msd_data['time'], msd_data['xymsd'], label="XYMSD")
    ax.plot(msd_data['time'], msd_data['xzmsd'], label="XZMSD")
    ax.plot(msd_data['time'], msd_data['yzmsd'], label="YZMSD")
    ax.plot(msd_data['time'], msd_data['xmsd'], label="XMSD")
    ax.plot(msd_data['time'], msd_data['ymsd'], label="YMSD")
    ax.plot(msd_data['time'], msd_data['zmsd'], label="ZMSD")
    ax.tick_params(labelsize=12)
    ax.set_xlabel("Time (ps)", fontsize=15)
    ax.set_ylabel("MSD ($\AA$)", fontsize=15)
    plt.legend()
    if output:
        plt.savefig(output, dpi=600)
    plt.show()
    plt.close()


def volume_plot(x, y, xlab="Timestep (ps)", ylab="System Volume ($\AA$)",
                output=None, set_style="default", palette="tab10",
                figsize=None):
    '''Plots the system volume vs timestep.

    Parameters
    ----------
    x : array like
        Timesteps
    y : array like
        Volume
    xlab : str (optional)
        X label
    ylab : str (optional)
        Y label
    output : str (optional)
        Output filename
    set_style : str (optional)
        Plot style
    palette : str (optional)
        Color palette
    figsize : bool (optional)
        Size of plot
    '''
    line_plot(x, y, xlab, ylab, output, set_style, palette,
              figsize)


def electric_field_plot(x, y, xlab="Coordinate ($\AA$)",
                        ylab="Electric Field", output=None,
                        set_style="default", palette="tab10",
                        figsize=None):
    '''Plots the electric field of a system.

    Parameters
    ----------
    x : array like
         X axis values - Coordinates of bins
    y : array like
        Y axis values - Electric field
    xlab : str (optional)
        X label
    ylab : str (optional)
        Y label
    output : str (optional)
        Output filename
    set_style : str (optional)
        Plot style
    palette : str (optional)
        Color palette
    figsize : bool (optional)
        Size of plot
    '''
    line_plot(x, y, xlab, ylab, output, set_style, palette,
              figsize)


def electrostatic_potential_plot(x, y, xlab="Coordinate ($\AA$)",
                                 ylab="Electrostatic Potential", output=None,
                                 set_style="default", palette="tab10",
                                 figsize=None):
    '''Plots the electrostatic potential of a system.

    Parameters
    ----------
    x : array like
         X axis values - Coordinates of bins
    y : array like
        Y axis values - Electrostatic potential
    xlab : str (optional)
        X label
    ylab : str (optional)
        Y label
    output : str (optional)
        Output filename
    set_style : str (optional)
        Plot style
    palette : str (optional)
        Color palette
    figsize : bool (optional)
        Size of plot
    '''
    line_plot(x, y, xlab, ylab, output, set_style, palette,
              figsize)


def one_dimensional_charge_density_plot(x, y, xlab="Coordinate ($\AA$)",
                                        ylab="Charge Density", output=None,
                                        set_style="default", palette="tab10",
                                        figsize=None):
    '''Plots the charge density of a system in one dimension.

    Parameters
    ----------
    x : array like
         X axis values - Coordinates of bins
    y : array like
        Y axis values - Charge density
    xlab : str (optional)
        X label
    ylab : str (optional)
        Y label
    output : str (optional)
        Output filename
    set_style : str (optional)
        Plot style
    palette : str (optional)
        Color palette
    figsize : bool (optional)
        Size of plot
    '''
    line_plot(x, y, xlab, ylab, output, set_style, palette,
              figsize)


def one_dimensional_density_plot(x, y, data_labels,
                                 xlab="X Coordinate ($\AA$)",
                                 ylab="Particle Density", output=None,
                                 set_style="default", palette="tab10",
                                 figsize=None):
    '''Plots the number density for a list of species.

    Parameters
    ----------
    x : list
        X axis values - Coordinates of bins
    y : list
        Y axis values - Number dnesity
    xlab : str (optional)
        X label
    ylab : str (optional)
        Y label
    output : str (optional)
        Output filename
    set_style : str (optional)
        Plot style
    palette : str (optional)
        Color palette
    figsize : bool (optional)
        Size of plot
    '''
    sns.palette = palette
    plt.style.use(set_style)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    for i in range(len(x)):
        ax.plot(x[i], y[i], label=data_labels[i])
    ax.set_xlabel(xlab, fontsize=13)
    ax.set_ylabel(ylab, fontsize=13)
    ax.tick_params(labelsize=12)
    plt.legend()
    plt.tight_layout()
    if output:
        plt.savefig(output, dpi=600)
    plt.show()
    plt.close()


def two_dimensional_charge_density_plot(x, y, z, xlab="X Coordinate ($\AA$)",
                                        ylab="Y Coordinate ($\AA$)",
                                        output=None,
                                        set_style="default", palette="seismic",
                                        figsize=None, colorbar=True):
    '''Plots the charge density of a system in two dimensions.

    Parameters
    ----------
    x : array like
        X axis - Coordinates
    y : array like
        Y axis - Coordinates
    z : array like
        Grid of charge densities
    xlab : str (optional)
        X label
    ylab : str (optional)
        Y label
    output : str (optional)
        Output filename
    set_style : str (optional)
        Plot style
    palette : str (optional)
        Color palette
    figsize : bool (optional)
        Size of plot
    colorbar : bool (optional)
        Colorbar on or off
    '''
    sns.palette = palette
    plt.style.use(set_style)
    fig = plt.figure(figsize=figsize)

    ax = fig.add_subplot(111)
    CM = ax.contourf(x, y, z, 50, cmap=palette)
    ax.set_xlabel(xlab, fontsize=15)
    ax.set_ylabel(ylab, fontsize=15)
    ax.tick_params(labelsize=12)
    if colorbar:
        cbar = fig.colorbar(CM)
        cbar.set_label('Charge Density', labelpad=-40, y=1.07, rotation=0)
    if output:
        plt.savefig(output, dpi=600)
    plt.show()
    plt.close()


def two_dimensional_density_plot(x, y, z, xlab="X Coordinate ($\AA$)",
                                 ylab="Y Coordinate ($\AA$)", output=None,
                                 set_style="default", palette="gray",
                                 figsize=None, colorbar=True):
    '''Plots the number density of atoms in a system in two dimensions.

    Parameters
    ----------
    x : array like
        X axis - Coordinates
    y : array like
        Y axis - Coordinates
    z : array like
        Grid of number densities
    xlab : str (optional)
        X label
    ylab : str (optional)
        Y label
    output : str (optional)
        Output filename
    set_style : str (optional)
        Plot style
    palette : str (optional)
        Color palette
    figsize : bool (optional)
        Size of plot
    colorbar : bool (optional)
        Colorbar on or off
    '''
    sns.palette = palette
    plt.style.use(set_style)
    fig = plt.figure(figsize=figsize)

    ax = fig.add_subplot(111)
    CM = ax.contourf(x, y, z, cmap=palette)
    ax.set_xlabel(xlab, fontsize=15)
    ax.set_ylabel(ylab, fontsize=15)
    ax.tick_params(labelsize=12)
    if colorbar:
        cbar = fig.colorbar(CM)
        cbar.set_label('Particle Density', labelpad=-40, y=1.07, rotation=0)
    plt.tight_layout()
    if output:
        plt.savefig(output, dpi=600)
    plt.show()
    plt.close()


def combined_density_plot(x, y, z, y2, xlab="X Coordinate ($\AA$)",
                          ylab="Y Coordinate ($\AA$)",
                          y2_lab="Particle Density",
                          output=None, set_style="default", palette="gray",
                          figsize=None):
    '''Plots the number density of atoms in a system in two dimensions
    and overlays the one dimensional plot.

    Parameters
    ----------
    x : array like
        X axis - Coordinates
    y : array like
        Y axis - Coordinates
    z : array like
        Grid of number densities
    y2 : array like
        Number density in one dimension
    xlab : str (optional)
        X label
    ylab : str (optional)
        Y label
    output : str (optional)
        Output filename
    set_style : str (optional)
        Plot style
    palette : str (optional)
        Color palette
    figsize : bool (optional)
        Size of plot
    '''
    sns.palette = palette
    plt.style.use(set_style)

    fig, ax1 = plt.subplots(figsize=figsize)
    ax2 = ax1.twinx()
    ax1.contourf(x, y, z, cmap=palette)
    ax1.set_xlabel(xlab, fontsize=15)
    ax1.set_ylabel(ylab, fontsize=15)
    ax1.set_xlim([np.amin(x), np.amax(x)])
    ax1.tick_params(labelsize=12)
    ax2.plot(x, y2)
    ax2.set_ylabel(y2_lab, fontsize=15)
    ax2.tick_params(labelsize=12)
    plt.tight_layout()

    if output:
        plt.savefig(output, dpi=600)
    plt.show()
