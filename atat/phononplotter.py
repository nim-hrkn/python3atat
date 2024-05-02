import os
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import numpy as np


def _gen_fvib(filepath, T0=0, T1=2000, dT=10, eps=1e-5):
    """
    Generate a DataFrame containing free energy values over a range of temperatures.

    Reads free energy data from a specified file and generates a DataFrame that pairs these
    energies with a series of temperature values starting from T0 to T1 in steps of dT.

    Parameters:
    - filepath (str): Path to the CSV file containing free energy data.
    - T0 (int): Start temperature for the data generation. Default is 0.
    - T1 (int): End temperature for the data generation. Default is 2000.
    - dT (int): Step size between temperatures. Default is 10.
    - eps (float): Small number to ensure inclusion of T1. Default is 1e-5.

    Returns:
    - pd.DataFrame: DataFrame containing two columns: "T" for temperature and "free_energy".
    """    
    df_fvib = pd.read_csv(filepath, index_col=None, header=None)
    df_fvib = df_fvib.rename(columns={0: "free_energy"})

    energy_list = []
    for energy in np.arange(T0, T1+eps, dT):
        energy_list.append(energy)
    df_T = pd.DataFrame(energy_list, columns=["T"])

    if df_fvib.shape[0] != df_T.shape[0]:
        print("df_fvib.shape=", df_fvib.shape, "df_T.shape[0]=", df_T.shape[0])
        raise RuntimeError("The size of fvib is different from the size of T.")

    df_fvib = pd.concat([df_fvib, df_T], axis=1)
    return df_fvib


class phononPlotter:
    """
    Class for plotting phonon-related properties from simulation data files.

    Provides methods to plot dispersion relations, density of states (DOS), and free energy
    as functions of temperature from data files located in a specified directory.

    parent_dir is the parent directory of "vol_0/".

    Attributes:
    - parent_dir (str): Path to the directory containing the data files.
    """    
    def __init__(self, parent_dir):
        """
        Initialize the phononPlotter class with a specific parent directory.

        Parameters:
        - parent_dir (str): Directory containing phonon data files.
        """        
        self.parent_dir = parent_dir

    def plot_dispersion(self, filename="vol_0/eigenfreq.out", filenamne_kpathinfo="kpath_info.pickle",
                        filename_png="eigenfreq.png"):
        """
        Plot the phonon dispersion relation.

        Loads phonon eigenfrequency data and the corresponding k-path information, then generates
        a plot saved as a PNG file.

        Parameters:
        - filename (str): Filename of the eigenfrequency data. Default is 'vol_0/eigenfreq.out'.
        - filenamne_kpathinfo (str): Filename for the k-path information. Default is 'kpath_info.pickle'.
        - filename_png (str): Filename for the output PNG plot. Default is 'eigenfreq.png'.
        """
        with open(os.path.join(self.parent_dir, filenamne_kpathinfo), "rb") as f:
            kpath = pickle.load(f)
        self.kpath = kpath

        filepath = os.path.join(self.parent_dir, filename)
        kpath = self.kpath
        _ = kpath.load_eigenfreq(filepath)  # kpath instance has data.
        kpath.gen_plot(filename_png)

    def plot_freeenergy(self, filename="fvib", filename_png="freeneergy.png"):
        """
        Plot the phonon dispersion relation.

        Parameters:
        - filename (str): Filename of the eigenfrequency data. Default is 'vol_0/eigenfreq.out'.
        - filenamne_kpathinfo (str): Filename for the k-path information. Default is 'kpath_info.pickle'.
        - filename_png (str): Filename for the output PNG plot. Default is 'eigenfreq.png'.
        """        
        filepath = os.path.join(self.parent_dir, filename)
        df_fvib_vasp = _gen_fvib(filepath)
        fig, ax = plt.subplots()
        df_fvib_vasp.plot(x="T", y="free_energy", ax=ax)
        fig.tight_layout()
        if filename_png is not None:
            fig.savefig(filename_png)

    def plot_dos(self, filename="vol_0/vdos.out", filename_png="dos.png"):
        """
        Plot the density of states (DOS).

        Reads vibration density of states data from a file and generates a plot,
        saving it as a PNG file if specified.

        Parameters:
        - filename (str): Filename of the DOS data. Default is 'vol_0/vdos.out'.
        - filename_png (str): Filename for the output PNG plot. Default is 'dos.png'.
        """     
        filepath = os.path.join(self.parent_dir, filename)
        vdos = np.loadtxt(filepath)
        fig, ax = plt.subplots()
        ax.plot(vdos[:, 0], vdos[:, 1])
        ylim = ax.get_ylim()
        ax.set_ylim([0, ylim[1]])
        fig.tight_layout()

        if filename_png is not None:
            fig.savefig(filename_png)


"""
Usage:

PARENT_DIR = "/home/user/work/Pt_fcc"

from atat import phononPlotter
plotter = phononPlotter(PARENT_DIR)
plotter.plot_dispersion()
plotter.plot_freeenergy()
plotter.plot_dos()
"""
