import seekpath
import spglib
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os


def _reciprocal_lattice_vectors(cell):
    """
    Calculate the reciprocal lattice vectors for a given real space lattice cell.

    Parameters:
    cell (tuple): Tuple of np.array representing the three real space lattice vectors (a1, a2, a3).

    Returns:
    np.array: A matrix where each row is a reciprocal lattice vector.
    """
    a1, a2, a3 = cell
    volume = np.dot(a1, np.cross(a2, a3))
    b1 = 2 * np.pi * np.cross(a2, a3) / volume
    b2 = 2 * np.pi * np.cross(a3, a1) / volume
    b3 = 2 * np.pi * np.cross(a1, a2) / volume
    return np.vstack([b1, b2, b3])


def reciprocal_distance(path_info, name1, name2):
    """
    Calculate the Euclidean distance between two points in reciprocal space.

    Parameters:
    path_info (dict): Dictionary containing 'reciprocal_point_coordinates'.
    name1 (str), name2 (str): Keys for the points to find the distance between.

    Returns:
    float: The distance between the two points.
    """
    point_coords = path_info.get('reciprocal_point_coodinates')
    rk1 = point_coords[name1]
    rk2 = point_coords[name2]
    d = np.linalg.norm(rk2-rk1)
    return d


def fractional_reciprocal_distance(path_info, name1, name2):
    """
    Calculate the Euclidean distance between two points using fractional reciprocal coordinates.

    Parameters:
    path_info (dict): Dictionary containing 'point_coords'.
    name1 (str), name2 (str): Keys for the points to find the distance between.

    Returns:
    float: The distance between the two points.
    """
    point_coords = path_info.get('point_coords')
    rk1 = point_coords[name1]
    rk2 = point_coords[name2]
    d = np.linalg.norm(rk2-rk1)
    return d


class Kpath:
    """
    Class to manage k-path generation and manipulation for band structure calculations.

    Attributes:
    atoms (ase.Atoms): Atoms object containing the atomic structure.
    """

    def __init__(self, atoms):
        """
        Initialize the Kpath object using an atomic structure.

        Parameters:
        atoms (ase.Atoms): Atoms object from which to generate the k-path.
        """
        self.atoms = atoms
        self.cell = (atoms.cell, atoms.get_scaled_positions(), atoms.get_atomic_numbers())
        self.primitive_cell = spglib.find_primitive(self.cell)
        self.cellinfo = spglib.get_symmetry_dataset(self.primitive_cell)
        self.path_info = seekpath.get_path(self.primitive_cell, True, "hpkot", 1.0e-3, 1.0e-3, 1.0)
        print("path_info", self.path_info)
        self.M = np.matrix(_reciprocal_lattice_vectors(self.path_info.get('conv_lattice')))

        reciprocal_point_coodinates = {}
        for name, fractional_k in self.path_info.get('point_coords').items():
            reciprocal_k = np.array(np.dot(self.M, np.array(fractional_k)))[0]
            reciprocal_point_coodinates.update({name: reciprocal_k.tolist()})
        self.path_info['reciprocal_point_coodinates'] = reciprocal_point_coodinates

    def reciprocal_distance(self, name1, name2):
        """
        Calculate the distance between two points in the reciprocal space defined by the path info of this Kpath object.

        Parameters:
        name1 (str), name2 (str): Keys for the points to find the distance between.

        Returns:
        float: The Euclidean distance between the two points in reciprocal space.
        """
        point_coords = self.path_info.get('reciprocal_point_coodinates')
        rk1 = np.array(point_coords[name1])
        rk2 = np.array(point_coords[name2])
        d = np.linalg.norm(rk2-rk1)
        return d

    def fractional_reciprocal_distance(self, name1, name2):
        """
        Calculate the distance between two points using fractional coordinates in the reciprocal space defined by the path info of this Kpath object.

        Parameters:
        name1 (str), name2 (str): Keys for the points to find the distance between.

        Returns:
        float: The Euclidean distance between the two points in fractional reciprocal space.
        """
        point_coords = self.path_info.get('point_coords')
        rk1 = np.array(point_coords[name1])
        rk2 = np.array(point_coords[name2])
        d = np.linalg.norm(rk2-rk1)
        return d

    def gen_kdiv(self, kpath_ndiv_min=10, distance='reciprocal'):
        """
        Generate a list of division numbers for each segment of the k-path based on the specified minimum divisions and distance metric.

        Parameters:
        ndiv_min (int): Minimum number of divisions for the smallest segment of the path.
        distance (str): Specifies the distance calculation method ('reciprocal' or 'fractional').

        Returns:
        list: Division numbers for each segment of the k-path.
        """
        self.kpath_ndiv_min = kpath_ndiv_min

        d_list = []
        for name1, name2 in self.path_info['path']:
            if distance == 'reciprocal':
                d = self.reciprocal_distance(name1, name2)
            else:
                d = self.fractional_reciprocal_distance(name1, name2)
            d_list.append(d)
        d_list = np.array(d_list)
        d_min = np.array(d_list).min()
        # 最小pathでndiv_minになるようにする。
        kdiv_list = list(map(int, d_list*kpath_ndiv_min/d_min))
        return kdiv_list

    def save_kpath(self, filename="kpath", path_division_min=50, format="atat"):
        """
        Save the generated k-path to a file.

        Parameters:
        ndiv_min (int): Minimum number of divisions for each segment of the path.
        filename (str): Filename to which the path should be saved.
        format (str): Format of the output file (e.g., 'atat').

        Effects:
        Writes the k-path data to a file in the specified format.
        """
        self.kpath_filename = filename
        self.kpath_division_min = path_division_min
        self.kpath_format = format
        point_coords = self.path_info.get('point_coords')
        self.path_info['path_division'] = self.gen_kdiv(path_division_min)
        kpath_list = []
        for kpath, kdiv in zip(self.path_info['path'], self.path_info['path_division']):
            name1 = kpath[0]
            name2 = kpath[1]
            k1 = point_coords[name1]
            k2 = point_coords[name2]
            line = [kdiv]
            line.extend(k1)
            line.extend(k2)
            line = list(map(str, line))
            kpath_list.append(" ".join(line))
            print(line, "  #", name1, name2)

        with open(filename, "w", encoding="utf-8") as f:
            f.write("\n".join(kpath_list)+"\n")

    def validate_unit(self, unit="THz"):
        unit_lower = unit.lower()
        if unit_lower not in ["thz", "ev","mev"]:
            raise ValueError(f'unknown unit={unit}')
        if unit_lower == "thz":
            unit_factor = 1e-12
        elif unit_lower == "ev":
            unit_factor = 4.135665538536E-15
        elif unit_lower == "mev":
            unit_factor = 4.135665538536E-12
        return unit_factor

    def gen_plot(self, filename=None, unit="THz"):
        unit_factor = self.validate_unit(unit)

        point_coords = self.path_info.get('reciprocal_point_coodinates')
        kstart = 0
        all_distance = []
        for kpath, kdiv in zip(self.path_info['path'], self.path_info['path_division']):
            name1 = kpath[0]
            name2 = kpath[1]
            rk1 = np.array(point_coords[name1])
            rk2 = np.array(point_coords[name2])
            path_distance = np.linalg.norm(rk2-rk1)

            # [k1,k2]をkdiv分割している。
            kx = np.linspace(0, 1, kdiv)*path_distance+kstart
            all_distance.append(kx)
            kstart += path_distance
        self.path_info['reciprocal_kpath'] = all_distance

        kstart = 0
        kticks = []
        for kpath, kdiv in zip(self.path_info['path'], self.path_info['path_division']):
            name1 = kpath[0]
            name2 = kpath[1]
            rk1 = np.array(point_coords[name1])
            rk2 = np.array(point_coords[name2])
            path_distance = np.linalg.norm(rk2-rk1)
            if len(kticks) > 0:
                if kticks[-1][1] != name1:
                    lastitem = kticks.pop()
                    kticks.append([kstart, ",".join([lastitem[1], name1])])
            kstart += path_distance
            kticks.append([kstart, name2])
        self.path_info['kpath_label'] = kticks

        fig, ax = plt.subplots()

        for k, freq in zip(self.path_info.get('reciprocal_kpath'), self.path_info.get('eigenfreq')):
            freq = np.array(freq)*unit_factor
            ax.plot(k, freq, "blue")

        xticks = []
        xticklabels = []
        for kpath, label in self.path_info['kpath_label']:
            xticks.append(kpath)
            if label == "GAMMA":
                label = r"$\Gamma$"
            xticklabels.append(label)
        lastx = kpath
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)

        # kpath <- the last k
        ylim = ax.get_ylim()
        ax.set_ylabel(unit)

        ax.vlines(xticks, ylim[0], ylim[1], "black")

        ax.set_xlim([0, lastx])
        ax.set_ylim(ylim)

        fig.tight_layout()
        if filename is not None:
            fig.savefig(filename)

    def load_eigenfreq(self, filename="vol_0/eigenfreq.out"):
        filepath = filename
        eigenfreq = np.loadtxt(filepath)
        start = 0
        all_eigenfreq = []
        path_division = self.path_info.get('path_division')
        for division in path_division:
            # print(start, division)
            all_eigenfreq.append(eigenfreq[start:start+division, :])
            start += division
        self.path_info['eigenfreq'] = all_eigenfreq
        return eigenfreq

    def dump(self, parent_dir, filename="kpath_info.pickle"):

        filepath = os.path.join(parent_dir, filename)
        with open(filepath, "wb") as f:
            pickle.dump(self, f)


"""
Usage:

PARENT_DIR = "/home/user/work/Pt_fcc"
filepath = os.path.join(PARENT_DIR, "opt.vasp")
print("load",filepath)
atoms = ase.io.read(filepath,format="vasp")
print(atoms)
kpath = Kpath(atoms)
filepath = os.path.join(PARENT_DIR,"kpath")
kpath.save_kpath(filename=filepath)
kpath.dump(PARENT_DIR) <--- to read later in phononPlotter
"""
