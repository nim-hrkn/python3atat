
"""
The original code is https: // github.com/jkitchin/python-atat
"""

from ase import Atom, Atoms
from numpy import array, dot, transpose
import numpy as np


def _reorder_atoms(atoms):
    """
    Reorder the Atoms object by atomic number.

    This function takes an Atoms object, identifies all unique atomic numbers,
    and reorders all atoms in the Atoms object based on these numbers. Atoms with the same atomic
    number are grouped together.

    Parameters:
    - atoms (Atoms): An ASE Atoms object containing a collection of atoms.

    Returns:
    - Atoms: A new ASE Atoms object with atoms reordered by atomic numbers.
    """
    numbers = np.unique(atoms.get_atomic_numbers())
    print(numbers)
    scaled_position_list = []
    atomic_number_list = []
    for number in numbers:
        for atom in atoms:
            if atom.number == number:
                scaled_position_list.append(atom.scaled_position)
                atomic_number_list.append(number)
    assert (len(atoms) == len(scaled_position_list))
    atoms = Atoms(atomic_number_list, cell=atoms.cell, scaled_positions=scaled_position_list)
    return atoms


def str2atoms(file='str.out'):
    """
    Convert a structure output file (str.out) into an ASE Atoms object.

    This function reads a file in a specific format where the first three lines
    define the global coordinate system vectors, the next three lines define the lattice
    vectors in terms of the global coordinate system, and the remaining lines detail
    the atomic positions and types.

    The str.out file starts with 3 lines defining the coordinate system
    [ax]  [ay]  [az]
    [bx]  [by]  [bz]
    [cx]  [cy]  [cz]

    then there are 3 lines defining the lattice vectors in terms of the
    coordinate system:
    [ua]  [ub]  [uc]
    [va]  [vb]  [vc]
    [va]  [wb]  [wc]

    where the actual lattice vectors are given by dot([u,v,w],[a,b,c])
    I think this defines the unit cell

    then, the rest of the lines correspond to atom positions:
    [atom a] [atom b] [atom c]  [atomtype]

    where the position is given by dot([[atom a] [atom b] [atom c],[a,b,c]    

    Parameters:
    - file (str): Path to the structure output file. Default is 'str.out'.

    Returns:
    - Atoms: An ASE Atoms object constructed from the data in the file.
    """
    with open(file, 'r') as f:
        lines = f.readlines()

    numlines = len(lines)

    # these are the cooridinate system vectors: A, B, C
    A = array([float(a) for a in lines[0].split()])
    B = array([float(a) for a in lines[1].split()])
    C = array([float(a) for a in lines[2].split()])

    GCS = array([A, B, C])

    # these are the lattice vectors: U, V, W
    # U = aA + bB + cC

    U = array([float(a) for a in lines[3].split()])
    V = array([float(a) for a in lines[4].split()])
    W = array([float(a) for a in lines[5].split()])

    unitcell = array([dot(transpose(GCS), U),
                      dot(transpose(GCS), V),
                      dot(transpose(GCS), W)])

    atoms = Atoms([], cell=unitcell)

    for n in range(6, numlines):
        # each line is [a,b,c] in terms of the global coordinates
        fields = lines[n].split()

        pos = array([float(a) for a in fields[0:3]])
        pos = dot(transpose(GCS), pos)

        type = fields[-1]

        atoms.append(Atom(type, pos, tag=1))

    # nmake sure all atoms are in the cell
    spos = atoms.get_scaled_positions()

    spos_wrapped = []
    for i, pos in enumerate(spos):
        pos = pos % [1, 1, 1]

        truth = abs(pos - 1) < 1e-4
        if truth.any():
            pos[truth] = 0.0
        spos_wrapped.append(pos)

    atoms.set_scaled_positions(spos_wrapped)
    atoms.set_pbc(True)
    return _reorder_atoms(atoms)


def atoms2str(atoms, strout='str.out'):
    """
    Write an ASE Atoms object to a structure output file (str.out).

    This function creates a file that represents the atomic structure in a format
    that starts with the cell vectors, followed by default unit vectors, and then
    each atom's scaled position and atomic symbol.

    Parameters:
    - atoms (Atoms): An ASE Atoms object to be written to file.
    - strout (str): The output file path where the structure will be saved. Default is 'str.out'.

    Returns:
    - None
    """
    unitcell = atoms.get_cell()

    with open(strout, 'w') as f:
        f.write('%1.5f %1.5f %1.5f\n' % tuple(unitcell[0]))
        f.write('%1.5f %1.5f %1.5f\n' % tuple(unitcell[1]))
        f.write('%1.5f %1.5f %1.5f\n' % tuple(unitcell[2]))
        f.write('1.0 0.0 0.0\n')  # It should be.
        f.write('0.0 1.0 0.0\n')
        f.write('0.0 0.0 1.0\n')
        for atom in atoms:
            f.write('%1.5f %1.5f %1.5f %s\n' % tuple(tuple(atom.scaled_position)+(atom.symbol,)))
