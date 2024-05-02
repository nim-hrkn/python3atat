import argparse
from ase.io import read
from atat import atoms2str


def main():
    parser = argparse.ArgumentParser(description='Convert POSCAR file to a structure file.')
    parser.add_argument('poscar', type=str, help='Path to the POSCAR file.')
    parser.add_argument('--strout', type=str, default='str.out', help='Output structure file name (default: str.out).')

    args = parser.parse_args()

    # Read atoms from the POSCAR file
    atoms = read(args.poscar, format='vasp')

    # Convert atoms to structure file
    atoms2str(atoms, args.strout)


if __name__ == "__main__":
    main()
