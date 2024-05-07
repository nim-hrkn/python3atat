import argparse
from atat import str2atoms


def main():
    parser = argparse.ArgumentParser(description='Convert POSCAR file to a structure file.')
    parser.add_argument('strout', type=str, default='str.out', help='Path to the str.out file.')
    parser.add_argument('--poscar', type=str, default='POSCAR', help='Output format is vasp poscar. (Default: POSCAR)')
    parser.add_argument('--cif', type=str, default=None, help='Output format is cif.')
    args = parser.parse_args()
    print(args)
    # Convert atoms to structure file
    atoms = str2atoms(args.strout)

    # Read atoms from the POSCAR file
    if args.cif is not None:
        atoms.write(args.cif, format='cif')
    elif args.poscar:
        # atoms.write(args.poscar, format='vasp', direct=True)
        atoms.write(args.poscar, format='vasp')


if __name__ == "__main__":
    main()
