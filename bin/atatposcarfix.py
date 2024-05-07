import argparse


def parse(poscarfile, output_poscarfile):

    with open(poscarfile, "r") as f:
        lines = f.readlines()

    newlines = [line.rstrip() for line in lines[:5]]

    # extract atomic species and their numbers
    positions = []
    species = []
    for line in lines[6:]:
        line = line.rstrip()
        positionandspecie = line.split()
        position = positionandspecie[:3]
        specie = positionandspecie[3]
        species.append(specie)
        positions.append(position)

    # and species and their numbers
    newlines.append(" ".join(species))
    newlines.append(" ".join(["1" for _ in species]))

    newlines.append(lines[5].rstrip())
    for position in positions:
        newlines.append(" ".join(position))

    with open(output_poscarfile, "w") as f:
        f.write("\n".join(newlines))


def main():

    parser = argparse.ArgumentParser(description='Convert POSCAR file to a structure file.')
    parser.add_argument('poscar', type=str, help='Path to the POSCAR file.')
    parser.add_argument('--output', type=str, default='POSCAR', help='Output POSCAR filepath (default: POSCAR).')

    args = parser.parse_args()
    poscarfile = args.poscar
    output_poscarfile = args.output
    parse(poscarfile, output_poscarfile)


if __name__ == "__main__":
    main()
