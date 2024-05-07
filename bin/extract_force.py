import os
import json
import numpy as np
import argparse


def extract_force(subdir=".", optjson="opt.json"):
    filepath_opt = os.path.join(subdir, optjson)
    with open(filepath_opt) as f:
        result = json.load(f)
    print(result.keys())
    filepath = os.path.join(subdir, "energy")
    energy = str(result["total_energy"])
    with open(filepath, "w") as f:
        f.write(energy)

    filepath = os.path.join(subdir, "force.out")
    forces = result["forces"]
    with open(filepath, "w") as f:
        for force in forces:
            force_str = " ".join(map(str, force))
            f.write(force_str+"\n")

    # Matlantis stress is eV/A^3
    # 1 eV/Angstrom3 = 160.21766208 GPa
    # 10 kbar = GPa
    filepath = os.path.join(subdir, "stress.out")
    stress = np.array(result["stress"])
    print(stress)
    stress = stress*160.21766208*10  # Pa -> kbar
    print(stress)
    # I'm not sure, but probably this order of stress.
    stressmatrix = [[stress[0], stress[3], stress[5]],
                    [stress[3], stress[1], stress[4]],
                    [stress[5], stress[4], stress[2]]]
    with open(filepath, "w") as f:
        for stress in stressmatrix:
            stress_str = " ".join(map(str, stress))
            f.write(stress_str+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert opt.json file to energy and stress file.')
    parser.add_argument('optjson', type=str, default='opt.json', help='opt.json file name (default: opt.json).')

    args = parser.parse_args()

    extract_force(subdir=".", optjson=args.optjson)
