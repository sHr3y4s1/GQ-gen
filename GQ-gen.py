# Import packages
import numpy as np
from Bio.PDB import PDBIO,  Atom, Residue, Chain, Model, Structure

# Generate a single guanine base
def generate_guanine_base(center, rot):
    base = []
    base_coords = {
    "P":[1.654,8.443,2.132],
    "O1P":[1.699,9.714,2.889],
    "O2P":[0.676,8.381,1.024],
    "C5\'":[4.031,7.381,2.43],
    "O5\'":[3.115,8.094,1.579],
    "C4\'":[4.736,6.297,1.639],
    "O4\'":[3.996,5.042,1.636],
    "C3\'":[4.962,6.595,0.157],
    "O3\'":[6.22,6.057,-0.234],
    "C2\'":[3.858,5.823,-0.566],
    "C1\'":[3.923,4.569,0.3],
    "N1":[2.243,-0.242,0.309],
    "C2":[3.526,0.248,0.456],
    "N2":[4.486,-0.664,0.628],
    "N3":[3.821,1.546,0.433],
    "C4":[2.719,2.32,0.249],
    "C5":[1.411,1.925,0.094],
    "C6":[1.096,0.542,0.119],
    "O6":[0,0,0],
    "N7":[0.58,3.032,-0.07],
    "C8":[1.404,4.051,-0.01],
    "N9":[2.721,3.697,0.183]
    }
# base coordinates centered at 0
    keys = list(base_coords.keys())
    coords = np.array([base_coords[k] for k in keys])
    new_coords = []
    for i in range(len(coords)):
        new_coords.append(np.dot(rot, coords[i]) + center)
    alt_base_coords = {k: list(c) for k, c in zip(keys, new_coords)}

    for i, (atom_name, coord) in enumerate(alt_base_coords.items()):
        atom = Atom.Atom(
            name=atom_name,
            coord=coord,
            bfactor=1.0,
            occupancy=1.0,
            altloc=' ',
            fullname=atom_name.rjust(4),
            serial_number=1 + i,
            element=atom_name[0] if atom_name[0] != 'O' else 'O'  # O3', O5' etc.
        )
        base.append(atom)
    return base

# Generate GQ tetrad
def generate_tetrad(center=np.array([0.0, 0.0, 0.0])):
    tetrad_atoms = []
    center = [[1.7,1.7,0],[-1.7,1.7,0],[-1.7,-1.7,0],[1.7,-1.7,0]]  # define center
    angles = [45, 135, 225, 315]  # angle in degrees

    for i in range(len(angles)):
        theta = np.radians(angles[i])
        rot = np.array([
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta),  np.cos(theta), 0],
            [0,              0,             1]
        ])
        guanine_atoms = generate_guanine_base(center[i], rot)

        # Create a Residue object for each guanine
        res = Residue.Residue((' ', i + 1, ' '), 'G', '')
        for atom in guanine_atoms:
            res.add(atom)
        tetrad_atoms.append(res)

    return tetrad_atoms

# Save GQ tetrad to PDB file
def save_tetrad_to_pdb(tetrad_atoms, filename="tetrad.pdb"):
    structure = Structure.Structure('GQ')
    model = Model.Model(0)
    chain = Chain.Chain('A')

    for res in tetrad_atoms:
        chain.add(res)
    model.add(chain)
    structure.add(model)

    io = PDBIO()
    io.set_structure(structure)
    io.save(filename)
    print(f"Tetrad saved to {filename}")

# Call tetrad generating function
if __name__ == "__main__":
    tetrad = generate_tetrad()
    save_tetrad_to_pdb(tetrad)
