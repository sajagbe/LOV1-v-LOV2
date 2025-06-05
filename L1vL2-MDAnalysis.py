import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align

# Load the multi-frame PDB file
u = mda.Universe("trajectory.pdb", "trajectory.pdb")

# Select protein and CA atoms
protein = u.select_atoms("protein")
ca_atoms = protein.select_atoms("name CA")

# Align all frames to the first frame using Cα atoms
align.AlignTraj(u, u, select="protein and name CA", in_memory=True).run()

# Number of residues
n_residues = len(protein.residues)

# Initialize an array for CA coordinates per frame
all_coords = np.zeros((len(u.trajectory), n_residues, 3))

# Collect CA coordinates for each frame
for i, ts in enumerate(u.trajectory):
    all_coords[i] = ca_atoms.positions

# Calculate per-residue RMSD (i.e., fluctuation of each Cα over time)
rmsd_per_residue = np.sqrt(((all_coords - all_coords.mean(axis=0))**2).sum(axis=2).mean(axis=0))

# # ##################################################################################
# Or alternatively Add In normalization
# Normalize the RMSD values to the range [0, 1]
# rmsd_min = rmsd_per_residue.min()
# rmsd_max = rmsd_per_residue.max()
# rmsd_normalized = (rmsd_per_residue - rmsd_min) / (rmsd_max - rmsd_min)



# Assign the RMSD value to the B-factor of all atoms in each corresponding residue
u.trajectory[0]  # Go to the first frame
for residue, rmsd in zip(protein.residues, rmsd_per_residue):
    for atom in residue.atoms:
        atom.tempfactor = rmsd

# ####################################################################################
# # To Assign the normalized RMSD to B-factors of all atoms in each residue
# u.trajectory[0]  # Go to the first frame
# for residue, norm_rmsd in zip(protein.residues, rmsd_normalized):
#     for atom in residue.atoms:
#         atom.tempfactor = norm_rmsd       
# 
# #################################################################################### 

# Write the structure with updated B-factors
with mda.Writer("rmsd_per_residue.pdb", multiframe=False) as w:
    w.write(protein)


# # Write the structure with updated B-factors
# with mda.Writer("normalized_rmsd_per_residue.pdb", multiframe=False) as w:
#     w.write(protein)
