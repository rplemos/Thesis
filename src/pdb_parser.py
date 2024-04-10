from classes import Protein, Chain, Residue, Atom
import numpy as np
import sys
import os

    
stacking = {
    'PHE':[11, 'CG','CZ'],
    'TYR':[12, 'CG','CZ'],
    'TRP':[14, 'CD2','CE2'],
    'HIS':[10, 'CG','ND1','CE1','NE2','CD2']
}
    
# not used! instead, use sysfunctions.cl_parse()
def open_pdbs():
    """
    Gets PDB files from command line parameters and handles errors.
    
    Returns:
        list: The PDB file names provided as command line arguments.
    
    Raises:
        ValueError: If the file type is not ".pdb" or if incorrect usage of the script is detected.
        FileNotFoundError: If any of the provided PDB files do not exist.
    """
    try:
        # if len(sys.argv) < 3:
        #     raise ValueError("Usage: python3 parse.py <pdb_file1> <pdb_file2> ...")

        for pdb_file_arg in sys.argv[1:]:
            if not pdb_file_arg.lower().endswith(".pdb"):
                raise ValueError(f"Invalid file type. Please provide a .pdb file: {pdb_file_arg}")

            if not os.path.exists(pdb_file_arg):
                raise FileNotFoundError(f"File not found: {pdb_file_arg}")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1) # error exiting
        
    pdb_files = []
    for file in sys.argv[1:]:
        pdb_files.append(file)
        
    return pdb_files

def parse_pdb(pdb_files):
    """
    Parses the PDB files and constructs Protein, Chain, Residue, and Atom objects.

    Args:
        pdb_files (list): List of PDB file names to parse.

    Returns:
        list: List of Protein objects representing the parsed proteins.
    """

    proteins = []
    for file in pdb_files:
        current_protein = Protein()
        current_chain = None
        current_residue = None

        with open(file) as f:
            for line in f:
                line = line.strip()

                if line.startswith("HEADER"):
                    current_protein.id = line[62:]
                elif line.startswith("TITLE"):
                    current_protein.set_title(line[10:])
                elif line.startswith("ATOM"):
                    chain_id = line[21]
                    resnum = int(line[22:26])
                    resname = line[17:20]

                    if resname == "HIE" or resname == "HID":  # alternative names for protonated histidines
                        resname = "HIS"                        

                    if current_chain is None or current_chain.id != chain_id:  # new chain
                        current_chain = Chain()
                        current_chain.id = chain_id
                        current_protein.chains.append(current_chain)

                    if current_residue is None or current_residue.resnum != resnum:  # new residue
                        current_residue = Residue()
                        current_residue.resnum = resnum
                        current_residue.resname = resname
                        current_residue.chain = current_chain  # Set the parent chain
                        current_chain.residues.append(current_residue)
                                                
                    atom = Atom()  # each line is a new atom anyway
                    
                    atomname = line[12:16].replace(" ", "")

                    if atomname == "OXT": # OXT is the C-terminal Oxygen atom. However, it exhibits the same properties of any Oxygen
                        atomname = atomname.replace("OXT","O")
                                      
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    
                    atom.set_atom_info(atomname, x, y, z)
                    atom.residue = current_residue  # Set the parent residue
                    current_residue.atoms.append(atom)
                    
                    # checks if the residue is aromatic and if the atoms are complete (all populated)
                    if current_residue.resname in stacking and len(current_residue.atoms) == stacking[current_residue.resname][0]: 
                        ring_atoms = np.array([[atom.x, atom.y, atom.z] for atom in current_residue.atoms[5:]]) # ignores [N, CA, C, O] and [RNG] atoms
                        centroid_atom = centroid(current_residue)
                        centroid_atom2 = centroid2(current_residue, ring_atoms)
                        current_residue.atoms.append(centroid_atom2)
                        current_residue.ring = True
                        # print("two atoms:", current_residue.resname, current_residue.resnum, centroid_atom.x, centroid_atom.y, centroid_atom.z)
                        # print("all atoms:", current_residue.resname, current_residue.resnum, centroid_atom2.x, centroid_atom2.y, centroid_atom2.z)
                        # print()

                        normal_vector = calc_normal_vector(ring_atoms)
                        current_residue.normal_vector = normal_vector

                elif line.startswith("END"):  # finishes and resets everything for the new protein
                    proteins.append(current_protein)
                    current_protein = None
                    current_chain = None
                    current_residue = None
    return proteins

def centroid(residue):
    size = len(stacking[residue.resname])-1
    sum_x, sum_y, sum_z = 0, 0, 0
    
    for atom in residue.atoms:
        if atom.atomname in stacking[residue.resname]:
            sum_x += atom.x
            sum_y += atom.y
            sum_z += atom.z

    centroid_x = sum_x / size
    centroid_y = sum_y / size
    centroid_z = sum_z / size
    
    centroid_atom = Atom()
    centroid_atom.set_atom_info("RNG", centroid_x, centroid_y, centroid_z)
    centroid_atom.residue = residue
        
    return centroid_atom

def centroid2(residue, ring_atoms):
    centroid = np.mean(ring_atoms, axis = 0)
    centroid_atom = Atom()
    centroid_atom.set_atom_info("RNG", centroid[0], centroid[1], centroid[2])
    centroid_atom.residue = residue
    return centroid_atom

def calc_normal_vector(ring_atoms):
    centroid = np.mean(ring_atoms, axis = 0) # axis=0 -> mean through the columns
    centered_ring_atoms = ring_atoms - centroid # normalizes to origin
    
    # Use singular value decomposition (SVD) to calculate the plane
    _, _, vh = np.linalg.svd(centered_ring_atoms) # vh = V^T
    normal_vector = vh[2]  # The normal vector is the last row of the V^T matrix
    
    return normal_vector