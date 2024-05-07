from classes import Protein, Chain, Residue, Atom
from numpy import mean, array
from numpy.linalg import svd
  

stacking = {
    'HIS':[10, 'CG','ND1','CE1','NE2','CD2'],
    'PHE':[11, 'CG','CD1','CE1','CZ','CE2','CD2'],
    'TRP':[14, 'CG','CD1','NE1','CE2','CZ2','CH2','CZ3','CE3','CD2'],
    'TYR':[12, 'CG','CD1','CE1','CZ','CE2','CD2'],
}


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
                        residues = []
                        current_chain = Chain(chain_id, residues)
                        current_protein.chains.append(current_chain)

                    if current_residue is None or current_residue.resnum != resnum:  # new residue
                        atoms = []
                        current_residue = Residue(resnum, resname, atoms, current_chain, False, None)
                        current_chain.residues.append(current_residue)
                                                                    
                    atomname = line[12:16].replace(" ", "")
                    if atomname == "OXT": # OXT is the C-terminal Oxygen atom. However, it exhibits the same properties of any Oxygen
                        atomname = atomname.replace("OXT","O")                
                    
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    occupancy = float(line[55:60])
                    
                    if occupancy >= 0: # ignores low quality atoms (arbitrary value!)
                        atom = Atom(atomname, x, y, z, occupancy, current_residue) # creates atom
                        current_residue.atoms.append(atom)
                    else:
                        print(atomname, occupancy, current_residue.resname, current_residue.resnum, current_chain.id, current_protein.id)
                    
                    # CHECKING FOR AROMATICS
                    if current_residue.resname in stacking:
                        allowed = stacking[current_residue.resname][1:]
                        all_atoms_have_occupancy_one = all(atom.occupancy == 1 for atom in current_residue.atoms if atom.atomname in allowed)
                        
                        # if ring has only one conformation and the residue is complete (all atoms populated)
                        if all_atoms_have_occupancy_one and len(current_residue.atoms) == stacking[current_residue.resname][0]:
                            ring_atoms = array([[atom.x, atom.y, atom.z] for atom in current_residue.atoms[5:]]) # ignores [N, CA, C, O] and [RNG] atoms
                            centroid_atom = centroid(current_residue, ring_atoms)
                            current_residue.atoms.append(centroid_atom)
                            current_residue.ring = True

                            normal_vector = calc_normal_vector(ring_atoms)
                            current_residue.normal_vector = normal_vector

                elif line.startswith("END"):  # finishes and resets everything for the new protein
                    
                    # Handling cases where there is no ID
                    if current_protein.id is None:
                        id = str(file).split("/")[-1]
                        id = id.split(".")[0]
                        current_protein.id = id  
                                         
                    proteins.append(current_protein)

                    current_protein = None
                    current_chain = None
                    current_residue = None
                
    return proteins


def centroid(residue, ring_atoms):
    centroid = mean(ring_atoms, axis = 0)
    centroid_atom = Atom("RNG", centroid[0], centroid[1], centroid[2], 1, residue)
    
    return centroid_atom


def calc_normal_vector(ring_atoms):
    centroid = mean(ring_atoms, axis = 0) # axis=0 -> mean through the columns
    centered_ring_atoms = ring_atoms - centroid # normalizes to origin
    
    # Use singular value decomposition (SVD) to calculate the plane
    _, _, vh = svd(centered_ring_atoms) # vh = V^T
    normal_vector = vh[2]  # The normal vector is the last row of the V^T matrix
    
    return normal_vector