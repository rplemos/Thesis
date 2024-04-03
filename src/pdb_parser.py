from classes import Protein, Chain, Residue, Atom
import sys
import os

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
                    current_protein.set_id(line[62:])
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
                        current_chain.set_id(chain_id)
                        current_protein.chains.append(current_chain)

                    if current_residue is None or current_residue.resnum != resnum:  # new residue
                        current_residue = Residue()
                        current_residue.set_res_info(resnum, resname)
                        current_residue.chain = current_chain  # Set the parent chain
                        current_chain.residues.append(current_residue)

                    atom = Atom()  # each line is a new atom anyway
                    
                    atomname = line[12:16].replace(" ", "")
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    
                    if atomname == "OXT": # OXT is the C-terminal Oxygen atom. However, it exhibits the same properties of any Oxygen
                        atomname = atomname.replace("OXT","O")
                        
                    atom.set_atom_info(atomname, x, y, z)
                    atom.residue = current_residue  # Set the parent residue
                    current_residue.atoms.append(atom)

                elif line.startswith("END"):  # finishes and resets everything for the new protein
                    proteins.append(current_protein)
                    current_protein = None
                    current_chain = None
                    current_residue = None
    return proteins
