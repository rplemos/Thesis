import pdb_parser
import superimposers
import sysfunctions
import contacts_fast

def main():
    mode, pdb_files, ref_pdb, atoms_to_be_aligned, rmsd = sysfunctions.cl_parse()
    if mode == "BioPython": 
        file_list = superimposers.biopython_superimpose(pdb_files, ref_pdb, atoms_to_be_aligned, rmsd)
    elif mode == "TMAlign":
        file_list = superimposers.tmalign_superimpose(pdb_files, ref_pdb, rmsd)
    
    parsed_proteins = pdb_parser.parse_pdb(file_list)

    distances = contacts_fast.fast_contacts(parsed_proteins[0], parsed_proteins[1])
    contacts_fast.show_contacts(distances)
    
    # for protein in parsed_proteins:
    #     for chain in protein.get_chains():
    #         for residue in chain.get_residues():
    #             for atom1 in residue.get_atoms():
    #                 for atom2 in residue.get_atoms():
    #                     if atom1 != atom2:
    #                         print(atom1.atomname, atom2.atomname)
    #                 #contacts.euclid(atom)
    #                 #info_list = [protein.id, chain.id, residue.resnum, residue.resname, atom.atomname, atom.x, atom.y, atom.z]
    #                 #print(info_list)
    #     print(protein.id, protein.title)
    
    # a = {}
    
    # for protein in parsed_proteins:
    #     for chain in protein.get_chains():
    #         for residue in chain.get_residues():
    #             for atom in residue.get_atoms():
    #                 a[f"{residue.resname}:{atom.atomname}"] = [atom.x, atom.y, atom.z]
                    
    # for key,value in a.items():
    #     print(key, value)            
    
if __name__ == "__main__":
    main()