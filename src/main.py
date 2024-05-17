import os
import sys
from timeit import default_timer as timer
import pdb_parser
import superimposers
import sysfunctions
import contacts_fast

import contact_map_plot

def main():
    
    start = timer()

    mode, pdb_files, ref_pdb, atoms_to_be_aligned, rmsd, avd_cutoff, plot, fast = sysfunctions.cl_parse()
    
    try:
        if not os.path.exists(ref_pdb):
            raise FileNotFoundError(f"Reference file not found: {ref_pdb}")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1) # error exiting
        
    if mode == "BioPython": 
        file_list = superimposers.biopython_superimpose(pdb_files, ref_pdb, atoms_to_be_aligned, rmsd)
    elif mode == "TMAlign":
        file_list = superimposers.tmalign_superimpose(pdb_files, ref_pdb, rmsd)
    elif mode == "Single":
        file_list = [ref_pdb]
        
    parsed_proteins = pdb_parser.parse_pdb(file_list)
    
    ref_contacts = None
    ref_protein = parsed_proteins[0].id
    for protein in parsed_proteins:
        print(f"Detecting contacts for {protein.id}")
        contacts = contacts_fast.fast_contacts(protein, fast)
        contacts_fast.show_contacts(contacts)
        
        chain_residues, total_size = protein.full_count()
        
        # don't run both at the same time (they return the same thing, the first just plots as well)
        if plot:
            contact_map_plot.plot_matrix(contacts, chain_residues, total_size)
        #contact_map_plot.contact_matrix(contacts, chain_residues, total_size)
        
        if ref_contacts is None:
            ref_contacts = contacts
        else:
            match_list, average_avd, contact_matches = contacts_fast.avd(ref_contacts, contacts, avd_cutoff)
            if match_list is not None:
                print(f"Average AVD for {ref_protein} and {protein.id} (old): {average_avd}\nNumber of contact matches found: {contact_matches}\n")
            else:
                print(f"No contact matches found between {ref_protein} and {protein.id}.\nTry increasing the cutoff value.\n")     
        print("-------------------------------------\n")

    # if match_list:
    #     sorted_match_list = sorted(match_list, key=lambda x: x.avd)
    #     for match in sorted_match_list:
    #         if match.d3d4:
    #             print(match.avd, match.contact1.print_values(), match.contact2.print_values(), match.d3d4)
    #         else:
    #             print("\t", match.avd, match.contact1.print_values(), match.contact2.print_values(), match.d3d4)

    # if match_list:
    #     cont = 0
    #     for match in match_list:
    #         if match.d3d4:
    #             print("\t",match.avd, match.contact1, match.contact2, match.d3d4)
    #             cont += 1
    #         else:
    #             print(match.avd, match.contact1.get_values(), match.contact2.get_values(), match.d3d4)
    # print(cont)
    
    end = timer()
    print(f"Total time elapsed: {end - start}\n")

if __name__ == "__main__":
    main()
