import os
import sys
from timeit import default_timer as timer
import pdb_parser
import superimposers
import sysfunctions
import contacts_fast

def main():
    
    start = timer()

    mode, pdb_files, ref_pdb, atoms_to_be_aligned, rmsd, avd_cutoff = sysfunctions.cl_parse()
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

    parsed_proteins = pdb_parser.parse_pdb(file_list)
    
    ref_distance = None
    ref_protein = parsed_proteins[0].id
    for protein in parsed_proteins:
        print(f"Detecting contacts for {protein.id} against {protein.id}")
        distances = contacts_fast.fast_contacts(protein, protein)
        contacts_fast.show_contacts(distances)
        if ref_distance is None:
            ref_distance = distances
        else:
            avd_list_new, average_avd_new, contact_matches_new = contacts_fast.new_avd(ref_distance, distances, avd_cutoff)
            avd_list, average_avd, contact_matches = contacts_fast.avd(ref_distance, distances, avd_cutoff, avd_list_new)
            if avd_list is not None:
                print(f"Average AVD for {ref_protein} and {protein.id} (new): {average_avd_new}\nNumber of contact matches found: {contact_matches_new}\n")
                print(f"Average AVD for {ref_protein} and {protein.id} (old): {average_avd}\nNumber of contact matches found: {contact_matches}\n")
            else:
                print(f"No contact matches found between {ref_protein} and {protein.id}.\nTry increasing the cutoff value.\n")     
        print("-------------------------------------\n")

    if avd_list_new:
        sorted_list = sorted(avd_list_new, key=lambda x: x[0])
        for contact in sorted_list:
            pass
            #print(contact)

    end = timer()
    print(f"Total time elapsed: {end - start}\n")

if __name__ == "__main__":
    main()
