import os
import sys
from timeit import default_timer as timer
import pdb_parser
import superimposers
import sysfunctions
import contacts_fast
import contact_map_plot

import concurrent.futures
import multiprocessing
import pickle
import resource

def main():

    start = timer()

    mode, pdb_files, ref_pdb, atoms_to_be_aligned, rmsd, avd_cutoff, plot, fast, show_contacts = sysfunctions.cl_parse()

    try:
        if not os.path.exists(ref_pdb):
            raise FileNotFoundError(f"Reference file not found: {ref_pdb}")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1) # error exiting
    
    if mode == "BioPython": 
        file_list = superimposers.biopython_superimpose(pdb_files, ref_pdb, atoms_to_be_aligned, rmsd)
    elif mode == "TMAlign":
        file_list = superimposers.prepare_tmalign_superimpose(pdb_files, ref_pdb, rmsd)
    elif mode == "Single":
        file_list = [ref_pdb]
    elif mode == "Benchmark":
        file_list = pdb_files
    
    if mode != "Benchmark" and mode != "Single":
        ref_protein, ref_contacts, ref_time = process_file(ref_pdb, fast)
        print(f"Reference file processing time: {ref_time}")
        
    manager = multiprocessing.Manager()
    maximum_distances = manager.dict()

    progress_list = manager.list()
    lock = manager.Lock()
    batch_size = 100
    batch_results = []
    batch_number = 0
    
    total_files = len(file_list)
            
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        future_to_file = {executor.submit(process_file, file_path, fast, maximum_distances, lock, progress_list): file_path for file_path in file_list}
        
        for i, future in enumerate(concurrent.futures.as_completed(future_to_file)):
            
            try:
                protein, contacts, time, maximum_distances = future.result()
                
                #batch_results.append(future.result())

                # Save and clear results every batch_size files
                # if (i + 1) % batch_size == 0:
                #     save_results(batch_results, batch_number)
                #     batch_results.clear()  # Clear the list to free memory
                #     gc.collect()  # Trigger garbage collection to free memory
                #     batch_number += 1
                                    
                if protein.title == "DNA/RNA":
                    print(f"File '{protein.id}' is not a Protein! Ignoring File\n")
                    print("-----------------------------------------------------\n")
                    continue
                
                # print(f"Detecting contacts for {protein.id}")
                # chains = [chain.id for chain in protein.chains]
                # print(f"Number of chains: {len(chains)}")
                # print(f"Chains to be analyzed: {chains}")
                # print(f"Protein size:\n\tFull (includes gaps): {protein.full_count()[1]}\n\tTrue (number of residues in the PDB file): {protein.true_count()}\n")
                # print(f"Number of contacts: {len(contacts)}")
                # print(f"Contact Detection - Time elapsed: {time}\n")
                
                print_memory_usage()
                chains = [chain.id for chain in protein.chains]
                print(protein.id, len(chains), chains)
                print(protein.full_count()[1], protein.true_count(), len(contacts), f"{time:.4f}")
                
                if show_contacts:
                    contacts_fast.show_contacts(contacts)
                    
                chain_residues, total_size = protein.full_count()
                true_size = protein.true_count()
                
                if plot:
                    contact_map_plot.plot_matrix(contacts, chain_residues, total_size)

                if mode != "Benchmark" and mode != "Single":            
                    match_list, average_avd, contact_matches = contacts_fast.avd(ref_contacts, contacts, avd_cutoff)
                    if match_list is not None:
                        print(f"Average AVD for {ref_protein.id} and {protein.id} (old): {average_avd}\nNumber of contact matches found: {contact_matches}\n")
                    else:
                        print(f"No contact matches found between {ref_protein.id} and {protein.id}.\nTry increasing the cutoff value.\n")     
                        
                print("-----------------------------------------------------\n")
                
                with lock:
                    processed_files = len(progress_list)
                    #print(f"Progress: {processed_files}/{total_files}")
                    print(processed_files,"/",total_files)
                
                del future_to_file
                                
            except Exception as e:
                #pass
                print(f"Error: {e}")
    
    end = timer()
    print(f"Total time elapsed: {end - start}\n")
    #print(maximum_distances)       

def process_file(file_path, fast, maximum_distances, lock, progress_list):

    # Update progress list
    with lock:
        progress_list.append(1)
    try:
        if file_path.endswith(".pdb"):
            parsed_data = pdb_parser.parse_pdb(file_path)
        else:
            parsed_data = pdb_parser.parse_pdbx(file_path)
        contacts, time, maximum_distances = contacts_fast.fast_contacts(parsed_data, fast, maximum_distances)
        return parsed_data, contacts, time, maximum_distances
    
    except KeyError as e:
        # Handle KeyError or any other exception
        return (file_path, None, e)  # Return tuple with file_path, None result, and exception
    
def save_results(batch_results, batch_number):
    # Save batch results to disk (e.g., as a file or database)
    with open(f"batch_{batch_number}.pkl", "wb") as f:
        pickle.dump(batch_results, f)

def print_memory_usage():
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000.0  # Convert to MB
    print(f"Memory usage: {usage} MB")

    # print("[ Top 10 memory-consuming lines ]")
    # for stat in top_stats[:10]:
    #     print(stat)
    
    

    # # ascending_distances = sorted(maximum_distances.items(), key=lambda x:x[1][0]) 
    # # for item in ascending_distances:
    # #     print(item)
    # # print(len(ascending_distances))    

    # # if match_list:
    # #     sorted_match_list = sorted(match_list, key=lambda x: x.avd)
    # #     for match in sorted_match_list:
    # #         if match.d3d4:
    # #             print(match.avd, match.contact1.print_values(), match.contact2.print_values(), match.d3d4)
    # #         else:
    # #             print("\t", match.avd, match.contact1.print_values(), match.contact2.print_values(), match.d3d4)
    # # if match_list:
    # #     cont = 0
    # #     for match in match_list:
    # #         if match.d3d4:
    # #             print("\t",match.avd, match.contact1, match.contact2, match.d3d4)
    # #             cont += 1
    # #         else:
    # #             print(match.avd, match.contact1.get_values(), match.contact2.get_values(), match.d3d4)
    # # print(cont)
    
    # end = timer()
    # print(f"Total time elapsed: {end - start}\n")
    
if __name__ == "__main__":
    main()