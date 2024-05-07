import os
from timeit import default_timer as timer
import pdb_parser
import contacts_fast

def main():
    
    start = timer()
    
    folder_path = "./cs1/"
    files_in_folder = [f"{folder_path}{f}" for f in os.listdir(folder_path) if f.endswith(".pdb")]
    
    for file in files_in_folder:
        parsed_protein = pdb_parser.parse_pdb(file)
        print(f"Detecting contacts for {parsed_protein.id}")
        contacts_fast.fast_contacts(parsed_protein)

    end = timer()
    print(f"Total time elapsed: {end - start}\n")

if __name__ == "__main__":
    main()
