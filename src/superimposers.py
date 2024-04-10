
import sysfunctions
import os
import subprocess
import re
import glob
try: 
    import Bio.PDB
except ImportError:
    print("Bio.PDB is not installed!")
    

def tmalign_superimpose(pdb_files, ref_pdb, rmsd):
    """
    Performs superimposition using TMAlign. This is the default mode to superimpose structures.
    Only structures that pass the RMSD threshold parameter are saved.
    
    Args:
        pdb_files (list): List of paths to PDB files to be aligned.
        ref_pdb (str): Path to the reference PDB file.
        rmsd (float): RMSD threshold for accepting alignments.

    Returns:
        .pdb files of the aligned proteins that pass the RMSD threshold. These files are stored in the 'tmalignoutputs' folder.
        list: List of paths to the aligned PDB files.
    """
    
    print("---TMAlign Mode Selected---\n")
    sysfunctions.TMAlign_check()  # check the status of the TMAlign executable
    
    executable = "./TMAlign"
    
    ref_name = os.path.basename(ref_pdb)    
    outputs = [ref_pdb] # starts the output with the reference (samples will be added after their processing)
    output_dir = "tmalignoutputs"
    
    if not os.path.exists(output_dir): # creates output folder if it doesn't exist
        os.mkdir(output_dir)
    print(f"Comparing to: {ref_name}")
    for sample in pdb_files:
        sample_name = os.path.basename(sample)
        if sample_name == ref_name: # exclude same protein comparison
            continue
        output_name = f"{output_dir}/{sample_name.split('.')[0]}_aligned"
        
        # runs TMAlign (current sample vs. reference protein)
        # suppresses the verbosity and captures it to process the RMSD values
        result = subprocess.run([executable, sample, ref_pdb, "-o", output_name], capture_output=True, text=True)
        
        # uses RE to get the RMSD value from the suppressed output
        rmsd_value = re.findall(r'RMSD\s*=\s*([\d.]+)', result.stdout)
        
        # if the RMSD is greater than the cutoff, deletes the generated files
        # this has room to improve! is there a way to just not create the files, instead of deleting them?
        if float(rmsd_value[0]) > rmsd: 
            print(f"RMSD value for file '{sample_name}' too high: {rmsd_value[0]}")
            files_to_delete = glob.glob(os.path.join(output_name + '*'))
            for file in files_to_delete:
                os.remove(file)
        else: # RMSD is lower than the cutoff
            outputs.append(f"{output_name}_rotate.pdb")
            print(f"RMSD value for file '{sample_name}': {rmsd_value[0]}.\t PDB File '{output_name}_rotate.pdb' created!")
            files_to_delete = glob.glob(os.path.join(output_name + '*'))
            # keeps only the .pdb files
            files_to_delete = [file for file in files_to_delete if not file.endswith(".pdb")]
            for file in files_to_delete:
                os.remove(file)
    print("\n-------------------------------------\n")
    return outputs
    

def biopython_superimpose(pdb_files, ref_pdb, atoms_to_be_aligned, rmsd):
    """
    Performs superimposition using BioPython.
    Only structures that pass the RMSD threshold parameter are saved.
    
    Args:
        pdb_files (list): List of paths to PDB files to be aligned.
        ref_pdb (str): Path to the reference PDB file.
        atoms_to_be_aligned (range): range of atoms defined as the start and end id parameters. 
        rmsd (float): RMSD threshold for accepting alignments.

    Returns:
        .pdb files of the aligned proteins that pass the RMSD threshold.
        list: List of paths to the aligned PDB files.
    """
    print("---BioPython Mode Selected---\n")   
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)

    ref_structure = pdb_parser.get_structure("reference", ref_pdb)
    ref_name = os.path.basename(ref_pdb)
    outputs = [ref_name]
    for sample in pdb_files:
        sample_structure = pdb_parser.get_structure("sample", sample)
        sample_name = os.path.basename(sample)

        # Use the first model in the pdb-files for alignment
        ref_model = ref_structure[0]
        sample_model = sample_structure[0]

        # Make a list of the atoms (in the structures) you wish to align.
        # In this case we use CA atoms whose index is in the specified range
        ref_atoms = []
        sample_atoms = []

        for ref_chain, sample_chain in zip(ref_model, sample_model):
            for ref_res, sample_res in zip(ref_chain, sample_chain):
                # Check if residue number ( .get_id() ) is in the lists
                if (ref_res.get_id()[1] in atoms_to_be_aligned and sample_res.get_id()[1] in atoms_to_be_aligned):
                    # Append CA atom to lists
                    ref_atoms.append(ref_res["CA"])
                    sample_atoms.append(sample_res["CA"])

        # Now we initiate the superimposer:
        super_imposer = Bio.PDB.Superimposer()
        super_imposer.set_atoms(ref_atoms, sample_atoms)
        super_imposer.apply(sample_model.get_atoms())

        print(f"RMSD between {ref_name} and {sample_name}: {super_imposer.rms}")  # RMSD
        output_name = f"{sample_name.split('.')[0]}_aligned.pdb"

        if super_imposer.rms < rmsd:
            io = Bio.PDB.PDBIO()
            io.set_structure(sample_structure)
            outputs.append(output_name)
            io.save(output_name)
            print(f"\tFile saved: '{output_name}'!")
        else:
            print(f"\tRMSD too high for file: '{output_name}'!")
    return outputs


def lovoalign_superimpose():
    # https://github.com/m3g/lovoalign/
    pass
