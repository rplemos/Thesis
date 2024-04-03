import os
import sys
import subprocess
import argparse

def cl_parse():
    """
    Parses command-line arguments for the program execution.

    Raises:
        ValueError: If an invalid mode is provided, start or end arguments are not positive integers, start argument 
            is greater than or equal to the end argument, or the RMSD value is less than 0.
        argparse.ArgumentError: If there is an error parsing command-line arguments.
        Exception: If an unexpected error occurs during argument parsing.

    Returns:
        tuple: A tuple containing the mode of alignment (TMAlign or BioPython), a list of PDB files to be processed, 
            the reference PDB file, a range of atom IDs to be aligned, and the RMSD cutoff value.
    """

    try:
        parser = argparse.ArgumentParser(description='Process PDB files and align atoms.')
        parser.add_argument('-mode', required=False, default='TMAlign', help='Select Program to Run the Alignments. Default = TMAlign')
        parser.add_argument('-pdb', nargs='+', required=True, type=validate_file, help='List of PDB files (at least one required)')
        parser.add_argument('-ref', required=True, type=validate_file, help='Reference PDB file')
        parser.add_argument('-s', '--start_id', type=int, required=False, default=0, help='Start ID for atom alignment. Only for BioPython mode. Default = 0')
        parser.add_argument('-e', '--end_id', type=int, required=False, default=100, help='End ID for atom alignment. Only for BioPython mode. Default = 100')
        parser.add_argument('-rmsd', type=float, required=False, default=999, help='Sets the RMSD cutoff value. Only for TMAlign mode')
        parser.add_argument('-avd', '--avd_cutoff', type=float, required=False, default=1.0, help='Sets the AVD cutoff value. Default = 1.0')

        args = parser.parse_args()

        mode = args.mode
        modes = ["TMAlign", "BioPython"]
        pdb_files = args.pdb
        ref_pdb = args.ref
        start_id = args.start_id
        end_id = args.end_id
        rmsd_value = args.rmsd
        avd_cutoff = args.avd_cutoff
        
        if mode not in modes:
            raise ValueError("Invalid Mode!")

        if start_id < 0 or end_id < 0:
            raise ValueError("Start and end arguments must be positive integers.")

        if start_id >= end_id:
            raise ValueError("Start argument must be lower than end argument.")
        
        if rmsd_value < 0:
            raise ValueError("RMSD value must be higher than 0.")

        # add option to end_id to be the length of the reference file
        # if end_id > 100:
        #     raise ValueError("End argument must be lower than 100.")
        
        atoms_to_be_aligned = range(start_id, end_id + 1)
        
    except argparse.ArgumentError as e:
        print(f"Argument Error: {str(e)}")
        sys.exit(1)

    except ValueError as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")
        sys.exit(1)
        
    return mode, pdb_files, ref_pdb, atoms_to_be_aligned, rmsd_value, avd_cutoff
        
def validate_file(value):
    """
    Validate if the files passed as the 'pdb_files' and 'reference_pdb' are .pdb.

    Args:
        value (str): The file path to be validated.

    Raises:
        argparse.ArgumentTypeError: Raised if the file is not a valid PDB file.

    Returns:
        str: The validated file path.
    """
    if not value.endswith('.pdb'):
        raise argparse.ArgumentTypeError(f"{value} is not a valid PDB file. File must end with '.pdb'")
    return value

def TMAlign_check():
    """
    Check the availability of TMAlign executable and compile it if necessary.

    This function checks if the TMAlign executable is present in the current directory. If it is not found, 
    it attempts to compile it from the TMAlign.cpp source file. If compilation is successful, the TMAlign 
    executable is created in the current directory.

    Raises:
        SystemExit: If TMAlign is not present in the folder or compilation fails.
    """
    print("Checking TMAlign...")
    tmalign_exe = "TMAlign"
    tmalign_cpp = "TMAlign.cpp"

    if os.path.exists(tmalign_exe):
        print ("TMAlign is ready to go!\n")
    elif os.path.exists(tmalign_cpp):
        print("Compiling TMAlign...")
        try:
            compile_command = ["g++", "-static", "-O3", "-ffast-math", "-lm", "-o", tmalign_exe, tmalign_cpp]
            subprocess.run(compile_command, check=True)
            print ("Compilation Successful!\n")
        except subprocess.CalledProcessError as e:
            print (f"Compilation failed with error code {e.returncode}.")
            sys.exit(1)
    else:
        print ("TMAlign not present in the folder!")
        sys.exit(1)
