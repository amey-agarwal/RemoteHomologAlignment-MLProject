import os, sys, warnings, argparse, subprocess, glob
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio import SeqIO
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import PDBIO
################################################
# Ensure these executables are in your system's PATH
# or provide the full absolute path.
foldseek = "foldseek"
mafft = "mafft"
clustalw2 = "clustalw2"

# The Foldseek matrix file MUST be in the same directory as this script
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
fs_mat_file = os.path.join(SCRIPT_DIR, "fs_mat.dat")
################################################

def validate_input_directory(pdb_dir):
    """
    Checks if the provided directory exists and contains at least two
    structure files (.pdb or .cif).
    ... (rest of function is unchanged) ...
    """
    print(f"--- 1. Validating Input Directory ---")

    # 1. Check if the directory exists and is a directory
    if not os.path.isdir(pdb_dir):
        print(f"Error: Directory not found or is not a directory: {pdb_dir}")
        sys.exit(1) # Exit with a non-zero status code

    # 2. Find all .pdb and .cif files
    structure_files = []
    for filename in os.listdir(pdb_dir):
        # Ignore hidden files like on mac like .DS_Store
        if not filename.startswith('.'):
            if filename.endswith(".pdb") or filename.endswith(".cif"):
                structure_files.append(os.path.join(pdb_dir, filename))

    # 3. Check if we have at least 2 files
    file_count = len(structure_files)
    if file_count < 2:
        print(f"Error: Found only {file_count} structure files (.pdb or .cif) in '{pdb_dir}'.")
        print("At least 2 structure files are required to build an alignment.")
        sys.exit(1)

    # 4. Success
    print(f"  Validation successful: Found {file_count} structure files.")
    print("-------------------------------------\n")
    return structure_files

def check_chain_counts(structure_file_list):
    """
    Iterates through a list of structure files and warns if any contain
    more than one chain in their first model.
    """
    print("--- 2. Checking Chain Counts ---")
    
    # Suppress Biopython' construct warnings 
    warnings.simplefilter('ignore', PDBConstructionWarning) 
    
    pdb_parser = PDBParser(QUIET=True)
    cif_parser = MMCIFParser(QUIET=True)
    
    found_multi_chain = False

    for file_path in structure_file_list:
        filename = os.path.basename(file_path)
        structure = None
        
        try:
            if file_path.endswith(".pdb"):
                structure = pdb_parser.get_structure(filename, file_path)
            elif file_path.endswith(".cif"):
                structure = cif_parser.get_structure(filename, file_path)
            
            if structure:
                if len(list(structure.get_models())) == 0:
                    print(f"  Warning: No models found in '{filename}'. Skipping chain check.")
                    continue
                
                model = structure[0]
                chains = list(model.get_chains())
                
                if len(chains) > 1:
                    chain_ids = [chain.id for chain in chains]
                    print(f"  Warning: '{filename}' contains {len(chains)} chains (IDs: {', '.join(chain_ids)}).")
                elif len(chains) == 0:
                     print(f"  Warning: No chains found in '{filename}'. This file may cause errors.")
            
        except Exception as e:
            print(f"  Warning: Could not parse '{filename}' to check chains. Error: {e}")
    
    if not found_multi_chain:
        print("  All structure files appear to have one chain. Nicely done.")
    
    print("----------------------------------\n")

def get_3di_fasta(pdb_directory,output_file):
    """
    Runs Foldseek to create a 3Di sequence FASTA file from a directory
    of PDB/CIF files. Cleans up temporary files.
    """
    print("--- 3. Running Foldseek to Get 3Di ---")
    
    # Define file/db names
    temp_db_name = "FSDB_temp"
    output_fasta = output_file
    
    abs_pdb_dir = os.path.abspath(pdb_directory)

    try:
        # 1. Create Foldseek database
        print(f"  Creating Foldseek DB... (using '{abs_pdb_dir}')")
        cmd_createdb = [foldseek, "createdb", abs_pdb_dir, temp_db_name]
        subprocess.run(
            cmd_createdb, capture_output=True, text=True, check=True, encoding='utf-8'
        )

        # 2. Create 3Di lookup database
        print("  Creating 3Di lookup DB...")
        cmd_lndb = [foldseek, "lndb", f"{temp_db_name}_h", f"{temp_db_name}_ss_h"]
        subprocess.run(
            cmd_lndb, capture_output=True, text=True, check=True, encoding='utf-8'
        )

        # 3. Convert 3Di sequences to a FASTA file
        print(f"  Converting 3Di sequences to '{output_fasta}'...")
        cmd_convert = [foldseek, "convert2fasta", f"{temp_db_name}_ss", output_fasta]
        subprocess.run(
            cmd_convert, capture_output=True, text=True, check=True, encoding='utf-8'
        )

    except subprocess.CalledProcessError as e:
        print(f"\nError: Foldseek command failed ('{' '.join(e.cmd)}').")
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        print("Please ensure 'foldseek' is installed and accessible in your system's PATH.")
        sys.exit(1)
    except FileNotFoundError:
        print(f"\nError: 'foldseek' command not found.")
        print("Please ensure 'foldseek' is installed and accessible in your system's PATH.")
        sys.exit(1)
    finally:
        # 4. Clean up all temporary FSDB files
        print("  Cleaning up temporary DB files (FSDB_temp*)...")
        for f in glob.glob(f"{temp_db_name}*"):
            try:
                os.remove(f)
            except OSError as e:
                print(f"Warning: could not remove temp file {f}: {e}")

    # 5. Final check for the output FASTA file
    if not os.path.exists(output_fasta) or os.path.getsize(output_fasta) == 0:
        print(f"Error: Foldseek ran, but the output file '{output_fasta}' is missing or empty.")
        sys.exit(1)

    print(f"  Successfully created '{output_fasta}'.")
    print("----------------------------------------\n")
    return output_fasta

def run_alignment(fasta_file, aligner_method):
    """
    Runs MAFFT or ClustalW2 on the 3Di FASTA file using the Foldseek matrix.

    Args:
        fasta_file (str): Path to the input 3Di.fa file.
        aligner_method (str): 'mafft' or 'clustalw'.

    Returns:
        str: Path to the output alignment file ('3Di.aln').
        
    Raises:
        SystemExit: If the aligner command fails.
    """
    print(f"--- 4. Running Alignment via {aligner_method.upper()} ---")
    output_alignment = "pdb_to_3di_3Di.fa"
    
    cmd = []
    
    try:
        if aligner_method == 'mafft':
            # MAFFT command: mafft --clustalout --aamatrix <matrix> <input> > <output>
            # We must run this via shell=True to handle the '>' redirection
            cmd_string = f"{mafft} --aamatrix {fs_mat_file} --globalpair --maxiterate 1000 {fasta_file} > {output_alignment}"
            print(f"  Running command: {cmd_string}")
            subprocess.run(
                cmd_string, shell=True, check=True, capture_output=True, text=True, encoding='utf-8'
            )
        
        else: # Default to clustalw
            # ClustalW command: clustalw2 -infile=<input> -outfile=<output> -pwmatrix=<matrix> -matrix=<matrix>
            cmd = [
                clustalw2,
                f"-infile={fasta_file}",
                f"-outfile={output_alignment}",
                f"-pwmatrix={fs_mat_file}",
                f"-matrix={fs_mat_file}"
            ]
            print(f"  Running command: {' '.join(cmd)}")
            subprocess.run(
                cmd, capture_output=True, text=True, check=True, encoding='utf-8'
            )

    except subprocess.CalledProcessError as e:
        # Handle errors from both shell=True and shell=False commands
        cmd_str = e.cmd if isinstance(e.cmd, str) else ' '.join(e.cmd)
        print(f"\nError: Aligner command failed ('{cmd_str}').")
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        print(f"Please ensure '{aligner_method}' is installed and accessible in your system's PATH.")
        sys.exit(1)
    except FileNotFoundError:
        print(f"\nError: '{aligner_method}' command not found.")
        print(f"Please ensure '{aligner_method}' is installed and accessible in your system's PATH.")
        sys.exit(1)

    # Final check for the output alignment file
    if not os.path.exists(output_alignment) or os.path.getsize(output_alignment) == 0:
        print(f"Error: {aligner_method} ran, but the output file '{output_alignment}' is missing or empty.")
        sys.exit(1)

    print(f"  Successfully created '{output_alignment}'.")
    print("----------------------------------------------\n")
    return output_alignment
    

    
#####################################################################################
# confidence stuff
def read_pwmatrix(pwmatrix_file):
    """
    Reads and parses the PWMATRIX file into a dictionary.

    :param pwmatrix_file: Path to the PWMATRIX file.
    :return: Dictionary of substitution scores {('A', 'C'): -3, ...}
    """
    pwmatrix = {}
    try:
        with open(pwmatrix_file, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: Could not find matrix file: {pwmatrix_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading matrix file: {e}")
        sys.exit(1)

    # Extract column headers (Residue labels)
    headers = lines[0].split()
    
    # Parse substitution values
    for line in lines[1:]:  # Skip the first row (header)
        values = line.split()
        if not values: # Skip empty lines
            continue
        row_residue = values[0]  # The first column is the row residue

        # Map each column residue with its corresponding score
        for col_index, col_residue in enumerate(headers):
            try:
                pwmatrix[(row_residue, col_residue)] = int(values[col_index + 1])
            except (IndexError, ValueError):
                print(f"Warning: Skipping malformed line in matrix file: {line.strip()}")
                continue
    return pwmatrix

def compute_confidence_tracks(msa_matrix, pwmatrix_file, valid_columns, block_threshold=2.0):
    """
    Computes per-column and global confidence scores from an MSA matrix.
    Ported from the Structome-AlignViewer backend.
    
    Returns:
        per_column_track_data (list): Saguaro-formatted list of dicts (for console)
        block_data (list): Saguaro-formatted blocks (unused)
        global_confidence_score (float): The final score, averaged from valid columns
        normalized_scores (list): The raw list of per-column normalized scores
    """
    pwmatrix = read_pwmatrix(pwmatrix_file)
    num_seqs = len(msa_matrix)
    alignment_length = len(msa_matrix[0])

    # Compute per-column confidence scores
    per_column_scores = []
    for col in range(alignment_length):
        col_pairs = [(msa_matrix[row][col], msa_matrix[row + 1][col]) for row in range(num_seqs - 1)]
        col_score = sum(pwmatrix.get(pair, -3.0) for pair in col_pairs if "-" not in pair)
        per_column_scores.append(col_score)

    # Normalize scores (Min-Max Scaling to 0-1)
    min_score, max_score = min(per_column_scores), max(per_column_scores)
    normalized_scores = [(s - min_score) / (max_score - min_score) if max_score != min_score else 0 for s in per_column_scores]

    # Convert per-column scores into JSON-like list (for Saguaro/console)
    per_column_track_data = [{"begin": col + 1, "end": col + 1, "value": round(score, 3)} for col, score in enumerate(normalized_scores)]

    # Compute confidence blocks (originally for Saguaro)
    block_data = []
    if valid_columns:
        start, prev_score = valid_columns[0], normalized_scores[valid_columns[0]]
        for i in range(1, len(valid_columns)):
            col = valid_columns[i]
            col_score = normalized_scores[col]
            if abs(col_score - prev_score) > block_threshold:
                block_data.append({"begin": start + 1, "end": col, "value": round(prev_score, 3)})
                start = col
            prev_score = col_score
        block_data.append({"begin": start + 1, "end": valid_columns[-1] + 1, "value": round(prev_score, 3)})
    
    # Calculate global score *only from valid (non-gappy) columns*
    global_confidence_score = sum(normalized_scores[col] for col in valid_columns) / len(valid_columns) if valid_columns else 0
    
    # --- MODIFIED: Return the raw normalized_scores list as well ---
    return per_column_track_data, block_data, round(global_confidence_score, 3), normalized_scores

def calculate_confidence_scores(input_file, mat_file, threshold=0.5):
    """
    Main analysis function (adapted from trim_msa_by_gblocks).
    - Reads the alignment.
    - Identifies valid/dropped columns (>= 50% non-gap).
    - Calls compute_confidence_tracks to get all scores.
    - Writes the per-column confidence file.
    - Creates the alignment-to-residue lookup dictionary.
    """
    print("  Reading alignment for scoring...")
    try:
        records = list(SeqIO.parse(input_file, "clustal"))
    except Exception as e:
        print(f"Error: Could not parse alignment file '{input_file}': {e}")
        sys.exit(1)
        
    if not records:
        print(f"[ERROR] No sequences found in alignment: {input_file}")
        sys.exit(1)

    msa_matrix = [list(str(record.seq)) for record in records]
    num_seqs = len(records)
    alignment_length = len(msa_matrix[0])

    # Identify valid columns (>=50% non-gap)
    valid_columns = []
    dropped_columns = []
    for col in range(alignment_length):
        col_data = [msa_matrix[row][col] for row in range(num_seqs)]
        non_gap_count = sum(1 for char in col_data if char != "-")
        if non_gap_count / num_seqs >= threshold:
            valid_columns.append(col)
        else:
            dropped_columns.append(col)
            
    if not valid_columns:
        print("[ERROR] No columns met the 50% threshold. All columns are gaps. Cannot proceed.")
        sys.exit(1)
    
    # --- Call Confidence Calculation ---
    print("  Calculating confidence scores...")
    pcc, _, global_score, norm_scores = compute_confidence_tracks(
        msa_matrix, mat_file, valid_columns
    )
    
    # --- Write the per-column confidence file ---
    output_conf_file = "confidence_per_column.txt"
    print(f"  Writing per-column confidence to '{output_conf_file}'...")
    try:
        with open(output_conf_file, "w") as f_out:
            f_out.write("Column_Index\tColumn_Chars\tNormalized_Confidence\n")
            for col_idx in range(alignment_length):
                index_str = str(col_idx + 1)
                col_chars = "".join([msa_matrix[row][col_idx] for row in range(num_seqs)])
                col_score = norm_scores[col_idx]
                f_out.write(f"{index_str}\t{col_chars}\t{col_score:.4f}\n")
                
    except Exception as e:
        print(f"Error: Failed to write confidence file '{output_conf_file}': {e}")
        sys.exit(1)

    print(f"  Successfully saved confidence data: {output_conf_file}")

    # --- Generate lookup dict for B-factor mapping ---
    print("  Generating alignment-to-residue mapping...")
    msa_dict = {
        os.path.basename(rec.id).replace(".cif", "").replace(".pdb", ""): str(rec.seq) 
        for rec in records
    }
    lookup_dict = generate_alignment_lookup(msa_dict)

    # Return `valid_columns` instead of the undefined `valid_cols`
    return (
        global_score, alignment_length, valid_columns, dropped_columns, 
        output_conf_file, lookup_dict, norm_scores
    )
    
def z_score_scop(confidence_value, aligner='clustalw2'):
    """
    Compute Z-score for SCOP AvgConfidence.
    Stats are based on the aligner used.
    """
    if aligner == 'mafft':
        # MAFFT-derived statistics
        mean = 0.4802
        std = 0.0628
    else:
        # Original ClustalW statistics 
        mean = 0.4932
        std = 0.0585
    return (confidence_value - mean) / std

def z_score_cath(confidence_value, aligner='clustalw2'):
    """
    Compute Z-score for CATH AvgConfidence.
    Stats are based on the aligner used.
    """
    if aligner == 'mafft':
        # MAFFT-derived statistics
        mean = 0.4783
        std = 0.0686
    else:
        # Original ClustalW statistics 
        mean = 0.5006
        std = 0.0656
    return (confidence_value - mean) / std
    
def clean_alignment_file(input_file, output_file):
    """
    Cleans a Clustal alignment file for downstream use.
    - Replaces consensus markers ('*', ':', '.') with spaces.
    - Keeps the Clustal header, all sequence lines, and all blank lines
      to preserve the file's block structure.

    Args:
        input_file (str): Path to the raw alignment file.
        output_file (str): Path to write the cleaned alignment.

    Returns:
        str: The path to the cleaned output file.
        
    Raises:
        SystemExit: If cleaning fails.
    """
    print(f"  Cleaning alignment file: {input_file} -> {output_file}")
    try:
        with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
            for line in f_in:
                # Check if it's a consensus line
                # These lines start with a space but are not just a blank line
                if line.startswith(" ") and line.strip() != "":
                    # It's a consensus line. Replace markers with spaces.
                    clean_line = line.replace("*", " ").replace(":", " ").replace(".", " ")
                    f_out.write(clean_line)
                else:
                    # It's a header, a sequence line, or a blank line between blocks.
                    # Keep it as-is to preserve the file structure.
                    f_out.write(line)

    except Exception as e:
        print(f"\nError: Failed to clean alignment file: {e}")
        sys.exit(1)
        
    return output_file
    
#########################################################################################################
def generate_alignment_lookup(msa_dict):
    """
    Given an MSA dictionary, create an index lookup for each sequence
    that maps alignment column index to a 1-based residue index.
    
    Ported from the Structome-AlignViewer backend.

    :param msa_dict: Dictionary { "seq_id": "aligned_sequence" }
    :return: Dictionary { "seq_id": [index_lookup_array] }
    """
    lookup_dict = {}

    for seq_id, aligned_seq in msa_dict.items():
        lookup_array = []
        original_index = 0  # Tracks the 1-based position in the original sequence

        for char in aligned_seq:
            if char == "-":
                lookup_array.append(None)  # Gaps get None
            else:
                original_index += 1  # Increment original sequence index
                lookup_array.append(original_index)  # Store the mapped index

        lookup_dict[seq_id] = lookup_array

    return lookup_dict

def save_confidence_to_bfactor(pdb_directory, lookup_dict, norm_scores):
    """
    For each structure, maps per-column confidence scores to residue
    B-factors and saves a new structure file.
    
    Adapted from Structome-AlignViewer backend.

    Args:
        pdb_directory (str): The directory containing original .pdb/.cif files.
        lookup_dict (dict): The map of {seq_id: [aln_col_to_res_idx_map]}.
        norm_scores (list): The raw list of per-column normalized confidence scores.
    """
    print("--- 8. Saving PDBs with B-Factor Confidence ---")
    
    output_dir = "PDBs_with_confidence"
    os.makedirs(output_dir, exist_ok=True)
    print(f"  Saving new files to: '{output_dir}'")

    pdb_parser = PDBParser(QUIET=True)
    cif_parser = MMCIFParser(QUIET=True)
    io = PDBIO()
    
    warnings.simplefilter('ignore', PDBConstructionWarning)

    for structure_id, position_map in lookup_dict.items():
        # Find the original file (could be .pdb or .cif)
        base_filename = os.path.splitext(structure_id)[0]
        # We must find the *original* file from the input dir
        orig_file_path = None
        for ext in [".pdb", ".cif"]:
            test_path = os.path.join(pdb_directory, f"{base_filename}{ext}")
            if os.path.exists(test_path):
                orig_file_path = test_path
                break
        
        if not orig_file_path:
            print(f"  Warning: Could not find original structure file for '{structure_id}'. Skipping.")
            continue
        
        # Parse the structure
        try:
            if orig_file_path.endswith(".pdb"):
                structure = pdb_parser.get_structure(structure_id, orig_file_path)
            else:
                structure = cif_parser.get_structure(structure_id, orig_file_path)
        except Exception as e:
            print(f"  Warning: Could not parse '{orig_file_path}'. Skipping. Error: {e}")
            continue

        # --- Create a map of {1-based_residue_index: confidence} ---
        confidence_per_residue = {}
        for aln_index, res_1_based_idx in enumerate(position_map):
            if res_1_based_idx is not None:  # Not a gap
                confidence = norm_scores[aln_index]
                confidence_per_residue[res_1_based_idx] = confidence
        
        # --- Apply scores to B-factor field ---
        for model in structure:
            for chain in model:
                res_counter = 0 # This is the 1-based residue index
                for residue in chain:
                    if not is_aa(residue):
                        continue
                    
                    res_counter += 1 # Increment for each valid AA residue
                    
                    # Check if this residue (by its 1-based index) has a score
                    if res_counter in confidence_per_residue:
                        conf_score = confidence_per_residue[res_counter]
                        for atom in residue:
                            atom.set_bfactor(conf_score * 100) # Scale 0-1 to 0-100
                    else:
                        # This residue was gapped in the alignment
                        for atom in residue:
                            atom.set_bfactor(0.0)
        
        # Save the new PDB file
        output_file = os.path.join(output_dir, f"{base_filename}_conf.pdb")
        io.set_structure(structure)
        io.save(output_file)

    print(f"  Successfully saved new PDB files.")
    print("------------------------------------------------\n")
    return output_dir    
    
def save_alignment_index_to_bfactor(pdb_directory, lookup_dict):
    """
    For each structure, maps the 1-based alignment column index to the
    residue B-factors and saves a new structure file.

    Args:
        pdb_directory (str): The directory containing original .pdb/.cif files.
        lookup_dict (dict): The map of {seq_id: [aln_col_to_res_idx_map]}.
    """
    print("--- 9. Saving PDBs with Alignment-Index B-Factors ---")
    
    output_dir = "PDBs_with_aln_index"
    os.makedirs(output_dir, exist_ok=True)
    print(f"  Saving new files to: '{output_dir}'")

    pdb_parser = PDBParser(QUIET=True)
    cif_parser = MMCIFParser(QUIET=True)
    io = PDBIO()
    
    warnings.simplefilter('ignore', PDBConstructionWarning)

    for structure_id, position_map in lookup_dict.items():
        # Find the original file
        base_filename = os.path.splitext(structure_id)[0]
        orig_file_path = None
        for ext in [".pdb", ".cif"]:
            test_path = os.path.join(pdb_directory, f"{base_filename}{ext}")
            if os.path.exists(test_path):
                orig_file_path = test_path
                break
        
        if not orig_file_path:
            print(f"  Warning: Could not find original structure file for '{structure_id}'. Skipping.")
            continue
        
        # Parse the structure
        try:
            if orig_file_path.endswith(".pdb"):
                structure = pdb_parser.get_structure(structure_id, orig_file_path)
            else:
                structure = cif_parser.get_structure(structure_id, orig_file_path)
        except Exception as e:
            print(f"  Warning: Could not parse '{orig_file_path}'. Skipping. Error: {e}")
            continue

        # --- Create a map of {1-based_residue_index: 1-based_alignment_index} ---
        index_per_residue = {}
        for aln_index_0_based, res_1_based_idx in enumerate(position_map):
            if res_1_based_idx is not None:  # Not a gap
                alignment_index_1_based = aln_index_0_based + 1
                index_per_residue[res_1_based_idx] = alignment_index_1_based
        
        # --- Apply scores to B-factor field ---
        for model in structure:
            for chain in model:
                res_counter = 0 # This is our 1-based residue index
                for residue in chain:
                    if not is_aa(residue):
                        continue
                    
                    res_counter += 1 # Increment for each valid AA residue
                    
                    # Check if this residue (by its 1-based index) has an index
                    if res_counter in index_per_residue:
                        aln_idx = index_per_residue[res_counter]
                        for atom in residue:
                            atom.set_bfactor(aln_idx) # Save the 1-based alignment index
                    else:
                        # This residue was gapped in the alignment
                        for atom in residue:
                            atom.set_bfactor(0.0) # Gapped columns get 0
        
        # Save the new PDB file
        output_file = os.path.join(output_dir, f"{base_filename}_aln_idx.pdb")
        io.set_structure(structure)
        io.save(output_file)

    print(f"  Successfully saved new PDB files.")
    print("---------------------------------------------------\n")
    return output_dir
#########################################################################################################
def main(output_base_dir=""):
    """
    Main pipeline execution function.
    """
    parser = argparse.ArgumentParser(
        description="Standalone pipeline for structure-aware alignment and confidence scoring."
    )
    parser.add_argument(
        "pdb_directory", 
        type=str,
        help="The directory containing your .pdb or .cif structure files."
    )
    parser.add_argument(
        "--method",
        type=str,
        choices=['mafft', 'clustalw2'],
        default='clustalw2',
        help="The alignment method to use. (default: clustalw2)"
    )
    args = parser.parse_args()

    # --- Check for essential matrix file ---
    if not os.path.exists(fs_mat_file):
        print(f"Error: Required matrix file not found: {fs_mat_file}")
        print(f"Please place 'fs_mat.dat' in the same directory as this script.")
        sys.exit(1)
    
    # --- Step 1: Validate Input ---
    structure_file_list = validate_input_directory(args.pdb_directory)
    
    # --- Step 2: Check Chains ---
    check_chain_counts(structure_file_list)

    # --- Step 3: Run Foldseek ---
    fasta_3di_file = get_3di_fasta(args.pdb_directory)

    # --- Step 4: Run Aligner ---
    raw_alignment_file = run_alignment(fasta_3di_file, args.method)

    # --- Step 5: Clean Alignment File ---
    print("--- 5. Cleaning Alignment File ---")
    cleaned_alignment_file = "pdb_to_3di_3Di_clean.aln"
    clean_alignment_file(raw_alignment_file, cleaned_alignment_file)
    print("------------------------------------\n")

    # --- Step 6: Calculate Confidence Scores ---
    print("--- 6. Calculating Confidence Scores ---")
    
    (
        global_score, total_cols, valid_cols_list, dropped_cols_list, 
        conf_file, lookup_dict, norm_scores
    ) = calculate_confidence_scores(
        cleaned_alignment_file, 
        fs_mat_file
    )
    print("-------------------------------------------\n")

    # --- Step 7: Get Z-Score Assessment ---
    print("--- 7. Calculating Z-Score Assessment ---")
    
    zs = z_score_scop(global_score, args.method)
    zc = z_score_cath(global_score, args.method)
    avg_z = (zs + zc) / 2
    
    if avg_z >= 2.0: assessment = "Statistically Rare Alignment (high)"
    elif avg_z >= 1.0: assessment = "Above Average, Less Common"
    elif avg_z > -1.0: assessment = "Common, Expected Range"
    elif avg_z > -2.0: assessment = "Below Average, Less Common"
    else: assessment = "Statistically Rare Alignment (low)"
    
    print("---------------------------------------\n")
    
    # --- Step 8: Save Confidence B-Factor PDBs ---
    bfactor_dir = save_confidence_to_bfactor(
        args.pdb_directory, lookup_dict, norm_scores
    )

    # --- Step 9: Save Alignment-Index B-Factor PDBs ---
    index_bfactor_dir = save_alignment_index_to_bfactor(
        args.pdb_directory, lookup_dict
    )

    # --- Final Summary Report (MODIFIED) ---
    valid_cols_count = len(valid_cols_list)
    dropped_cols_count = len(dropped_cols_list)

    print("========= Results Summary =========")
    print(f"  Alignment Method:     {args.method}")
    print(f"  Global Confidence Score: {global_score:.4f} (Avg. of valid cols)")
    print(f"  Average Z-Score (v-v): {avg_z:.4f}")
    print(f"  Assessment:             {assessment}")
    print("---")
    print(f"  Alignment Columns:    {total_cols} (Total)")
    print(f"  Valid Columns (>=50%): {valid_cols_count} ({valid_cols_count/total_cols*100:.1f}%)")
    print(f"  Gappy Columns (<50%): {dropped_cols_count} ({dropped_cols_count/total_cols*100:.1f}%)")
    print("---")
    print(f"  Output Files:")
    print(f"  - 3Di Sequences:        {fasta_3di_file}")
    print(f"  - 3Di Alignment (clean): {cleaned_alignment_file}")
    print(f"  - Per-Column Data:    {conf_file}")
    print(f"  - PDBs w/ Confidence:   {bfactor_dir}/")
    print(f"  - PDBs w/ Aln. Index:   {index_bfactor_dir}/")
    print(f"  - (Raw Alignment:      {raw_alignment_file})")
    print("=================================")

# if __name__ == "__main__":
#     main()