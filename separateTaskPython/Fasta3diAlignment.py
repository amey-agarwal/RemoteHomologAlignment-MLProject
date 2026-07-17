import renamedWorkflow as w
import os

def threediFastaAlign(fasta_3di_file = "reserveFastaFiles/uL24+bL19/3Di.fa",
                       alignment_method = "mafft",
                       aligned_filename=""):
        SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
        
        fs_mat_file = os.path.join(SCRIPT_DIR, "fs_mat.dat")
        
        # --- Step 4: Run Aligner ---
        raw_alignment_file = w.run_alignment(fasta_3di_file, alignment_method)

        # --- Step 5: Clean Alignment File ---
        print("--- 5. Cleaning Alignment File ---")
        cleaned_alignment_file = ""
        if not aligned_filename:
                cleaned_alignment_file = "reserveFastaFiles/uL24+bL19/3Di_aligned.fa"
        else:
                cleaned_alignment_file = aligned_filename
        w.clean_alignment_file(raw_alignment_file, cleaned_alignment_file)
        print("------------------------------------\n")

# threediFastaAlign()
# --- Step 6: Calculate Confidence Scores ---
# print("--- 6. Calculating Confidence Scores ---")

# (
#     global_score, total_cols, valid_cols_list, dropped_cols_list, 
#     conf_file, lookup_dict, norm_scores
# ) = w.calculate_confidence_scores(
#     cleaned_alignment_file, 
#     fs_mat_file
# )
# print("-------------------------------------------\n")

# # --- Step 7: Get Z-Score Assessment ---
# print("--- 7. Calculating Z-Score Assessment ---")

# zs = w.z_score_scop(global_score, alignment_method)
# zc = w.z_score_cath(global_score, alignment_method)
# avg_z = (zs + zc) / 2

# if avg_z >= 2.0: assessment = "Statistically Rare Alignment (high)"
# elif avg_z >= 1.0: assessment = "Above Average, Less Common"
# elif avg_z > -1.0: assessment = "Common, Expected Range"
# elif avg_z > -2.0: assessment = "Below Average, Less Common"
# else: assessment = "Statistically Rare Alignment (low)"

# print("---------------------------------------\n")

# # --- Step 8: Save Confidence B-Factor PDBs ---
# bfactor_dir = w.save_confidence_to_bfactor(
#     alignment_method, lookup_dict, norm_scores
# )

# # --- Step 9: Save Alignment-Index B-Factor PDBs ---
# index_bfactor_dir = w.save_alignment_index_to_bfactor(
#     alignment_method, lookup_dict
# )

# # --- Final Summary Report (MODIFIED) ---
# valid_cols_count = len(valid_cols_list)
# dropped_cols_count = len(dropped_cols_list)

# print("========= Results Summary =========")
# print(f"  Alignment Method:     {alignment_method}")
# print(f"  Global Confidence Score: {global_score:.4f} (Avg. of valid cols)")
# print(f"  Average Z-Score (v-v): {avg_z:.4f}")
# print(f"  Assessment:             {assessment}")
# print("---")
# print(f"  Alignment Columns:    {total_cols} (Total)")
# print(f"  Valid Columns (>=50%): {valid_cols_count} ({valid_cols_count/total_cols*100:.1f}%)")
# print(f"  Gappy Columns (<50%): {dropped_cols_count} ({dropped_cols_count/total_cols*100:.1f}%)")
# print("---")
# print(f"  Output Files:")
# print(f"  - 3Di Sequences:        {fasta_3di_file}")
# print(f"  - 3Di Alignment (clean): {cleaned_alignment_file}")
# print(f"  - Per-Column Data:    {conf_file}")
# print(f"  - PDBs w/ Confidence:   {bfactor_dir}/")
# print(f"  - PDBs w/ Aln. Index:   {index_bfactor_dir}/")
# print(f"  - (Raw Alignment:      {raw_alignment_file})")
# print("=================================")