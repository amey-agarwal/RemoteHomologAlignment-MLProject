import os
import shutil
from transferPdbFromzip import extract_pdbs_from_nested_zip, extract_pdb_files
from structTo3diSeq import PDBto3DiSeqs
from Fasta3diAlignment import threediFastaAlign
from threedi_to_protein import threedi_to_protein
from trimPDBnames import newfixPdb
import renamedWorkflow as workflowCode
from trimFastaNamesandSeqs import combine_fasta_files

def checkFileExists(file_path):
    if not os.path.exists(file_path):
        return False
    return True

class HelperFunctions:
    def __init__(self):
        pass

    def PDBzip2PDBfolder(self,src_zip,dst_dir,PDBzipType,protein_name):
        if PDBzipType == "nested":
            extract_pdbs_from_nested_zip(src_zip,dst_dir)
        elif PDBzipType == "root":
            extract_pdb_files(src_zip,dst_dir,protein_name)
        else:
            print("PDBzipType not given")
            exit(0)
    
    def PDBfolderNamefix(self,pdb_folder,protein_name,fixType):
        if fixType=='fix AF2 names':
            newfixPdb(pdb_folder,protein_name)
        else:
            print("No fixes for PDB files")

    def PDBto3di(self,pdb_folder,output_file):
        PDBto3DiSeqs(pdb_folder,output_file)
    
    def align3diFasta(self,fasta_file, alignment_method, output_file):
        threediFastaAlign(fasta_file, alignment_method, output_file)

    def threedi2protein(self, trimmed_fasta, aligned_3di_file, output_file):
        threedi_to_protein(trimmed_fasta, aligned_3di_file, output_file)

def PDBstoOrigSeqAlignments(zip_files, seq_folder_name, protein_domain, protein, orig_trimmed_seq_files=[], orig_trimmed_seqs_combined_filename="", PDBzipType=None, PDBfixType=None, eval=False):
    """
    zip_files: [], all zip files to be placed in PDB folder
    protein: [], all proteins aligned
    orig_trimmed_seq_files: [], all orig trimmed seq files - SH3_uL24_trimmed.fa
    PDBfixType: 'fix AF2 names'
    """
    fasta_folder = f"reserveFastafiles/{seq_folder_name}"
    pdb_folder = f"reservePDBfilegroups/{seq_folder_name}"

    # cp ~/Desktop/Assignments/winter/sharedAlignments/ALL-CODED-ALIGNMENTS/SH3_uL02_trimmed.fa reservePDBFiles/{seq_folder_name}
    dst = fasta_folder
    if not orig_trimmed_seq_files:
        # src = f"~/Desktop/Assignments/winter/sharedFiles/ALL-SH3-OB-PERMUTED/{protein_domain}_{protein}_trimmed.fa"
        print("trimmed seq file not provided")
        return
    else:
        for src in orig_trimmed_seq_files:
            if not os.path.exists(src):
                print(f"Original trimmed seqs alignment file not found at {src}\nNo trimmed seq file given")
                return
    os.makedirs(fasta_folder,exist_ok=True)
    os.makedirs(pdb_folder,exist_ok=True)

    for src in orig_trimmed_seq_files:
        shutil.copy(src,dst)

    # PDB zip -> PDB folder
    if zip_files:
        for pdb_files_zip in zip_files:
            if not os.path.exists(pdb_files_zip):
                print(f"{pdb_files_zip} not found")
                exit(0)
        for i in range(len(zip_files)):
                print(zip_files[i])
                # protein name included
                HelperFunctions().PDBzip2PDBfolder(zip_files[i],pdb_folder,PDBzipType,protein[i]); print("PDB zip -> PDB folder done")
    else:
        print("No Zip files given.Proceeding without zip files")
            
    # folder name fix
    # if PDBfixType:
    #     HelperFunctions().PDBfolderNamefix(pdb_folder,protein,PDBfixType); print("PDB folder name fix done")

    # python workflow.py reservePDBfiles/{seq_folder_name} --method mafft
    # mv 3Di.fa 3Di_clean.aln reserveFastafiles/{seq_folder_name}
    # workflowCode.main(f"reserveFastafiles/{seq_folder_name}")
    protein_combined_name = '+'.join(protein)

    # PDB folder -> 3di seqs
    threedi_seqs_unaligned_file = fasta_folder + f"/3di_{protein_domain}_{protein_combined_name}.fa"
    HelperFunctions().PDBto3di(pdb_folder, threedi_seqs_unaligned_file);  print("PDB folder -> 3di seqs done")

    # 3di seqs -> 3di aligned seqs
    threedi_seqs_aligned_file = fasta_folder + f"/3di_aligned_{protein_domain}_{protein_combined_name}.fa"
    HelperFunctions().align3diFasta(threedi_seqs_unaligned_file, "mafft", threedi_seqs_aligned_file);  print("3di seqs -> 3di aligned seqs done")
    
    # 3di aligned seqs -> protein residue
    orig_aligned_output_file = fasta_folder + f"/orig_aligned_{protein_domain}_{protein_combined_name}"
    #orig_trimmed_seq_file = orig_trimmed_seq_files[0]
    #if len(orig_trimmed_seq_files)>1:
    # to rename the files as <protein>_<filename>
    combine_fasta_files(orig_trimmed_seq_files,fasta_folder+f"/{orig_trimmed_seqs_combined_filename}.fa",protein)
    orig_trimmed_seq_file = fasta_folder+f"/{orig_trimmed_seqs_combined_filename}.fa"

    print(orig_trimmed_seq_file,threedi_seqs_aligned_file,orig_aligned_output_file)

    HelperFunctions().threedi2protein(orig_trimmed_seq_file,
                                      threedi_seqs_aligned_file,
                                      orig_aligned_output_file);  print("3di aligned seqs -> protein residue done")

    if not eval:
        print("Alignment not eval")
        return
    # FIX ALL NAMES -> CREATE MAPPING OF OLD NAMES TO NEW NAMES AND STORE IT
    # NAME FIXING | python separateTaskPython/trimFastaNamesandSeqs.py
    # NAME FIXING | python separateTaskPython/3di_to_protein.py

#PDBstoOrigSeqAlignments('SH3_uL02_trimmed.fa','uL02_AF2_files_bundle.zip','uL02_allSeqs_AF2Pred')

def different3DiFastaAlignments(files_orig:list,aligned_files_3di:list,final_folder_name:str,final_aligned_filename:str) -> None:
    """Used when there are 2 or more folders with their own alignment files, and these folders need to be aligned into single alignment file

    files_orig:['SH3_uL24.fa','SH3_bL19.fa']
    aligned_files_3di:['3di_aligned_SH3_uL24.fa','3di_aligned_SH3_bL19.fa']
    final_folder_name:"CompositeEverySeqPred_uL24+bL19_aligned.fa"
    final_aligned_filename:"aligned_SH3_uL24+bL19.fa"
    """
    if len(files_orig) < 2 or len(aligned_files_3di) < 2:
        print("need more than 1 file input")
        return
    if len(files_orig)!=len(aligned_files_3di):
        print("file length not same")
        return

    fasta_folder = f'reserveFastaFiles/{final_folder_name}'
    os.makedirs(fasta_folder,exist_ok=True)

    # copy a base file into the destniation folder
    src = files_orig[0]
    if not checkFileExists(src): print(f"{src} not found"); return

    protein_orig_dst = fasta_folder + "/combined_orig_seqs.fa" # base file
    shutil.copy(src,protein_orig_dst)

    # copy a base file into the destniation folder
    src = aligned_files_3di[0]
    if not checkFileExists(src): print(f"{src} not found"); return

    protein_3di_dst = fasta_folder + "/combined_3di_seqs.fa"
    shutil.copy(src,protein_3di_dst)

    # copy and append all the files into single file
    for i in range(1,len(files_orig)):
        protein_file_orig = files_orig[i]
        protein_file_3di = aligned_files_3di[i]

        if not checkFileExists(protein_file_orig) and checkFileExists(protein_file_3di):
            print(f"{protein_file_orig} or {protein_file_3di} not found")
            continue

        with open(protein_file_orig, "r") as source, open(protein_orig_dst, "a") as dest:
            dest.write(source.read())
        
        with open(protein_file_3di, "r") as source, open(protein_3di_dst, "a") as dest:
            dest.write(source.read())

    # run code for alignment
    aligned_3di = fasta_folder + f"/3di_{final_aligned_filename}.fa"
    HelperFunctions().align3diFasta(protein_3di_dst, "mafft", aligned_3di)

    # map back to protein residue from 3di
    output_file = fasta_folder + f"/orig_{final_aligned_filename}" # don't write fa
    HelperFunctions().threedi2protein(protein_orig_dst, aligned_3di, output_file)

    return

def alignDifferent3DiFastaAlignments(base_dir, file_list_alignment, aligned_folder, aligned_filename):
    """
    file_list_alignment:[], list of folders that contain the files for alignment
    aligned_folder:"", final aligned folder
    aligned_filename:"", filename of aligned file iin aligned folder
    """
    if len(file_list_alignment) < 2:
        print("need more than 1 file input")
        return
        
    fasta_folder = base_dir + "/" + aligned_folder
    protein = file_list_alignment[0].split('_')[0]
    
    os.makedirs(fasta_folder,exist_ok=True)

    # copy the first folder file to be aligned into the destination folder as a base file
    src = base_dir + "/" + file_list_alignment[0] + "/" + f"SH3_{protein}_trimmed.fa"
    if not checkFileExists(src): print(f"{src} not found"); return
    protein_orig_dst = fasta_folder + "/" + "orig_" + aligned_filename + ".fa"
    shutil.copy(src,protein_orig_dst)

    src = base_dir + "/" + file_list_alignment[0] + "/" + f"3di_aligned_SH3_{protein}.fa"
    if not checkFileExists(src): print(f"{src} not found"); return
    protein_3di_dst = fasta_folder + "/" + "3di_" + aligned_filename + ".fa"
    shutil.copy(src,protein_3di_dst)

    for i in range(1,len(file_list_alignment)):
        protein = file_list_alignment[i].split('_')[0]
        protein_file_orig = base_dir + "/" + file_list_alignment[i] + "/" + f"SH3_{protein}_trimmed.fa"
        protein_file_3di = base_dir + "/" + file_list_alignment[i] + "/" + f"3di_aligned_SH3_{protein}.fa"

        if not checkFileExists(protein_file_orig) and checkFileExists(protein_file_3di):
            print(f"{protein_file_orig} or {protein_file_3di} not found")
            continue

        with open(protein_file_orig, "r") as source, open(protein_orig_dst, "a") as dest:
            dest.write(source.read())
        
        with open(protein_file_3di, "r") as source, open(protein_3di_dst, "a") as dest:
            dest.write(source.read())

    # run code for alignment
    aligned_3di = fasta_folder + f"/3di_aligned_{aligned_filename}.fa"
    HelperFunctions().align3diFasta(protein_3di_dst, "mafft", aligned_3di)

    # map back to protein residue from 3di
    output_file = fasta_folder + f"/orig_aligned_{aligned_filename}" # don't write fa
    HelperFunctions().threedi2protein(protein_orig_dst, aligned_3di, output_file)

def Progressive_alignDifferent3DiFastaAlignments(base_dir, file_list_alignment, aligned_folder, aligned_filename):

    if len(file_list_alignment) < 3:
        print("need more than 2 file input")
        return
        
    fasta_folder = base_dir + "/" + aligned_folder
    protein = file_list_alignment[0].split('_')[0]
    
    os.makedirs(fasta_folder,exist_ok=True)

    src = base_dir + "/" + file_list_alignment[0] + "/" + f"SH3_{protein}_trimmed.fa"
    if not checkFileExists(src): print(f"{src} not found"); return
    protein_orig_dst = fasta_folder + "/" + "orig_" + aligned_filename + ".fa"
    shutil.copy(src,protein_orig_dst)

    src = base_dir + "/" + file_list_alignment[0] + "/" + f"3di_aligned_SH3_{protein}.fa"
    if not checkFileExists(src): print(f"{src} not found"); return
    protein_3di_dst = fasta_folder + "/" + "3di_" + aligned_filename + ".fa"
    shutil.copy(src,protein_3di_dst)

    print(f"taken {protein}")

    for i in range(1,len(file_list_alignment)):
        protein = file_list_alignment[i].split('_')[0]
        protein_file_orig = base_dir + "/" + file_list_alignment[i] + "/" + f"SH3_{protein}_trimmed.fa"
        protein_file_3di = base_dir + "/" + file_list_alignment[i] + "/" + f"3di_aligned_SH3_{protein}.fa"

        if not checkFileExists(protein_file_orig) and checkFileExists(protein_file_3di):
            print(f"{protein_file_orig} or {protein_file_3di} not found")
            continue

        with open(protein_file_orig, "r") as source, open(protein_orig_dst, "a") as dest:
            dest.write(source.read())
        
        with open(protein_file_3di, "r") as source, open(protein_3di_dst, "a") as dest:
            dest.write(source.read())

        # run code for alignment
        #aligned_3di = fasta_folder + f"/3di_aligned_{aligned_filename}.fa"
        HelperFunctions().align3diFasta(protein_3di_dst, "mafft", protein_3di_dst)

        # map back to protein residue from 3di
        output_file = fasta_folder + f"/orig_aligned_{aligned_filename}" # don't write fa
        HelperFunctions().threedi2protein(protein_orig_dst, protein_3di_dst, output_file)

        print(f"aligned {protein}")
    # if len(file_list_alignment) < 3:
    #     print("need more than 2 files input")
    #     return
    
    # first_two_files = [file_list_alignment[0],file_list_alignment[1]]
    # alignDifferent3DiFastaAlignments(base_dir, first_two_files, aligned_folder, aligned_filename)
        
    # fasta_folder = base_dir + "/" + aligned_folder
    # protein = file_list_alignment[0].split('_')[0]

    # protein_orig_dst = fasta_folder + f"/orig_{aligned_filename}.fa"
    # protein_3di_dst = fasta_folder + f"/3di_aligned_{aligned_filename}.fa"

    # for i in range(2,len(file_list_alignment)):
    #     progressive_list = [fasta_folder, file_list_alignment[i]]
        
    #     protein = file_list_alignment[i].split('_')[0]
    #     protein_file_orig = base_dir + "/" + file_list_alignment[i] + "/" + f"SH3_{protein}_trimmed.fa"
    #     protein_file_3di = base_dir + "/" + file_list_alignment[i] + "/" + f"3di_aligned_SH3_{protein}.fa"

    #     if not checkFileExists(protein_file_orig) and checkFileExists(protein_file_3di):
    #         print(f"{protein_file_orig} or {protein_file_3di} not found")
    #         continue

    #     with open(protein_file_orig, "r") as source, open(protein_orig_dst, "a") as dest:
    #         dest.write(source.read())
        
    #     with open(protein_file_3di, "r") as source, open(protein_3di_dst, "a") as dest:
    #         dest.write(source.read())

    # # run code for alignment
    # aligned_3di = fasta_folder + f"/3di_aligned_{aligned_filename}.fa"
    # HelperFunctions().align3diFasta(protein_3di_dst, "mafft", aligned_3di)

    # # map back to protein residue from 3di
    # output_file = fasta_folder + f"/orig_aligned_{aligned_filename}" # don't write fa
    # HelperFunctions().threedi2protein(protein_orig_dst, aligned_3di, output_file)

def main():
    PDBstoOrigSeqAlignments(zip_file="reservePDBfilegroups/rawpredictionzipFiles-af2/localSeqAlignment/uL24_AF2_files_bundle.zip",
                            seq_folder_name="uL24_allSeqs_AF2Pred",
                            protein_domain="SH3",
                            protein="uL24",
                            orig_trimmed_seq_file="/Users/ameyagarwal/Desktop/Assignments/winter/sharedFiles/ALL-SH3-OB-PERMUTED/SH3_uL24_trimmed.fa")

    files_to_be_aligned = ['uL02_allSeqs_AF2Pred','uL24_allSeqs_AF2Pred', 'bL19_allSeqs_AF2Pred', 
                           'aL14_allSeqs_AF2Pred','aL21_allSeqs_AF2Pred','aS04_allSeqs_AF2Pred', 
                           'bEFP_allSeqs_AF2Pred','aIF5a_allSeqs_AF2Pred']
    base_dir="reserveFastaFiles"
    aligned_folder_name = "all_SH3"
    aligned_filename = "aligned_all_SH3"
    alignDifferent3DiFastaAlignments(base_dir, files_to_be_aligned, aligned_folder_name, 
                                     aligned_filename)
    
    # progrressively align approach
    # aligned_folder_name = "Progressive_allSH3exceptuL24"
    # aligned_filename = "Progressive_aligned_allSH3exceptuL24"
    # Progressive_alignDifferent3DiFastaAlignments(base_dir, files_to_be_aligned, aligned_folder_name, 
    #                                  aligned_filename)

# main()

# rm -rf reserveFastaFiles/new
# rm -rf reservePDBfilegroups/new
# PDBstoOrigSeqAlignments(['reservePDBfilegroups/rawpredictionzipFiles-af2/AFDatabase/SH3_uL24_PDBs.zip'],
#                          'CompositeEverySeqPred_uL24_alignment',
#                          'SH3',['uL24'],orig_trimmed_seq_files=['/Users/ameyagarwal/Downloads/SH3_uL24.fas'],
#                                                           orig_trimmed_seqs_combined_filename="uL24_orig_seqs_fullFasta",
#                                                           PDBzipType="root",
#                                                           PDBfixType=None)

# PDBstoOrigSeqAlignments(['reservePDBfilegroups/rawpredictionzipFiles-af2/AFDatabase/SH3_bL19_PDBs.zip'],
#                          'CompositeEverySeqPred_bL19_alignment',
#                          'SH3',['bL19'],orig_trimmed_seq_files=['/Users/ameyagarwal/Downloads/SH3_bL19.fas'],
#                                                           orig_trimmed_seqs_combined_filename="bL19_orig_seqs_fullFasta",
#                                                           PDBzipType="root",
#                                                           PDBfixType=None)

PDBstoOrigSeqAlignments(['reservePDBfilegroups/rawpredictionzipFiles-af2/AFDatabase/SH3_uL24_PDBs.zip',
                         'reservePDBfilegroups/rawpredictionzipFiles-af2/AFDatabase/SH3_bL19_PDBs.zip'],
                         'CompositeEverySeqPred_uL24+bL19_alignment',
                         'SH3',['uL24','bL19'],orig_trimmed_seq_files=['/Users/ameyagarwal/Downloads/SH3_uL24.fas',
                                                                       '/Users/ameyagarwal/Downloads/SH3_bL19.fas'],
                                                          orig_trimmed_seqs_combined_filename="orig_seqs_uL24+bL19_fullFasta",
                                                          PDBzipType="root",
                                                          PDBfixType=None)


# different3DiFastaAlignments(['reserveFastaFiles/CompositeEverySeqPred_uL24_alignment/orig_aligned_SH3_uL24.fa',
#                              'reserveFastaFiles/CompositeEverySeqPred_bL19_alignment/orig_aligned_SH3_bL19.fa'],
#                              ['reserveFastaFiles/CompositeEverySeqPred_uL24_alignment/3di_aligned_SH3_uL24.fa',
#                               'reserveFastaFiles/CompositeEverySeqPred_bL19_alignment/3di_aligned_SH3_bL19.fa'],
#                               "CompositeEverySeqPred_bL19_uL24_alignment",
#                               "CompositeEverySeqPred_bL19_uL24")