
# fixed_len_3Di_name="SH3_uL02b_10022"
# fasta_file = "reserveFastaFiles/uL02_allSeqs_AF2Pred/SH3_uL02_trimmed.fa"
# name_new_file="reserveFastaFiles/uL02_allSeqs_AF2Pred/namesTrimmed_gapsRemoved_SH3_uL02_trimmed" #don't add .fa

# def readFastatoHashmap(path,N):
#   name_seqs={}
#   contig,name=[],""
#   for line in open(path):
#     line = line.strip()
#     if not line:
#       continue
#     if line[0] == '>':
#       if contig:
#         name_seqs[name] = ''.join(contig).replace('-','')
#         contig=[]
#       name = line[1:][:N]
#       name_seqs[name] = ""
#     else:
#       contig.append(line)
#   name_seqs[name] = ''.join(contig).replace('-','')
#   return name_seqs

# def putAllArraySeqstoFile(names,aligned_family,filename="output"):
#   """
#   Input : names : [seq1,seq2] , aligned_family : [seq1,seq2]
#   """
#   with open(f'{filename}.fa', 'w') as f:
#     for i in range(len(aligned_family)):
#         f.write(f">{names[i]}\n")
#         f.write(f"{aligned_family[i]}\n")
#   print(f"File written {filename}.fa fasta")
#   return

def combine_fasta_files(fasta_files, output_file,protein_names):
    """
    Combine multiple FASTA files into a single FASTA file.

    Parameters
    ----------
    fasta_files : list[str]
        List of FASTA file paths.

    output_file : str
        Output FASTA filename.

    Example
    -------
    combine_fasta_files(
        [
            "SH3_bL19.fas",
            "SH3_uL24.fas",
            "SH3_aS04.fas"
        ],
        "SH3_all.fas"
    )
    """

    total_sequences = 0
    with open(output_file, "w") as outfile:
        for i in range(len(fasta_files)):
            fasta = fasta_files[i]
            protein = protein_names[i]
            print(f"Appending {fasta}")
            with open(fasta, "r") as infile:
                for line in infile:
                    if line.startswith(">"):
                        total_sequences += 1
                        line = f">{protein}_{line[line.index('>')+1:]}"
                        outfile.write(line)
                    else:
                        line = line.strip().replace('-','')
                        outfile.write(line)
                        outfile.write("\n")
                # Ensure exactly one newline between files
                if not line.endswith("\n"):
                    outfile.write("\n")

    print(f"\nCreated: {output_file}")
    print(f"Total sequences: {total_sequences}")

# n = len(fixed_len_3Di_name)
# fastaFile = readFastatoHashmap(fasta_file,n)
# putAllArraySeqstoFile(list(fastaFile.keys()),list(fastaFile.values()),name_new_file)