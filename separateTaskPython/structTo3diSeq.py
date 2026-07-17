import renamedWorkflow as w

def readFasta(path,trimDashes):
  name_seqs={}
  contig,name=[],""
  for line in open(path):
    line = line.strip()
    if not line:
      continue
    if line[0] == '>':
      if contig:
        name_seqs[name] = ''.join(contig)
        if trimDashes:
          name_seqs[name] = name_seqs[name].replace('-','')
        contig=[]
      name = line[1:]
      name_seqs[name] = ""
    else:
      contig.append(line)
  name_seqs[name] = ''.join(contig)
  if trimDashes:
      name_seqs[name] = name_seqs[name].replace('-','')
  return name_seqs

def putAllArraySeqstoFile(names,aligned_family,filename="output.fa"):
  """
  Input : names : [seq1,seq2] , aligned_family : [seq1,seq2]
  """
  with open(f'{filename}', 'w') as f:
    for i in range(len(aligned_family)):
        f.write(f">{names[i]}\n")
        f.write(f"{aligned_family[i]}\n")
  print(f"File written {filename} fasta")
  return

def PDBto3DiSeqs(single_pdb_directory="reservePDBfilegroups/uL02_allSeqs_AF2Pred",
                 output_file="reserveFastaFiles/uL02_allSeqs_AF2Pred/3di_uL02.fa"):
    w.check_chain_counts(single_pdb_directory)
    
    fasta_3di_file = w.get_3di_fasta(single_pdb_directory,output_file)
    
    fasta_dict = readFasta(output_file,False)
    fasta_names = list(fasta_dict.keys())
    
    flag = 0
    for i in range(len(fasta_names)):
       if '[forward_slash_character]' in fasta_names[i]:
         flag = 1
         fasta_names[i] = fasta_names[i].replace('[forward_slash_character]','/')
    
    if flag == 1:
       print("\n\nIdentified sequence names with '/' character, renamed to [forward_slash_character]\n\n")

    putAllArraySeqstoFile(fasta_names, list(fasta_dict.values()),output_file)

    print("fasta 3di file created")