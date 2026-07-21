import re
# orig_seqs_file="reserveFastaFiles/bL19+uL02/SH3_bL19+uL02.fa"
# threedi_seqs_file="reserveFastaFiles/bL19+uL02/3di_aligned_aligned_bL19+uL02.fa"

# assuming the names of the sequences are exactly the same
# assuming all the names are short

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

def addDashes(seq_name,protein_seq,threedi_seq):
   protein_seq = list(protein_seq)
   threedi_seq = list(threedi_seq)
   
   dashPointer, protein_pointer = len(threedi_seq), len(protein_seq)

   if len(''.join(threedi_seq).replace('-',''))!=len(''.join(protein_seq)):
     print("===============================================================")
     print(f"{seq_name} found Protein sequence and 3Di not same length")
     print(f"{''.join(threedi_seq).replace('-','')}\n{''.join(protein_seq)}")
     print("Skipping this protein")
     print("===============================================================")
     return ''
   
   while dashPointer>0:
    dashPointer-=1
    if threedi_seq[dashPointer] == '-':
      protein_seq = protein_seq[:protein_pointer] + ['-'] + protein_seq[protein_pointer:]
    else:
      protein_pointer-=1
   
   return ''.join(protein_seq) 

pattern = r'^(\S+)\s+(\S+)'

def read3DialnFile(threedi_file):
  name_and_seq = {}
  firstLine = False
  with open(threedi_file,'r') as f:
    for line in f:
      line = line.strip()
      
      if not firstLine:
        firstLine = True
        continue
      
      match = re.search(pattern, line)

      if match:
        name=match.group(1)
        sequence=match.group(2)
        
        if name in name_and_seq:
          name_and_seq[name] = ''.join(list(name_and_seq[name]) + list(sequence))
        else:
          name_and_seq[name] = sequence
  return name_and_seq

def putAllArraySeqstoFile(names,aligned_family,filename="output"):
  """
  Input : names : [seq1,seq2] , aligned_family : [seq1,seq2]
  """
  with open(f'{filename}.fa', 'w') as f:
    for i in range(len(aligned_family)):
        f.write(f">{names[i]}\n")
        f.write(f"{aligned_family[i]}\n")
  print(f"File written {filename} fasta")
  return

def threedi_to_protein(orig_seqs_file,threedi_seqs_file,orig_aligned_filename):
  orig_seqs_hashmap=readFasta(orig_seqs_file,True)
  threedi_seqs_hashmap=readFasta(threedi_seqs_file,False)
  new_orig_seqs_hashmap={}
  new_threedi_seqs_hashmap={}
  #threedi_seqs_hashmap=read3DialnFile(threedi_seqs_file)

  for seq_name in threedi_seqs_hashmap:
    if "ALPHAFOLD MONOMER" in seq_name:
      new_seq_name = seq_name[:seq_name.index("ALPHAFOLD MONOMER")-1]
      new_threedi_seqs_hashmap[new_seq_name] = threedi_seqs_hashmap[seq_name]
    else:
      new_threedi_seqs_hashmap[seq_name] = threedi_seqs_hashmap[seq_name]

  for seq_name in orig_seqs_hashmap:
    # threedi_seqs_hashmap_name = threedi_seqs_hashmap[k]
    if seq_name in new_threedi_seqs_hashmap:
      print(f"found {seq_name} common in orig_seqs and threedi_seqs")
      if '[forward_slash_character]' in new_threedi_seqs_hashmap:
        new_threedi_seqs_hashmap[seq_name] = new_threedi_seqs_hashmap[seq_name].replace('[forward_slash_character]','/')
      if '[forward_slash_character]' in orig_seqs_hashmap[seq_name]:
        orig_seqs_hashmap[seq_name] = orig_seqs_hashmap[seq_name].replace('[forward_slash_character]','/')


      aligned_protein_seq = addDashes(seq_name,orig_seqs_hashmap[seq_name],new_threedi_seqs_hashmap[seq_name])
      if len(aligned_protein_seq):
        new_orig_seqs_hashmap[seq_name] = aligned_protein_seq

  seqs_names=list(new_orig_seqs_hashmap.keys())
  orig_seqs_aligned=list(new_orig_seqs_hashmap.values())
  threedi_seqs_aligned=list(new_threedi_seqs_hashmap.values())

  new_orig_seqs_aligned_name = ""
  if not orig_aligned_filename:
    new_orig_seqs_aligned_name = 'reserveFastaFiles/orig_aligned_custom'
  else:
    new_orig_seqs_aligned_name = orig_aligned_filename
  putAllArraySeqstoFile(seqs_names, orig_seqs_aligned, new_orig_seqs_aligned_name)
#putAllArraySeqstoFile(seqs_names,threedi_seqs_aligned,'reserveFastaFiles/bL19_allSeqs_AF2Pred/3di_aligned_bL19.fa')

# threedi_to_protein(orig_seqs_file="reserveFastaFiles/uL24+bL19/orig.fa",
#                    threedi_seqs_file="reserveFastaFiles/uL24+bL19/3Di_aligned.fa",
#                    orig_aligned_filename="reserveFastaFiles/uL24+bL19/orig_aligned_uL24+bL19.fa")

# threedi_to_protein(orig_seqs_file="")