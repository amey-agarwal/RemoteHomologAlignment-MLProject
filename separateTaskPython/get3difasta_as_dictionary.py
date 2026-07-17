def readFasta(path):
  name_seqs={}
  contig,name=[],""
  for line in open(path):
    line = line.strip()
    if not line:
      continue
    if line[0] == '>':
      if contig:
        name_seqs[name] = ''.join(contig)
        contig=[]
      name = line[1:]
      name_seqs[name] = ""
    else:
      contig.append(line)
  name_seqs[name] = ''.join(contig)
  return name_seqs

#print(readFasta('reserveFastaFiles/seqTrimAfterFoldSeek/3Di.fa'))