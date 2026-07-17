import os
# import re
# input_directory="reservePDBfilegroups/aL14_new"

# # SH3_bL19_1148_SYNYG_243346556180_d9ca6.result_SH3_bL19_1148_SYNYG_243346556180_d9ca6_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb
# # SH3_bL19_1148_SYNYG_243346556180
# # SH3_bL19_1148_SYNYG_24-33,46-55,61-80

# # separate code to fix bL19
# def readFastatoHashmap(path):
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
#       name = line[1:]
#       name_seqs[name] = ""
#     else:
#       contig.append(line)
#   name_seqs[name] = ''.join(contig).replace('-','')
#   return name_seqs

# fasta_dict = readFastatoHashmap("reserveFastaFiles/uL02_allSeqs_AF2Pred/SH3_uL02_trimmed.fa")

# # def fixbL19(fasta_file="SH3_bL19_trimmed.fa"):
# #     fasta_dict = readFastatoHashmap(fasta_file="SH3_bL19_trimmed.fa")

# def fixPdbNames(input_dir):
#     allpdbnames = []
#     allfilenames = [x for x in os.listdir(input_dir) if x.endswith('.pdb')]
#     for filename in allfilenames:
#         underScoreSplitArray = filename.split('_') # assuming SH3_[name] or OB_[name]
#         subName = underScoreSplitArray[0]+"_"+underScoreSplitArray[1]
#         findName = ".result_"+subName

#         trimName = filename[:filename.index(findName)] # trim filename to just before .result in the name
#         nameOfPdb = trimName[:trimName.rindex('_')]

#         subNameofPDB = nameOfPdb[:nameOfPdb.rindex('_')]
#         # if 'SH3_uL02b_62928' in subNameofPDB:
#         #    print("=======================",subNameofPDB)
#         if bool(re.search(r'\d', subNameofPDB.split('_')[-1])):
#            #print(subNameofPDB)
#            first_index = next((idx for idx, char in enumerate(subNameofPDB.split('_')[-1]) if char.isdigit()), None)
#            subNameofPDB = '_'.join(subNameofPDB.split('_')[:-1]) + "_" + subNameofPDB.split('_')[-1][:first_index]
#         for properFileName in fasta_dict.keys():
#            if '/' in properFileName:
#               if subNameofPDB in properFileName[:properFileName.index('/')]:
#                 nameOfPdb = properFileName
#                 break
#            else:
#               if subNameofPDB in properFileName[:properFileName.rindex('_')]:
#                 nameOfPdb = properFileName
#                 break
#         allpdbnames.append(nameOfPdb)
#     for i in range(len(allfilenames)):
#         filepath = f"{input_directory}/{allfilenames[i]}"
#         if not os.path.exists(filepath):
#             print(f"{filepath} file not found")
#             exit(0)
#     for pdbname in allpdbnames:
#        if pdbname not in fasta_dict.keys():
#           print(pdbname)
#     for i in range(len(allpdbnames)):
#        if '/' in allpdbnames[i]:
#           allpdbnames[i] = allpdbnames[i][:allpdbnames[i].index('/')] + '[forward_slash_character]' + allpdbnames[i][allpdbnames[i].index('/')+1:]
#           #print(allpdbnames[i])
#     for i in range(len(allpdbnames)):
#        if '/' in allpdbnames[i] or '\\' in allpdbnames[i]:
#           print(allpdbnames[i])
#     for i in range(len(allfilenames)):
#         src_filepath = f"{input_directory}/{allfilenames[i]}"
#         dst_filepath = f"{input_directory}/{allpdbnames[i]}.pdb"
#         os.rename(src_filepath,dst_filepath)

# fixPdbNames(input_directory)

def newfixPdb(input_dir,protein_name):
    allpdbnames = []
    allfilenames = [x for x in os.listdir(input_dir) if x.endswith('.pdb')]

    for filename in allfilenames:
        fixed_filename = filename[:filename.index('_unrelaxed')]
        fixed_filename = filename[:filename.rindex(f'_{protein_name}')]

        allpdbnames.append(fixed_filename)

    for i in range(len(allpdbnames)):
       if '/' in allpdbnames[i]:
          allpdbnames[i] = allpdbnames[i][:allpdbnames[i].index('/')] + '[forward_slash_character]' + allpdbnames[i][allpdbnames[i].index('/')+1:]
          #print(allpdbnames[i])
    for i in range(len(allpdbnames)):
       if '/' in allpdbnames[i] or '\\' in allpdbnames[i]:
          print(allpdbnames[i])
          print("pdb files not renamed")
          exit(0)

    for i in range(len(allfilenames)):
        src_filepath = f"{input_dir}/{allfilenames[i]}"
        dst_filepath = f"{input_dir}/{allpdbnames[i]}.pdb"
        os.rename(src_filepath,dst_filepath)
    print("pdb files renamed")


# newfixPdb(input_directory, "aL14")