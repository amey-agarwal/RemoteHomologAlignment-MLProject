# Remote Homolog Alignment - ML Project
~ Amey Agarwal, Biswajit Banerjee, Anton Petrov
> Licensed under Williams Research Group, School of Biological Sciences, Georgia Tech, Atlanta GA [[iCOOL Lab](https://sites.gatech.edu/icool/)]

<p>This project is aimed at understanding the origins of the genetic code of the cell. It includes a pipeline with alphafold structure prediction, 3Di sequence encodings and alignments giving an output of aligned remote homologous protein.</p>

## Directory structure

```
├── README.MD
├── projectCmds.sh              : project start-up codes

# Alphafold 2 protein structure prediction
├── RunAlphaFold2.ipynb         : Alphafold 2 structure prediction Jupyter Notebook

# PDB files - experimental and predicted
├── reservePDBfilegroups        : reference PDB files
│   ├── folder/

# Fasta File Outputs
├── reserveFastaFiles           
│   ├── folder/
│   │   ├── 3di_aligned_*.fa    : 3di sequences aligned
│   │   ├── 3di_*.fa
│   │   ├── orig_aligned_*.fa   : original sequences aligned wrt 3di sequence alignment
│   │   └── *_trimmed.fa        : original sequences

# Structome-AlignViewer
├── workflow_v0.2.py            : structome module
├── fs_mat.dat                  : structome alignment matrix
├── readme.md.txt               : structome README

# Python code for pipelines
├── separateTaskPython          : helper functions
│   ├── runPipelines.py         : primary code run

# Environments
├── requirements.txt            : Python dependencies to run project
└── environment.yml             : Conda environment dependencies to run project

├── notes.txt                   : experimentation notes
```

## Getting Started 

<p>Project is run using python libraries and conda environment for Biopython and Bioinformatics tools</p>

- projectCmds.sh demonstrates the commands used for the run
- take a protein sequence and predict its structure using RunAlphaFold2.ipynb Jupyter Notebook
- Structome-AlignViewer code can be used to encode protein structures as 1D seqeunces - 3Di sequences
- Structome-AlignViewer code can be used to align 3Di sequences

<p>Demo Structome-AlignViewer run</p>

- PDB files => 3Di sequence => aligned 3Di sequence

```
python workflow_v0.2.py reservePDBfilegroups/<pick-folder-with-PDB-files> --method mafft
```

<p>Demo Pipeline run</p>

- PDB files => 3Di sequence => aligned 3Di sequence => Original protein sequences aligned

```
python separateTaskPython/runPipelines.py
```

## References
- Alphafold Prediction on Colab [Notebook Link](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) ~ Martin Steinegger and Sergey Ovchinnikov
- Structome-AlignViewer [Website Link](https://biosig.lab.uq.edu.au/structome_alignviewer/)
