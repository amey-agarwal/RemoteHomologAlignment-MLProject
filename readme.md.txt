# Structome‑AlignViewer — Local Walkthrough (for private, scalable runs)

This page is a step‑by‑step guide for running the standalone CLI on your own machine or cluster. It complements the online app by removing size, queueing, and privacy constraints.

Why run locally instead of the online app?

* Limit per job: The hosted app is capped at ~50 structures per submission.
* Shared server load: At busy times, runs can slow down.
* Input idiosyncrasies: There may be data edge‑cases we haven’t trapped yet.
* Confidentiality: Some users prefer not to upload proprietary structures.

The Bitbucket release + this walkthrough let you scale and validate the approach on your own infrastructure and datasets.

What you’ll do in this walkthrough

* Set up a small, reproducible environment.
* Prepare a directory of PDB/mmCIF structures (ideally single‑chain per file).
* Run the CLI to:
    - convert structures → Foldseek 3Di FASTA,
    - align 3Di sequences (MAFFT or ClustalW2) using fs_mat.dat,
    - compute per‑column confidence + global score + Z‑scores, and
    - write new PDBs with B‑factors set to confidence or alignment index.
* Inspect the results and validate quality.

Prerequisites

* Python 3.9+
* Biopython (installed via requirements.txt)
* External tools on your PATH:
    - foldseek + 3Di substitution matrix file "fs_mat.dat"
    - mafft and clustalw2
    - The programs should be available on command-line and the matrix file should be in the same location as the workflow.py script.

Tip: You can override tool paths by making the change at the top of the script.

Example dataset layout

```
my_structs/
├─ 1abc_A.pdb
├─ 2def_A.pdb
├─ 3ghi_A.cif
└─ ... (≥2 files total)
```
Recommendation: Ensure each file contains one polymer chain. If not, results will not be reliable.

# Run the CLI (ClustalW2 or MAFFT)

* python workflow_v0.2.py my_structs --method clustalw2
* python workflow_v0.2.py my_structs --method mafft

Expected outputs in the working directory:

  - 3Di.fa – 3Di sequences from Foldseek
  - 3Di.aln – raw multiple alignment
  - 3Di_clean.aln – consensus lines cleaned (*:. → spaces)
  - confidence_per_column.txt – per‑column confidence report
  - PDBs_with_confidence/ – PDBs with B‑factor = normalized confidence × 100
  - PDBs_with_aln_index/ – PDBs with B‑factor = alignment column index

Interpreting results

1) Console summary

You’ll see the global confidence, Z‑scores (SCOP/CATH baselines), the fraction of valid (non‑gappy) columns, and output paths.

  - Global Confidence: Average of normalized per‑column scores across columns with ≥50% non‑gaps.
  - Z‑scores: Relative to SCOP/CATH baselines for the chosen aligner.

2) Visual inspection

    Open any *_conf.pdb in your viewer of choice and color by B‑factor to reveal high/low‑confidence regions.

3) Column‑level table

    Inspect confidence_per_column.txt to spot low‑confidence columns, unusual residue mixes, or alignment artifacts.

Troubleshooting (common)

  - foldseek/clustalw2/mafft: command not found → Install/add to PATH.
  - **fs_mat.dat missing
  - output files not produced because input structures did not have single chain per structure.
