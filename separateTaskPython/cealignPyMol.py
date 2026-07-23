import os
import glob
from collections import OrderedDict

import pymol
from pymol import cmd


class CEAlignFastaBuilder:
    """
    Uses PyMOL's CEAlign implementation to create a
    structure-guided multiple sequence alignment.
    """

    def __init__(
        self,
        pdb_dir,
        output_fasta,
        output_aligned_dir=None,
        reference=None,
    ):
        self.pdb_dir = pdb_dir
        self.output_fasta = output_fasta
        self.output_aligned_dir = output_aligned_dir
        self.reference = reference

        self.pdb_files = []
        self.objects = OrderedDict()
        self.display_names = OrderedDict()

        pymol.finish_launching(["pymol", "-cq"])

    ##############################################################
    # Utility functions
    ##############################################################

    def discover_pdbs(self):
        self.pdb_files = sorted(
            glob.glob(os.path.join(self.pdb_dir, "*.pdb"))
        )

        if len(self.pdb_files) == 0:
            raise RuntimeError("No PDB files found.")

        if self.reference is None:
            self.reference = self.pdb_files[0]

        for i, pdb in enumerate(self.pdb_files):
            obj = f"mol{i}"

            cmd.load(pdb, obj)

            self.objects[pdb] = obj

            self.display_names[obj] = os.path.splitext(
                os.path.basename(pdb)
            )[0]

    ##############################################################
    # Loading
    ##############################################################

    def load_structures(self):
        cmd.delete("all")

        for pdb in self.pdb_files:
            obj = os.path.splitext(os.path.basename(pdb))[0]
            cmd.load(pdb, obj)
            self.objects[pdb] = obj

    ##############################################################
    # Sequence extraction
    ##############################################################

    def get_sequence(self, obj):
        """
        Returns:
            residues : list of residue numbers
            sequence : one-letter AA sequence
        """

        model = cmd.get_model(f"{obj} and name CA")

        residues = []
        seq = []

        aa = {
            "ALA":"A","ARG":"R","ASN":"N","ASP":"D",
            "CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
            "HIS":"H","ILE":"I","LEU":"L","LYS":"K",
            "MET":"M","PHE":"F","PRO":"P","SER":"S",
            "THR":"T","TRP":"W","TYR":"Y","VAL":"V"
        }

        seen = set()

        for atom in model.atom:

            key = (atom.chain, atom.resi)

            if key in seen:
                continue

            seen.add(key)

            residues.append(atom.resi)
            seq.append(aa.get(atom.resn, "X"))

        return residues, seq

    ##############################################################
    # Alignment
    ##############################################################

    def align_structure(self, mobile_obj, target_obj):
        """
        Align mobile to target.

        Returns
        -------
        dict returned by patched cealign()

        {
            "alignment_length": ...,
            "RMSD": ...,
            "rotation_matrix": ...,
            "mapping": [...]
        }
        """

        return cmd.cealign(target_obj, mobile_obj)
        # result = cmd.cealign(target_obj, mobile_obj)

        # return result

    ##############################################################
    # Residue mapping
    ##############################################################
    # def residue_mapping(self, target_obj, mobile_obj):
    #     """
    #     Build residue mapping using the alignment object.

    #     Returns
    #     -------
    #     list of tuples

    #     (target_resi, mobile_resi)
    #     """

    #     aln_name = "__tmp_alignment__"

    #     cmd.delete(aln_name)

    #     cmd.align(
    #         mobile_obj,
    #         target_obj,
    #         object=aln_name,
    #         cycles=0,
    #         transform=0,
    #     )

    #     mapping = []

    #     idx = cmd.index(aln_name)

    #     atoms = []

    #     for obj, index in idx:

    #         atom = cmd.get_model(f"{obj} and index {index}").atom[0]
    #         atoms.append((obj, atom.resi))

    #     for i in range(0, len(atoms), 2):
    #         mapping.append((atoms[i][1], atoms[i + 1][1]))

    #     cmd.delete(aln_name)

    #     return mapping

    def residue_lookup(self, obj):

        AA = {
            "ALA":"A","ARG":"R","ASN":"N","ASP":"D",
            "CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
            "HIS":"H","ILE":"I","LEU":"L","LYS":"K",
            "MET":"M","PHE":"F","PRO":"P","SER":"S",
            "THR":"T","TRP":"W","TYR":"Y","VAL":"V"
            }

        model = cmd.get_model(f"{obj} and guide")

        lookup = {}

        seen = set()

        for atom in model.atom:

            key = (atom.chain, atom.resi)

            if key in seen:
                continue

            seen.add(key)

            lookup[key] = AA.get(atom.resn, "X")

        return lookup
    ##############################################################
    # FASTA construction
    ##############################################################

    def build_alignment(self):

        ref_obj = self.objects[self.reference]

        ref_lookup = self.residue_lookup(ref_obj)

        ref_keys = list(ref_lookup.keys())

        aligned = OrderedDict()

        aligned[ref_obj] = list(ref_lookup.values())

        for pdb in self.pdb_files:

            if pdb == self.reference:
                continue

            obj = self.objects[pdb]

            result = self.align_structure(obj, ref_obj)

            seq = ["-"] * len(ref_keys)

            mobile_lookup = self.residue_lookup(obj)

            ref_index = {
                key:i
                for i,key in enumerate(ref_keys)
            }

            for pair in result["mapping"]:

                rk = (
                    pair["target_chain"],
                    pair["target_resi"]
                )

                mk = (
                    pair["mobile_chain"],
                    pair["mobile_resi"]
                )

                if rk not in ref_index:
                    continue

                if mk not in mobile_lookup:
                    continue

                seq[
                    ref_index[rk]
                ] = mobile_lookup[mk]

            aligned[obj] = seq

            if self.output_aligned_dir:

                os.makedirs(
                    self.output_aligned_dir,
                    exist_ok=True
                )

                cmd.save(
                    os.path.join(
                        self.output_aligned_dir,
                        obj + ".pdb"
                    ),
                    obj
                )

        return aligned

    # def build_alignment(self):

    #     ref_obj = self.objects[self.reference]

    #     ref_resi, ref_seq = self.get_sequence(ref_obj)

    #     aligned_sequences = OrderedDict()

    #     aligned_sequences[ref_obj] = list(ref_seq)

    #     ref_lookup = {
    #         r: i for i, r in enumerate(ref_resi)
    #     }

    #     for pdb in self.pdb_files:

    #         if pdb == self.reference:
    #             continue

    #         obj = self.objects[pdb]

    #         self.align_structure(obj, ref_obj)

    #         mobile_resi, mobile_seq = self.get_sequence(obj)

    #         mobile_lookup = {
    #             r: aa
    #             for r, aa in zip(mobile_resi, mobile_seq)
    #         }

    #         aligned = ["-"] * len(ref_seq)

    #         mapping = self.residue_mapping(ref_obj, obj)

    #         for ref_r, mob_r in mapping:

    #             if ref_r not in ref_lookup:
    #                 continue

    #             if mob_r not in mobile_lookup:
    #                 continue

    #             aligned[ref_lookup[ref_r]] = mobile_lookup[mob_r]

    #         aligned_sequences[obj] = aligned

    #         if self.output_aligned_dir is not None:

    #             os.makedirs(self.output_aligned_dir, exist_ok=True)

    #             out = os.path.join(
    #                 self.output_aligned_dir,
    #                 obj + "_aligned.pdb",
    #             )

    #             cmd.save(out, obj)

    #     return aligned_sequences

    ##############################################################
    # FASTA writer
    ##############################################################

    def write_fasta(self, aligned_sequences):
        with open(self.output_fasta, "w") as f:
            for name, seq in aligned_sequences.items():
                f.write(f">{name}\n")
                f.write("".join(seq))
                f.write("\n")

    ##############################################################
    # Main
    ##############################################################

    def run(self):
        self.discover_pdbs()
        self.load_structures()
        msa = self.build_alignment()
        self.write_fasta(msa)
        cmd.quit()

##############################################################
# Main
##############################################################

def main():

    builder = CEAlignFastaBuilder(
        pdb_dir="reservePDBfilegroups/CompositeEverySeqPred_uL24+bL19_alignment",
        output_fasta="uL24+bL19_alignment.fa",
        output_aligned_dir="cealign_alignment",
        reference=None,
    )

    builder.run()


if __name__ == "__main__":
    main()