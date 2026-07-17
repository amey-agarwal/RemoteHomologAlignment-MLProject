from itertools import combinations


def read_fasta_alignment(filename):
    """
    Reads an aligned FASTA file.

    Returns
    -------
    names : list[str]
    seqs  : list[str]
    """
    names = []
    seqs = []

    with open(filename) as f:
        name = None
        seq = []

        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    names.append(name)
                    seqs.append("".join(seq))

                name = line[1:]
                seq = []
            else:
                seq.append(line)

        if name is not None:
            names.append(name)
            seqs.append("".join(seq))

    return names, seqs


def build_alignment_maps(names, aligned):
    """
    For every sequence and every alignment column,
    determine which original residue occupies that column.

    Example

    A-CD-E

    returns

    column : residue index

      0 -> 0
      1 -> None
      2 -> 1
      3 -> 2
      4 -> None
      5 -> 3

    """
    maps = {}

    for name, seq in zip(names, aligned):

        residue_index = 0
        mapping = []

        for aa in seq:

            if aa == "-":
                mapping.append(None)
            else:
                mapping.append(residue_index)
                residue_index += 1

        maps[name] = mapping

    return maps


def aligned_pairs(names, aligned):
    """
    Produce every residue pair that is aligned.

    Each residue is identified by

        (sequence_name, residue_index)

    NOT by amino acid.
    """

    maps = build_alignment_maps(names, aligned)

    L = len(aligned[0])

    pairs = set()

    for col in range(L):

        residues = []

        for name in names:

            idx = maps[name][col]

            if idx is not None:
                residues.append((name, idx))

        # all pairwise combinations
        for a, b in combinations(residues, 2):

            if a > b:
                a, b = b, a

            pairs.add((a, b))

    return pairs


def exact_column_sets(names, aligned):
    """
    Each column is represented as the set of residue identities.

    This ignores the column index completely.
    """

    maps = build_alignment_maps(names, aligned)

    columns = []

    L = len(aligned[0])

    for col in range(L):

        column = frozenset(
            (name, maps[name][col])
            for name in names
            if maps[name][col] is not None
        )

        if len(column) >= 2:
            columns.append(column)

    return columns


def evaluate_alignment(gt_file, pred_file):

    gt_names, gt = read_fasta_alignment(gt_file)
    pr_names, pr = read_fasta_alignment(pred_file)

    if gt_names != pr_names:
        raise ValueError("Sequence order/names differ.")

    gt_pairs = aligned_pairs(gt_names, gt)
    pr_pairs = aligned_pairs(pr_names, pr)

    TP = len(gt_pairs & pr_pairs)
    FP = len(pr_pairs - gt_pairs)
    FN = len(gt_pairs - pr_pairs)

    precision = TP / (TP + FP) if TP + FP else 0
    recall = TP / (TP + FN) if TP + FN else 0

    if precision + recall == 0:
        f1 = 0
    else:
        f1 = 2 * precision * recall / (precision + recall)

    gt_cols = exact_column_sets(gt_names, gt)
    pr_cols = set(exact_column_sets(pr_names, pr))

    correct_cols = sum(c in pr_cols for c in gt_cols)

    column_score = correct_cols / len(gt_cols)

    return {
        "TP": TP,
        "FP": FP,
        "FN": FN,
        "Precision": precision,
        "Recall": recall,
        "F1": f1,
        "Exact Column Score": column_score,
    }


if __name__ == "__main__":

    gt = "reserveFastaFiles/ancestral_seqs/orig_expected_align.fa"
    pred = "reserveFastaFiles/seqTrimAfterFoldSeek/orig_seqs_aligned-seqchopAfterFoldseek.fa"

    print(gt.split('/')[-1])
    print(pred)

    results = evaluate_alignment(gt, pred)

    print("\nAlignment Evaluation")
    print("--------------------")

    for k, v in results.items():

        if isinstance(v, float):
            print(f"{k:20s}: {v:.4f}")
        else:
            print(f"{k:20s}: {v}")