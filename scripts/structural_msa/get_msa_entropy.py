#!/usr/bin/env python

# import packages
from Bio import AlignIO
import pandas as pd
import numpy as np

# read in structural msa
struct_msa = list(AlignIO.read("../../data/structural_msa/result.afasta", "fasta"))
human_align_seq = [s for s in struct_msa if s.name.startswith("0_relaxed_")][0]
struct_msa_len = len(human_align_seq.seq)

# function to compute shannon entropy
def get_entropy(aligned_residues):
    unique_aligned_residues = set(aligned_residues)
    N = len(aligned_residues)
    entropy = 0
    for res in unique_aligned_residues:
        p = aligned_residues.count(res)/N
        entropy += p * np.log2(p)
    return -entropy

# gather data
count = 1
residue_pos_list, entropy_list = [], []
for i in range(struct_msa_len):
    if human_align_seq.seq[i] != "-":
        aligned_residues = [seq.seq[i] for seq in struct_msa]
        entropy = get_entropy(aligned_residues)
        residue_pos_list.append(count)
        entropy_list.append(entropy)
        count += 1

# write data to file
df = pd.DataFrame()
df["Human_L1_ORF2_Residue_Index"] = residue_pos_list
df["Structural_MSA_Shannon_Entropy"] = entropy_list
df.to_csv("../../results/structural_msa/struct_msa_shannon_entropy.tsv", sep="\t", index=False)
