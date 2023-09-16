#!/usr/bin/env python

# import packages
import pandas as pd
import numpy as np
from sklearn import manifold
import argparse

# input arguments
parser = argparse.ArgumentParser()
parser.add_argument("chain_type", help="chain type", choices=["RT_chains", "bioassembly_chains"])
ref_types = ["ORF2-full", "ORF2-N", "ORF2-tower", "ORF2-RT", "ORF2-fingers", "ORF2-palm", "ORF2-thumb", "ORF2-wrist", "ORF2-C"]
parser.add_argument("reference_pdb", help="reference pdb", 
                    choices = ref_types)
parser.add_argument("distance_metric", help="distance metric", 
                    choices = ["Normalized_Perplexity"])
args = parser.parse_args()
chain_type = args.chain_type
reference_pdb = args.reference_pdb
distance_metric = args.distance_metric

# read in data
df = pd.read_csv(f"../../data/perplexity/{chain_type}_mmligner_stats.tsv", sep="\t")
df["Perplexity"] = 2**(-df["Compression"])

exclude_orf2 = [x for x in ref_types if x != reference_pdb]
df = df[(~df["PDB_1"].str.fullmatch("|".join(exclude_orf2))) & (~df["PDB_2"].str.fullmatch("|".join(exclude_orf2)))]

df = df[(~df["PDB_1"].str.contains("class")) & (~df["PDB_2"].str.contains("class"))]
df = df[(~df["PDB_1"].str.contains("crystal")) & (~df["PDB_2"].str.contains("crystal"))]

# compute information-theoretic values
df["Normalized_Compression_Distance"] = (df["I(A & <S,T>)"] - df["NULL(T)"])/df["NULL(S)"]
df["Normalized_Perplexity"] = 2**(df["Normalized_Compression_Distance"])

# make matrix
unique_names = list(set(list(df["PDB_1"].unique()) + list(df["PDB_2"].unique())))
mat = []
for name_1 in unique_names:
    row = []
    for name_2 in unique_names:
        if name_1 != name_2:
            val_1 = df[(df["PDB_1"] == name_1) & (df["PDB_2"] == name_2)][distance_metric].to_numpy()
            val_2 = df[(df["PDB_1"] == name_2) & (df["PDB_2"] == name_1)][distance_metric].to_numpy()
            for val in [val_1, val_2]:
                if val.size > 0:
                    val = val[0]
                    row.append(val)
                    break
        else:
            row.append(0)
    mat.append(row)

dist_mat = np.array(mat)

mean_dist, std_dist = np.mean(dist_mat), np.std(dist_mat)
dist_mat = (dist_mat - mean_dist)/std_dist

# 2D distance cluster
mds_model = manifold.MDS(n_components = 2, random_state = 123, 
            dissimilarity = 'precomputed',
            normalized_stress='auto', n_jobs=-1, n_init=1000, max_iter=1000)

mds_coords = mds_model.fit_transform(dist_mat)

# save coordinates to file
df_share = pd.DataFrame()
df_share["X_coord"] = mds_coords[:,0]
df_share["Y_coord"] = mds_coords[:,1]
df_share["Names"] = unique_names
df_share.to_csv(f"../../results/perplexity/{chain_type}_{reference_pdb}_{distance_metric}_coords.tsv", sep="\t", index=False)
