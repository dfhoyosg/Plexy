#!/usr/bin/bash

# go across groups of protein structures
for chain_type in "bioassembly_chains" "RT_chains"
do
    # full ORF2 or RT domain
    for ref in "ORF2-full" "ORF2-RT"
    do
        # run script
        python make_dist_mat.py $chain_type $ref "Normalized_Perplexity"
    done
done
