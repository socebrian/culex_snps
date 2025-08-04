#!/bin/bash
# Run IQ-TREE
echo "Starting IQ-TREE with model selection and 1000 bootstraps..."
/home/scamison/miniconda3/envs/trees/bin/iqtree2 -s  aligned_mitoc_seq_pipiens1234.fasta -m MFP -b 1000 -nt 8 -seed 12345 -pre mitoc_pipiens1234_tree
echo "Finished."
echo "Tree with bootstrap support: pREFIX.treefile"
