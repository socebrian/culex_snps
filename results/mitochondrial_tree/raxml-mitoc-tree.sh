#!/bin/bash

# Input file (aligned fasta)
INPUT="aligned_mitoc_seq_pipiens1234.fasta"

# Prefixes for outputs
BOOTSTRAP_NAME="bootstraps"
BEST_TREE_NAME="best_tree"
FINAL_NAME="final_tree"

# Model
MODEL="GTRGAMMA"

# Random seeds
SEED_START=12345
SEED_BOOT=54321

# Bootstrap replicates
BS_REPS=1000

echo "Step 1: Generating bootstrap trees..."
raxmlHPC -s "$INPUT" -m "$MODEL" -p $SEED_START -b $SEED_BOOT -# $BS_REPS -n "$BOOTSTRAP_NAME"

echo "Step 2: Searching for best ML tree..."
raxmlHPC -s "$INPUT" -m "$MODEL" -p $SEED_START -n "$BEST_TREE_NAME"

echo "Step 3: Mapping bootstrap support onto best tree..."
raxmlHPC -f b -m "$MODEL" -t "RAxML_bestTree.$BEST_TREE_NAME" -z "RAxML_bootstrap.$BOOTSTRAP_NAME" -n "$FINAL_NAME"

echo "All done. Final tree with bootstrap support: RAxML_bipartitions.$FINAL_NAME"
