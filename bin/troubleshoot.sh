#!/bin/bash

base_dir=$1

# Compile all test datasets
SCRIPT_DIR=$(dirname "$0")
echo $SCRIPT_DIR
end=$2

echo $base_dir
ls $base_dir/ghost_and_sampled_msprime_scripts/training_set 
ls $base_dir/ghost_only_msprime_scripts/training_set 
ls $base_dir/ils_only_msprime_scripts/training_set 
ls $base_dir/sampled_only_msprime_scripts/training_set
# Call the R script with the compiled training set directories
Rscript $SCRIPT_DIR/train_ghostbuster.R $base_dir/ghost_and_sampled_msprime_scripts/training_set $base_dir/ghost_only_msprime_scripts/training_set $base_dir/ils_only_msprime_scripts/training_set $base_dir/sampled_only_msprime_scripts/training_set

