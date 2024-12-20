#!/bin/bash
# Exit on error and debug output
set -e
set -x

# Check for base directory argument
if [ -z "$1" ]; then
    echo "Error: No base directory provided."
    exit 1
fi

base_dir=$1

# Compile all test datasets
SCRIPT_DIR=$(dirname "$0")
mkdir -p "$base_dir/ghost_and_sampled_msprime_scripts/training_set"
mkdir -p "$base_dir/ghost_only_msprime_scripts/training_set"
mkdir -p "$base_dir/ils_only_msprime_scripts/training_set"
mkdir -p "$base_dir/sampled_only_msprime_scripts/training_set"

end=$2
for i in $(seq 1 "$end")
do

	mv "$base_dir/ghost_and_sampled_msprime_scripts/simulation_$i/summary_stats.csv" \
   		"$base_dir/ghost_and_sampled_msprime_scripts/training_set/sim_"$i".csv"
	mv "$base_dir/ghost_only_msprime_scripts/simulation_$i/summary_stats.csv" \
		"$base_dir/ghost_only_msprime_scripts/training_set/sim_"$i".csv"
	mv "$base_dir/ils_only_msprime_scripts/simulation_$i/summary_stats.csv" \
   		"$base_dir/ils_only_msprime_scripts/training_set/sim_"$i".csv"
	mv "$base_dir/sampled_only_msprime_scripts/simulation_$i/summary_stats.csv" \
   		"$base_dir/sampled_only_msprime_scripts/training_set/sim_"$i".csv"
done


# Call the R script with the compiled training set directories
Rscript $SCRIPT_DIR/train_ghostbuster.R $base_dir/ghost_and_sampled_msprime_scripts/training_set $base_dir/sampled_only_msprime_scripts/training_set $base_dir/ils_only_msprime_scripts/training_set $base_dir/ghost_only_msprime_scripts/training_set

mv $base_dir/ghost_only_msprime_scripts/training_set/*pdf $base_dir
mv $base_dir/ghost_only_msprime_scripts/training_set/*RData $base_dir
