#!/bin/bash

# Export the command for GNU Parallel
export iqtree_cmd='iqtree -s {} -m TEST'

# Use GNU Parallel to run the command on all files in the specified directory
find "$1" -type f | parallel -j "$(nproc)" $iqtree_cmd


cat $1/*treefile > $2/calculated_trees.txt
