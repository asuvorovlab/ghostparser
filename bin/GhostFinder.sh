#!/bin/bash

# Initialize variables
out_taxa=""
ghostbuster_output=""
input_trees=""
output_file=""

# Function to display usage
usage() {
  echo "Usage: $0 --out_taxa <value> --ghostbuster_output <value> --input_trees <value> --output_file <value>"
  echo
  echo "Options:"
  echo "  --out_taxa <value>             File with outgroup taxa, one per line."
  echo "  --ghostbuster_output <value>   GhostBuster output file."
  echo "  --input_trees <value>          File containing the gene trees, one per line."
  echo "  --output_file <value>          Path to the final output file."
  exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --out_taxa)
      out_taxa="$2"
      shift 2
      ;;
    --ghostbuster_output)
      ghostbuster_output="$2"
      shift 2
      ;;
    --input_trees)
      input_trees="$2"
      shift 2
      ;;
    --output_file)
      output_file="$2"
      shift 2
      ;;   
    *)
      echo "Error: Unknown option $key"
      usage
      ;;
  esac
done

# Check if all required arguments are provided
if [[ -z "$out_taxa" || -z "$ghostbuster_output" || -z "$input_trees" || -z "$output_file" ]]; then
  echo "Error: Missing required arguments."
  usage
fi

# Generate taxa files from ghostbuster_output.
awk -F'\t' '$2==" Evidence of unsampled introgression" {print $1}' "$ghostbuster_output" | cut -d" " -f2 | sort | uniq > GhostFinder_A_taxa.txt 
awk -F'\t' '$2==" Evidence of unsampled introgression" {print $1}' "$ghostbuster_output" | cut -d" " -f3 | sort | uniq > GhostFinder_B_taxa.txt 
awk -F'\t' '$2==" Evidence of unsampled introgression" {print $1}' "$ghostbuster_output" | cut -d" " -f1 | sort | uniq > additional_taxa_to_exclude.txt

# Build GhostFinder_C_taxa.txt by excluding taxa found in the other files.
# (Note: Adjust the filename condition if out_taxa is not literally "Out.txt".)
awk 'FNR==NR { ex[$1]; next }
     FILENAME=="additional_taxa_to_exclude.txt" { ex[$1]; next }
     FILENAME=="GhostFinder_A_taxa.txt" { ex[$1]; next }
     FILENAME=="Out.txt" { ex[$1]; next }
     { gsub(/[(),;]/, " ");
       for(i=1; i<=NF; i++) {
         sub(/:.*/, "", $i);  # Remove branch lengths
         if($i != "" && !($i in ex) && $i !~ /^[0-9]+$/)
           print $i
       }
     }' GhostFinder_B_taxa.txt additional_taxa_to_exclude.txt GhostFinder_A_taxa.txt "$out_taxa" "$input_trees" | sort -u > GhostFinder_C_taxa.txt

# Run GhostBuster.sh from the same directory as this script.
python "$(dirname "$0")/ghostbuster.py" --out_taxa "$out_taxa" --A_taxa GhostFinder_A_taxa.txt --B_taxa GhostFinder_B_taxa.txt --C_taxa GhostFinder_C_taxa.txt --input_trees "$input_trees" --output_file Full_GhostFinder_output.txt

# Extract lines with " sampled" from the GhostBuster output.
grep " sampled" Full_GhostFinder_output.txt > "$output_file"

