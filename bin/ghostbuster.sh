#!/bin/bash

# Initialize variables
out_taxa=""
A_taxa=""
B_taxa=""
C_taxa=""
input_trees=""
output_file=""

# Function to display usage
usage() {
  echo "Usage: $0 --out_taxa <value> --A_taxa <value> --B_taxa <value> --C_taxa <value> --input_trees <value> --output_file <value>"
  echo
  echo "Options:"
  echo "  --out_taxa <value>     File containing each of the outgroup taxa, one per line."
  echo "  --A_taxa <value>       File containing each of the A taxa, one per line."
  echo "  --B_taxa <value>       File containing each of the B taxa, one per line."
  echo "  --C_taxa <value>       File containing each of the C taxa, one per line."
  echo "  --input_trees <value>  File containing the gene trees, one per line."
  echo "  --output_file <value>  User specified path to an output file."
  exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --out_taxa)
      out_taxa="$2"
      shift # Skip argument name
      shift # Skip argument value
      ;;
    --A_taxa)
      A_taxa="$2"
      shift
      shift
      ;;
    --B_taxa)
      B_taxa="$2"
      shift
      shift
      ;;
    --C_taxa)
      C_taxa="$2"
      shift
      shift
      ;;
    --input_trees)
      input_trees="$2"
      shift
      shift
      ;;
   --output_file)
      output_file="$2"
      shift
      shift
      ;;   
    *)
      echo "Error: Unknown option $key"
      usage
      ;;
  esac
done

# Check if all required arguments are provided
if [[ -z $out_taxa || -z $A_taxa || -z $B_taxa || -z $C_taxa || -z $input_trees || -z $output_file ]]; then
  echo "Error: Missing required arguments."
  usage
fi

#Check to see if all of the taxa are valid

# Array to hold taxa that are missing in the input trees
missing_taxa=()

# Array containing all taxa file variables
taxa_files=("$out_taxa" "$A_taxa" "$B_taxa" "$C_taxa")

# Loop through each taxa file
for taxa_file in "${taxa_files[@]}"; do
  # Read each taxon in the file
  while IFS= read -r taxon || [ -n "$taxon" ]; do
    # Skip empty lines
    if [[ -z "$taxon" ]]; then
      continue
    fi
    # Check if the taxon appears in any of the input trees.
    # The -w option makes grep match whole words.
    if ! grep -q -w "$taxon" "$input_trees"; then
      missing_taxa+=("$taxon")
    fi
  done < "$taxa_file"
done

# If any missing taxa were found, print an error message and exit.
if [ ${#missing_taxa[@]} -gt 0 ]; then
  echo "Error: The following taxa were not found in any input tree:"
  for taxon in "${missing_taxa[@]}"; do
    echo " - $taxon"
  done
  exit 1
fi


# Display parsed arguments
echo "Parsed arguments:"
echo "  out_taxa: $out_taxa"
echo "  A_taxa: $A_taxa"
echo "  B_taxa: $B_taxa"
echo "  C_taxa: $C_taxa"
echo "  input_trees: $input_trees"


#basically need to tell it all ABC and OUT taxa

echo -e "Triplet Tested\tConclusion\tDCT result\tBLT result\tIBL result\tSister_introggression" > $output_file
# Generate every possible combination of lines from the four files
#Add an if then statement, if this is not provided, just provide every triplet from the root node
while read a; do
    while read b; do
        while read c; do
            while read d; do
                # Ensure no value is repeated in the combination
                if [[ "$a" != "$b" && "$a" != "$c" && "$a" != "$d" && "$b" != "$c" && "$b" != "$d" && "$c" != "$d" ]]; then
                    echo "$a $b $c $d"
                fi
            done < $out_taxa
        done < $A_taxa
    done < $B_taxa
done < $C_taxa > possible_combinations.txt
awk '{
  delete seen; dup=0;
  for(i=1; i<=NF; i++){
    if($i in seen){ dup=1; break }
    seen[$i]=1
  }
  if(!dup) print
}' possible_combinations.txt > tmp && mv tmp possible_combinations.txt


#order is out, c, b, a


while read line; do
    out=$(echo $line | cut -d" " -f4)
    c=$(echo $line | cut -d" " -f1)
    b=$(echo $line | cut -d" " -f2)
    a=$(echo $line | cut -d" " -f3)
    nw_prune -v $input_trees $out $a $b $c | sed "s/$a/A/g" | sed "s/$b/B/g" | sed "s/$c/C/g" | sed "s/$out/Out/g" > triplet.txt

    # Get summary stats
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
    python $SCRIPT_DIR/process_trees.py triplet.txt "$input_trees.tree_stats.txt"
#    Sis_introgression=$(Rscript $SCRIPT_DIR/introgression_between_sister_pairs.r $input_trees.tree_stats.txt)
    #python $SCRIPT_DIR/pairwise_topology.py "$input_trees.tree_stats.txt" "$input_trees.means.csv" > "$input_trees.statistical_comparisons.txt"
    bl_stats=$(Rscript $SCRIPT_DIR/branch_stats2.0.r $input_trees.tree_stats.txt)

    # If results are significantly different
    result=$(echo $bl_stats | cut -d"~" -f1 | grep "are significantly different" | wc -l)
    chisq=$(echo $bl_stats | cut -d"~" -f1)
    out_pval=$(echo $bl_stats | cut -d"~" -f2 | cut -d":" -f4)
    out_higher=$(echo $bl_stats | cut -d"~" -f2 | cut -d";" -f1 | cut -d":" -f2 | cut -d";" -f1)
    blt=$(echo $bl_stats | cut -d"~" -f2)
    if [ "$result" -gt 0 ]; then
        # Convert scientific notation to decimal using awk
        out_pval_decimal=$(echo "$out_pval" | awk '{printf "%f", $1}')

        if [[ $(echo "$out_pval_decimal < 0.05" | bc -l 2>/dev/null) -eq 1 && "$out_higher" == AC ]]; then
		echo -e "$line\t Evidence of unsampled introgression\t $chisq\t $blt" >> $output_file
        else
                echo -e "$line\t Evidence of sampled introgression\t $chisq\t $blt" >> $output_file
        fi
    else
	    echo -e "$line\t No evidence of introgression\t $chisq\t $blt" >> $output_file
fi
done < possible_combinations.txt


#Parse through the output file, determine all the taxa that received unsampled introgression
#Go through and replace the "C" taxa, do any of them hit as sampled?

#rm possible_combinations.txt
#rm triplet.txt
#rm $input_trees.tree_stats.txt
#rm $input_trees.chisq_results.txt
