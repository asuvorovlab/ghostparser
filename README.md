# Ghostbuster Pipeline

The Ghostbuster Pipeline is a tool designed to analyze gene trees for evidence of introgression between taxa. The pipeline processes multiple taxa files and gene trees, validates the taxa, generates all possible combinations, and runs downstream analyses using a mix of Bash, Python, and R scripts.

## Features

- **Taxa Validation:** Ensures that each taxon listed in your taxa files is present in the input gene trees.
- **Combination Generation:** Creates all valid combinations from provided taxa files.
- **Tree Pruning:** Uses `nw_prune` to extract relevant portions of the gene trees.
- **Statistical Analysis:** Runs Python and R scripts to generate summary statistics and test for introgression.
- **Customizable Output:** Results are written to a user-specified output file.

## Requirements

- **Git** – For version control.
- **Bash** – To run the main pipeline script.
- **nw_prune** – A tool for pruning phylogenetic trees.  
  (Installation instructions can be found on the [nw_prune GitHub page](https://github.com/josephryan/nw_prune) or your package manager.)
- **Python 3** – For running the `process_trees.py` script.
- **R** – For running R scripts (`introgression_between_sister_pairs.r` and `branch_stats2.0.r`).
- **Conda** – To easily install all the required dependencies from the `environment.yml` file.

## Installation

### 1. Clone the Repository

Clone the repository to your local machine using either SSH or HTTPS:

Using SSH (recommended if you have a passkey set up):
```bash
git clone git@github.com:asuvorovlab/ghostbuster.git

git clone https://github.com/asuvorovlab/ghostbuster.git
```
This repository includes an environment.yml file that lists all necessary dependencies. Follow these steps to set up the environment:

Ensure Conda is installed.
You can download Miniconda if you don't have it.
#### 2. Download the dependencies using conda
```bash
conda env create -f environment.yml
```

## Usage
The pipeline is executed using a the bash script bin/ghostbuster.sh. The script requires several arguments that specify input files and an output file.

Command-Line Arguments
--out_taxa: Path to a file containing each of the outgroup taxa (one per line).
--A_taxa: Path to a file containing each of the A taxa (one per line).
--B_taxa: Path to a file containing each of the B taxa (one per line).
--C_taxa: Path to a file containing each of the C taxa (one per line).
--input_trees: Path to a file containing the gene trees (one per line).
--output_file: Path where the output will be written.
Example
Run the pipeline with your files as follows:

```bash
bin/ghostbuster.sh --out_taxa data/out_taxa.txt --A_taxa data/A_taxa.txt --B_taxa data/B_taxa.txt --C_taxa data/C_taxa.txt --input_trees data/input_trees.txt --output_file results/output.txt
```
What the Script Does
Argument Parsing & Validation:
Checks that all required arguments are provided and validates that each taxon in the taxa files appears in the gene trees.

Combination Generation:
Generates all possible combinations of taxa (ensuring no repetition) and writes them to possible_combinations.txt.

Tree Processing:
For each combination, the script prunes the input gene trees using nw_prune and substitutes taxa names. It then runs:

A Python script (process_trees.py) to summarize tree statistics.
Two R scripts (introgression_between_sister_pairs.r and branch_stats2.0.r) to assess introgression and branch statistics.
Output:
The results, including evidence for sampled or unsampled introgression, are written to the user-specified output file.

Additional Notes
Script Location:
Ensure that all required scripts (pipeline_script.sh, process_trees.py, introgression_between_sister_pairs.r, and branch_stats2.0.r) are located in the repository’s root or correctly referenced in the script.

Temporary Files:
The script generates temporary files (e.g., possible_combinations.txt and triplet.txt). Modify or uncomment the cleanup commands at the bottom of the script as needed.

Troubleshooting:
If you encounter issues, ensure that:

The Conda environment is activated.
All dependencies in environment.yml are correctly installed.
File paths provided as arguments are correct.
