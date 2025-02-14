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
### What the Script Does

1. **Input Parsing & Taxa Validation:**  
   - The pipeline begins by parsing command-line arguments to obtain file paths for the outgroup taxa (`--out_taxa`), three ingroup taxa (`--A_taxa`, `--B_taxa`, and `--C_taxa`), the gene trees (`--input_trees`), and the output file (`--output_file`).
   - It then validates that every taxon listed in the taxa files is present in the gene trees. If any taxon is missing, the pipeline prints an error message and exits.

2. **Combination Generation:**  
   - The script generates every possible combination by selecting one taxon from each of the four taxa files, ensuring no taxon is repeated within a combination.
   - These combinations are saved to a temporary file (`possible_combinations.txt`) for further processing.

3. **Tree Pruning & Standardization:**  
   - For each valid combination, the pipeline uses the `nw_prune` tool to extract a subtree from the gene trees that contains the selected taxa.
   - The pruned tree is saved to a temporary file (`triplet.txt`), and the taxon names are standardized (e.g., replaced with labels like `A`, `B`, `C`, and `Out`) to simplify subsequent analysis.

4. **Tree Statistics Calculation (Python):**  
   - The Python script `process_trees.py` processes the pruned trees from `triplet.txt`:
     - It reads each Newick-formatted tree, ensuring that all required taxa (`A`, `B`, `C`, and `Out`) are present.
     - The script uses the ETE3 library to root the tree at the `Out` taxon and prunes out any extraneous taxa (e.g., a "Ghost" taxon).
     - It identifies sister relationships among taxa `A`, `B`, and `C` by checking which pair shares the same parent node.
     - For each sister pair, it calculates half the distance between them (as a measure of their branch separation).
     - It then computes additional branch metrics:
       - The distance from the non-sister taxon to the most recent common ancestor (MRCA) of `A`, `B`, and `C`.
       - The internal branch length from the MRCA of `A`, `B`, and `C` to the MRCA of the sister pair.
       - The internal branch length from the MRCA of all four taxa (`A`, `B`, `C`, and `Out`) to the MRCA of `A`, `B`, and `C`.
     - The computed values are written to a tab-separated values (TSV) file with the columns: **Topology**, **Distance**, **Out**, **Internal**, and **Internal_ABCO**.

5. **Branch-Level Statistical Analysis (R):**  
   - The R script `branch_stats2.0.r` then processes the TSV file generated by the Python script:
     - It reads the TSV file and ensures that a **Topology** column is present.
     - The script counts the occurrences of each topology (`AC`, `BC`, and `AB`).
     - It performs a two-proportion z-test comparing the proportions of `AC` versus `BC` topologies:
       - If the p-value is less than 0.01, the difference is considered statistically significant.
       - Otherwise, it is deemed not significant.
     - In addition, the script conducts a Kolmogorov–Smirnov (KS) test on the **Out** branch length data for topologies `AC` versus `AB`:
       - It calculates median branch lengths for these topologies and identifies which one has the higher median.
       - The KS test yields a D-statistic and p-value that are reported alongside the proportion test result.
     - Finally, the R script prints a combined result (separated by a tilde `~`) that summarizes both the proportion test and the KS test outcomes.

6. **Final Output:**  
   - For each taxa combination processed, the pipeline appends a summary of the computed branch metrics and statistical test results to the user-specified output file.

