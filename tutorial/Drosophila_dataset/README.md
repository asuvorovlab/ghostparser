# Drosophila Tutorial

This tutorial tests for ghost introgression in clade nine of the fly genus Drosophila. 

## Overview 
(From Tolman et al. in preparation)
Ghostbuster operates under the theory that gene trees with the topology A,C|B would have a longer outgroup branch length than genes supporting the topology A,B|C if ghost introgression has occurred(Fig. 1). The user provides, as input, paths to five separate files, one file each containing the putative “A”, “B”, “C”, and “Out” taxa (used to root the triplet) with one taxa per line, and a fifth file containing all of the gene trees. Within this framework, the A and B taxa should be sister, C taxa are hypothesized to have introgressed with A taxa, and B taxa are the hypothesized recipients of ghost introgression. The GhostBuster pipeline tests all possible combinations of A,B,C and Out taxa. For each triplet, the pipeline  first employs the Discordant Count Test (DCT) (Suvorov, Kim, et al. 2022; Suvorov, Scornavacca, et al. 2022) to test for any signal of introgression. This implementation of the DCT uses a two proportion Z-test to determine if the counts of the two discordant topologies for a given triplet are significantly different. If these counts are not significantly different (p > .01), Ghostbuster concludes that there is no evidence of introgression in the species triplet. The stringent cutoff (p < .01) is used to reduce the occurrence of false positives.

GhostBuster then uses ete3 v3.1.3 (Huerta-Cepas et al. 2016) and newick_utils v1.6 (Junier 2024) to extract branch length statistics from each gene tree. It then employs the ks test, implemented using the rstatix v0.7.2 package (Kassambara 2023) in R v4.4.2 (R Core Team 2021), to determine if the distribution of outgroup branch lengths of gene trees with the topology A,C|B (Fig. 1) significantly differ from the distribution of outgroup branch lengths of genes supporting the concordant topology  A,B|C. If the DCT test is significant (p < .01) and the distribution of A,C|B outgroup branch lengths is significantly larger than the distribution of A,B|C outgroup branch lengths, the “B” taxa is considered to have been a recipient of unsampled or ghost introgression with an older time T (Fig. 1). If the DCT test is passed, but the latter condition is not met, this is considered evidence of sampled introgression between the “A” and “C” taxa. Regardless of the conclusion, GhostBuster outputs the supported hypothesis and the results of the DCT and outgroup branch tests for each triplet tested.

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

3. **Tree Pruning & Standardization:**  
   - For each valid combination, the pipeline uses the `nw_prune` tool to extract a subtree from the gene trees that contains the selected taxa.
   - The pruned tree is saved to a temporary file (`triplet.txt`), and the taxon names are standardized (e.g., replaced with labels like `A`, `B`, `C`, and `Out`) to simplify subsequent analysis.

4. **Tree Statistics Calculation (Python):**  
   - The Python script `process_trees.py` processes the pruned trees from `triplet.txt`:
     - It reads each Newick-formatted tree, ensuring that all required taxa (`A`, `B`, `C`, and `Out`) are present.
     - The script uses the ETE3 library to root the tree at the `Out` taxon and prunes out any extraneous taxa (e.g., a "Ghost" taxon).
     - It identifies sister relationships among taxa `A`, `B`, and `C` by checking which pair shares the same parent node.
     - For each sister pair, it computes the distance from the non-sister taxon to the most recent common ancestor (MRCA) of `A`, `B`, and `C`.

5. **Branch-Length Statistical Analysis (R):**  
   - The R script `branch_stats2.0.r` then processes the TSV file generated by the Python script:
     - It reads the tree statistics output file.
     - Counts the occurrences of each topology (`AB|C`, `AC|B`, and `BC|A`).
     - It performs a two-proportion z-test comparing the proportions of `AC` versus `BC` topologies:
       - If the p-value is less than 0.01, the difference is considered statistically significant, and evidence of introgression.
       - Otherwise, it is deemed not significant.
     - The script conducts a Kolmogorov–Smirnov (KS) test on the **A** branch length data for topologies `AC|B` versus `AB|C`:
       - It calculates median branch lengths for these topologies and identifies which one has the higher median.
       - The KS test yields a D-statistic and p-value that are reported alongside the proportion test result. If the outgroup branch length for `AC|B` has a higher median than `AB|C`, and the distributions are significantly different (p < .05)
         this is considered evidence of unsampled or "ghost introgression" between an unsampled taxon and "B." If this is not the case, this is considered to be a case of introgression between the "A" and "C" taxa.
       - Finally, the R script prints a combined result  that summarizes both the proportion test and the KS test outcomes.

6. **Final Output:**  
   - For each taxa combination processed, the pipeline appends a summary of the computed branch metrics and statistical test results to the user-specified output file.
  
## Theoretical rationale
![Figure 1](https://github.com/e-tolman/images/blob/main/Figure_1.png?raw=true)
Theoretical framework of GhostBuster, with the expected coalescence of the “outgroup” taxa demonstrated. Outgroup branch lengths are shown from empirical systems identified with BPP as evolving without introgression  (Jaltomata: J. auriculata, J. yungayensis, J. biflora), under sampled introgression (Jaltomata: J. repandidentata, J. procumbens, J. darcyana), and under ghost introgression (Thuja: T. standishii, T. sutchuenensis, T. plicata). 

For a rooted species triplet of the topology AB|C, a majority of gene trees are expected to display the “concordant” topology. The gene trees displaying the two “discordant” topologies of AC|B and BC|A are expected to occur at similar frequencies in the absence of introgression (Suvorov, Kim, et al. 2022; Suvorov, Scornavacca, et al. 2022). Both sampled and ghost introgression can inflate the counts of one discordant topology over the other. If introgression occurs between C and A (Fig. 1), a higher proportion of the gene trees will support the topology AC|B than BC|A. Likewise, if B receives gene flow from a ghost lineage that predates the divergence of C from A and B, the introgressed loci will “slip” through this ghost lineage to the older coalescence point, making A and C appear as sister, elevating the count of gene trees that support AC|B (Fig. 1). In the case of an introgressed loci from C to A, B would be expected to coalesce with A and C no earlier than the divergence between A and B at time T (Fig. 1). For an introgressed loci from a ghost lineage to B, B would be expected to coalesce with A and C, at time T, no earlier than the divergence between the ghost lineage, and the triplet A,B,C (Fig. 1). Thus, ghost introgression should be identifiable by an older time T on average as compared to the coalescence time T in triplets uninfluenced by introgression (Fig. 1). Within this framework the branch length of B in the triplet A,C|B and C in the triplet A,B|C would be expected to be influenced by the timing of coalescence, and an indicator of ghost introgression.
