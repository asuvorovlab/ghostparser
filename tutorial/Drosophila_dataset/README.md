# Drosophila Tutorial

This tutorial tests for ghost introgression in clade nine of the fly genus _Drosophila_. Suvorov et al. (2021) found evidence of introgression between the clade comprised of _D. arawakana_, _D. dunni_, _D. funebris_ etc. and _D.immigrans_. Here we test if this is in fact introgression between these two lineages, or introgression involving _D. immigrans_, _D._ albomicans_, _D. neonasuta_,etc. and a "ghost clade." In this scenario, _D. arawakana_, _D. dunni_, _D. funebris_ etc. is the A clade, hypothesized to have introgressed with the C clade, _D. immigrans_. The clade sister to _D. immigrans_ is the B clade, hypothesized to have introgressed with a ghost lineage.

#Add a picture.

## Step one: run ghostbuster
Ghostbuster operates under the theory that gene trees with the topology A,C|B would have a longer outgroup branch length than genes supporting the topology A,B|C if ghost introgression has occurred(Fig. 1). The user provides, as input, paths to five separate files, one file each containing the putative “A”, “B”, “C”, and “Out” taxa (used to root the triplet) with one taxa per line, and a fifth file containing all of the gene trees. Within this framework, the A and B taxa should be sister, C taxa are hypothesized to have introgressed with A taxa, and B taxa are the hypothesized recipients of ghost introgression. The GhostBuster pipeline tests all possible combinations of A,B,C and Out taxa. For each triplet, the pipeline  first employs the Discordant Count Test (DCT) (Suvorov, Kim, et al. 2022; Suvorov, Scornavacca, et al. 2022) to test for any signal of introgression. This implementation of the DCT uses a two proportion Z-test to determine if the counts of the two discordant topologies for a given triplet are significantly different. If these counts are not significantly different (p > .01), Ghostbuster concludes that there is no evidence of introgression in the species triplet. The stringent cutoff (p < .01) is used to reduce the occurrence of false positives.

GhostBuster then uses ete3 v3.1.3 (Huerta-Cepas et al. 2016) and newick_utils v1.6 (Junier 2024) to extract branch length statistics from each gene tree. It then employs the ks test, implemented using the rstatix v0.7.2 package (Kassambara 2023) in R v4.4.2 (R Core Team 2021), to determine if the distribution of outgroup branch lengths of gene trees with the topology A,C|B (Fig. 1) significantly differ from the distribution of outgroup branch lengths of genes supporting the concordant topology  A,B|C. If the DCT test is significant (p < .01) and the distribution of A,C|B outgroup branch lengths is significantly larger than the distribution of A,B|C outgroup branch lengths, the “B” taxa is considered to have been a recipient of unsampled or ghost introgression with an older time T (Fig. 1). If the DCT test is passed, but the latter condition is not met, this is considered evidence of sampled introgression between the “A” and “C” taxa. Regardless of the conclusion, GhostBuster outputs the supported hypothesis and the results of the DCT and outgroup branch tests for each triplet tested.

## Step two: consider the ghostbuster output

- **Git** – For version control.
- **Bash** – To run the main pipeline script.
- **nw_prune** – A tool for pruning phylogenetic trees.  
  (Installation instructions can be found on the [nw_prune GitHub page](https://github.com/josephryan/nw_prune) or your package manager.)
- **Python 3** – For running the `process_trees.py` script.
- **R** – For running R scripts (`introgression_between_sister_pairs.r` and `branch_stats2.0.r`).
- **Conda** – To easily install all the required dependencies from the `environment.yml` file.

## Step 3: Run ghostfinder to search for the ghost lineage in the drosophila tree

```bash
bin/ghostbuster.sh --out_taxa data/out_taxa.txt --A_taxa data/A_taxa.txt --B_taxa data/B_taxa.txt --C_taxa data/C_taxa.txt --input_trees data/input_trees.txt --output_file results/output.txt
```

## Step 4: Look at the ghostfinder output
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
## Step 4: Look at the ghostfinder output

