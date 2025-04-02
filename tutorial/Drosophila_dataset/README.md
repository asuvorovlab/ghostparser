# _Drosophila_ ghostbuster Tutorial

This tutorial tests for ghost introgression in clade nine of the fly genus _Drosophila_. Suvorov et al. (2021) found evidence of introgression between the clade comprised of _D. arawakana_, _D. dunni_, _D. funebris_ etc. and _D.immigrans_. Here we test if this is in fact introgression between these two lineages, or introgression involving _D. immigrans_, _D._ albomicans_, _D. neonasuta_,etc. and a "ghost clade." In this scenario, _D. arawakana_, _D. dunni_, _D. funebris_ etc. is the A clade, hypothesized to have introgressed with the C clade, _D. immigrans_. The clade sister to _D. immigrans_ is the B clade, hypothesized to have introgressed with a ghost lineage. 

#Add a picture.

## Step one: run ghostbuster
To run ghostbuster, we need to first prepare several input files.   

The first four files specify the taxa hypothesized to be the A,B,C taxa, as well as the outgroup taxa used to root triplets. They contain the name of one taxa per line, as they appear in the gene trees.   

In this example, the files are called "A.txt", "B.txt", "C.txt", and "Out.txt" respectively.

```bash
bin/ghostbuster.sh --out_taxa data/out_taxa.txt --A_taxa data/A_taxa.txt --B_taxa data/B_taxa.txt --C_taxa data/C_taxa.txt --input_trees data/input_trees.txt --output_file results/output.txt
```

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

