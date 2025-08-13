# _Drosophila_ ghostparser Tutorial

This tutorial tests for ghost introgression in clade seven of the fly genus _Drosophila_. Suvorov et al. (2022) found evidence of introgression between _D. hydei_ and a clade comprised of _D. arizonae_, _D. mojavensis_, and _D. navojoa_ using the Branch Length Test (BLT) and Discordant Count Test (DCT). Here we recreate the "ghostparser" run to test if this is in fact introgression between these two lineages, or introgression involving _D. repleta_ and an unsampled "ghost" lineage. 
![Figure 1](https://github.com/e-tolman/images/blob/main/Drosophila_github_figure.png?raw=true)


## Step one: Install ghostparser

### 1. Clone the Repository

Clone the repository to your local machine using either SSH or HTTPS:

Using SSH (recommended if you have a passkey set up):
```bash
git clone git@github.com:asuvorovlab/ghostparser.git

git clone https://github.com/asuvorovlab/ghostparser.git
```
This repository includes an environment.yml file that lists all necessary dependencies. Follow these steps to set up the environment:

Ensure Conda is installed.
You can download Miniconda if you don't have it.

#### 2. Download the dependencies using conda
```bash
conda env create -f environment.yml
```

## Step two: Run ghostparser
To run ghostparser, we need to first prepare several input files. The file "drosophila.trees" contains the gene trees for the _Drosophila_ dataset, and "species.tree" is the maximum likelihood species tree, both from Suvorov et al. (2022). "triplets.csv" contains every possible triplet combination we can use to test our two competing hypotheses. Each triplet contains _D. hydei_ which is the only representative of one of the clades hypothesized to have been involved in a sampled introgression event; one of _D. arizonae_, _D. mojavensis_, or _D. navojoa_, which are all in the second clade hypothesized to have been involved in a sampled introgression event; and _D. repleta_ which is the hypothesized ghost recipient. If the triplets argument is not provided, ghostparser will automatically test every possible triplet in the species tree. "Out.txt" contains the outgroup taxon which will be used to root every tested triplet.

With these input files prepared, ghostparser can be executed like so:

```bash
python ~/ghostparser_final/bin/ghostparser.py --out_taxa Out.txt --input_trees drosophila.trees --species_tree species.tree --triplets triplets.txt --output_file ghostparser_output.txt
```

By default, ghostparser will utilize all available threads. This can be changed using the optional --threads option. Using eight threads, this execution should take 1-2 minutes.

## Step three: Examine the ghostparser output
Following the completion the ghostparser, the file "ghostparser_output.txt" will have been created. The file is printed below:

```bash
Taxon_C Taxon_B Taxon_A Taxon_out       DCT_statistic   DCT_p_value     AB_count        BC_count        AC_count        THT_statistic   THT_p_value     AB_median       BC_median       Test_conclusion
D_hydei D_arizonae      D_repleta       Anopheles_gambiae       361.2245        1.5239e-80      1115    953     364     0.0657  2.3577e-02      0.082073        0.085841        Evidence of unsampled introgression
D_hydei D_mojavensis    D_repleta       Anopheles_gambiae       368.7658        3.4747e-82      1163    969     368     0.0703  1.0813e-02      0.083483        0.086902        Evidence of unsampled introgression
D_hydei D_navojoa       D_repleta       Anopheles_gambiae       366.4096        1.1323e-81      1168    980     377     0.0656  2.0400e-02      0.083401        0.086669        Evidence of unsampled introgression
D_hydei D_seriema       D_repleta       Anopheles_gambiae       324.3658        1.6216e-72      1079    905     358     0.0736  9.6875e-03      0.084029        0.087615        Evidence of unsampled introgression
```

The first column, Taxon_C is the "C" taxon, Taxon_B is the "B" taxon in the triplet, and Taxon_A is "A" taxon in the triplet. Taxon_C is the outgroup in the triplet. Taxon_B appears as sister to Taxon_C more frequently than Taxon_A does in the gene tree set. DCT_statistic is the Z-score from the two-sample Z-test used test if the count of discordant topology one (B,C)A is significantly greater than the count of discordant topology two (A,C)B in the gene trees. DCT_p_value is the p-value from this test, ghostparser considers any value below 0.01 to be significant. If the p-value is larger than 0.01, Test_conclusion will be "No evidence of introgression." AB_count, BC_count, and AC_count are the counts of gene trees reflecting the triplet topologies (A,B)C, (B,C)A and (A,C)B respectively. THT_statistic is the test statistic of a KS test comparing the distributions of triplet tree heights of concordant gene trees ((A,B)C) and discordant-one gene trees ((B,C)A). THT_p_value is the p-value from this test, with 0.05 used as the cutoff for significance. Assuming the Z-test is passed, Test conclusion will be "Evidence of sampled introgression, likely not involving inflow" if the THT test is not significant, "Evidence of sampled introgression, possibly involving inflow" if the distribution of tree heights of discordant-one gene trees is significantly _smaller_ than the distribution of tree heights of concordant gene trees, and "Evidence of unsampled introgression" if the distribution of tree heights of discordant-one gene trees is significantly _larger_ than the distribution of tree heights of concordant gene trees



The ghostparser results for all tested triplets can be summarized by counting the unique values in the fourteenth column of the output file.

```bash
cut -f14 tutorial_output.txt | grep -v "Conclusion" | sort | uniq -c
```

## Step four: Run ghostfinder to search for the ghost lineage in the drosophila tree
Once ghostparser is complete, ghostfinder can take the output of the ghostparser run, consider the taxa hypothesized as recipients of recipients of ghost introgression, and determine if there is evidence of sampled introgression involving the putative ghosts 

```bash
python ghostparser/bin/ghostfinder.py --out_taxa Out.txt --ghostparser_output tutorial_output.txt --input_trees drosophila.trees --output_file ghostfinder_output.txt
```
As the ghostfinder pipeline is testing many more triplets than the first ghostparser run, it should take between 15 and 30 minutes to run on 8 threads.

## Step five: Examine at the ghostfinder output
The ghostfinder output file will either contain ghostparser results, listing lines where another taxa shows in the tree shows evidence of introgressing with the putative ghost recipient, or will display the line "No putative ghost lineages found." In this case, there is not evidence that the putative ghost recipients introgressed with another lineage in sampled _Drosophila_. The mystery continues!

## Citations
Suvorov A, Kim BY, Wang J, Armstrong EE, Peede D, D’Agostino ERR, et al. Widespread introgression across a phylogeny of 155 Drosophila genomes. Curr Biol. 2022 Jan 10;32(1):111-123.e5. 


