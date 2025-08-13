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

The first column, Taxon_C is the "C" taxon, Taxon_B is the "B" taxon in the triplet, and Taxon_A is "A" taxon in the triplet. Taxon_C is the outgroup in the triplet. Taxon_B appears as sister to Taxon_C more frequently than Taxon_A does in the gene tree set. DCT_statistic is the Z-score from the two-sample Z-test used test if the count of discordant topology one (B,C)A is significantly greater than the count of discordant topology two (A,C)B in the gene trees. DCT_p_value is the p-value from this test, ghostparser considers any value below 0.01 to be significant. If the p-value is larger than 0.01, Test_conclusion will be "No evidence of introgression." AB_count, BC_count, and AC_count are the counts of gene trees reflecting the triplet topologies (A,B)C, (B,C)A and (A,C)B respectively. THT_statistic is the test statistic of a KS test comparing the distributions of triplet tree heights of concordant gene trees ((A,B)C) and discordant-one gene trees ((B,C)A). THT_p_value is the p-value from this test, with 0.05 used as the cutoff for significance. Assuming the Z-test is passed, Test conclusion will be "Evidence of sampled introgression, likely not involving inflow" if the THT test is not significant, "Evidence of sampled introgression, possibly involving inflow" if the distribution of tree heights of discordant-one gene trees is significantly _smaller_ than the distribution of tree heights of concordant gene trees, and "Evidence of unsampled introgression" if the distribution of tree heights of discordant-one gene trees is significantly _larger_ than the distribution of tree heights of concordant gene trees. Here, we find strong evidence that _D. repleta_ is likely the recipient of ghost introgression!



The ghostparser results for all tested triplets can be summarized by counting the unique values in the fourteenth column of the output file.

```bash
cut -f14 tutorial_output.txt | grep -v "Conclusion" | sort | uniq -c
```

## Step four: Run ghostfinder to search for the ghost lineage in the drosophila tree
Once ghostparser is complete, ghostfinder can take the output of the ghostparser run, consider the taxa hypothesized as recipients of recipients of ghost introgression, and determine if there is evidence of sampled introgression involving the putative ghost recipients in the species tree. It takes as input, the path to the ghostparser output file you would like to test, paths to the species tree gene trees used to run ghostparser, a path to the output file, and (optionally) the number of threads to use.  Ghostfinder takes each pairing of “A” and “B”, where the “A” taxon is the hypothesized recipient of Ghost introgression, and the outgroup in discordant topology one. It then tests all possible “C” taxa from the species tree to determine if any “C” taxa will flip the assignments of the sister taxa in the triplet. In instances where the putative ghost recipient is reassigned as the “B” taxon, ghostfinder uses GhostParser to determine which of these “C” taxa show evidence of inflow introgression with the putative ghost, and flag these introgression donors as putative ghost lineages.

Ghostfinder can be run like so:


```bash
python ~/ghostparser_final/bin/ghostfinder.py --ghostparser_output ghostparser_output.txt --species_tree species.tree --gene_trees drosophila.trees --output_file ghostfinder_output.txt --threads 8
```
As the ghostfinder pipeline is testing many more triplets than the first ghostparser run, it should take between 15 and 30 minutes to run on 8 threads. 

The output file is the output of the ghostparser run, tested on triplets where the "C" taxon could possible be the "ghost" lineage. The first few lines of this output file are copied below. Here, each "C" taxon is the possible "ghost" lineage, and the "Test_conclusion" tells us whether there is evidence that this possible "ghost" lineage introgressed with the ghost introgression recipient (now Taxon_B).

```bash
Taxon_C Taxon_B Taxon_A Taxon_out       DCT_statistic   DCT_p_value     AB_count        BC_count        AC_count        THT_statistic   THT_p_value     AB_median       BC_median       Test_conclusion
D_affinis       D_repleta       D_arizonae      Anopheles_gambiae       0.2231  6.3666e-01      2184    10      8       0.4617  2.8691e-02      0.357083        0.250138        No evidence of introgression
D_americana     D_repleta       D_arizonae      Anopheles_gambiae       0.4034  5.2535e-01      2344    22      18      0.3946  2.2537e-03      0.149138        0.214234        No evidence of introgression
D_ananassae     D_repleta       D_arizonae      Anopheles_gambiae       0.2508  6.1649e-01      2384    9       7       0.5488  9.0257e-03      0.424612        0.267989        No evidence of introgression
D_asahinai      D_repleta       D_arizonae      Anopheles_gambiae       0.5313  4.6607e-01      2399    10      7       0.4595  2.9820e-02      0.390431        0.266370        No evidence of introgression
D_athabasca     D_repleta       D_arizonae      Anopheles_gambiae       0.0590  8.0803e-01      2346    9       8       0.5358  1.1624e-02      0.357962        0.239299        No evidence of introgression
D_aurauria      D_repleta       D_arizonae      Anopheles_gambiae       1.0033  3.1650e-01      2389    10      6       0.4447  3.8944e-02      0.388156        0.268297        No evidence of introgression
D_azteca        D_repleta       D_arizonae      Anopheles_gambiae       0.2508  6.1649e-01      2388    9       7       0.5281  1.3471e-02      0.356711        0.238022        No evidence of introgression
D_bakoue        D_repleta       D_arizonae      Anopheles_gambiae       1.0034  3.1650e-01      2369    10      6       0.4641  2.7430e-02      0.404875        0.275614        No evidence of introgression
D_bifasciata    D_repleta       D_arizonae      Anopheles_gambiae       0.2230  6.3673e-01      2424    10      8       0.4465  3.7691e-02      0.348984        0.248426        No evidence of introgression
D_bipectinata   D_repleta       D_arizonae      Anopheles_gambiae       1.0033  3.1651e-01      2393    10      6       0.4713  2.3977e-02      0.423895        0.292385        No evidence of introgression
D_birchii       D_repleta       D_arizonae      Anopheles_gambiae       0.5313  4.6608e-01      2430    10      7       0.4329  4.7838e-02      0.408202        0.299717        No evidence of introgression
D_bocqueti      D_repleta       D_arizonae      Anopheles_gambiae       0.5313  4.6606e-01      2364    10      7       0.4638  2.7577e-02      0.409935        0.275414        No evidence of introgression
D_bromeliae     D_repleta       D_arizonae      Anopheles_gambiae       0.0698  7.9162e-01      2359    30      28      0.342   1.9569e-03      0.182749        0.146744        No evidence of introgression
D_bunnanda      D_repleta       D_arizonae      Anopheles_gambiae       1.0033  3.1651e-01      2412    10      6       0.4524  3.3933e-02      0.404414        0.277463        No evidence of introgression
D_burlai        D_repleta       D_arizonae      Anopheles_gambiae       1.1462  2.8434e-01      2352    9       5       0.5264  1.3914e-02      0.405312        0.262112        No evidence of introgression
D_cyrtoloma     D_repleta       D_arizonae      Anopheles_gambiae       0.3622  5.4727e-01      2001    14      11      0.4689  4.4205e-03      0.205342        0.149819        No evidence of introgression
D_elegans       D_repleta       D_arizonae      Anopheles_gambiae       0.5313  4.6608e-01      2429    10      7       0.4378  4.3963e-02      0.397830        0.276320        No evidence of introgression
D_ercepeae      D_repleta       D_arizonae      Anopheles_gambiae       0.2231  6.3672e-01      2399    10      8       0.4754  2.2203e-02      0.416983        0.289683        No evidence of introgression
D_erecta        D_repleta       D_arizonae      Anopheles_gambiae       1.0034  3.1648e-01      2328    10      6       0.4686  2.5249e-02      0.428392        0.291743        No evidence of introgression
D_funebris      D_repleta       D_arizonae      Anopheles_gambiae       0.0402  8.4108e-01      2421    13      12      0.395   3.5345e-02      0.233241        0.309983        No evidence of introgression
D_fuyamai       D_repleta       D_arizonae      Anopheles_gambiae       0.5313  4.6607e-01      2405    10      7       0.4403  4.2060e-02      0.404704        0.299363        No evidence of introgression
```

We do not identify any other lineages in our species tree that have introgressed with the ghost introgression recipient. Thus, we can conclude that the "ghost" lineage is likely an unsampled sister lineage to the clade we originally tested in ghostparser.

## Citations
Suvorov A, Kim BY, Wang J, Armstrong EE, Peede D, D’Agostino ERR, et al. Widespread introgression across a phylogeny of 155 Drosophila genomes. Curr Biol. 2022 Jan 10;32(1):111-123.e5. 


