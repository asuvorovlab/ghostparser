# GhostParser: A highly scalable phylogenomic approach for the identification of ghost introgression 

This respository contains supplementary information for "GhostParser: A highly scalable phylogenomic approach for the identification of ghost introgression" and instructions for running GhostParser, and its follow up analysis, GhostFinder. 

## GhostParser Overview (from Tolman and Suvorov (2025))

The GhostParser pipeline represents a python script which uses a combination of python and R languages. In order to launch GhostParser analysis the user provides the following: (i) an input file with gene trees including branch lengths, one per line in newick format, (ii) a text file that contains a single outgroup taxon  to root the triplet, and (iii) the rooted species tree topology in newick format. The user may also provide an optional file with the specific triplets of interest for testing. If this optional file is provided, GhostParser will only analyze each user-specified triplet. Otherwise, it will, by default, analyze every possible triplet from the species tree (excluding the triplets that contain an outgroup taxon). GhostParser  will classify each of the three taxa in the triplet as the A, B, or C taxon. For any species triplet, the A and B taxa are always defined as sister in concordant gene tree , whereas the B taxon is sister to C in the most common discordant topology.
GhostFinder is an optional follow-up analysis to GhostParser. It uses the output from GhostParser—which identifies recipients of ghost introgression—to search for the potential “ghost” lineage in the species tree. It identifies, from the GhostParser output file, the triplets where evidence of unsampled introgression has been found. It takes each pairing of A and  B, where the A taxon is the hypothesized to be a putative recipient of ghost introgression, and the outgroup in Gdis1 . It then tests all possible C taxa from the species tree to determine if any C taxa will change the assignments of the sister taxa in the triplet. In instances where the putative ghost recipient is reassigned as the B taxon, GhostFinder uses GhostParser to determine which of these C taxa show evidence of inflow introgression with the putative ghost, and flag these introgression donors as putative ghost lineages. 


## bin
Contains the ghostparser.py and ghostfinder.py scripts, as well as the scripts "branch_stats2.0.r", and "process_trees.py" which are necessary to run ghosotparser.py.

## empirical_datasets
Contains the ghostparser runs one the _Jaltomata_, _Thuja_, _Drosophila_, and _Heliconius_ datasets. The simulated gene trees, alignments, and GhostParser output are all available on FigShare (10.6084/m9.figshare.29921969).

## tutorial

The README.rmd in this directories provides all instructions and files necessary for running ghostparser on a _Drosophila_ datasets, as well as the GhostParser output for this run. Refer to this README to run GhostParser on your dataset.

