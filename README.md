# Complex structural variation, phylogeny, and disease associations of the mucin pangenome

This is a dedicated repo for custom scripts used to analyze complex structural variation in the mucin genes in the Human Pangenome Reference Consortium and Human Genome Structural Variation Consortium haplotypes. 

## Description

There are five custom scripts: 

### Snakefile_mapNHP
--> A Snakemake pipeline that measures VNTR-containing exon sizes in non-human primate (and other) genome assemblies using a flank-anchor approach. For each mucin locus, the pipeline extracts upstream and downstream flanking sequences from the T2T-CHM13v2 reference, automatically trims softmasked (repetitive) bases from the inner edges to produce clean anchors, and discards anchors shorter than a minimum length. Anchors are mapped to each target assembly with minimap2 (allowing secondary alignments), and best anchor pairs are selected per locus using a per-mucin maximum span constraint to prevent paralog mis-pairing in clustered mucin regions. The pipeline reports exon size, locus coordinates, and a status flag (including a special case for very small exons where the two anchors overlap) in a single TSV.

### vntr_motif_finding.R
--> Detects and tiles mucin VNTR repeat units within a protein sequence given a user-supplied canonical repeat motif. The analyze_repeats() function uses Biostrings::matchPattern (with tolerated mismatches against both the single and doubled canonical repeat) to locate all repeat hits, optionally splits them into separate blocks based on a maximum allowed gap, and then determines the best reading-frame phase by rotating the canonical motif and scoring matches. The repeat region is tiled into full and partial repeat units, and for each unit the script reports absolute coordinates, mismatch count, percent identity, percent variability, and the repeat sequence — enabling downstream analysis of VNTR motif diversity across haplotypes.

### infer_haplogroups.R
--> Infers haplogroups from a phylogenetic tree annotated with VNTR length metadata. The script computes pairwise patristic distances between tips, performs hierarchical clustering, and scans across a range of tree-cut heights to find the partitioning that maximizes the ANOVA R² between cluster assignment and VNTR length, while restricting clusters to those supported by a bootstrap threshold (default ≥70) and a minimum cluster size. Orphan tips dropped by the bootstrap filter are reassigned to the nearest surviving tip's cluster by patristic distance. The final clustering is renumbered sequentially and visualized as a tip-colored tree.

### number_genomes_95concordance_locityper.R
--> A command-line Rscript that estimates how many reference genomes are needed to reach 95% leave-one-out genotyping concordance for each locus, based on subsampling experiments at varying reference panel sizes. For every (locus, metric) combination, it fits a 3-parameter logistic saturation model — concordance(n) = A / (1 + exp(-k·(n − n0))) — using minpack.lm::nlsLM, then solves analytically for the panel size at which the fitted curve crosses the threshold. Each locus is classified as already achieved, achievable (with the extrapolated n95), or plateau (asymptote below threshold). Outputs include a per-fit results TSV, a dense prediction grid for plotting, and a faceted PDF showing observed points, fitted curves, and the 95% threshold line.

### CF_association_testing_allmodels.R
--> A command-line Rscript that tests whether mucin VNTR variations associated with the continuous KNoRMA phenotype in the Cystic Fibrosis Foundation's WGS dataset. each model adjusts for sex and ancestry PCs (1-10) computed previously using somalier. The script runs four models:
1.) Haplogroup dosage: linear regression of KNoRMA on the number of copies (0, 1, or 2) of each haplogroup, fit separately for each mucin and haplogroup.
2.) VNTR length: linear regression of KNoRMA on z-scored summed VNTR length (the sum of the two haplotype lengths, standardized per mucin).
3.) Binned median-split length dosage: linear regression of KNoRMA on the number of haplotypes (0, 1, or 2) falling in each median-split length bin (short or long).
4.) Binned length-extremes dosage: linear regression of KNoRMA on the number of haplotypes (0, 1, or 2) falling in the bottom 20%, middle 60%, or top 20% of the per-mucin length distribution.




## Required tools/packages

### infer_haplogroups.R (R)
ape — for cophenetic.phylo() and getMRCA()
ggtree — for tree plotting and node data extraction
dplyr — for %>%, filter(), select(), mutate(), left_join(), n_distinct(), pull()

### vntr_motif_finding.R (R / Bioconductor)
Biostrings (Bioconductor) — for AAString() and matchPattern()

### Snakefile_mapNHP (Python + command-line tools)
Python packages: pandas, plus stdlib subprocess, re, collections
Workflow manager: snakemake
Command-line tools (must be on $PATH): samtools (for faidx) and minimap2

### number_genomes_95concordance_locityper.R (R)
optparse — command-line argument parsing
dplyr — data manipulation
tidyr — long-format reshaping
ggplot2 — diagnostic plot
minpack.lm — provides nlsLM() for the Levenberg–Marquardt logistic fit

### CF_association_testing_allmodels.R
dplyr - data manipulation
broom - tidy() and glance() to extract coefficients and model fit from lm()


## Help/Contact

For questions, contact Lizzie Plender at plendere@uw.edu

