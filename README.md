# MIRNA-RIP

### MIRNA-RIP: The miRNA Regulatory Interaction Predictor

## About MIRNA-RIP

MIRNA-RIP is a software package for R (www.r-project.org) to predict miRNAs as regulators of potential target genes using matched gene- and miRNA expression profiles. MIRNARIP implements two model types using Mixed Integer Linear Programming: simple linear regression and piecewise linear regression. In both model types, specific parameters are optimized to estimate the target gene expression using the expression of the corresponding miRNA. Target gene predictions from both models are summarized in a combined model. More information about installation and usage including an example data set can be found in the manual.

## Downloads

### MIRNA-RIP 1.0

    MIRNA-RIP R package

    MIRNA-RIP manual
    start_MIRNARIP.R

### Example data

    Gene expression data, taken from [1]

    miRNA expression data, taken from [1]

    miRNA – target gene interactions, taken from [2]

    miRNA mature2pre-miRNA_ID mapping, extracted from [3]

    miRNA experimental2pre-miRNA_ID mapping, extracted from [3]

    miRNA candidates

Please cite: Ast et al. MiR-192, miR-200c and miR-17 are fibroblast-mediated inhibitors of colorectal cancer invasion. Oncotarget (2018).

### References

[1] The Cancer Genome Atlas Network.
Comprehensive Molecular Characterization of Human Colon and Rectal Cancer.
Nature. 2012;487(7407):330-337.
[2] Karagkouni et al.
DIANA-TarBase v8: a decade-long collection of experimentally supported miRNA–gene interactions.
Nucleic Acids Res. 2018 Jan 4; 46(Database issue): D239–D245.
[3] Kozomara and Griffiths-Jones
miRBase: annotating high confidence microRNAs using deep sequencing data.
Nucleic Acids Research, Volume 42, Issue D1, 1 January 2014, Pages D68–D73
