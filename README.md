# RMCL
Kate Zielinski,
Enright Laboratory, 
Dept. Of Pathology,
University of Cambridge.

## Introduction
RMCL is a set of R utilities for expression-based clustering of gene-expression data from large-scale microarray or NGS studies.

RMCL utilises Markov Clustering (MCL) as implemented by Stijn van Dongen (micans.org). A Gene-expression matrix is transformed into an expression correlation matrix using a correlation function (e.g. Pearson correlation). This produces a correlation matrix where for each gene or transcript it's expression correlation across all samples is computed against *all* other genes in the original matrix.

This correlation matrix is transformed into a weighted network where genes are connected if they exhibit a positive correlation above a specific value. Markov clustering is then performed on this network to detect clusters of genes with shared expression properties.

Finally the clusters obtained are passed back to R to allow the user to expore gene-expression and Gene-Ontology characteristics of the gene-expression clusters.

## Steps Involved

1. The user provides a gene-expression matrix or dataframe in R (rows are genes, columns are samples). This data is usually normalised.
2. RMCL takes this matrix as input and produces an mcl formatted temporary text file
3. RMCL launches *mcxarray* from the MCL package to perform large-scale memory-efficient correlation analysis.
4. A correlation cut-off is applied to produce a network, typically this is set between 0.7 and 0.95 depending on sample number.
5. The *mcl* clustering is applied and takes as its main parameter the *Inflation* value, usually set between 2 and 5.
6. Finally clusters are returned to R as an R dataframe where each row represents the cluster number for each row in the input matrix.


## Usage


