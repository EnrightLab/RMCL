# RMCL
Kate Zielinski (kez22@cam.ac.uk),
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

Given an input gene-exprsession matrix such as this:
```
> head(input_matrix)
                X1029479_RD_CORE_S30 X1337985_CM_MEM_S28 X1343208_RD_WASH_S59 X1440723_RD_WASH_S58 X1466854_RD_WASH_S9   .....
ENSG00000000003              0.00000          170.366754             75.73753            213.88457            88.79548   .....
ENSG00000000005              0.00000            0.000000              0.00000              0.00000             0.00000   .....
ENSG00000000419            105.34711           99.050439             88.20067            145.62353           111.71044   .....
ENSG00000000457             60.19835          134.708597            117.92046            167.23953           111.13756   .....
ENSG00000000460             62.34829           26.743618             42.18293             38.68125            55.56878   .....
ENSG00000000938           1337.26333            8.914539            160.10338            122.86986           173.00796   .....
```

You perform the clustering analysis as follows:

```
mcl_cluster(input_matrix, inflation=3.0, corr="pearson", threshold=0.7)
```

In this case we run correlation analysis, filtering and MCL clustering with the following parameters: 
a) ```threshold=0.7``` - correlation cutoff at >0.7 - Genes with higher than this value are kept in the weighted expression network
b) ```corr="pearson"``` - We choose Pearson correlation, other options are (Pearson, Cosine, Jaccard).
c) ```inflation=3.0``` - Use 3.0 as the Inflation parameter (granularity setting) for MCL, its main parameter

The output is as follows:

```
[1] "Clustering:32675" "Clustering:65"   
[1] "Inflation:3"
[1] "Correlation Function:pearson"
[1] "Processing Matrix using mcxarray to generate correlations"
[mcxarray] read table with 65 rows and 32675 columns
[mcxarray] 2123875 entries in table
___ [mcxarray] constant data for label <ENSG00000183169> - no pearson
___ [mcxarray] constant data for label <ENSG00000199595> - no pearson
___ [mcxarray] constant data for label <ENSG00000214514> - no pearson
___ [mcxarray] constant data for label <ENSG00000224107> - no pearson
___ [mcxarray] constant data for label <ENSG00000228842> - no pearson
___ [mcxarray] constant data for label <ENSG00000248444> - no pearson
___ [mcxarray] constant data for label <ENSG00000268823> - no pearson
___ [mcxarray] constant data for label <ENSG00000270066> - no pearson
0000000000000000000
[mclIO] writing </var/folders/mp/16qzjn013vx0yj1gv6hwq_rc0000gn/T//RtmpBrL2J2mcl_matrix_corr.mx>
.......................................
[mclIO] wrote native interchange 32675x32675 matrix with 3430923 entries to stream </var/folders/mp/16qzjn013vx0yj1gv6hwq_rc0000gn/T//RtmpBrL2J2mcl_matrix_corr.mx>
[1] "Running Markov Clustering with MCL on correlation matrix"
[mclIO] reading </var/folders/mp/16qzjn013vx0yj1gv6hwq_rc0000gn/T//RtmpBrL2J2mcl_matrix_corr.mx>
.......................................
[mclIO] read native interchange 32675x32675 matrix with 3430923 entries
[mcl] pid 97338
 ite -------------------  chaos  time hom(avg,lo,hi) m-ie m-ex i-ex fmv
  1  ................... 130.57  2.29 4902.96/0.00/159226336.00 7.18 4.05 4.05  45
  2  ................... 219.12 12.86 0.82/0.06/48.73 11.19 0.76 3.07  71
  3  ................... 148.96  5.72 0.79/0.02/10.69 6.26 0.39 1.21  71
  4  ...................  91.03  1.41 0.94/0.06/60.33 2.63 0.33 0.40  53
  5  ...................  29.79  0.42 0.92/0.09/39.12 1.41 0.38 0.15  24
  6  ...................  12.63  0.17 0.86/0.16/14.38 1.08 0.51 0.08  11
  7  ...................   6.93  0.04 0.94/0.33/5.19 1.02 0.53 0.04   0
  8  ...................   3.13  0.02 0.98/0.23/1.69 1.00 0.43 0.02   0
  9  ...................   2.15  0.01 0.99/0.73/1.70 1.00 0.76 0.01   0
 10  ...................   1.12  0.01 1.00/0.76/1.47 1.00 0.89 0.01   0
 11  ...................   1.92  0.01 1.00/0.56/1.00 1.00 0.99 0.01   0
 12  ...................   1.74  0.01 1.00/0.76/1.00 1.00 0.99 0.01   0
 13  ...................   0.04  0.01 1.00/0.96/1.00 1.00 1.00 0.01   0
 14  ...................   0.07  0.01 1.00/1.00/1.00 1.00 1.00 0.01   0
 15  ...................   0.22  0.01 1.00/0.96/1.00 1.00 1.00 0.01   0
 16  ...................   0.71  0.01 0.99/0.72/1.00 1.00 1.00 0.01   0
 17  ...................   1.89  0.01 0.99/0.36/1.00 1.00 1.00 0.01   0
 18  ...................   0.93  0.01 0.99/0.64/1.00 1.00 1.00 0.01   0
 19  ...................   0.96  0.01 1.00/0.50/1.00 1.00 0.88 0.01   0
 20  ...................   0.19  0.01 1.00/0.92/1.00 1.00 0.97 0.01   0
 21  ...................   0.00  0.01 1.00/1.00/1.00 1.00 0.99 0.01   0
[mcl] jury pruning marks: <91,94,98>, out of 100
[mcl] jury pruning synopsis: <92.6 or scrumptious> (cf -scheme, -do log)
[mcl] output is in /var/folders/mp/16qzjn013vx0yj1gv6hwq_rc0000gn/T//RtmpBrL2J2mcl_matrix.clusters
[mcl] 12067 clusters found
[mcl] output is in /var/folders/mp/16qzjn013vx0yj1gv6hwq_rc0000gn/T//RtmpBrL2J2mcl_matrix.clusters

Please cite:
    Stijn van Dongen, Graph Clustering by Flow Simulation.  PhD thesis,
    University of Utrecht, May 2000.
       (  http://www.library.uu.nl/digiarchief/dip/diss/1895620/full.pdf
       or  http://micans.org/mcl/lit/svdthesis.pdf.gz)
OR
    Stijn van Dongen, A cluster algorithm for graphs. Technical
    Report INS-R0010, National Research Institute for Mathematics
    and Computer Science in the Netherlands, Amsterdam, May 2000.
       (  http://www.cwi.nl/ftp/CWIreports/INS/INS-R0010.ps.Z
       or  http://micans.org/mcl/lit/INS-R0010.ps.Z)

[1] "Clustering Complete!"
```


## Output Results

The data returned is an R Dataframe with Genes as Rows and Cluster assignments for each gene in a *Cluster* named column. Cluster IDs are numeric with the largest cluster starting at 1.

The object looks like this:
```
                Cluster
ENSG00000000003     587
ENSG00000000005      60
ENSG00000000419      18
ENSG00000000457    2666
ENSG00000000460      55
ENSG00000000938       6
ENSG00000000971     418
ENSG00000001036     588
ENSG00000001084      33
ENSG00000001167    2667
ENSG00000001460       1
ENSG00000001461     101
ENSG00000001497    1417
ENSG00000001561    2668
ENSG00000001617       6
ENSG00000001626       3
ENSG00000001629    2669
ENSG00000001630    2670
ENSG00000001631     892
ENSG00000002016    2671
```

## Accessory Functions

The package provides toolkits to explore the shared co-expression within clusters and also to perform Gene Ontology Analysis via the *clusterProfiler* package.

## Expression Analysis
