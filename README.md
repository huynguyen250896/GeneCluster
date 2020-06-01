# GeneCluster v0.1.0
#### I. Introduction
- The package geneCor is built to serve as a support tool for the paper "...".
- A previous study defined subtype-specific genes are the ones mutated predominantly in the samples assigned to one single subtype than in the other subtypes [1]. Subsequently, those genes are features that reflect the difference between subgroups of heterogeneous cancers [2]. To computationally detect subtype-specific genes, we built the R package GeneCluster from the idea of the reference paper [3]. In brief, given a gene from a list of genes of interest, it will be specifically distributed to either of the identified subgroups based on the mean values (e.g., CNA changes, MET changes, and expression levels). Then, a gene was considered as a subtype-specific one if P-value <= 0.05 (one-way ANOVA test).


#### IV. Implementation
Use the following command to install directly from GitHub;
```sh
devtools::install_github("huynguyen250896/GeneCluster")
```
Call the library;
```sh
library(GeneCluster)
```
running example:
```sh
SubtypeSpecificGene(omics = df, cluster = groups)
```

