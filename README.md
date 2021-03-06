# GeneCluster v0.1.0
#### I. Introduction
---
The package GeneCluster is built to serve as a support tool for the paper "*[Multi-omics analysis detects novel prognostic subgroups of breast cancer](https://www.frontiersin.org/articles/10.3389/fgene.2020.574661/full?utm_source=F-NTF&utm_medium=EMLX&utm_campaign=PRD_FEOPS_20170000_ARTICLE#F5)*". </br>
A previous study defined subtype-specific genes are the ones mutated predominantly in the samples assigned to one single subtype than in the other subtypes [1]. Subsequently, those genes are features that reflect the difference between subgroups of heterogeneous cancers [1, 2]. To computationally detect subtype-specific genes, we built the R package GeneCluster from the idea of the reference paper [3]. In brief, given a gene from a list of genes of interest, it will be specifically distributed to either of predictive subgroups based on the mean values (e.g., CNA changes, MET changes, and expression levels). Then, a gene was considered as a subtype-specific one if P-value <= 0.05 (one-way ANOVA test). </br>

*[1] Cyll, K., et al., Tumour heterogeneity poses a significant challenge to cancer biomarker research. British journal of cancer, 2017. 117(3): p. 367-375.*

*[2] Alizadeh, A.A., et al., Toward understanding and exploiting tumor heterogeneity. Nature medicine, 2015. 21(8): p. 846-853.*

*[3] Shen, R., et al., Integrative Subtype Discovery in Glioblastoma Using iCluster. PLOS ONE, 2012. 7(4): p. e35236.*

#### II. Understanding the tool
---
The following are parameters provided by GeneCluster:
- omics: data.frame or matrix. The first input data includes its rows are samples and its columns are genes.

- cluster: Predictive subgroups correspond with samples after running a clustering tool for your own data (e.g., k-means, hclust,...).

- adjustedP: logical. Whether we should adjust the gained P-values (One-way ANOVA test) using the Benjamini-Hochberg procedure. Default is `adjustedP = T`.

Please see [data_n_code](https://github.com/huynguyen250896/GeneCluster/tree/master/data_n_code) to grasp a data format required by GeneCluster and its usage well.

#### III. Pipeline
---
![Figure](https://imgur.com/rMuwnLh.png)
**Figure:** Pipeline of the package GeneCluster.

#### IV. Implementation
---
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
SubtypeSpecificGene(omics = exp, cluster = groups)
```

#### IV. What's new
---
- 2021-01-26 : The function now can adjust gained P-values using the Benjamini-Hochberg procedure.

#### V. Citation 
---
Please kindly cite the following paper (and Star this Github repository if you find this tool of interest) if you use the tool in this repo: </br>
```sh
Author: Nguyen, Quang-Huy
Nguyen, Hung
Nguyen, Tin
Le, Duc-Hau
Year: 2020
Title: Multi-omics analysis detects novel prognostic subgroups of breast cancer
Journal: Frontiers in Genetics
Type of Article: ORIGINAL RESEARCH
DOI: 10.3389/fgene.2020.574661
```
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) <huynguyen96.dnu AT gmail DOT com> for any questions about the code and results.

