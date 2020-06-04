# GeneCluster v0.1.0
#### I. Introduction
- The package geneCor is built to serve as a support tool for the paper "...".
- A previous study defined subtype-specific genes are the ones mutated predominantly in the samples assigned to one single subtype than in the other subtypes [1]. Subsequently, those genes are features that reflect the difference between subgroups of heterogeneous cancers [1, 2]. To computationally detect subtype-specific genes, we built the R package GeneCluster from the idea of the reference paper [3]. In brief, given a gene from a list of genes of interest, it will be specifically distributed to either of the identified subgroups based on the mean values (e.g., CNA changes, MET changes, and expression levels). Then, a gene was considered as a subtype-specific one if P-value <= 0.05 (one-way ANOVA test).

#### IV. Reference
[1] Cyll, K., et al., Tumour heterogeneity poses a significant challenge to cancer biomarker research. British journal of cancer, 2017. 117(3): p. 367-375.
[2] Alizadeh, A.A., et al., Toward understanding and exploiting tumor heterogeneity. Nature medicine, 2015. 21(8): p. 846-853.
[3] Shen, R., et al., Integrative Subtype Discovery in Glioblastoma Using iCluster. PLOS ONE, 2012. 7(4): p. e35236.

#### V. Implementation
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

#### VI.Citation 
Please kindly cite the two repositories if you use the code, datasets or any results in this repo: </br>
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3686388.svg)](https://doi.org/10.5281/zenodo.3686388)
```sh
@software{quang_huy_nguyen_2020_3686388,
  author       = {Nguyen, Quang-Huy},
  title        = {huynguyen250896/computeQ: computeQ v 0.1.0},
  month        = feb,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {0.1.0},
  doi          = {10.5281/zenodo.3686388},
  url          = {https://doi.org/10.5281/zenodo.3686388}
}
```
</br> And </br>
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3829388.svg)](https://doi.org/10.5281/zenodo.3829388)
```sh
@software{nguyen_quang_huy_2020_3829388,
  author       = {Nguyen, Quang-Huy},
  title        = {huynguyen250896/geneSA: GeneSA v 0.1.0},
  month        = may,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {10.5281/zenodo.3829388},
  url          = {https://doi.org/10.5281/zenodo.3829388}
}
```
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) <huynguyen96.dnu AT gmail DOT com> for any questions about the code and results.

