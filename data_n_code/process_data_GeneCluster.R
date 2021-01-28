#Call library
devtools::install_github("huynguyen250896/GeneCluster")
library(GeneCluster)

#load example data
miRNA = read.table("data_expression_miRNA.txt", sep = "\t", check.names = F, header = T, row.names = 1)

dim(miRNA)
# data is of dimension 398 miRNA genes x 300 subjects

#Hierarchical clustering of subjects
hc = hclust(dist(t(miRNA)))
memb <- cutree(hc, k = 3) #Suppose we want to have three separate clusters of patients

#Predictive subgroups
head(memb)
# TCGA-BH-A0BD-01 TCGA-AN-A0FJ-01 TCGA-A8-A06Q-01 TCGA-A8-A09D-01 TCGA-A1-A0SM-01 TCGA-B6-A0IE-01 
#       1               1               1               1               2               1 

# Run GeneCluster 
GeneCluster::SubtypeSpecificGene(omics = t(miRNA), cluster = memb)
#                GeneCluster       X1       X2       X3      P.value      Q.value
# hsa.miR.34a              3 6.185874 4.947188 6.255715 1.364334e-10 5.430049e-08
# hsa.miR.29b.2.           1 4.858117 4.166922 3.300045 5.898720e-10 1.173845e-07
# hsa.miR.708              1 4.625680 3.526900 4.337173 1.408336e-09 1.868393e-07
# hsa.miR.339.5p           3 3.158102 1.789311 3.587808 1.703861e-09 1.695342e-07
# hsa.miR.193b             1 6.163570 5.002492 6.052717 5.308790e-09 4.225797e-07
# hsa.miR.26b              3 9.174499 8.206727 9.221825 6.971175e-09 4.624213e-07
# 
# NOTE: 
# *the results shown above are incomplete.
# *subtype_specific_gene.csv placed in your current working directory.
# *Please check to identify which gene is specifically assigned to which subgroup.


