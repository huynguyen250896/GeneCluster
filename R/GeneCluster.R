#' @title Tool for Identification of subtype-specific genes
#'
#' @description Given a gene from a list of genes of interest, it will be specifically distributed to either of the identified subgroups based on its mean values (e.g., CNA changes, MET changes, and expression levels). Then, a gene was considered as a subtype-specific one if P-value <= 0.05 (one-way ANOVA test). For further information on requirements as well as how to implement this tool, please visit my Github repository: https://github.com/huynguyen250896/GeneCluster
#'
#' @param omics,cluster
#'
#' @return NULL
#'
#' @examples SubtypeSpecificGene(omics = df, cluster = groups)
#'
#' @export

# Compute mean
SubtypeSpecificGene = function(omics = NULL, cluster = NULL){
  #library
  library(maps)
  library(purrr)
  library(tidyr)
  library(tidyverse)
  library(dplyr)

  #subtype
  group = data.frame(cluster)
  data = omics %>% data.frame()
  data = data[rownames(group),]
  data$sub = group$cluster

  #implementation
  Mean <- round(aggregate(as.matrix(data[,-length(data)]) ~ data[,length(data)], FUN = mean),3)
  Mean=t(Mean)
  colnames(Mean)<-Mean[1,]
  Mean<-Mean[-1,]
  #compute p-value
  lis_out <- data[,-length(data)] %>% purrr::map(~ aov(.x ~ data[,length(data)] )) %>% purrr::map(~summary(.))
  Pvalue_lis<- lis_out %>% purrr::map(~.x[[1]][5] %>% unlist)

  Pvalue_out <- do.call(rbind, Pvalue_lis) %>% as.data.frame()
  colnames(Pvalue_out)[1] = "P.Value"
  #merge mean and p-value
  mean_overall = data.frame(Mean,Pvalue_out[1])
  #create sub matrix with adj.p.val<0.05 ~ GOOD statistically significant difference
  sig = subset(mean_overall, P.Value <0.05)

  #find specific group of each gene
  sig1 = sig %>%
    mutate_all(as.numeric) %>%
    mutate(row = row_number()) %>%
    pivot_longer(-row) %>%
    group_by(row) %>%
    mutate(max_value = max(value),
           group_number = str_extract(name,"[:digit:]") %>% as.numeric(),
           GeneCluster = if_else(value == max_value ,group_number,NA_real_)) %>%
    fill(GeneCluster,.direction = c("updown")) %>%
    select(-group_number,-max_value) %>%
    pivot_wider(names_from = name,values_from = value)
  #tidy to a dataframe
  sig1 <- as.data.frame(sig1)
  sig1 = sig1[,-1]
  rownames(sig1) = rownames(sig);

  #write the results as the csv file
  write.table(sig1,"specific_subtype_gene.txt",sep = "\t", quote = FALSE)

  #return the results
  print(head(sig1))

  #warning
  print(writeLines("\nNOTE: \n*the results, which are being shown above, are as incomplete ones.\n*specific_subtype_gene.csv placed in your current working directory.\n*Please check to identify which genes are specifically assigned to which subgroups."))
}



