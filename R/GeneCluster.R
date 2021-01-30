#' @title GeneCluster: Identification of subtype-specific genes
#'
#' @description A previous study defined subtype-specific genes are the ones mutated predominantly in the samples assigned to one single subtype than in 
#' the other subtypes [1]. Subsequently, those genes are features that reflect the difference between subgroups of heterogeneous cancers [1, 2]. To 
#' computationally detect subtype-specific genes, we built the R package GeneCluster from the idea of the reference paper [3]. In brief, given a gene from 
#' a list of genes of interest, it will be specifically distributed to either of the identified subgroups based on the mean values (e.g., CNA changes, MET 
#' changes, and expression levels). Then, a gene was considered as a subtype-specific one if P-value <= 0.05 (one-way ANOVA test). 
#'
#' @docType package
#'
#' @author Quang-Huy Nguyen
#'
#' @usage SubtypeSpecificGene(omics, cluster, adjustedP = T)
#'
#' @param omics data.frame or matrix. \code{omics} includes its rows are samples and its columns are genomic features.
#'
#' @param cluster numeric. Predictive subgroups correspond with samples after running a clustering tool for your data (e.g., k-means, hclust,...).
#'
#' @param adjustedP logical. Whether we should adjust the gained P-values (ANOVA test) using the Benjamini-Hochberg procedure. Default is \code{adjustedP = T}
#'
#' @return NULL
#'
#' @examples SubtypeSpecificGene(omics = df, cluster = groups)
#'
#' @export

SubtypeSpecificGene = function(omics = NULL, cluster = NULL, adjustedP = T){
  #library
  library(maps)
  library(purrr)
  library(tidyr)
  library(tidyverse)
  library(dplyr)
  
  #Errors
  if(missing(omics)){
    stop("Error: input data is missing. \n")
  }
  
  if(missing(cluster)){
    stop("Error: predictive subgroups correspond with samples are missing. \n")
  }
  
  if(nrow(omics) != length(cluster)){
    stop("Error: Please make sure the samples in omics and clusters are the same and in exactly the same order. \n")
  }
  
  #define the computeQ function may adjust gained P-values following Benjamini-Hochberg FDR
  computeQ <- function(x)
  {
    (x$P.Value*nrow(x))/(x$rank)
  }
  
  #subtype
  group = data.frame(cluster)
  data = omics %>% data.frame()
  data = data[rownames(group),]
  data$sub = group$cluster
  
  #implementation
  Mean <- aggregate(as.matrix(data[,-length(data)]) ~ data[,length(data)], FUN = mean)
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
  
  if(adjustedP == T | TRUE){
    mean_overall = mean_overall[order(mean_overall$P.Value),] #order rows following p-value
    mean_overall$rank = rank(mean_overall$P.Value)
    mean_overall$Q.value = computeQ(mean_overall)
    #create sub matrix with adj.P.Value =< 0.05 and Q.value =< 0.05 
    sig = subset(mean_overall, P.Value <= 0.05 & Q.value <= 0.05)
    sig = sig %>%
      dplyr::select(-rank)
    
    #find specific group of each gene
    sig1 = sig %>%
      dplyr::mutate_all(as.numeric) %>%
      dplyr::mutate(row = dplyr::row_number()) %>%
      tidyr::pivot_longer(-row)
    
    sig2 = sig1[!(sig1$name == "P.Value" | sig1$name == "Q.value"),] #temporarily remove P.Value and Q.Value
    sig2 = sig2 %>% 
      dplyr::group_by(row) %>%
      dplyr::mutate(max_value = max(value),
                    group_number = stringr::str_extract(name,"[:digit:]") %>% as.numeric(),
                    GeneCluster = dplyr::if_else(value == max_value ,group_number,NA_real_)) %>%
      tidyr::fill(GeneCluster,.direction = c("updown")) %>%
      dplyr::select(-group_number,-max_value) %>%
      tidyr::pivot_wider(names_from = name,values_from = value)
    P.value = sig$P.Value
    Q.value = sig$Q.value
    sig2 = data.frame(sig2, P.value, Q.value)
    
    #tidy to a dataframe
    sig2 <- as.data.frame(sig2)
    sig2 = sig2[,-1] #remove 'row' column
    rownames(sig2) = rownames(sig)
    
    #replace 'X' by 'Subgroup'
    names(sig2) = gsub("X","mean in Subgroup ", names(sig2))
    
    #write the results as the csv file
    write.table(sig2,"subtype_specific_gene.txt",sep = "\t", quote = FALSE)
    
    #return the results
    print(head(sig2))

    #warning
    print(writeLines("\nNOTE: \n*the results shown above are incomplete.\n*subtype_specific_gene.csv placed in your current working directory.\n*Please check to identify which gene is specifically assigned to which subgroup."))
  
    
  } else{
    #create sub matrix with P.Value =< 0.05 
    sig = subset(mean_overall, P.Value <= 0.05)
  
  #find specific group of each gene
  sig1 = sig %>%
    dplyr::mutate_all(as.numeric) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    tidyr::pivot_longer(-row)

  sig2 = sig1[!(sig1$name == "P.Value"),] #temporarily remove P.Value
  sig2 = sig2 %>% 
    dplyr::group_by(row) %>%
    dplyr::mutate(max_value = max(value),
                  group_number = stringr::str_extract(name,"[:digit:]") %>% as.numeric(),
                  GeneCluster = dplyr::if_else(value == max_value ,group_number,NA_real_)) %>%
    tidyr::fill(GeneCluster,.direction = c("updown")) %>%
    dplyr::select(-group_number,-max_value) %>%
    tidyr::pivot_wider(names_from = name,values_from = value)
  
  P.value = sig$P.Value
  sig2 = data.frame(sig2,P.value)

  #tidy to a dataframe
  sig2 <- as.data.frame(sig2)
  sig2 = sig2[,-1] #remove 'row' column
  rownames(sig2) = rownames(sig)
  
  #replace 'X' by 'Subgroup'
  names(sig2) = gsub("X","mean in Subgroup ", names(sig2))
  
  #write the results as the csv file
  write.table(sig2,"subtype_specific_gene.txt",sep = "\t", quote = FALSE)
  
  #return the results
  print(head(sig2))
  
  #warning
  print(writeLines("\nNOTE: \n*the results shown above are incomplete.\n*subtype_specific_gene.csv placed in your current working directory.\n*Please check to identify which gene is specifically assigned to which subgroup."))}
}
