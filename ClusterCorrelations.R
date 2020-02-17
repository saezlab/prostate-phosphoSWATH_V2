rm(list=ls());cat('\014');if(length(dev.list()>0)){dev.off()}

library(OmnipathR)
library(dplyr)
library(readr)
library(igraph)
library(bc3net)

setwd("/home/alvaldeolivas/Documents/GitHub/Saezlab/prostate-phosphoSWATH_V2/")
ResultsLinearModel <- read_tsv("Data/limma_model_results_20190301.tsv")

LinearModelData_df <- ResultsLinearModel %>% 
    dplyr::filter(!is.na(residues_str)) %>% 
    dplyr::mutate(residues_str = strsplit(residues_str, "_")) %>% 
    tidyr::unnest(residues_str) %>% 
    dplyr::mutate(GeneSymbol_Residue = paste(GeneSymbol, residues_str, sep="_")) 

covariate_vars <- c("Intercept","ms_day180409", "ms_day180410", "ms_day180412", 
                    "ms_day180414", "Culture_batch", "fraction_missing", "LNCaP")

## From linear model to Matrix suitable to run Viper.
LinearModelData_df_Clean <-     
    LinearModelData_df %>% 
    dplyr::select(statistic,GeneSymbol_Residue, term, p.value) %>%
    dplyr::filter(!(term %in% covariate_vars)) %>%
    dplyr::filter(GeneSymbol_Residue != "NA_SNA") %>%
    dplyr::filter(GeneSymbol_Residue != "NA_TNA") %>%
    dplyr::group_by(GeneSymbol_Residue, term) %>%
    dplyr::filter(p.value == min(p.value)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>% 
    dplyr::select(-p.value)  

## Most variable genes (Phosposites)
MostVariableSites <- LinearModelData_df_Clean %>% 
    dplyr::group_by(GeneSymbol_Residue) %>%
    dplyr::mutate(variance = var(statistic))  %>%
    dplyr::arrange(desc(variance)) %>%
    dplyr::distinct(GeneSymbol_Residue, .keep_all = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::top_n(500) %>%
    dplyr::select(GeneSymbol_Residue)
    
MatrixStatistic <- 
    dplyr::semi_join(LinearModelData_df_Clean,MostVariableSites) %>%
    tidyr::pivot_wider(names_from = term, values_from = statistic)  %>%
    tibble::column_to_rownames(var = "GeneSymbol_Residue") %>% 
    as.matrix()

corMatrix <- 
    cor(t(MatrixStatistic), use = "pairwise.complete.obs",method = c("pearson"))

corGraph_AdjMatrix <- 
    apply(corMatrix, 2, function(x){ ifelse(x >= 0.7,1,0)})
AnticorGraph_AdjMatrix <- 
    apply(corMatrix, 2, function(x){ ifelse(x <= -0.7,1,0)})    

### We transform from tha adjacency matrix to a graph and we apply clustering
### methods:
corGraph_igraph <- graph_from_adjacency_matrix(corGraph_AdjMatrix, 
    mode = c("undirected"), diag = FALSE)

## We can check the components of the graph:

clusters(corGraph_igraph)

## So, I get the largest connected component: 

corGraph_igraph <- getgcc(corGraph_igraph)


## I can apply recursive Loavien in the clusters larger than 40

## We apply clustering on the networks.
# 2.- Perform of iterative Louvain Clustering.

IterativeLouvain <- function(x,OriginalNetwork,minSize,maxSize) {
  a <- 0 
  n <- length(x)
  newList <- list()
  for (i in seq(n)){
    currentGenes <- x[[i]]
    if (length(currentGenes) >= minSize){
      if(length(currentGenes) <= maxSize){
        a <- a + 1 
        newList[a]  <-x[i]
      } else {
        Subnetwork <-induced_subgraph(OriginalNetwork, currentGenes)
        Modules_rec <- cluster_louvain(Subnetwork, weights = E(Subnetwork)$weight)  
        m <- length(Modules_rec)
        for (j in seq(m)){
          if (length(Modules_rec[[j]]) >= minSize ){
            a <- a + 1 
            newList[[a]] <- Modules_rec[[j]]
          }
        }  
      } 
    }
  }
  return(newList)
}


## Sizes of the modules we want to consider. 
smallest_size <- 5
largest_size <- 40

## First Louvain clustering to obtain the first modules

Modules <- cluster_louvain(corGraph_igraph)

iter <- 0

## Iterative clusters to split the first modules into smaller clusters.

while (max(sapply(Modules[],length)) >= largest_size && iter <= 100) {
  iter <-iter +1 
  print (iter)
  NewModules <- IterativeLouvain(Modules,corGraph_igraph,smallest_size,largest_size)
  Modules <- NewModules  
}

