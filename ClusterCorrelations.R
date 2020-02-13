rm(list=ls());cat('\014');if(length(dev.list()>0)){dev.off()}

library(OmnipathR)
library(dplyr)
library(readr)


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

corGraph_AdjMatrix <- apply(corMatrix, 2, function(x){ ifelse(x >= 0.5,1,0)})
AnticorGraph_AdjMatrix <- apply(corMatrix, 2, function(x){ ifelse(x <= -0.5,1,0)})    


