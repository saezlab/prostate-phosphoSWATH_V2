rm(list=ls());cat('\014');if(length(dev.list()>0)){dev.off()}

library(CARNIVAL)
library(OmnipathR)
library(dplyr)
library(readr)
library(caret)
library(bc3net)

range <- function(x){ (x - min(x))/(max(x)-min(x)) * (1 - (-1)) + -1 }

setwd("/home/alvaldeolivas/Documents/GitHub/Saezlab/prostate-phosphoSWATH_V2/")

#### Generate the Network

OmnipathInteractions <- import_Omnipath_Interactions() %>%
    dplyr::filter(consensus_stimulation != 0 | consensus_inhibition != 0)  %>%  
    dplyr::mutate(sign = if_else(consensus_stimulation==1,1,-1)) %>%
    dplyr::select(source_genesymbol, sign,  target_genesymbol) %>%
    dplyr::rename(source ="source_genesymbol", target ="target_genesymbol") 
    
PDTs_df <- 
    read_tsv("Data/41587_2019_391_MOESM5_ESM.csv",col_names = FALSE,  skip = 1) %>%
    dplyr::rename(kinases = "X1", substrates = "X2", kinaseFamily = "X3", sign = "X4") %>%
    dplyr::mutate(substrates = gsub("\\(", "_", substrates)) %>%
    dplyr::mutate(substrates = gsub("\\)", "",substrates)) %>%
    dplyr::select(-kinaseFamily) %>%
    dplyr::distinct()  

KSN_PDTs <- as.data.frame(PDTs_df) %>%
    dplyr::rename(source="kinases", target="substrates") %>%
    dplyr::select(source, sign, target)

CarnivalNetwork <- dplyr::bind_rows(OmnipathInteractions, KSN_PDTs) %>%
    dplyr::distinct()

CarnivalNetwork$target <- gsub("[/]","_",CarnivalNetwork$target)
CarnivalNetwork$target <- gsub("[space]","_",CarnivalNetwork$target)
CarnivalNetwork$source <- gsub("[-]", "_", CarnivalNetwork$source)
CarnivalNetwork$target <- gsub("[-]", "_", CarnivalNetwork$target)
CarnivalNetwork$source <- gsub(pattern = "/", replacement = "_", x = CarnivalNetwork$source, fixed = TRUE)
CarnivalNetwork$source <- gsub(pattern = " ", replacement = "_", x = CarnivalNetwork$source, fixed = TRUE)
CarnivalNetwork$target <- gsub(pattern = "/", replacement = "_", x = CarnivalNetwork$target, fixed = TRUE)
CarnivalNetwork$target <- gsub(pattern = " ", replacement = "_", x = CarnivalNetwork$target, fixed = TRUE)

### We are going to keep only the largest connected component of our network

CarnivalNetwork_igraph <- 
    graph_from_data_frame(CarnivalNetwork[,c(1,3)], directed = TRUE) %>% 
    getgcc() %>%
    igraph::as_data_frame()  

CarnivalNetwork_gcc <- dplyr::semi_join(CarnivalNetwork,CarnivalNetwork_igraph,
    by = c("source" = "from", "target" = "to"))

# CarnivalNetwork$source <- gsub("[-+{},;() ]","___",CarnivalNetwork$source)
# CarnivalNetwork$target <- gsub("[-+{},;() ]","___",CarnivalNetwork$target)

###### Linear model to take the t-values and use them as the dowstream targets
## for CARNIVAL, the measurements objects

ResultsLinearModel <- read_tsv("Data/limma_model_results_20190301.tsv")

LinearModelData_df <- ResultsLinearModel %>% 
    dplyr::filter(!is.na(residues_str)) %>% 
    dplyr::mutate(residues_str = strsplit(residues_str, "_")) %>% 
    tidyr::unnest(residues_str) %>% 
    dplyr::mutate(GeneSymbol_Residue = paste(GeneSymbol, residues_str, sep="_")) 

LinearModelData_LNCaP_noInhib_t2_EGF <- LinearModelData_df %>%
    dplyr::filter(term == "LNCaP_noInhib_t2_EGF") %>%
    dplyr::filter(p.value < 0.05)  %>%
    dplyr::select(GeneSymbol_Residue, statistic, p.value)  %>% 
    dplyr::group_by(GeneSymbol_Residue) %>%
    dplyr::filter(p.value == min(p.value))  %>%
    dplyr::ungroup() %>%
    dplyr::select(-p.value)  %>%
    dplyr::filter(GeneSymbol_Residue %in% CarnivalNetwork_gcc$target) %>%
    tibble::column_to_rownames(var = "GeneSymbol_Residue") %>%
    as.data.frame() 
    
# rownames(LinearModelData_LNCaP_noInhib_t2_EGF) <- 
#    gsub("[-+{},;() ]","___",rownames(LinearModelData_LNCaP_noInhib_t2_EGF))

###### NES Scores. I am going to use these scores as Progeny scores to drive the network

Kin_activity_PDTs <- as.data.frame(readRDS("Results/Kin_activity_PDTs.rds")) 

Kin_activity_LNCaP_noInhib_t2_EGF <- Kin_activity_PDTs  %>%
    dplyr::select(LNCaP_noInhib_t2_EGF) 

## We have to scale the NES between 1 and 0. 
Kin_activity_LNCaP_noInhib_t2_EGF$LNCaP_noInhib_t2_EGF <- 
    as.data.frame(range(Kin_activity_LNCaP_noInhib_t2_EGF$LNCaP_noInhib_t2_EGF))

# rownames(Kin_activity_LNCaP_noInhib_t2_EGF) <- 
#     gsub("[-+{},;() ]","___",rownames(Kin_activity_LNCaP_noInhib_t2_EGF))

##### And I have a perturbation EGF
inputObj <- data.frame(EGF = 1)

##  We run CARNIVAL for our particular condition
# counter_CL <- detectCores() - 2


CarnivalResults <-runCARNIVAL(
    solverPath="/opt/ibm/ILOG/CPLEX_Studio129/cplex/bin/x86-64_linux/cplex",
    netObj=CarnivalNetwork_gcc,
    measObj=as.data.frame(t(LinearModelData_LNCaP_noInhib_t2_EGF)),
    inputObj = inputObj,
    DOTfig=TRUE, 
    dir_name="Results",
    weightObj=as.data.frame(t(Kin_activity_LNCaP_noInhib_t2_EGF)),
    # nodeID = 'gene',
    timelimit = 7200,
    solver = "cplex")
            # progenyMembers = progenyMembers_mice,
            # parallelIdx1 = counter_CL)
saveRDS(CarnivalResults, file = "Results/CarnivalResults.rds")