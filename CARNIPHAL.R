library(CARNIVAL)
library(OmnipathR)
library(dplyr)
library(readr)

vignette("CARNIVAL-vignette")

setwd("/home/alvaldeolivas/Documents/GitHub/Saezlab/prostate-phosphoSWATH_V2/")

OmnipathInteractions <- import_Omnipath_Interactions() %>%
    dplyr::filter(is_stimulation != 0 | is_inhibition != 0)  %>%  
    dplyr::mutate(sign = if_else(is_stimulation==1,1,-1)) %>%
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

CarnivalNetwork <- unique(as.data.frame(rbind(OmnipathInteractions,KSN_PDTs)))


counter_CL <- detectCores() - 2
##  We run CARNIVAL for the different populations 
CarnivalResults <-
runCARNIVAL(solverPath="/opt/ibm/ILOG/CPLEX_Studio129/cplex/bin/x86-64_linux/cplex",
            netObj=InputCarnival,
            measObj=t(ViperScores_clus3),
            dir_name="Results_CARNIVAL_clus3",
            weightObj=t(ProgenyResults_clus3),
            nodeID = 'gene',
            timelimit = 7200,
            solver = "cplex",
            progenyMembers = progenyMembers_mice,
            parallelIdx1 = counter_CL)