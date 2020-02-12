rm(list=ls());cat('\014');if(length(dev.list()>0)){dev.off()}

library(readr)
library(piano)
library(dplyr)
library(omicToolsTest)
library(ggplot2)

CarnivalResults <- 
    readRDS("/home/alvaldeolivas/Documents/GitHub/Saezlab/prostate-phosphoSWATH_V2/Results/CarnivalResults.rds")


sif_carni <- as.data.frame(CarnivalResults$weightedSIF, 
                           stringsAsFactors = FALSE)
    
colnames(sif_carni) <- c("source", "sign", "target", "Weight")
sucesses <- unique(c(gsub("_.*","",sif_carni$source),gsub("_.*","",sif_carni$target)))


CARNI_att_res <- as.data.frame(CarnivalResults$nodesAttributes, 
    stringsAsFactors = FALSE)

bg <- unique(gsub("_.*","",CARNI_att_res$Node))                   

pathways <- gmt_to_csv("/home/alvaldeolivas/Documents/GitHub/Saezlab/prostate-phosphoSWATH_V2/Data/c2.cp.v7.0.symbols.gmt")
sig_pathways <- runGSAhyper(sucesses, universe = bg, gsc = loadGSC(pathways))
sig_pathways_df <- as.data.frame(sig_pathways$resTab)


Kin_activity_PDTs <- as.data.frame(readRDS("/home/alvaldeolivas/Documents/GitHub/Saezlab/prostate-phosphoSWATH_V2/Results/Kin_activity_PDTs.rds")) 

Kin_activity_LNCaP_noInhib_t2_EGF <- Kin_activity_PDTs  %>%
    dplyr::select(LNCaP_noInhib_t2_EGF) 

kinases <- as.data.frame(Kin_activity_LNCaP_noInhib_t2_EGF)

sig_pathways_df$sign <- unlist(lapply(row.names(sig_pathways_df), function(x, kinases, pathways){
  return(mean(kinases[row.names(kinases) %in% pathways[pathways$term == x,1],1]))
},kinases = kinases, pathways = pathways))


sig_pathways_df <- sig_pathways_df[!is.nan(sig_pathways_df$sign),]

PathwaysSelect <- sig_pathways_df[,c(1,2,7)]
PathwaysSign <- PathwaysSelect[PathwaysSelect$`Adjusted p-value` <= 0.00001,]
PathwaysSign$pathway <- row.names(PathwaysSign) 
PathwaysSign <- PathwaysSign[,c(4,1,2,3)]
colnames(PathwaysSign) <- c("pathway","pvalue","AdjPvalu","sign")
PathwaysSign$pathway <- as.factor(PathwaysSign$pathway )


Interesting_pathways <- c("KEGG_REGULATION_OF_ACTIN_CYTOSKELETON",
                          "KEGG_PROSTATE_CANCER",
                          "PID_AR_NONGENOMIC_PATHWAY")

p <- ggplot(PathwaysSign, aes(x = reorder(pathway, pvalue), 
                              y = -log10(pvalue) , fill=sign)) + 
    geom_col() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, 
          # size = ifelse(PathwaysSign$pathway %in% Interesting_pathways, 8, 6), 
          colour = ifelse(levels(reorder(PathwaysSign$pathway, PathwaysSign$pvalue)) %in% Interesting_pathways, "red", "grey40"),
          face = ifelse(levels(reorder(PathwaysSign$pathway, PathwaysSign$pvalue)) %in% Interesting_pathways, "bold", "plain")),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    xlab("")
   #  scale_alpha(range = c(0.4, 0.8))
    
    #  geom_bar(stat = "identity", alpha = PathwaysSign$sign) + 
    # t
    # theme(legend.position = "none",
    #       axis.title.x=element_blank(),
    #      axis.title.y=element_blank(),
    #       axis.text.x=element_blank(),
    #      axis.ticks.x=element_blank(),
    #      plot.margin = NULL) 
    # scale_y_continuous(expand = c(0,0))

