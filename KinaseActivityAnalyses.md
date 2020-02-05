Prostate phosphoSWATH: Kinase Activity Analyses
================
Alberto Valdeolivas: <alberto.valdeolivas@bioquant.uni-heidelberg.de>;
Date:
23/01/2020

## Abstract

Study the phospho-proteome of two prostate cancer cell lines upon
perturbation with a combination of different ligands and inhibitors.
This vignette is focused on the application of Kinase Activity
Estimation Methods.

## Source Overlaps

In this section, we compare the overlap between the phophosites measured
in our study and those already available in **Omnipath**
(<http://omnipathdb.org/>),  
**KEA2** (<https://www.maayanlab.net/KEA2/>) and the Putative Downstream
Targets  
**PDTs** article (<https://www.nature.com/articles/s41587-019-0391-9>).

We first load the required libraries:

``` r
library(readr)
library(dplyr)
library(OmnipathR)
library(VennDiagram)
library(viper)
library(tidyr)
library(pheatmap)
```

We read here the dataset containing the results of our linear model
performed in the script. We extract the ID of the unique phosphosites to
compare with the prior knowledge resources.

``` r
ResultsLinearModel <- read_tsv("Data/limma_model_results_20190301.tsv")

LinearModelData_df <- ResultsLinearModel %>% 
    dplyr::filter(!is.na(residues_str)) %>% 
    dplyr::mutate(residues_str = strsplit(residues_str, "_")) %>% 
    tidyr::unnest(residues_str) %>% 
    dplyr::mutate(GeneSymbol_Residue = paste(GeneSymbol, residues_str, sep="_")) 

PhophositesData <- LinearModelData_df %>%
    dplyr::pull(GeneSymbol_Residue) %>% 
    unique()
length(PhophositesData)
## [1] 2235
```

We now use the **OmnipathR** package to retrieve from the Omnipath
webserver post-translational modifications. In particular, we select
phosphorylation and dephosphorylation events and get the associated
phosphosites.

``` r
Omnipath_df <- import_Omnipath_PTMS() %>%
    dplyr::filter(modification %in% c("dephosphorylation","phosphorylation")) %>%
    dplyr::mutate(GeneResidue = paste(substrate_genesymbol,(paste0(residue_type,
        residue_offset)), sep = "_")) 

PhophositesOmnipath <- Omnipath_df %>%
    dplyr::pull(GeneResidue)  %>%
    unique()
length(PhophositesOmnipath)
## [1] 13535
```

In order to increase our coverage, we also read the curated list of
phosphosites available in <https://www.maayanlab.net/KEA2>.

``` r
KEA2_database_path <- "Data/kinase-substrate_phospho-site_level_set_library.tsv"

## This dataset contains a totally different format. We created the following
## function to convert their enzyme-substrates interactions to the same format

processFile = function(filepath) {
    kinases <- character()
    substrates <- character()
    con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    currentLine <- unlist(strsplit(line,"\t"))
    currentKinase <- currentLine[1] 
    Allsubstrates <- currentLine[3:length(currentLine)]
    
    for (j in seq(Allsubstrates)){
        kinases <- c(kinases, currentKinase)
        substrates <-  c(substrates, Allsubstrates[j])
    }
    
  }
  
  close(con)
  df <- data.frame(kinases = kinases, substrates = substrates, 
    stringsAsFactors = FALSE)
  return(df)
}

KEA2_df <- processFile(KEA2_database_path)

PhosphositesKEA2 <- unique(KEA2_df$substrates)
length(PhosphositesKEA2)
## [1] 7857
```

In addition, we included the putative downstream targets (PDTs) reported
in Hijazi et al (2020)
(<https://www.nature.com/articles/s41587-019-0391-9>)

``` r
## We read the datasets from Hijazi et al 2020. Supplementary Dataset 3
PDTs_df <- 
    read_tsv("Data/41587_2019_391_MOESM5_ESM.csv",col_names = FALSE,  skip = 1) %>%
    dplyr::rename(kinases = "X1", substrates = "X2", kinaseFamily = "X3", sign = "X4") %>%
    dplyr::mutate(substrates = gsub("\\(", "_", substrates)) %>%
    dplyr::mutate(substrates = gsub("\\)", "",substrates)) %>%
    dplyr::select(-kinaseFamily) %>%
    dplyr::distinct()
    
PhophositesPDTs <- PDTs_df %>%
    dplyr::pull(substrates)  %>%
    unique()
length(PhophositesPDTs)
## [1] 6195
```

We now check the overlap in the number of phosphosites between the
different sources:

``` r
## Individual Sources
area1 <- length(PhophositesData)
area2 <- length(PhophositesOmnipath)
area3 <- length(PhosphositesKEA2)
area4 <- length(PhophositesPDTs)

## Pairwise Overlap between Sources
n12 <- length(dplyr::intersect(PhophositesData,PhophositesOmnipath))
n13 <- length(dplyr::intersect(PhophositesData,PhosphositesKEA2))
n14 <- length(dplyr::intersect(PhophositesData,PhophositesPDTs))
n23 <- length(dplyr::intersect(PhophositesOmnipath,PhosphositesKEA2))
n24 <- length(dplyr::intersect(PhophositesOmnipath,PhophositesPDTs))
n34 <- length(dplyr::intersect(PhosphositesKEA2,PhophositesPDTs))

## Triple overlaps
n123 <- length(dplyr::intersect(dplyr::intersect(PhophositesData, 
        PhophositesOmnipath), PhosphositesKEA2))
n124 <- length(dplyr::intersect(dplyr::intersect(PhophositesData, 
        PhophositesOmnipath), PhophositesPDTs))
n134 <- length(dplyr::intersect(dplyr::intersect(PhophositesData, 
        PhosphositesKEA2), PhophositesPDTs))
n234 <- length(dplyr::intersect(dplyr::intersect(PhophositesOmnipath, 
        PhosphositesKEA2), PhophositesPDTs))
    
## Total overlap    
n1234 <- length(dplyr::intersect(dplyr::intersect(dplyr::intersect(PhophositesData, 
        PhophositesOmnipath), PhosphositesKEA2),PhophositesPDTs))
    
Venn_plot <- draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24,
    n34, n123, n124, n134, n234, n1234, 
    category = c("Data Phosphosites", "Omnipath Phosphosites",
    "KEA2 Phosphosites", "PDTs Phosphosites"), 
    lty = rep("blank", 4), fill = c("light blue", "red","orange","grey"), 
    alpha = rep(0.25, 4), euler.d = TRUE, scaled=TRUE,
    rotation.degree = 0, reverse=TRUE, cex=1.25,  
    cat.dist = rep(0.075, 4), cat.cex = 1.25)
grid.draw(Venn_plot)
```

![](KinaseActivityAnalyses_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Kinase Activity Estimation

We are now going to run the Kinase Activity estimation (KAE) using viper
based on the kinase substrate networks (KSN) provided by the different
data sources alone and combined. We want to check out wether the
increase in the phosphosite coverage leads to more interpretable
results.

Next, we define here the general function to generate viper regulons
from KSN:

``` r
#'\code{df_to_viper_regulon}
#'
#'This function is designed to generate a ready to use regulon object for viper
#'from a 3 column dataframe representation of a target set collection.
#'
#'@param df a dataframe of n*3 dimension. The first column corresponds the targets,
#'and the second column indicates which regulon does each target belongs to.
#'The third column corresponds to the weight and sign of the interaction between
#'a regulon and its targets.
#'
#'@return a list where each element is a regulon in the viper format. This list
#'is ready to be used as a regulon set in viper.
df_to_viper_regulon <- function(df)
{
  names(df) <- c("feature","pathway","sign")
  df <- df[complete.cases(df),]
  pathway_regulon <- list(0)
  i <- 1
  for(pathway in unique(df$pathway))
  {
    pathway_feature_list <- list(0)
    features <- df[df$pathway == pathway, 3]
    names(features) <- df[df$pathway == pathway, 1]
    pathway_feature_list[[1]] <- features
    pathway_feature_list[[2]] <- rep(1,length(features))
    names(pathway_feature_list) <- c("tfmode","likelihood")
    pathway_regulon[[i]] <- pathway_feature_list
    i <- i+1
  }
  names(pathway_regulon) <- unique(df$pathway)
  return(pathway_regulon)
}
```

We also transform the results of the linear model intro a matrix with
the proper format to run Viper.

``` r
## create sample annotation for using colors with the heatmaps
non_annotated_vars <- c("Intercept", "LNCaP", "ms_day180409", "ms_day180410", 
    "ms_day180412", "ms_day180414", "Culture_batch", "fraction_missing")
sample_annotation <- ResultsLinearModel %>% 
    select(term) %>% 
    distinct()
sample_annotation <- sample_annotation %>% 
    filter(!(term %in% non_annotated_vars)) %>% 
    separate(term, into=c("cell_line", "inhibitor", "time", "ligand"),
             sep="_", remove=FALSE)
sample_annotation <- sample_annotation %>% 
    mutate_if(is.character, as.factor) %>% 
    as.data.frame()
rownames(sample_annotation) <- sample_annotation$term
sample_annotation <- sample_annotation %>% select(-term)

## Specify colors
ann_colors = list(
  cell_line = c(abl = "#E69F00", LNCaP = "#56B4E9"),
  inhibitor = c(noInhib = "#56B4E9", iPI3K = "#E69F00", iMEK = "#009E73"),
  time = c(t0 = "#56B4E9", t1 = "#E69F00", t2 = "#009E73"),
  ligand = c(noLigand = "#56B4E9", EGF = "#E69F00", DHT = "#009E73")
)

covariate_vars <- c("Intercept","ms_day180409", "ms_day180410", "ms_day180412", 
                    "ms_day180414", "Culture_batch", "fraction_missing")

## From linear model to Matrix suitable to run Viper.
MatrixStatistic <-     
    LinearModelData_df %>% 
    dplyr::select(statistic,GeneSymbol_Residue, term) %>%
    dplyr::filter(!(term %in% covariate_vars)) %>%
    dplyr::group_by(GeneSymbol_Residue, term) %>%
    dplyr::filter(statistic == max(abs(statistic))) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()  %>%
    tidyr::pivot_wider(names_from = term, values_from = statistic)  %>%
    tibble::column_to_rownames(var = "GeneSymbol_Residue") %>% 
    as.matrix()
```

### KAE using Omnipath KSN

We generate the KSN only taking into account the information available
in **Omnipath**.

``` r
Omnipath_df_reduced <- Omnipath_df %>%
  dplyr::select(enzyme_genesymbol,GeneResidue) %>%
  dplyr::rename(kinases = "enzyme_genesymbol", substrates = "GeneResidue")

KSN_Omnipath <- Omnipath_df_reduced %>%
  dplyr::distinct() %>%
  dplyr::mutate(sign = 1)
nrow(KSN_Omnipath)
## [1] 32877
```

We generate the regulons and run viper to estimate kinase activity using
the Omnipath KSN and the t-statistic from the limma linear model. We
consider only those kinases regulating at least five phosphosites in our
analysis.

``` r
KSN_Omnipath_regulon <- df_to_viper_regulon(KSN_Omnipath[,c(2,1,3)])

Kin_activity_Omnipath <- viper(MatrixStatistic, regulon = KSN_Omnipath_regulon, 
   minsize = 5, adaptive.size = TRUE, eset.filter = TRUE, verbose = FALSE)
```

We finally display the activity of the kinases for the different
conditions under study:

``` r
pheatmap(Kin_activity_Omnipath, cluster_rows=TRUE, cluster_cols=TRUE,
    annotation_col=sample_annotation, annotation_colors = ann_colors, 
    fontsize=8,  fontsize_row = 8,  treeheight_col = 1.5,  border_color = NA,
    treeheight_row = 1.5)
```

![](KinaseActivityAnalyses_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### KAE using KEA2 KSN

We generate the KSN only taking into account the information available
in **KEA2** (<https://www.maayanlab.net/KEA2/>).

``` r
KSN_KEA2 <- KEA2_df %>%
  dplyr::distinct() %>%
  dplyr::mutate(sign = 1)
nrow(KSN_KEA2)
## [1] 13889
```

We generate the regulons and run viper to estimate kinase activity using
the KEA2 KSN and the t-statistic from the limma linear model. We
consider only those kinases regulating at least five phosphosites in our
analysis.

``` r
KSN_KEA2_regulon <- df_to_viper_regulon(KSN_KEA2[,c(2,1,3)])

Kin_activity_KEA2 <- viper(MatrixStatistic, regulon = KSN_KEA2_regulon, 
   minsize = 5, adaptive.size = TRUE, eset.filter = TRUE, verbose = FALSE)
```

We finally display the activity of the kinases for the different
conditions under study:

``` r
pheatmap(Kin_activity_KEA2, cluster_rows=TRUE, cluster_cols=TRUE,
    annotation_col=sample_annotation, annotation_colors = ann_colors, 
    fontsize=8,  fontsize_row = 8,  treeheight_col = 1.5,  border_color = NA,
    treeheight_row = 1.5)
```

![](KinaseActivityAnalyses_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

### KAE using PDTs KSN

We generate the KSN only taking into account the information reported in
the **PDTs** study (Hijazi et al. 2020).

``` r
KSN_PDTs <- as.data.frame(PDTs_df)
nrow(KSN_PDTs)
## [1] 19326
```

We generate the regulons and run viper to estimate kinase activity using
the PDTs KSN and the t-statistic from the limma linear model. We
consider only those kinases regulating at least five phosphosites in our
analysis.

``` r
KSN_PDTs_regulon <- df_to_viper_regulon(KSN_PDTs[,c(2,1,3)])

Kin_activity_PDTs <- viper(MatrixStatistic, regulon = KSN_PDTs_regulon, 
   minsize = 5, adaptive.size = TRUE, eset.filter = TRUE, verbose = FALSE)
```

We finally display the activity of the kinases for the different
conditions under study:

``` r
pheatmap(Kin_activity_PDTs, cluster_rows=TRUE, cluster_cols=TRUE,
    annotation_col=sample_annotation, annotation_colors = ann_colors, 
    fontsize=8,  fontsize_row = 8,  treeheight_col = 1.5,  border_color = NA,
    treeheight_row = 1.5)
```

![](KinaseActivityAnalyses_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

### KAE merging all the resources to build KSN

We finally merge all the resources (**Omnipath** + **KEA2** + **PDTs**)
to generate a new Kinase regulatory network and we run Viper again.
Note: If finally use this information, we have to take into account
(discuss) that PDTs also include inderect interactions (downstream
estimations).

``` r
KSN_merged_Allsources <- dplyr::bind_rows(KSN_Omnipath, KSN_KEA2, KSN_PDTs) %>%
  dplyr::distinct() %>%
  dplyr::mutate(sign = 1)
nrow(KSN_merged_Allsources)
## [1] 59978
```

We generate the regulons and run viper to estimate kinase activity using
all the sources to build a KSN and the t-statistic from the limma linear
model. We consider only those kinases regulating at least five
phosphosites in our analysis.

``` r
KSN_Allsources_regulon <- df_to_viper_regulon(KSN_merged_Allsources[,c(2,1,3)])

Kin_activity_Allsources <- viper(MatrixStatistic, regulon = KSN_Allsources_regulon, 
   minsize = 5, adaptive.size = TRUE, eset.filter = TRUE, verbose = FALSE)
```

We finally display the activity of the kinases for the different
conditions under study:

``` r
pheatmap(Kin_activity_Allsources, cluster_rows=TRUE, cluster_cols=TRUE,
    annotation_col=sample_annotation, annotation_colors = ann_colors, 
    fontsize=8,  fontsize_row = 8,  treeheight_col = 1.5,  border_color = NA,
    treeheight_row = 1.5)
```

![](KinaseActivityAnalyses_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

## References

D Turei, T Korcsmaros and J Saez-Rodriguez (2016) OmniPath: guidelines
and gateway for literature-curated signaling pathway resources. Nature
Methods 13 (12) PMID: 27898060

Hijazi, M., Smith, R., Rajeeve, V. et al. Reconstructing kinase network
topologies from phosphoproteomics data reveals cancer-associated
rewiring. Nat Biotechnol (2020).
<https://doi.org/10.1038/s41587-019-0391-9>