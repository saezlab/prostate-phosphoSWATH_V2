# Prostate PhosphoSWATH V2

Study the phospho-proteome of two prostate cancer cell lines upon perturbation 
with a combination of different ligands and inhibitors.

## Analysis descibed in this repository

These are the analysis described in this repository:

- **Kinase activity estimation**: First, we evaluate the number of phosphosites
detected in our experiment that are already known in the databases/literature. 
We then estimate the activity of the kinases onthe different experimental 
conditions under study using Viper and different kinase substate networks. 

<https://github.com/saezlab/prostate-phosphoSWATH_V2/blob/master/KinaseActivityAnalyses.md>

- **CarniPHal**: we run CARNIVAL with a particular input setup that allows us
to generate a causal network aiming at explaining the mechanism underlying 
the phosphosite activity in a given condition. (Aurelien's idea) 

<https://github.com/saezlab/prostate-phosphoSWATH_V2/blob/master/RunCARNIPHAL.md>

- **Enrichment of CarniPHal results**: we use the outputs of the two previous
scripts to perform an enrichment analysis in order to detect the most active
pathways in a given condition.

<https://github.com/saezlab/prostate-phosphoSWATH_V2/blob/master/EnrichmentCARNIPHAL_Results.md>

- **Clustering of Phosphosites**: we cluster the most variable phosphosite based
on the correlation between their t-values. We aim at finding cluster of 
phosphosites whose expression follows similar trends along time or for different
ligands or inhibitors. 

<https://github.com/saezlab/prostate-phosphoSWATH_V2/blob/master/ClusterCorrelations.md>