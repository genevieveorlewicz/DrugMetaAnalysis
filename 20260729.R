getwd()
list.files()
mydata <- read.csv("~/Desktop/Neuropsych Internship/R/DE_Results_GoodAnnotation.csv")

#Example pipeline for aligning our results across datasets:
#Genevieve
#July 27 2026 --> July 29 2026

############

#Goals:
#Each dataset has differential expression results from a slightly different list of genes
#Depending on the exact tissue dissected, the sensitivity of the transcriptional profiling platform, the representation on the transcriptional profiling platform (for microarray), and the experimental conditions
#The differential expression results from different datasets will also be in a slightly different order
#We want to align these results so that the differential expression results from each dataset are columns, with each row representing a different gene

############
library(dplyr)
library(plyr)

if (!requireNamespace("plyr", quietly = TRUE)) {
  install.packages("plyr")
}

0

############

#Download the necessary functions from our Github repository into your working directory:
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDEResults_2024/Function_AligningDEResults.R

############

#Reading in the functions:

source("Function_AligningDEResults.R")

###########

#Aligning the mouse datasets with each other:

#Example Usage:

ListOfMouseDEResults <-list (DEResults_GSE43410, DEResults_GSE18751, DEResults_GSE63748, DEResults_GSE9134)

str(ListOfMouseDEResults[[3]])

AligningMouseDatasets(ListOfMouseDEResults)

str(Mouse_MetaAnalysis_FoldChanges)

ListOfRatDEResults <-list (DEResults_GSE169665, DEResults_GSE66349)

AligningRatDatasets(ListOfRatDEResults)

str(ListOfRatDEResults[[3]])

str(Rat_MetaAnalysis_FoldChanges)


#################

#Code for aligning the rat and mice results:

#This code isn't nicely functionalized yet
#It also assumes that there are mouse datasets
#It will break if there are only rat datasets - this needs to be fixed.

################

#First: What are gene orthologs?

# Homology refers to biological features including genes and their products that are descended from a feature present in a common ancestor.

# Homologous genes become separated in evolution in two different ways: separation of two populations with the ancestral gene into two species or gene duplication of the ancestral gene within a lineage:

### Genes separated by speciation are called orthologs.
### Genes separated by gene duplication events are called paralogs.

#This definition came from NCBI (https://www.nlm.nih.gov/ncbi/workshops/2023-08_BLAST_evol/ortho_para.html)

#We have the ortholog database that we downloaded from Jackson Labs on April 25, 2024
#This database was trimmed and formatted using the code "FormattingRatMouseOrthologDatabase_20240425.R"

MouseVsRat_NCBI_Entrez<-read.csv("HOM_MouseVsRat_EntrezEnsemblAgree_NoMultimapped_20260511.csv", header=TRUE, stringsAsFactors = FALSE, row.names=1, colClasses=c("character", "character", "character"))

str(MouseVsRat_NCBI_Entrez)

#We want to join this ortholog database to our mouse results (Log2FC and SV):

Mouse_MetaAnalysis_FoldChanges_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_FoldChanges, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_FoldChanges_wOrthologs)
#'data.frame':	28920 obs. of  30 variables:

Mouse_MetaAnalysis_SV_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_SV, by="Mouse_EntrezGene.ID", type="full")

str(Mouse_MetaAnalysis_SV_wOrthologs)
#'data.frame':	28920 obs. of  30 variables:


#*If there are rat datasets*, we then want to join our mouse Log2FC and SV results to the rat results using the ortholog information:
MetaAnalysis_FoldChanges<-join(Mouse_MetaAnalysis_FoldChanges_wOrthologs, Rat_MetaAnalysis_FoldChanges, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_FoldChanges)
#'data.frame':	36792 obs. of  32 variables:

MetaAnalysis_SV<-join(Mouse_MetaAnalysis_SV_wOrthologs, Rat_MetaAnalysis_SV, by="Rat_EntrezGene.ID", type="full")
str(MetaAnalysis_SV)
#'data.frame':	36792 obs. of  32 variables:


#*If there aren't any rat datasets*, we just rename the dataframes so that our downstream code works:
MetaAnalysis_FoldChanges<-Mouse_MetaAnalysis_FoldChanges_wOrthologs
str(MetaAnalysis_FoldChanges)

MetaAnalysis_SV<-Mouse_MetaAnalysis_SV_wOrthologs
str(MetaAnalysis_SV)


###############################

#Not all of these have annotation...
#Those would be the genes that have Entrez IDs that aren't 1:1 with Ensembl
#This matters more this year because we will be joining across datasets with different annotation

sum(is.na(MetaAnalysis_FoldChanges$Mouse_ENSEMBLGene.ID) & is.na(MetaAnalysis_FoldChanges$Rat_Ensembl))
#[1] 14097

MetaAnalysis_FoldChanges<-MetaAnalysis_FoldChanges[(is.na(MetaAnalysis_FoldChanges$Mouse_ENSEMBLGene.ID) & is.na(MetaAnalysis_FoldChanges$Rat_Ensembl))==FALSE,]

dim(MetaAnalysis_FoldChanges)
#[1] 22695    32

MetaAnalysis_SV<-MetaAnalysis_SV[(is.na(MetaAnalysis_SV$Mouse_ENSEMBLGene.ID) & is.na(MetaAnalysis_SV$Rat_Ensembl))==FALSE,]

dim(MetaAnalysis_SV)
#[1] 22695    32

#For simplicity's sake for labeling later charts, let's make a combo annotation with both mouse and rat gene symbol:
MetaAnalysis_FoldChanges$MouseRat_GeneSymbol<-paste(MetaAnalysis_FoldChanges$Mouse_Symbol, MetaAnalysis_FoldChanges$Rat_Symbol, sep="_")

MetaAnalysis_SV$MouseVsRat_GeneSymbol<-paste(MetaAnalysis_SV$Mouse_Symbol, MetaAnalysis_SV$Rat_Symbol, sep="_")


#######################

#We should probably spend some time renaming the comparisons included in our MetaAnalysis_FoldChanges and MetaAnalysis_SV objects now...

colnames(MetaAnalysis_FoldChanges)
# [1] "DB.Class.Key"                                          
# [2] "Mouse_Common.Organism.Name"                            
# [3] "Mouse_NCBI.Taxon.ID"                                   
# [4] "Mouse_Symbol"                                          
# [5] "Mouse_EntrezGene.ID"                                   
# [6] "Mouse_Mouse.MGI.ID"                                    
# [7] "Mouse_HGNC.ID"                                         
# [8] "Mouse_OMIM.Gene.ID"                                    
# [9] "Mouse_Genetic.Location"                                
# [10] "Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38."
# [11] "Mouse_Name"                                            
# [12] "Mouse_Synonyms"                                        
# [13] "Mouse_ENSEMBLGene.ID"                                  
# [14] "Ensembl_Entrez"                                        
# [15] "Rat_Common.Organism.Name"                              
# [16] "Rat_NCBI.Taxon.ID"                                     
# [17] "Rat_Symbol"                                            
# [18] "Rat_EntrezGene.ID"                                     
# [19] "Rat_Mouse.MGI.ID"                                      
# [20] "Rat_HGNC.ID"                                           
# [21] "Rat_OMIM.Gene.ID"                                      
# [22] "Rat_Genetic.Location"                                  
# [23] "Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38."  
# [24] "Rat_Name"                                              
# [25] "Rat_Synonyms"                                          
# [26] "Rat_Ensembl"                                           
# [27] "Rat_DBSymbol"                                          
# [28] "Rat_DBName"                                            
# [29] "GSE111212_exercise"                                    
# [30] "GSE123582_F1..sedentary...F1..sedentary."              
# [31] "GSE270831_social.isolation"                            
# [32] "GSE299436_exercise"                                    
# [33] "MouseRat_GeneSymbol" 


write.csv(colnames(MetaAnalysis_FoldChanges), "Colnames_MetaAnalysisFoldChanges.csv")




###################

#Comparing Log2FC across datasets

#Simple scatterplot... not so promising:

#Example scatter plot comparing two datasets:


#Note - many people prefer to plot these relationships using RRHOs (Rank rank hypergeometric overlap plots)
#I like using both.
#The code for the RRHOs is a little complicated, but I'm happy to share if folks are interested.

#Here's code for looking at the correlation of all of our log2FC results with all of our other log2FC results
#This is called a correlation matrix:

cor(
  as.matrix(MetaAnalysis_FoldChanges[,29:34]),
  use="pairwise.complete.obs",
  method="spearman"
)
#There isn't much similarity across conditions here (outside of comparisons within the same experiment)

#An illustration of the correlation matrix using a hierarchically clustered heatmap, although somewhat pathetic:
heatmap(
  cor(
    as.matrix(MetaAnalysis_FoldChanges[,29:34]),
    use="pairwise.complete.obs",
    method="spearman"
  )
)


str(MetaAnalysis_FoldChanges[,29:34])












################
#Example code showing the pipeline for running a basic meta-analysis of Log2FC and sampling variance values using our previously generated objects MetaAnalysis_FoldChanges & MetaAnalysis_SV
#Megan Hagenauer
#Original version: July 25 2024
#In response to reviewers' comments, this function has been updated to include heterogeneity statistics, publication bias statistics, and robustness statistics
#Updated version: March 10, 2026
#Updated again: July 28, 2026 (to make the code more generalizable for the full 2026 cohort)

######################

#Grabbing input data and setting the working directory:

######################

#Installing and loading relevant code packages:

if (!require("metafor", quietly = TRUE)){
  install.packages("metafor")
}

library(metafor)

if (!require("plyr", quietly = TRUE)){
  install.packages("plyr")
}

library(plyr)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("multtest")                                 

library(multtest)

######################

#Download the two necessary functions from our Github repository:

#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_RunBasicMetaAnalysis.R

#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_FalseDiscoveryCorrection.R

######################

a
#Read in the functions:

source("Function_RunBasicMetaAnalysis.R")

source("Function_FalseDiscoveryCorrection.R")

######################

#Example Usage:

#Figure out which columns contain differential expression output:
colnames(MetaAnalysis_FoldChanges)

#Then make a variable that is a vector of with the column numbers containing the DE output:
#Example code for making a variable representing the columns 
Columns_DE<-c(29:34)

#Then make a variable that is the shorthand annotation that we will use to label the genes in our output, e.g. MouseRat_GeneSymbol:
Column_GeneName<-35

NumberOfComparisons=6
CutOffForNAs=2
#I have 6 comparisons
#2 NA is too many

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV, Columns_DE, Column_GeneName)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...

#Take a peek at the output:

str(metaOutput)
head(metaOutput)
tail(metaOutput)

write.csv(metaOutput, "metaOutput_wHeterogeneityPubBiasRobustMeasures.csv")
write.csv(MetaAnalysis_Annotation, "MetaAnalysis_Annotation_for_metaOutput_wHeterogeneityPubBiasRobustMeasures.csv")

colnames(metaOutput)

write.csv(influence_dfbs, "influence_dfbs.csv")
write.csv(influence_cookd, "influence_cookd.csv" )
write.csv(influence_TF, "influence_TF.csv")

###############

#FDR correction:

FalseDiscoveryCorrection(metaOutput, MetaAnalysis_Annotation)

