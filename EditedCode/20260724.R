####### Reading in Gemma
#Example code for importing Gemma's differential expression results for our screened datasets
#Genevieve Orlewicz
#2026-07-22 --> conducted on the 24th

#This is the code if we already narrowed our result sets down to only the ones that we want
ResultSet_contrasts<-read.csv("ResultSets_Screened.csv", header=TRUE, stringsAsFactors = FALSE )


#########
save.image("~/Desktop/Neuropsych Internship/R/Reading in Gemma.RData")
setwd("~/Desktop/Neuropsych Internship/R")

install.packages("tidyverse")

library(tidyverse)

install.packages(dplyr)

library(dplyr)






##########
#Here is code if we want to subset our result sets by the "Include" column:
list.files()
ResultSets_toScreen<-read.csv("ResultSets_toScreen_TeamDrug - ResultSets_toScreen.csv", header=TRUE, stringsAsFactors = FALSE)

str(ResultSets_toScreen)

ResultSet_contrasts<-ResultSets_toScreen[ResultSets_toScreen$Include=="Y", ]

str(ResultSet_contrasts)

#Reading in the function:

DownloadingDEResults<-function(ResultSet_contrasts){
  
  #Some ResultSets have more than one statistical contrast, so they are present more than once in our data frame, e.g.:
  #ResultSet_contrasts$ResultSetIDs
  #[1] 553805 553805 553805 570552 556647
  
  #To pull down the statistical results, we'll only want the unique result set ids:
  UniqueResultSetIDs<-unique(ResultSet_contrasts$ResultSetIDs)
  
  print("These are the result sets that you identified as being of interest")
  print(UniqueResultSetIDs)
  #553805 570552 556647
  
  differentials <- UniqueResultSetIDs %>% lapply(function(x){
    #   # take the first and only element of the output. the function returns a list 
    #   # because single experiments may have multiple resultSets. Here we use the 
    #   # resultSet argument to directly access the results we need
    get_differential_expression_values(resultSet = x)[[1]]
  })
  
  str(differentials)
  #That code worked. Excellent!
  
  # # some datasets might not have all the advertised differential expression results
  # # calculated due to a variety of factors. here we remove the empty differentials
  missing_contrasts <- differentials %>% sapply(nrow) %>% {.==0}
  #[1] FALSE FALSE FALSE
  differentials <<- differentials[!missing_contrasts]
  UniqueResultSetIDs<<-UniqueResultSetIDs[!missing_contrasts]
  
  print("These are the result sets that had differential expression results:")
  print(UniqueResultSetIDs)
  
  print("Your differential expression results for each of your result sets are stored in the object named differentials. This object is structured as a list of data frames. Each element in the list represetns a result set, with the data frame containing the differential expression results")
  
  #Within any particular Result Set, there are likely to be some contrasts that we want and others that we don't want
  #For example, a result set might contain a variety of stress interventions
  #And maybe we only want the acute stress contrast results
  
  #We already identified which statistical contrasts we wanted during our screening:
  #This is the object with the specific contrast ids that we want:
  #ResultSet_contrasts$ContrastIDs
  
  #Which will be these columns within the listed dataframes of differential expression results:
  
  print("These are the columns for the effect sizes for our statistical contrasts of interest (Log(2) Fold Changes")
  Contrasts_Log2FC<<-paste("contrast_", ResultSet_contrasts$ContrastIDs, "_log2fc", sep="")
  
  print(Contrasts_Log2FC)
  #[1] "contrast_151618_log2fc" "contrast_151617_log2fc" "contrast_151619_log2fc"
  #[4] "contrast_186753_log2fc" "contrast_204289_log2fc"
  
  print("these are the columns for the T-statistics for our statistical contrasts of interest - we will use that information to derive the sampling variances")
  
  Contrasts_Tstat<<-paste("contrast_", ResultSet_contrasts$ContrastIDs, "_tstat", sep="")
  
  print(Contrasts_Tstat)
  # [1] "contrast_151618_tstat" "contrast_151617_tstat" "contrast_151619_tstat"
  # [4] "contrast_186753_tstat" "contrast_204289_tstat"  
  
}


library(gemma.R)
library(tidyr)

DownloadingDEResults(ResultSet_contrasts)



#################
SavingGemmaDEResults_forEachResultSet<-function(differentials, UniqueResultSetIDs, ResultSet_contrasts){
  
  for (i in c(1:length(differentials))){
    
    ThisResultSet<-UniqueResultSetIDs[i]
    
    #Pulling out the dataset name from our other data frame
    #For some reason I can find this in the Gemma ResultSet differential expression output
    #Since some datasets have multiple result sets, we just grab the dataset name from the first entry
    ThisDataSet<-ResultSet_contrasts$ExperimentID[ResultSet_contrasts$ResultSetIDs==ThisResultSet][1] 
    
    #Write out a data frame containing the differential expression output for the result set
    #And name it with the dataset id and result set id:
    write.csv(differentials[[i]], paste("DEResults", ThisDataSet, ThisResultSet, ".csv", sep="_"))
    
    rm(ThisDataSet, ThisResultSet)
  }
}


#You can input this function by running the code discussed above to create the function in your R environment. 

#Alternatively, you can download the script for the function from our Github site and save the file in your working directory:
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDEResults_2024/Function_SavingGemmaDEResults_forEachResultSet.R

#And then source it from your working directory:
source("Function_SavingGemmaDEResults_forEachResultSet.R")

#Example usage:

SavingGemmaDEResults_forEachResultSet(differentials, UniqueResultSetIDs, ResultSet_contrasts)


######################

#The fast and dirty version of processing the Gemma DE files
#Megan Hagenauer 
#July 23 2026

#Note: by looping this instead of running it individually for each dataset the output will...
#have columns in the Log2FC,Tstat, SE, and SV output with stupid names (either very long or uninterpretable)
#But we can fix those later...
#The loop may also crash on some datasets, in which case just jump to the next dataset (iteration)

#####################

#Install and load the necessary code packages

#should already be installed from previous steps:

library(gemma.R)

install.packages("tidyverse")

library(tidyverse)

install.packages("dplyr")

library(dplyr)

########################


#Next, download the functions from the 2024 BDA github repository and put them in your working directory
https://github.com/hagenaue/BrainDataAlchemy/tree/main/MetaAnalysis_GemmaDEResults_2024

#source the functions from their files:

source("Function_FilteringDEResults_GoodAnnotation.R")

source("Function_ExtractingDEResultsForContrasts.R")

source("Function_CollapsingDEResults_OneResultPerGene.R")

########################

#And then apply the functions to your differentials object:

YourWorkingDirectory<-getwd()

for(i in c(1:length(differentials)) ){
  
  print(i)
  
  DE_Results<-differentials[[i]]
  FilteringDEResults_GoodAnnotation(DE_Results)
  
  ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
  
  CollapsingDEResults_OneResultPerGene(GSE_ID, DE_Results_GoodAnnotation, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns)
  
  setwd(YourWorkingDirectory)
}

########################

#Save your workspace!  (under R session)
#Save your code! (under file)







#############################
########

#This is code that provides functions for aligning our differential expression results across datasets
#Megan Hagenauer July 25 2024

#Goals:
#Each dataset has differential expression results from a slightly different list of genes
#Depending on the exact tissue dissected, the sensitivity of the transcriptional profiling platform, the representation on the transcriptional profiling platform (for microarray), and the experimental conditions
#The differential expression results from different datasets will also be in a slightly different order
#We want to align these results so that the differential expression results from each dataset are columns, with each row representing a different gene

#Reading in the plyr code package:
library(plyr)

################

#A function for aligning all of our rat differential expression results from different datasets into a single data frame for Log2FCs and sampling variances (SVs):

AligningRatDatasets<-function(ListOfRatDEResults){
  
  #Making an empty list to hold our results:
  Rat_MetaAnalysis_FoldChange_Dfs<-list()
  
  #Looping over all of the rat differential expression results:
  for(i in c(1:length(ListOfRatDEResults))){
    
    #Placing each of the log2FC results for each dataset into a single list
    #Each element in the list is formatted so that the rownames are Rat Entrez Gene ID and then there are columns containing the Log2FC for the differential expression results:
    Rat_MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(Rat_EntrezGene.ID=row.names(ListOfRatDEResults[[i]][[1]]),ListOfRatDEResults[[i]][[1]], stringsAsFactors=FALSE)
  }
  
  #Letting the user know the structure of the set of Log2FC dataframes that we are starting out with:
  print("Rat_MetaAnalysis_FoldChange_Dfs:")
  print(str(Rat_MetaAnalysis_FoldChange_Dfs))
  
  #Running "join all" on the list to align all of the results by Entrez Gene ID and make them a single data frame:
  Rat_MetaAnalysis_FoldChanges<<-join_all(Rat_MetaAnalysis_FoldChange_Dfs, by="Rat_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  #Letting the user know the structure of the dataframe that we just created:
  print("Rat_MetaAnalysis_FoldChanges:")
  print(str(Rat_MetaAnalysis_FoldChanges))
  
  #Doing the same steps for the sampling variances:
  
  Rat_MetaAnalysis_SV_Dfs<-list()
  
  for(i in c(1:length(ListOfRatDEResults))){
    Rat_MetaAnalysis_SV_Dfs[[i]]<-data.frame(Rat_EntrezGene.ID=row.names(ListOfRatDEResults[[i]][[4]]),ListOfRatDEResults[[i]][[4]], stringsAsFactors=FALSE)
  }
  
  print("Rat_MetaAnalysis_SV_Dfs:")
  print(str(Rat_MetaAnalysis_SV_Dfs))
  
  Rat_MetaAnalysis_SV<<-join_all(Rat_MetaAnalysis_SV_Dfs, by="Rat_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("Rat_MetaAnalysis_SV:")
  print(str(Rat_MetaAnalysis_SV))
  
  #Cleaning up our environment to remove unneeded objects:
  rm(Rat_MetaAnalysis_SV_Dfs, Rat_MetaAnalysis_FoldChange_Dfs)
}

#Example Usage;

#ListOfRatDEResults<-list(DEResults_GSE205325)
# 
#AligningRatDatasets(ListOfRatDEResults)
# [1] "Rat_MetaAnalysis_FoldChange_Dfs:"
# List of 1
# $ :'data.frame':	17196 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID    : chr [1:17196] "24153" "24157" "24158" "24159" ...
# ..$ GSE205325_LPS_Chronic: num [1:17196] 0.3485 0.0288 -0.1887 -0.2126 0.1235 ...
# NULL
# [1] "Rat_MetaAnalysis_FoldChanges:"
# 'data.frame':	17196 obs. of  2 variables:
#   $ Rat_EntrezGene.ID    : chr  "24153" "24157" "24158" "24159" ...
# $ GSE205325_LPS_Chronic: num  0.3485 0.0288 -0.1887 -0.2126 0.1235 ...
# NULL
# [1] "Rat_MetaAnalysis_SV_Dfs:"
# List of 1
# $ :'data.frame':	17196 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID    : chr [1:17196] "24153" "24157" "24158" "24159" ...
# ..$ GSE205325_LPS_Chronic: num [1:17196] 0.0745 0.0229 0.0383 0.0257 0.0135 ...
# NULL
# [1] "Rat_MetaAnalysis_SV:"
# 'data.frame':	17196 obs. of  2 variables:
#   $ Rat_EntrezGene.ID    : chr  "24153" "24157" "24158" "24159" ...
# $ GSE205325_LPS_Chronic: num  0.0745 0.0229 0.0383 0.0257 0.0135 ...
# NULL


###########

#A function for aligning all of our mouse differential expression results from different datasets into a single data frame for Log2FCs and sampling variances (SVs):

#This function works the same way as the rat alignment function:

AligningMouseDatasets<-function(ListOfMouseDEResults){
  
  Mouse_MetaAnalysis_FoldChange_Dfs<-list()
  
  for(i in c(1:length(ListOfMouseDEResults))){
    Mouse_MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(Mouse_EntrezGene.ID=row.names(ListOfMouseDEResults[[i]][[1]]),ListOfMouseDEResults[[i]][[1]], stringsAsFactors=FALSE)
  }
  
  print("Mouse_MetaAnalysis_FoldChange_Dfs:")
  print(str(Mouse_MetaAnalysis_FoldChange_Dfs))
  
  Mouse_MetaAnalysis_FoldChanges<<-join_all(Mouse_MetaAnalysis_FoldChange_Dfs, by="Mouse_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("Mouse_MetaAnalysis_FoldChanges:")
  print(str(Mouse_MetaAnalysis_FoldChanges))
  
  Mouse_MetaAnalysis_SV_Dfs<-list()
  
  for(i in c(1:length(ListOfMouseDEResults))){
    Mouse_MetaAnalysis_SV_Dfs[[i]]<-data.frame(Mouse_EntrezGene.ID=row.names(ListOfMouseDEResults[[i]][[4]]),ListOfMouseDEResults[[i]][[4]], stringsAsFactors=FALSE)
  }
  
  print("Mouse_MetaAnalysis_SV_Dfs:")
  print(str(Mouse_MetaAnalysis_SV_Dfs))
  
  Mouse_MetaAnalysis_SV<<-join_all(Mouse_MetaAnalysis_SV_Dfs, by="Mouse_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("Mouse_MetaAnalysis_SV:")
  print(str(Mouse_MetaAnalysis_SV))
  
  rm(Mouse_MetaAnalysis_SV_Dfs, Mouse_MetaAnalysis_FoldChange_Dfs)
}

#Example Usage;

#ListOfMouseDEResults<-list(DEResults_GSE126678, DEResults_GSE181285)

#AligningMouseDatasets(ListOfMouseDEResults)

# [1] "Mouse_MetaAnalysis_FoldChange_Dfs:"
# List of 2
# $ :'data.frame':	21614 obs. of  4 variables:
#   ..$ Mouse_EntrezGene.ID              : chr [1:21614] "11287" "11298" "11302" "11303" ...
# ..$ GSE126678_LPS_Acute              : num [1:21614] 1.9397 0.0805 0.0595 0.0306 0.276 ...
# ..$ GSE126678_LPS_SubchronicPlusAcute: num [1:21614] 1.2967 -0.0472 -0.1459 0.1367 1.5651 ...
# ..$ GSE126678_LPS_Subchronic         : num [1:21614] 0.0582 0.203 -0.1144 0.1361 -0.0051 ...
# $ :'data.frame':	18563 obs. of  2 variables:
#   ..$ Mouse_EntrezGene.ID: chr [1:18563] "100008567" "100009600" "100012" "100017" ...
# ..$ GSE181285_LPS_Acute: num [1:18563] 0.0198 0.0225 0.0641 -0.0049 -0.0588 ...
# NULL
# [1] "Mouse_MetaAnalysis_FoldChanges:"
# 'data.frame':	24287 obs. of  5 variables:
#   $ Mouse_EntrezGene.ID              : chr  "11287" "11298" "11302" "11303" ...
# $ GSE126678_LPS_Acute              : num  1.9397 0.0805 0.0595 0.0306 0.276 ...
# $ GSE126678_LPS_SubchronicPlusAcute: num  1.2967 -0.0472 -0.1459 0.1367 1.5651 ...
# $ GSE126678_LPS_Subchronic         : num  0.0582 0.203 -0.1144 0.1361 -0.0051 ...
# $ GSE181285_LPS_Acute              : num  -0.042 -0.0368 -0.0534 0.1067 -0.5258 ...
# NULL
# [1] "Mouse_MetaAnalysis_SV_Dfs:"
# List of 2
# $ :'data.frame':	21614 obs. of  4 variables:
#   ..$ Mouse_EntrezGene.ID              : chr [1:21614] "11287" "11298" "11302" "11303" ...
# ..$ GSE126678_LPS_Acute              : num [1:21614] 0.62127 0.14737 0.00437 0.01624 1.34369 ...
# ..$ GSE126678_LPS_SubchronicPlusAcute: num [1:21614] 0.66434 0.14559 0.00438 0.01567 0.98359 ...
# ..$ GSE126678_LPS_Subchronic         : num [1:21614] 0.84004 0.14211 0.00439 0.01592 1.40671 ...
# $ :'data.frame':	18563 obs. of  2 variables:
#   ..$ Mouse_EntrezGene.ID: chr [1:18563] "100008567" "100009600" "100012" "100017" ...
# ..$ GSE181285_LPS_Acute: num [1:18563] 0.11456 0.01612 0.00329 0.00487 0.00719 ...
# NULL
# [1] "Mouse_MetaAnalysis_SV:"
# 'data.frame':	24287 obs. of  5 variables:
#   $ Mouse_EntrezGene.ID              : chr  "11287" "11298" "11302" "11303" ...
# $ GSE126678_LPS_Acute              : num  0.62127 0.14737 0.00437 0.01624 1.34369 ...
# $ GSE126678_LPS_SubchronicPlusAcute: num  0.66434 0.14559 0.00438 0.01567 0.98359 ...
# $ GSE126678_LPS_Subchronic         : num  0.84004 0.14211 0.00439 0.01592 1.40671 ...
# $ GSE181285_LPS_Acute              : num  0.00391 0.02738 0.00601 0.0101 0.03332 ...
# NULL









#####################
#Example pipeline for aligning our results across datasets:
#Megan Hagenauer
#July 27 2026

############

#Goals:
#Each dataset has differential expression results from a slightly different list of genes
#Depending on the exact tissue dissected, the sensitivity of the transcriptional profiling platform, the representation on the transcriptional profiling platform (for microarray), and the experimental conditions
#The differential expression results from different datasets will also be in a slightly different order
#We want to align these results so that the differential expression results from each dataset are columns, with each row representing a different gene

############

if (!requireNamespace("plyr", quietly = TRUE)) {
  install.packages("plyr")
}


getwd()

############

#Download the necessary functions from our Github repository into your working directory:
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDEResults_2024/Function_AligningDEResults.R

############

#Reading in the functions:

source("Function_AligningDEResults.R")

###########

#Aligning the mouse datasets with each other:

#Example Usage:

ListOfMouseDEResults<-list(DEResults_GSE111212_515554, DEResults_GSE123582_519225)

AligningMouseDatasets(ListOfMouseDEResults)

# [1] "Mouse_MetaAnalysis_FoldChange_Dfs:"
# List of 2
# $ :'data.frame':	21614 obs. of  4 variables:
#   ..$ Mouse_EntrezGene.ID              : chr [1:21614] "11287" "11298" "11302" "11303" ...
# ..$ GSE126678_LPS_Acute              : num [1:21614] 1.9397 0.0805 0.0595 0.0306 0.276 ...
# ..$ GSE126678_LPS_SubchronicPlusAcute: num [1:21614] 1.2967 -0.0472 -0.1459 0.1367 1.5651 ...
# ..$ GSE126678_LPS_Subchronic         : num [1:21614] 0.0582 0.203 -0.1144 0.1361 -0.0051 ...
# $ :'data.frame':	18563 obs. of  2 variables:
#   ..$ Mouse_EntrezGene.ID: chr [1:18563] "100008567" "100009600" "100012" "100017" ...
# ..$ GSE181285_LPS_Acute: num [1:18563] 0.0198 0.0225 0.0641 -0.0049 -0.0588 ...
# NULL
# [1] "Mouse_MetaAnalysis_FoldChanges:"
# 'data.frame':	24287 obs. of  5 variables:
#   $ Mouse_EntrezGene.ID              : chr  "11287" "11298" "11302" "11303" ...
# $ GSE126678_LPS_Acute              : num  1.9397 0.0805 0.0595 0.0306 0.276 ...
# $ GSE126678_LPS_SubchronicPlusAcute: num  1.2967 -0.0472 -0.1459 0.1367 1.5651 ...
# $ GSE126678_LPS_Subchronic         : num  0.0582 0.203 -0.1144 0.1361 -0.0051 ...
# $ GSE181285_LPS_Acute              : num  -0.042 -0.0368 -0.0534 0.1067 -0.5258 ...
# NULL
# [1] "Mouse_MetaAnalysis_SV_Dfs:"
# List of 2
# $ :'data.frame':	21614 obs. of  4 variables:
#   ..$ Mouse_EntrezGene.ID              : chr [1:21614] "11287" "11298" "11302" "11303" ...
# ..$ GSE126678_LPS_Acute              : num [1:21614] 0.62127 0.14737 0.00437 0.01624 1.34369 ...
# ..$ GSE126678_LPS_SubchronicPlusAcute: num [1:21614] 0.66434 0.14559 0.00438 0.01567 0.98359 ...
# ..$ GSE126678_LPS_Subchronic         : num [1:21614] 0.84004 0.14211 0.00439 0.01592 1.40671 ...
# $ :'data.frame':	18563 obs. of  2 variables:
#   ..$ Mouse_EntrezGene.ID: chr [1:18563] "100008567" "100009600" "100012" "100017" ...
# ..$ GSE181285_LPS_Acute: num [1:18563] 0.11456 0.01612 0.00329 0.00487 0.00719 ...
# NULL
# [1] "Mouse_MetaAnalysis_SV:"
# 'data.frame':	24287 obs. of  5 variables:
#   $ Mouse_EntrezGene.ID              : chr  "11287" "11298" "11302" "11303" ...
# $ GSE126678_LPS_Acute              : num  0.62127 0.14737 0.00437 0.01624 1.34369 ...
# $ GSE126678_LPS_SubchronicPlusAcute: num  0.66434 0.14559 0.00438 0.01567 0.98359 ...
# $ GSE126678_LPS_Subchronic         : num  0.84004 0.14211 0.00439 0.01592 1.40671 ...
# $ GSE181285_LPS_Acute              : num  0.00391 0.02738 0.00601 0.0101 0.03332 ...
# NULL

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

###################

#Comparing Log2FC across datasets

#Simple scatterplot... not so promising:

#Example scatter plot comparing two datasets:
plot(MetaAnalysis_FoldChanges$GSE299436_exercise~MetaAnalysis_FoldChanges$GSE111212_exercise)

#Note - many people prefer to plot these relationships using RRHOs (Rank rank hypergeometric overlap plots)
#I like using both.
#The code for the RRHOs is a little complicated, but I'm happy to share if folks are interested.

#Here's code for looking at the correlation of all of our log2FC results with all of our other log2FC results
#This is called a correlation matrix:

cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:28,33)]), use="pairwise.complete.obs", method="spearman")
#There isn't much similarity across conditions here (outside of comparisons within the same experiment)

#An illustration of the correlation matrix using a hierarchically clustered heatmap, although somewhat pathetic:
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:28,33)]), use="pairwise.complete.obs", method="spearman"))

