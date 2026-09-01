1 + 1
library("gemma.R", character.only = TRUE)
DownloadingDEResults <- function(ResultSet_contrasts){
  print("hello")
}

#Reading in Gemma_Updated
#From Bansi Bhimjiani 
#2026-07-27


#This is the code if we already narrowed our result sets down to only the ones that we want
ResultSet_contrasts<-read.csv("ResultSets_toScreen_TeamDrug - ResultSets_toScreen.csv", header=TRUE, stringsAsFactors = FALSE)


#Here is code if we want to subset our result sets by the "Include" column:
list.files()
ResultSets_toScreen<-read.csv("ResultSets_toScreen_TeamDrug - ResultSets_toScreen.csv", header=TRUE, stringsAsFactors = FALSE)

str(ResultSets_toScreen)

ResultSet_contrasts<-ResultSets_toScreen[ResultSets_toScreen$Include=="Y", ]

str(ResultSet_contrasts)

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!requireNamespace("gemma.R", quietly = TRUE)) {
  BiocManager::install("gemma.R")
}

library("gemma.R", character.only = TRUE)

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
  
  names(differentials)<-UniqueResultSetIDs
  
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

DownloadingDEResults(ResultSet_contrasts)

search()
library(dplyr)

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

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(tidyverse)

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

library(dplyr)

########################

#Next, download the functions from our repository and put them in your working directory:

#Function FilteringDEResults_GoodAnnotation:
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDEResults_2024/Function_FilteringDEResults_GoodAnnotation.R

#Function ExtractingDEResultsForContrasts (this one was updated/debugged for 2026):
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_ExtractingDEResultsForContrasts.R

#Function_CollapsingDEResults_OneResultPerGene:
#https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GEOData_2026/Function_CollapsingDEResults_OneResultPerGene.R

##########################

#Then source the functions from their files in your working directory:

source("Function_FilteringDEResults_GoodAnnotation.R")

source("Function_ExtractingDEResultsForContrasts.R")

source("Function_CollapsingDEResults_OneResultPerGene.R")

########################

#And then apply the functions to your differentials object:

YourWorkingDirectory<-getwd()

for(i in c(1:length(differentials)) ){
  
  print(i)
  
  DE_Results<-differentials[[i]]
  CurrentResultSet<-names(differentials)[i]
  
  FilteringDEResults_GoodAnnotation(DE_Results)
  
  ExtractingDEResultsForContrasts(DE_Results_GoodAnnotation, Contrasts_Log2FC, Contrasts_Tstat, ResultSet_contrasts)
  
  CollapsingDEResults_OneResultPerGene(
    GSE_ID,
    DE_Results_GoodAnnotation,
    ComparisonsOfInterest,
    NamesOfFoldChangeColumns,
    NamesOfTstatColumns
  )
  
}


length(differentials)
str(differentials)


args(ExtractingDEResultsForContrasts)





#####
search()
DownloadingDEResults(ResultSet_contrasts)
str(differentials)
ls()
exists("differentials")

########################

#Save your workspace!  (under R session)
#Save your code! (under file)

library(dplyr)
library("gemma.R", character.only = TRUE)
ResultSets_toScreen <- read.csv(
  "ResultSets_toScreen_TeamDrug - ResultSets_toScreen.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

ResultSet_contrasts <- ResultSets_toScreen[ResultSets_toScreen$Include=="Y", ]
DownloadingDEResults(ResultSet_contrasts)



differentials <- UniqueResultSetIDs %>% lapply(function(x){
  print(paste("Downloading ResultSet:", x))
  tryCatch(
    get_differential_expression_values(resultSet = x)[[1]],
    error = function(e) {
      print(paste("FAILED:", x))
      print(e$message)
      return(NULL)
    }
  )
})
