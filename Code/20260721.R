#This is some code adapted from steps 8.1 and 8.2 of the 2024 protocol:
#2024 protocol:https://www.protocols.io/view/brain-data-alchemy-project-meta-analysis-of-re-ana-yxmvmejkng3p/v1?step=8.2
#Megan Hagenauer, July 20, 2026

############################################

#This part is a little different:

#Starting with "GemmaExpressionInfo.csv"
#I made a column in that file called "Keep" and marked off Y or N, then resaved the file as: "GemmaExpressionInfo_Keep.csv"

GemmaExpressionInfo_Keep<-read.csv("GemmaExpressionInfo_Keep.csv", header=TRUE, stringsAsFactors = FALSE)

str(GemmaExpressionInfo_Keep)

ExperimentIDs<-GemmaExpressionInfo_Keep$ExperimentIDs[GemmaExpressionInfo_Keep$Keep=="Y"]


############################
#The code from step 8.2 in the 2024 protocol:

library(gemma.R)

GettingResultSetInfoForDatasets<-function(ExperimentIDs){
  
  #Making an empty data.frame to store results:
  
  ResultSets_toScreen<-data.frame(ExperimentID="NA",ResultSetIDs="NA", ContrastIDs="NA", ExperimentIDs="NA", FactorCategory="NA", ExperimentalFactors="NA", BaselineFactors="NA", Subsetted=FALSE, SubsetBy="NA")
  
  str(ResultSets_toScreen)
  # 'data.frame':	1 obs. of  9 variables:
  # $ ExperimentID       : chr "NA"
  # $ ResultSetIDs       : chr "NA"
  # $ ContrastIDs        : chr "NA"
  # $ ExperimentIDs      : chr "NA"
  # $ FactorCategory     : chr "NA"
  # $ ExperimentalFactors: chr "NA"
  # $ BaselineFactors    : chr "NA"
  # $ Subsetted          : logi FALSE
  # $ SubsetBy           : chr "NA"
  
  #We will then loop over each of the datasets:
  
  for(i in c(1:length(ExperimentIDs))){
    
    #For each dataset, we will use Gemma's API to access the experimental design info:
    Design<-gemma.R::get_dataset_differential_expression_analyses(ExperimentIDs[i])
    
    if(nrow(Design)>0){
      #Next, we'll make some empty vectors to store the experimental factor and baseline factor information for each result id for the dataset:
      ExperimentalFactors<-vector(mode="character", length(Design$result.ID))
      BaselineFactors<-vector(mode="character", length(Design$result.ID))
      
      #We will then loop over each of the result ids for the dataset:
      for(j in c(1:length(Design$result.ID))){
        
        #And grab the vector of experimental factors associated with that result id
        ExperimentalFactorVector<-Design$experimental.factors[[j]]$summary
        #And collapse that info down to a single entry that will fit in our data.frame
        ExperimentalFactors[j]<-paste(ExperimentalFactorVector, collapse="; ")
        
        #And then grab the vector of baseline/control/reference values associated with that result id
        BaselineFactorVector<-Design$baseline.factors[[j]]$summary
        #And collapse that info down to a single entry that will fit in our data.frame
        BaselineFactors[j]<-paste(BaselineFactorVector, collapse="; ")
      }
      
      #Some of the datasets are subsetted for the differential expression analyses
      #We will make an empty vector to store subset information for each result id
      SubsetBy<-vector(mode="character", length(Design$result.ID))
      
      #Then we will determine whether the dataset is subsetted:
      if(Design$isSubset[1]==TRUE){
        
        #If it is subsetted, we will loop over each result id for the dataset
        for (j in c(1:length(Design$result.ID))){
          
          #And grab the vector of subsetting information
          SubsetByVector<-Design$subsetFactor[[j]]$summary
          
          #And then collapse that information down to a single entry that will fit in our dataframe
          SubsetBy[j]<-paste(SubsetByVector, collapse="; ")
        }  
        
        #if the dataset wasn't subsetted for the differential expression analysis:
      }else{
        #We'll just make a vector of NA values to put in the "Subsetted by" column
        SubsetBy<-rep(NA, length((Design$result.ID)))
      }
      
      #Then we combine all of the information for all of the result sets for the dataset into a dataframe
      ResultSets_ForExperiment<-cbind.data.frame(ExperimentID=rep(ExperimentIDs[i],length(Design$result.ID)),ResultSetIDs=Design$result.ID, ContrastIDs=Design$contrast.ID, ExperimentIDs=Design$experiment.ID, FactorCategory=Design$factor.category, ExperimentalFactors, BaselineFactors, Subsetted=Design$isSubset, SubsetBy)
      
      #And add that information as rows to our data frame including the result set information for all datasets:
      ResultSets_toScreen<-rbind.data.frame(ResultSets_toScreen, ResultSets_ForExperiment)
      
      #Then clean up our space before looping to the next dataset:
      rm(ResultSets_ForExperiment, Design, ExperimentalFactors, BaselineFactors, SubsetBy)
      
    }else{
      rm(Design)
    }
    
  }
  
  #When we're done, we'll want to remove the initial (empty) row in our data.frame:
  ResultSets_toScreen<-ResultSets_toScreen[-1,]
  
  #We can make some empty vectors that we can use to store screening notes:
  Include<-vector(mode="character", length=nrow(ResultSets_toScreen))
  WrongBaseline<-vector(mode="character", length=nrow(ResultSets_toScreen))
  ResultsNotRegionSpecific<-vector(mode="character", length=nrow(ResultSets_toScreen))
  ReAnalyze<-vector(mode="character", length=nrow(ResultSets_toScreen))
  
  #And add them as columns to our dataframe:            
  ResultSets_toScreen<-cbind.data.frame(ResultSets_toScreen, Include, WrongBaseline, ResultsNotRegionSpecific, ReAnalyze)
  
  #And then write everything out as a .csv file that we can easily mark up in a spreadsheet program:
  write.csv(ResultSets_toScreen, "ResultSets_toScreen.csv")
  
  print("The Result Sets for your Datasets have been outputted into ResultSets_toScreen.csv")
  print(str(ResultSets_toScreen))
  
  #And clean up our environment:
  rm(Include, WrongBaseline, ResultsNotRegionSpecific, ReAnalyze)
}


#Applying the function:
GettingResultSetInfoForDatasets(ExperimentIDs)


###########
ResultSet_contrasts<-read.csv("ResultSets_Screened.csv", header=TRUE, stringsAsFactors = FALSE )




######
list.files()

ResultSets_toScreen <- read.csv(
  "ResultSets_toScreen_TeamDrug - ResultSets_toScreen.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

str(ResultSets_toScreen)

ResultSet_contrasts <- ResultSets_toScreen[
  ResultSets_toScreen$Include == "Y",
]

str(ResultSet_contrasts)



#################
getwd()
setwd("~/Desktop/Neuropsych Internship/R")
list.files()

library(gemma.R)
library(tidyr)








############
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






