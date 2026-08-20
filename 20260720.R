#################
#' gemma.R package: Access curated gene expression data and differential expression analyses
#'
#' This package contains wrappers and convenience functions for Gemma's RESTful API
#' that enables access to curated expression and differential expression data
#' from over 15,000 published studies (as of mid-2022). Gemma (https://gemma.msl.ubc.ca) is a web site, database and a set
#' of tools for the meta-analysis, re-use and sharing of transcriptomics data,
#' currently primarily targeted at the analysis of gene expression profiles.
#'
#' Most users will want to start with the high-level functions like \code{\link{get_dataset_object}}, \code{\link{get_differential_expression_values}} and \code{\link{get_platform_annotations}}
#' Additional lower-level methods are available that directly map to the Gemma RESTful API methods.
#'
#' For more information and detailed usage instructions check the
#' \href{https://pavlidislab.github.io/gemma.R/index.html}{README}, the
#' \href{https://pavlidislab.github.io/gemma.R/reference/index.html}{function reference}
#' and the \href{https://pavlidislab.github.io/gemma.R/articles/gemma.R.html}{vignette}.
#'
#' All software-related questions should be posted to the Bioconductor Support Site:
#' \url{https://support.bioconductor.org}
#'
#' @references
#'
#' Lim, N. et al., Curation of over 10 000 transcriptomic studies to enable
#' data reuse, Database, 2021. \url{https://doi.org/10.1093/database/baab006}
#'
#' @author Javier Castillo-Arnemann, Jordan Sicherman, Ogan Mancarci, Guillaume Poirier-Morency
#'
#' @name gemma.R
#'
#' @import data.table
#' @import bit64
#' @import digest
#' @importFrom magrittr %>%
#' @importFrom rlang .data
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") utils::globalVariables(c("."))


utils::globalVariables(c("platform.ID", "analysis.ID", "analysis.Threshold", "baseline.category", 
                         "baseline.categoryURI", "baseline.factorValue", "baseline.factorValueURI", 
                         "cf.Val", "cf.ValLongUri", "experiment.ID", "experiment.shortName", "factor.ID", "factorId", 
                         "factorValue", "genes.Analyzed", "id", "probes.Analyzed", "result.ID", 
                         "resultIds", "stats.DE", "stats.Down", "stats.Up", "subsetFactor.category", 
                         "subsetFactor.categoryURI", "subsetFactor.Enabled", "subsetFactor.factorValue", 
                         "subsetFactor.factorValueURI", "valueUri", "category", 
                         "categoryURI", "experimental.factorValue","value",
                         "contrast.ID","%$%",'baseline.factors',"ID","sample.ID","filter_genes","print.listable_pheatmap",
                         "corrected_pvalue"))






################# 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gemma.R")

BiocManager::install("https://github.com/PavlidisLab/gemma.R/blob/devel/R/gemma.R.R
")





################# DOWNLOADING THE ASSORTMENT
install.packages("plyr")

library(plyr)

install.packages("tidyverse")

library(tidyverse)





############## DOWNLOADING GEMMA PT. 1
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gemma.R")

library(gemma.R)



############### SETTING WORKING DIRECTORY
setwd()




############### DOWNLOADING GEMMA PT. 2
#We can download the processed expression data for any particular dataset using this code:
Expression<-gemma.R::get_dataset_processed_expression("GSE81672")
#In this case "GSE81672" is the dataset id
#For reference, this example is an RNA-Seq dataset

#Structure of the new Expression object:
str(Expression)
#Classes ‘data.table’ and 'data.frame':	35205 obs. of  103 variables:
#The first four columns are row metadata: Probe, GeneSymbol, GeneName, NCBIid
#The rest of the columns are gene expression values for each subject



####### GSE HIST, MIN, MAX, ETC
#If you want to visualize the distribution of gene expression data for the entire study, you will need to grab all of the columns that aren't row metadata (i.e., not rows 1-4) and force them into the format of a numeric matrix first: 

ExpressionMatrix<-as.matrix(Expression[,-c(1:4)])

#You can visualize the distribution of the gene expression data for the dataset using a histogram
#We can make this pretty by adding a title and x-axis label
#You can also change the color and scaling
hist(ExpressionMatrix, main="Histogram", xlab="Log2 Expression", col="purple", cex.axis=1.3, cex.lab=1.3)

#The range of the x-axis is the range of the log2 gene expression values
#The large spike on the left side of the histogram ("floor effect") are all of the genes that aren't truly expressed or have too low of expression to be measurable

#You can save the histogram using "export" in the Plots window

#We can also pull out numeric values summarizing the distribution, e.g.:
min(ExpressionMatrix, na.rm=TRUE)
#[1] -5.8601
median(ExpressionMatrix, na.rm=TRUE)
#[1] -2.1651
max(ExpressionMatrix, na.rm=TRUE)
#[1] 12.312

#Or to get more of an overview:
summary(ExpressionMatrix)




########## LOADING IN MY GSEs
ExperimentIDs <- c("GSE88736", "GSE20885", "GSE47457", "GSE141520", "GSE134935", "GSE111600", "GSE110344", "GSE261663")
ExperimentIDs



######### W/ GOOGLE SHEETS CSV, OVERWRITES THE EXPERIMENTIDS PART
"GSE_Accessions_QueryResults_TeamDrug - Sheet1.csv"

list.files()

ExperimentIDs<-read.csv("GSE_Accessions_QueryResults_TeamDrug - Sheet1.csv", header=TRUE, stringsAsFactors = FALSE)
str(ExperimentIDs)

ExperimentIDs<-ExperimentIDs[,1]
str(ExperimentIDs)

GemmaExpressionInfo <- data.frame(
  ExperimentIDs = ExperimentIDs,
  InGemma = character(length(ExperimentIDs)),
  MinExpression = numeric(length(ExperimentIDs)),
  MedianExpression = numeric(length(ExperimentIDs)),
  MaxExpression = numeric(length(ExperimentIDs))
)
str(GemmaExpressionInfo)




########## BASIC STRUCTURE FOR A LOOP + Test whether GEO Accession ID is in Gemma w/ the if function
for(i in c(1:length(ExperimentIDs))){

  print(ExperimentIDs[i])
  
  if(inherits(try(get_dataset_processed_expression(ExperimentIDs[i]), silent=TRUE), "try-error")){
    
    print("Error: Not in Gemma")
    GemmaExpressionInfo$InGemma[i]<-"N"
    
  }else{
    
  GemmaExpressionInfo$InGemma[i]<-"Y"

  Expression<-gemma.R::get_dataset_processed_expression(ExperimentIDs[i])
  
  ExpressionMatrix<-as.matrix(Expression[,-c(1:4)])
  
  pdf(paste(ExperimentIDs[i], "_Histogram.pdf", sep=""), height=4, width=4)
  
  hist(ExpressionMatrix, main="Histogram", xlab="Log2 Expression", col="blue", cex.axis=1.3, cex.lab=1.3)
  
  dev.off()
  
  print(min(ExpressionMatrix, na.rm=TRUE))
  
  GemmaExpressionInfo$MinExpression[i]<-min(ExpressionMatrix, na.rm=TRUE)
  
  print(median(ExpressionMatrix, na.rm=TRUE))
  
  GemmaExpressionInfo$MedianExpression[i]<-median(ExpressionMatrix, na.rm=TRUE)
  
  print(max(ExpressionMatrix, na.rm=TRUE))
  
  GemmaExpressionInfo$MaxExpression[i]<-max(ExpressionMatrix, na.rm=TRUE)

  rm(Expression, ExpressionMatrix)
  
  }
  
}




######
str(GemmaExpressionInfo)

write.csv(GemmaExpressionInfo, "GemmaExpressionInfo.csv")



