# This is the code used to run a search in GEO for transcriptional profiling datasets for our cocaine meta-analysis project
# Genevieve Orlewicz
# Search date: June 26, 2026

######### 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

Biocmanager::install("GEOquery")

library(GEOquery)

n


########
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")



####

QueryResults <- searchGEO(MyQueryTerms)

# This will show you an overview of the identified GEO Records:

str(QueryResults)



########
MyQueryTerms <- '(cocain*[All Fields] OR psychostimulant*[All Fields] OR stimulant*[All Fields]) AND ("accumbens"[All Fields] OR hippocamp*[All Fields] OR "dentate gyrus"[All Fields] OR CA1[All Fields] OR CA2[All Fields] OR CA3[All Fields] OR CA4[All Fields] OR "CA field"[All Fields] OR subiculum[All Fields] OR fimbria[All Fields] OR "cornu ammonis"[All Fields]) AND ("Mus musculus"[ORGN] OR "Rattus norvegicus"[ORGN]) AND ("Expression profiling by high throughput sequencing"[DataSet Type] OR "Expression profiling by array"[DataSet Type]) AND gse[Filter]'

QueryResults <- searchGEO(MyQueryTerms)

str(QueryResults)



