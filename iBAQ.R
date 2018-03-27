# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#iBAQ-calculation V1.1
#Graham Lab, USC
#03/26/2018
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#loading required packages
library(ggplot2)
library(gridExtra)
library(plyr)
library(RColorBrewer)
library(dplyr)
library(reshape2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#loading iBaq normalization factors library  
#the iBaq normalization factors are the number of tryptic peptides with length between 6 and 30, inclusive

normFactor=read.table(file = "20170822-Human-uniprot-all-reviewed-iBAQ-counter-only.txt", sep = "\t", header=T) # library contain the number of tryptic peptides with length between 6 and 30, inclusive

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#load PD peptideGroups result
#please strip all special characters #,(,),',$, and ect
#data should be median normalized (done in excel), KNN imputed (done in R),
#For tutorial, please load PD peptideGroup output "test-in.txt"

pepGroupAll <- read.table(file = "test-in.txt", sep = "\t", header=T)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#select your ther starting and ending position of sample columns

colnames(pepGroupAll)
## below are the column names for the example PD peptideGroup output "test-in.txt"
# [1] "Checked"                                           "Confidence"                                       
# [3] "Sequence"                                          "Modifications"                                    
# [5] "Modifications.in.Master.Proteins"                  "Qvality.PEP"                                      
# [7] "Qvality.q.value"                                   "Protein.Groups"                                   
# [9] "Proteins"                                          "PSMs"                                             
# [11] "Master.Protein.Accessions"                         "Positions.in.Master.Proteins"                     
# [13] "Missed.Cleavages"                                  "Theo..MH...Da."                                   
# [15] "Razor.Quan.Results"                                "Abundance..F1..Sample..D1..senescent..A"          
# [17] "Abundance..F3..Sample..D1..senescent..C"           "Abundance..F5..Sample..D1..young..D"              
# [19] "Abundance..F7..Sample..D1..young..H"               "Abundance..F2..Sample..D2..senescent..A"          
# [21] "Abundance..F4..Sample..D2..senescent..C"           "Abundance..F6..Sample..D2..young..D"              
# [23] "Abundance..F8..Sample..D2..young..H"               "Quan.Info"                                        
# [25] "Confidence..by.Search.Engine...Sequest.HT"         "Percolator.q.Value..by.Search.Engine...Sequest.HT"
# [27] "Percolator.PEP..by.Search.Engine...Sequest.HT"     "XCorr..by

### The quatitation values are from column 16 to 23, total of 8 samples for test-in.txt 
sampleStart=16
sampleEnd=23
samples <- as.factor(colnames(pepGroupAll[sampleStart:sampleEnd])) # variable "samples" contains all your sample name



##### Reconstruct data frame
#### seperate rows with ";"
# For peptides that are mapped to multiple proteins (A; B; C), split these rows into individual rows
# one for each protein mapping (A; B; C becomes three rows, one for A, one for B, one for C)

pepUnique<-pepGroupAll[ which(grepl(";",pepGroupAll$Master.Protein.Accessions)==FALSE),]
uniqueness1 <- as.data.frame(pepUnique[,1])
colnames(uniqueness1) <- c("Uniqueness")
uniqueness1[]<-1
pepUnique<-cbind(pepUnique,uniqueness1)


pepNotUnique<-pepGroupAll[ which(grepl(";",pepGroupAll$Master.Protein.Accessions)==TRUE),]

pepNotUniqueNew={}
for(k in 1:length(pepNotUnique)){
  temp<-pepNotUnique[k,]
  for (i in 1:length(strsplit(as.character(temp$Master.Protein.Accessions),";")[[1]])) {
  temp$Master.Protein.Accessions<-strsplit(as.character(pepNotUnique[k,]$Master.Protein.Accessions),";")[[1]][i]
  pepNotUniqueNew <- rbind(pepNotUniqueNew,temp)
}
}

uniqueness2 <- as.data.frame(pepNotUniqueNew[,1])
colnames(uniqueness2) <- c("Uniqueness")
uniqueness2[]<-0
pepNotUniqueNew<-cbind(pepNotUniqueNew,uniqueness2)

Peptides<-rbind(pepUnique,pepNotUniqueNew)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# #iBaq function

pepGroup_to_iBaq <- function(PeptidesDataFrame,sampleStartPos,SampleEndPos,sampleOrder,outputName) {
  ##extract quantitation values
  pepN<-cbind(PeptidesDataFrame$Master.Protein.Accessions,PeptidesDataFrame[,(sampleStartPos:SampleEndPos)],PeptidesDataFrame$Uniqueness)
  colnames(pepN) <- c("Master.Protein.Accessions", as.character(sampleOrder),"Uniqueness")
  
  ##convert to long format for using ddply()
  pepNlong <- melt(pepN,
                   # ID variables - all the variables to keep but not split apart on
                   id.vars=c("Master.Protein.Accessions", "Uniqueness"),
                   # The source columns
                   measure.vars=as.character(sampleOrder),
                   # Name of the destination column that will identify the original
                   # column that the measurement came from
                   variable.name="samples",
                   value.name="Intensity"
  )
  
  ##New: summing over peptides intensities
  pepSum<-ddply(pepNlong, .(Master.Protein.Accessions,samples), summarize,
                sum = round(sum(Intensity), 4),
                uniquePepNum = sum(Uniqueness),
                quantPepNum=n())
  ##convert to wide format
  pepSumWide <- dcast(pepSum, Master.Protein.Accessions+uniquePepNum+quantPepNum ~ samples, value.var="sum")
  ##change columns order
  pepSumWide[ , c("iBaqFactor")] <- NA
  pepSumWide <- pepSumWide[c("Master.Protein.Accessions","uniquePepNum","quantPepNum","iBaqFactor",as.character(sampleOrder))]
  
  
  ###export sum of peptides
  #write.table(pepSumWide, file = "outputName-sumPep.txt",
  #            sep = "\t", quote=F, row.names=F)
  
  iBaqProt={}
  for (i in 1:nrow(normFactor)){
    index <-  which(grepl(paste("^",as.character(normFactor[i,1]),"$",sep = ""),pepSumWide$Master.Protein.Accessions, ignore.case = TRUE))
    index2 <-  which(grepl(paste("^",as.character(normFactor[i,1]),"$",sep = ""),PeptidesDataFrame$Master.Protein.Accessions, ignore.case = TRUE))
    pepSumWide[index,4]<-as.numeric(normFactor[i,2])
    iBaqProt <- rbind(iBaqProt,cbind(pepSumWide[index,1:4],pepSumWide[index,5:(length(sampleOrder)+4)]/ (as.numeric(normFactor[i,2]))))
  }
  
  iBaqProt <- iBaqProt[c("Master.Protein.Accessions",as.character(sampleOrder),"uniquePepNum","quantPepNum","iBaqFactor")]
  
  
  write.table(iBaqProt, file = outputName,
              sep = "\t", quote=F, row.names=F)

}

#Call pepGroup_to_iBaq()
#Peptides: reconstructed peptideGroups
#sampleStart: the starting position of your sample column number
#sampleEnd: the end position of your sample column number
#samples: sample names vector
#"test-out-test2.txt": change output name
pepGroup_to_iBaq(Peptides,sampleStart,sampleEnd,samples,"test-out-test2.txt")


