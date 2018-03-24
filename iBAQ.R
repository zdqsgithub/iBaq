####PD-iBAQ-calculation
###  v1.1 updates: 
###   1. 1keep only the first protein accesion number when isoforms are presented

#loading required packages
library(ggplot2)
library(gridExtra)
library(plyr)
library(RColorBrewer)
library(dplyr)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#loading library contains the number of tryptic peptides with length between 6 and 30, inclusive

#df1=read.table(file = "20170822-Human-uniprot-all-reviewed-iBAQ-number-completed.txt", sep = "\t", header=T) # full library contains prot descritption and other information
df2=read.table(file = "20170822-Human-uniprot-all-reviewed-iBAQ-counter-only.txt", sep = "\t", header=T) # library contain the number of tryptic peptides with length between 6 and 30, inclusive

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#load PD result
#data should be PD peptideGroups median normalized (done in excel), KNN imputed (done in R)

PeptidesAll <- read.table(file = "test-in.txt", sep = "\t", header=T)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#step 1: select your columns
colnames(PeptidesAll)
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

### The quatitation values are from column 16 to 23, total of 8 samples
sampleStart=16
sampleEnd=23
samples <- as.factor(colnames(PeptidesAll[sampleStart:sampleEnd])) ##select your sample quantitation columns



##### Reconstruct data frame

#### seperate rows with ";"

#find position of prot contain ;
#which(grepl(";",Peptides$Master.Protein.Accessions)==TRUE)

pepUnique<-PeptidesAll[ which(grepl(";",Peptides$Master.Protein.Accessions)==FALSE),]
uniqueness1 <- as.data.frame(pepUnique[,1])
colnames(uniqueness1) <- c("Uniqueness")
uniqueness1[]<-1
pepUnique<-cbind(pepUnique,uniqueness1)


pepNotUnique<-PeptidesAll[ which(grepl(";",Peptides$Master.Protein.Accessions)==TRUE),]

pep1<-pepNotUnique[1,]
strsplit(as.character(pep1$Master.Protein.Accessions),";")
#strsplit(as.character(pep1$Master.Protein.Accessions),";")[[1]][1]
#length(strsplit(as.character(pep1$Master.Protein.Accessions),";")[[1]])

df4={}
for(k in 1:length(pepNotUnique)){
  temp<-pepNotUnique[k,]
  for (i in 1:length(strsplit(as.character(temp$Master.Protein.Accessions),";")[[1]])) {
  temp$Master.Protein.Accessions<-strsplit(as.character(pepNotUnique[k,]$Master.Protein.Accessions),";")[[1]][i]
  df4 <- rbind(df4,temp)
}
}

uniqueness2 <- as.data.frame(df4[,1])
colnames(uniqueness2) <- c("Uniqueness")
uniqueness2[]<-0
df4<-cbind(df4,uniqueness2)

Peptides<-rbind(pepUnique,df4)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#iBAQ algorithm
#may need to re-define column number for sample intensities, here there are 10 samples from col 30 to 39

unique.proteins <- unique(Peptides$Master.Protein.Accessions)
sum.Peptides <- data.frame(matrix(data = NA, nrow = length(unique.proteins), ncol = 1))
colnames(sum.Peptides) <- c("Protein")

# which(colnames(Peptides)=="Master.Protein.Accessions")

#summing over peptides intensities
for (i in 1:length(samples)){
  temp.df <- data.frame(matrix(data = NA, nrow = nrow(Peptides), ncol = 2))
  colnames(temp.df) <- c("Protein", "Intensity")
  temp.sample <- as.character(samples[i])
  for (j in 1:nrow(Peptides)){
    temp.df$Protein[j] = as.character(Peptides[j,which(colnames(Peptides)=="Master.Protein.Accessions")])
    temp.df$Intensity[j] = Peptides[j,(i+sampleStart-1)]##scaning over samples intensities columns
  }
  sum.temp.df <- ddply(temp.df, .(Protein), summarize,
                       sum.Peptides = sum(Intensity, na.rm = TRUE))
  colnames(sum.temp.df) <- c("Protein", as.character(samples[i]))
  sum.Peptides <- cbind(sum.Peptides, sum.temp.df)
}

sum.Peptides[,1] <- NULL
Proteins <- as.data.frame(sum.Peptides[,1])
colnames(Proteins) <- c("Master.Protein.Accessions")

sum.Peptides.df <- sum.Peptides[,seq(2, 2*length(samples), 2)]
sum.Peptides.df <- cbind(Proteins, sum.Peptides.df)


#dividing intensities by the sum of possible tryptic pep between 6 and 30 and observed miss cleavages
df3={}
for (i in 1:nrow(df2)){
  
  # example match > which(grepl("P09960",sum.Peptides.df$Master.Protein.Accessions, ignore.case = TRUE))
  index <-  which(grepl(paste("^",as.character(df2[i,1]),"$",sep = ""),sum.Peptides.df$Master.Protein.Accessions, ignore.case = TRUE))
  index2 <-  which(grepl(paste("^",as.character(df2[i,1]),"$",sep = ""),Peptides$Master.Protein.Accessions, ignore.case = TRUE))
  df3 <- rbind(df3,cbind(sum.Peptides.df[index,1],sum.Peptides.df[index,2:(length(samples)+1)] / (as.numeric(df2[i,2])+as.numeric(Peptides[index2,]$Missed.Cleavages))))
}


#output, please change the output file name as needed
write.table(df3, file = "test-out.txt",
            sep = "\t", quote=F, row.names=F)

