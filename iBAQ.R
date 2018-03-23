####PD-iBAQ-calculation
###  v3 updates: 
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

df1=read.table(file = "20170822-Human-uniprot-all-reviewed-iBAQ-number-completed.txt", sep = "\t", header=T) # full library contains prot descritption and other information
df2=read.table(file = "20170822-Human-uniprot-all-reviewed-iBAQ-counter-only.txt", sep = "\t", header=T) # library contain the number of tryptic peptides with length between 6 and 30, inclusive

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#load PD result
#data should be PD peptideGroups median normalized (done in excel), KNN imputed (done in R)

Peptides <- read.table(file = "112717-MCF10A-short-time-gal-wcl-JS-30ACN_PeptideGroups-KNNImputed-ANOVA.txt", sep = "\t", header=T)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#iBAQ algorithm
#may need to re-define column number for sample intensities, here there are 10 samples from col 30 to 39
samples <- as.factor(colnames(Peptides[17:34]))
unique.proteins <- unique(Peptides$Master.Protein.Accessions)
sum.Peptides <- data.frame(matrix(data = NA, nrow = length(unique.proteins), ncol = 1))
colnames(sum.Peptides) <- c("Protein")

#summing over peptides intensities
for (i in 1:length(samples)){
  temp.df <- data.frame(matrix(data = NA, nrow = nrow(Peptides), ncol = 2))
  colnames(temp.df) <- c("Protein", "Intensity")
  temp.sample <- as.character(samples[i])
  for (j in 1:nrow(Peptides)){
    temp.df$Protein[j] = as.character(Peptides[j,11])
    temp.df$Intensity[j] = Peptides[j,i+16]#scaning over samples intensities columns
  }
  sum.temp.df <- ddply(temp.df, .(Protein), summarize,
                       sum.Peptides = sum(Intensity, na.rm = TRUE))
  colnames(sum.temp.df) <- c("Protein", as.character(samples[i]))
  sum.Peptides <- cbind(sum.Peptides, sum.temp.df)
}
sum.Peptides[,1] <- NULL
Proteins <- as.data.frame(sum.Peptides[,1])
colnames(Proteins) <- c("Master.Protein.Accessions")

sum.Peptides.df <- sum.Peptides[,seq(2, 36, 2)]#may need to change column indexes
sum.Peptides.df <- cbind(Proteins, sum.Peptides.df)


#dividing intensities by the sum of possible tryptic pep between 6 and 30 and observed miss cleavages
df3={}
for (i in 1:nrow(df2)){
  
  # example match > which(grepl("P09960",sum.Peptides.df$Master.Protein.Accessions, ignore.case = TRUE))
  index <-  which(grepl(paste("^",as.character(df2[i,1]),"$",sep = ""),sum.Peptides.df$Master.Protein.Accessions, ignore.case = TRUE))
  index2 <-  which(grepl(paste("^",as.character(df2[i,1]),"$",sep = ""),Peptides$Master.Protein.Accessions, ignore.case = TRUE))
  df3 <- rbind(df3,cbind(sum.Peptides.df[index,1],sum.Peptides.df[index,2:19] / (as.numeric(df2[i,2])+as.numeric(Peptides[index2,]$Missed.Cleavages))))
}


#output, please change the output file name as needed
write.table(df3, file = "iBAQ-out-20180119.txt",
            sep = "\t", quote=F, row.names=F)
