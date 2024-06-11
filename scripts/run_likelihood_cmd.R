#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

print(paste(args[1], "is the counts file you would like to use."))
print(paste(args[2], "is the PATH to where you would like your output.File will be named by script for downstream input"))

#I will use this to get the observed ML estimate of s across my windows. This time run all windows in parallel.


##### Run Likelihood on counts data
## Viability
## updated 2/10/24 - JMC

library("dplyr")
library("data.table")

##Import functions
source('ML_functions.R')

#input haplotype counts
haplotype_counts <- fread(args[1], header = TRUE, check.names=FALSE)

win <- as.numeric(2) #win size
pxhat <- as.numeric(0.1317675) #assumed px   

#Likelihood calculations
#give windows new numbers?
#yes!

winlist <- as.character(unique(haplotype_counts$Window))
ml_result <- lapply(winlist, function(window) likelihood.out.v(haplotype_counts, window, 
                                                               pxhat = pxhat, shrinkage=FALSE))

ml_output <- do.call(rbind, ml_result)
#ml_output$Window <- as.numeric(ml_output$Window)
#ml_output[order(ml_output$Window),] -> ml_output
ml_output <- ml_output %>% arrange(LRT_p_value)

write.table(ml_output, file=paste(args[2], "ml_output_observed", sep = ""), quote = F, col.names = T, row.names = F)


