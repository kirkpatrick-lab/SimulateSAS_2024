#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(paste(args[1], "is the counts file you would like to use."))
print(paste(args[2], "is where you would like your output."))

#### Runs the SAS simulation code remotely with two inputs - a "counts file", 
#giving number of observations per haplotype per sex per window across the genome
#and an output location, which will have the simulated values of 's' and 'F'

#Important! You will also need the MAF for the SNPs in your working directory and 
#The r2 values for those SNPs (see below on lines 23-34)
#AND! The ml calculations for the observed data in your working directory. 
#This must be named ml_output_observed, but the upstream script should do this for you.

####### Simulate SAS - interpolation code from Jared 2/2/24. Included again here
#for troubleshooting purposes.

library(data.table)
library(tidyverse)

##Import functions
source('ML_functions.R')

#make sure these files are in your working directory.
mafs <- scan(file = "UKB_mafs_filtered.txt") #distribution of MAFs
r2_vals <- scan(file = "UKB_r2_filtered.txt") #distribution of r^2 values

ukb_r2_mafs <- fread("UKB_mafs_r2.txt", 
                     sep="\t", header = TRUE, check.names=FALSE) #pairwise table of MAFs and r^2

counts <- fread(args[1], header = T) #this is your haplotype counts across all windows you want
output_dest=args[2] #output destination

####################################################################################
## Functions for simulations

##simulate counts 

##simulate counts 
sim_counts <- function(num_F,num_M,p1,px,p0,ld,s, 
                       use_r=FALSE, 
                       sampling=FALSE,
                       sample_p=FALSE){
  
  d<-runif(1, min = 0, max = 1) #sample d, uniform between 0-1
  
  if (sample_p == TRUE) {
    
    #sample pair of SNPs
    param_samp <- ukb_r2_mafs[sample(nrow(ukb_r2_mafs), 1), ]
    
    #unknown site frequency
    px <- sample(mafs, 1) #sample from MAF dist
    
    #pairwise MAFs
    p1 <- param_samp$MAF_A  
    p0 <- param_samp$MAF_B
    
    q1 <- 1 - p1
    qx <- 1 - px
    q0 <- 1 - p0
    
    #ld from pair
    ld <- sqrt(as.numeric(param_samp$R2)) #in terms of correlation r
    
    D_10 <- D_conv(ld, p1, p0)
    
    I_D_10 <- ifelse(D_10 < 0, -1, 1) #indicator
    
    D_1x <- I_D_10 * (D_10^(1-d)) * ((p1*q1)^(d/2)) * ((p0*q0)^((d-1)/2)) * ((px*qx)^(0.5))
    D_x0 <- (D_10^d) * ((p1*q1)^(-d/2)) * ((p0*q0)^((1-d)/2)) * ((px*qx)^(0.5))
    
    
  } else {
    
    q1 <- 1 - p1
    qx <- 1 - px
    q0 <- 1 - p0
    
    D_10 <- ifelse(use_r == TRUE, D_conv(ld, p1, p0), ld)
    
    I_D_10 <- ifelse(D_10 < 0, -1, 1) #indicator
    
    D_1x <- I_D_10 * (D_10^(1-d)) * ((p1*q1)^(d/2)) * ((p0*q0)^((d-1)/2)) * ((px*qx)^(0.5))
    D_x0 <- (D_10^d) * ((p1*q1)^(-d/2)) * ((p0*q0)^((1-d)/2)) * ((px*qx)^(0.5))
    
    
  }
  
  #get hap freqs
  hap_terms<-hap_f_calcs(p0,p1,px,D_1x,D_x0,D_10,k=1)
  
  #save
  f0_0 <- hap_terms$f0_0
  f1_0 <- hap_terms$f1_0
  f2_0 <- hap_terms$f2_0
  f3_0 <- hap_terms$f3_0
  
  f0_1 <- hap_terms$f0_1
  f1_1 <- hap_terms$f1_1
  f2_1 <- hap_terms$f2_1
  f3_1 <- hap_terms$f3_1
  
  fem0 <- ((1-s*px)*f0_0) + ((1+s*qx)*f0_1)
  fem1 <- ((1-s*px)*f1_0) + ((1+s*qx)*f1_1)
  fem2 <- ((1-s*px)*f2_0) + ((1+s*qx)*f2_1)
  fem3 <- ((1-s*px)*f3_0) + ((1+s*qx)*f3_1)
  
  male0 <- ((1+s*px)*f0_0) + ((1-s*qx)*f0_1)
  male1 <- ((1+s*px)*f1_0) + ((1-s*qx)*f1_1)
  male2 <- ((1+s*px)*f2_0) + ((1-s*qx)*f2_1)
  male3 <- ((1+s*px)*f3_0) + ((1-s*qx)*f3_1)
  
  while (any(c(fem0, fem1, fem2, fem3, male0, male1, male2, male3) < 0)) {
    
    #sample pair of SNPs
    param_samp <- ukb_r2_mafs[sample(nrow(ukb_r2_mafs), 1), ]
    
    #unknown site frequency
    px <- sample(mafs, 1) #sample from MAF dist
    
    #pairwise MAFs
    p1 <- param_samp$MAF_A  
    p0 <- param_samp$MAF_B
    
    q1 <- 1 - p1
    qx <- 1 - px
    q0 <- 1 - p0
    
    #ld from pair
    ld <- sqrt(as.numeric(param_samp$R2)) #in terms of correlation r
    
    D_10 <- D_conv(ld, p1, p0)
    
    I_D_10 <- ifelse(D_10 < 0, -1, 1) #indicator
    
    D_1x <- I_D_10 * (D_10^(1-d)) * ((p1*q1)^(d/2)) * ((p0*q0)^((d-1)/2)) * ((px*qx)^(0.5))
    D_x0 <- (D_10^d) * ((p1*q1)^(-d/2)) * ((p0*q0)^((1-d)/2)) * ((px*qx)^(0.5))
    
    #hap terms
    hap_terms<-hap_f_calcs(p0,p1,px,D_1x,D_x0,D_10,k=1)
    
    f0_0 <- hap_terms$f0_0
    f1_0 <- hap_terms$f1_0
    f2_0 <- hap_terms$f2_0
    f3_0 <- hap_terms$f3_0
    
    f0_1 <- hap_terms$f0_1
    f1_1 <- hap_terms$f1_1
    f2_1 <- hap_terms$f2_1
    f3_1 <- hap_terms$f3_1
    
    fem0 <- ((1-s*px)*f0_0) + ((1+s*qx)*f0_1)
    fem1 <- ((1-s*px)*f1_0) + ((1+s*qx)*f1_1)
    fem2 <- ((1-s*px)*f2_0) + ((1+s*qx)*f2_1)
    fem3 <- ((1-s*px)*f3_0) + ((1+s*qx)*f3_1)
    
    male0 <- ((1+s*px)*f0_0) + ((1-s*qx)*f0_1)
    male1 <- ((1+s*px)*f1_0) + ((1-s*qx)*f1_1)
    male2 <- ((1+s*px)*f2_0) + ((1-s*qx)*f2_1)
    male3 <- ((1+s*px)*f3_0) + ((1-s*qx)*f3_1)
    
  }
  
  if (sampling==TRUE){
    
    prob_f <- c(fem0,fem1,fem2,fem3)
    prob_m <- c(male0,male1,male2,male3)
    
    nF <- rmultinom(1, size = num_F, prob = prob_f)
    nM <- rmultinom(1, size = num_M, prob = prob_m)
    
    f_haps <- c(nF[1], nF[2], nF[3], nF[4])
    m_haps <- c(nM[1], nM[2], nM[3], nM[4])
    
  } else {
    
    nF0 <- fem0*num_F
    nF1 <- fem1*num_F
    nF2 <- fem2*num_F
    nF3 <- fem3*num_F
    
    nM0 <- male0*num_M
    nM1 <- male1*num_M
    nM2 <- male2*num_M
    nM3 <- male3*num_M
    
    f_haps <- c(nF0, nF1, nF2, nF3)
    m_haps <- c(nM0, nM1, nM2, nM3)
  }
  
  
  h_counts <- list("females" = f_haps, 
                   "males" = m_haps,
                   "ps" = c(p1,px,p0),
                   "d" = d,
                   "ld_D10" = D_10,
                   "ld_r" = r_conv(D_10,p1,p0),
                   "freqs" = hap_terms)
  
  return(h_counts)
  
}

#Simulate SAS and estimate selection coefficient 

sim_SAS <- function(nf,nm,s,p1,p0,px,ld,d=0.5,use_r=FALSE,
                    sampling=FALSE,
                    sample_p=FALSE){ 
  
  hcnt <- sim_counts(num_F = nf, num_M = nm, s=s,
                     p1 = p1, px = px, p0 = p0,
                     ld=ld,
                     use_r = use_r, 
                     sampling = sampling, sample_p = sample_p)
  
  px_hat = 0.1317675 #median MAF for all imputed SNPs
  
  #counts
  fem_counts <- hcnt$females
  male_counts <- hcnt$males
  
  #get freqs
  freqs <- hap_freqs(fem_counts, male_counts, p = px_hat, shrinkage = FALSE)
  freqs0 <- freqs$frq_0
  freqs1 <- freqs$frq_1
  
  #Null ML
  ml_null<-lik.function.n.log(s=0,p=px_hat,fem_counts = fem_counts, male_counts = male_counts,
                              freqs0 = freqs0, freqs1 = freqs1)
  
  
  #optimize (dynamic intervals)
  ml <- optimize(lik.function.n.log, interval = c(-5,5), maximum = TRUE,
                 p = px_hat, fem_counts = fem_counts, male_counts = male_counts,
                 freqs0 = freqs0, freqs1 = freqs1)
  
  if (is.nan(ml$objective)) {
    max_s <- 0.5
    inp_s <- seq(-max_s,max_s,length=10000)
    vect.lik.func <- Vectorize(lik.function.n.log, vectorize.args = "s")
    int_t <- range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                 fem_counts = fem_counts, male_counts = male_counts,
                                                 freqs0 = freqs0, freqs1 = freqs1))])
    if (any(is.nan(int_t)) || any(is.infinite(int_t)) || int_t[1] >= 0 || int_t[2] < 0) {
      inp_s <- seq(-max_s, max_s, length=100000) # Increase length to 100000
      int_t <- range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                   fem_counts = fem_counts, male_counts = male_counts,
                                                   freqs0 = freqs0, freqs1 = freqs1))])
    }
    
    ml <- optimize(lik.function.n.log, interval = c(int_t[1],int_t[2]), maximum = TRUE,
                   p = px_hat, fem_counts = fem_counts, male_counts = male_counts,
                   freqs0 = freqs0, freqs1 = freqs1)
  }
  
  #calculate p-values (LRT)
  LRT_stat<- -2 * (as.numeric(ml$objective) - as.numeric(ml_null))
  p_val <- pchisq(abs(LRT_stat), df=1, lower.tail = FALSE)
  
  #out
  out<-list("maximum" = ml$maximum,
            "objective" = ml$objective,
            "null_objective" = ml_null,
            "pvalue" = p_val)
  
  return(out)
}



data <- counts %>% group_by(Window) %>%
  mutate(all_counts = male_counts + female_counts)
nM_byw = data %>% summarise(nM = sum(male_counts))
nF_byw = data %>% summarise(nF = sum(female_counts))
hap_freq <- data %>% mutate(freqF = female_counts/sum(female_counts), freqM = male_counts/sum(male_counts), freq_all = all_counts/sum(all_counts))

#set my priors here, because we're going to do all windows in parallel
snp_under_sel <- runif(1, 0, 1)

nruns = 1
x_b = runif(nruns, 0, 1) # just sample one S for whole run
s_prior = 10^(-(4*x_b + 1))

x_f = runif(nruns, 0, 1) #and, similarly, just sample one F for whole run 
f_prior = 10^(-4*x_f)


uq_windows <- data$Window %>% unique()
MLout <- data.frame()
start <- Sys.time()
for(i in 1:length(uq_windows)){ #remove head
    snp_under_sel <- runif(1, 0, 1)
    b_alt <- ifelse(f_prior < snp_under_sel, 0, s_prior) #set s to zero based on the value of F - this connects the strength of selection
    #to the proportion of SNPs under selection
    simulate1 <- sim_SAS(nf=nF_byw[nF_byw$Window ==uq_windows[i],]$nF, nm=nM_byw[nM_byw$Window==uq_windows[i],]$nM, s=b_alt, sample_p = TRUE, sampling=TRUE)
    MLout <- rbind(MLout, simulate1)
}
end <- Sys.time()
runtime = end - start #for you to check the time of your run, if you would like to estimate total runtime 
#print(runtime)

MLout$s_prior <- s_prior
MLout$F_prior <- f_prior
#Write out all values from the ML function, in case of downstream errors. (You won't have to run the simulation again)
write.table(MLout, file=paste(output_dest, ".backupMLout", sep=""), col.names = F, quote = F, row.names = F)

#arrange by p-value to get distribution
MLout <- MLout %>% arrange(pvalue)

#Read in the "observed" ML values from the real data
obs <- read.table("ml_output_observed", header=T) %>% select(Window, LRT_p_value) %>% unique()

#Make sure we have the same number of windows in our observed versus simulated data
if(nrow(obs)==nrow(MLout)){
  outdf <- data.frame(s_prior=s_prior, f_prior=f_prior, SSE=sum((MLout$pvalue - obs$LRT_p_value)^2))
}

if(nrow(obs)!=nrow(MLout)){
  print("Help me! I don't have the correct number of entries in my distribution!")
}

#write your outputs (just the simulated S, simulated F, and SSE between this distribution and the observed distribution. Should only be one line)
write.table(outdf, file=output_dest, col.names = F, quote = F, row.names = F)





