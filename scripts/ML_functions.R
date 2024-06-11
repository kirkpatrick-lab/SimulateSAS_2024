### Functions used in ML pipeline for UK Biobank data
### updated 4/29/24
### Jared Cole

###############################################

# D conversion (to get r)
D_conv <- function(r,p1,p2){
  D <- r*(sqrt(p1*(1-p1))*sqrt(p2*(1-p2)))
  return(D)
}

#r conversion (to get D)
r_conv <- function(D,p1,p2){
  r <- D/(sqrt(p1*(1-p1))*sqrt(p2*(1-p2)))
  return(r)
}

# zygote haplotype frequiencies for 3-sites
hap_f_calcs <- function(p0,p1,px,D_1x,D_x0,D_10,k){
  
  q0 <- 1-p0
  q1 <- 1-p1
  qx <- 1-px
  
  #with given value of k
  f0_0 <- q0*q1*qx + qx*D_10 + k*q0*D_1x + k*q1*D_x0
  f1_0 <- p0*q1*qx - qx*D_10 + k*p0*D_1x - k*q1*D_x0
  f2_0 <- p1*q0*qx - qx*D_10 - k*q0*D_1x + k*p1*D_x0
  f3_0 <- p0*p1*qx + qx*D_10 - k*p0*D_1x - k*p1*D_x0
  f0_1 <- px*q0*q1 + px*D_10 - k*q0*D_1x - k*q1*D_x0
  f1_1 <- p0*px*q1 - px*D_10 - k*p0*D_1x + k*q1*D_x0
  f2_1 <- p1*px*q0 - px*D_10 + k*q0*D_1x - k*p1*D_x0
  f3_1 <- p0*p1*px + px*D_10 + k*p0*D_1x + k*p1*D_x0
  
  fq <- list(f0_0 = f0_0, f1_0 = f1_0, f2_0 = f2_0, f3_0 = f3_0, 
             f0_1 = f0_1, f1_1 = f1_1, f2_1 = f2_1, f3_1 = f3_1)
  
  return(fq)
}

### Infer the haplotype frequencies in zygotes with 2 observable sites

# get zygote freqs (assumes target in center)
hap_freqs <- function(f,m,p, k=1, shrinkage=FALSE){
  
  #counts
  fcounts <- f
  mcounts <- m
  
  #sums
  NF <- sum(fcounts)
  NM <- sum(mcounts)
  
  #allele freqs
  pF0 <- (fcounts[2]+fcounts[4])/NF
  pF1 <- (fcounts[3]+fcounts[4])/NF
  
  pM0 <- (mcounts[2]+mcounts[4])/NM
  pM1 <- (mcounts[3]+mcounts[4])/NM
  
  qF0 <- 1-pF0
  qF1 <- 1-pF1
  
  qM0 <- 1-pM0
  qM1 <- 1-pM1
  
  p0 <- (pF0+pM0)/2
  p1 <- (pF1+pM1)/2
  px <- p
  
  q0 <- 1-p0
  q1 <- 1-p1
  
  qx <- 1-px
  
  ## LD (infer from two sites, averaged between sexes)
  f4 <- ((fcounts[4]/NF)+(mcounts[4]/NM))/2
  D_10 <- f4 - (p0*p1)
  
  I_D_10 <- ifelse(D_10 < 0, -1, 1) #indicator
  
  D_1x <- (abs(D_10)^(0.5))*I_D_10*(px*qx)^(0.5)*((p1*q1)/(p0*q0))^(1/4)
  D_x0 <- (abs(D_10)^(0.5))*(px*qx)^(0.5)*((p0*q0)/(p1*q1))^(1/4)
  
  r10 <- D_10/sqrt(p0*q0*p1*q1)
  r1x <- D_1x/sqrt(p1*q1*px*qx)
  rx0 <- D_x0/sqrt(p0*q0*px*qx)
  
  #starting haplotype freqs, determine k
  
  if (shrinkage==TRUE){
    
    k<-1
    
    repeat {
      hap_terms <- hap_f_calcs(p0,p1,px,D_1x,D_x0,D_10,k)
      
      k_values <- numeric() # store k
      
      for (term_name in names(hap_terms)) {
        term_value <- hap_terms[[term_name]]
        if (!is.na(term_value) && term_value < 0) {
          if (term_name == "f0_0") {
            k_values <- c(k_values, (-q0*q1*qx - qx*D_10)/(q0*D_1x + q1*D_x0))
          } else if (term_name == "f1_0") {
            k_values <- c(k_values, (-p0*q1*qx + qx*D_10)/(p0*D_1x - q1*D_x0))
          } else if (term_name == "f2_0") {
            k_values <- c(k_values, (-p1*q0*qx + qx*D_10)/(-q0*D_1x + p1*D_x0))
          } else if (term_name == "f3_0") {
            k_values <- c(k_values, (p0*p1*qx + qx*D_10)/(p0*D_1x + p1*D_x0))
          } else if (term_name == "f0_1") {
            k_values <- c(k_values, (px*q0*q1 + px*D_10)/(q0*D_1x + q1*D_x0))
          } else if (term_name == "f1_1") {
            k_values <- c(k_values, (p0*px*q1 - px*D_10)/(p0*D_1x - q1*D_x0))
          } else if (term_name == "f2_1") {
            k_values <- c(k_values, (p1*px*q0 - px*D_10)/(-q0*D_1x + p1*D_x0))
          } else if (term_name == "f3_1") {
            k_values <- c(k_values, (-p0*p1*px - px*D_10)/(p0*D_1x + p1*D_x0))
          }
        }
      }
      
      # If all terms are non-negative, we exit the loop
      if (!any(is.na(unlist(hap_terms))) && all(unlist(hap_terms) >= 0)) {
        break
      } else {
        # Update k to the minimum positive k value calculated
        k <- min(k_values[k_values > 0])
        k <- as.numeric(substr(as.character(k), 1, 8))
      }
    }
    
    hap_terms<-hap_f_calcs(p0,p1,px,D_1x,D_x0,D_10,k)
    
  } else {
    
    hap_terms<-hap_f_calcs(p0,p1,px,D_1x,D_x0,D_10,k)
    
  }
  
  hapfreqs_0 <- c(hap_terms$f0_0,hap_terms$f1_0,hap_terms$f2_0,hap_terms$f3_0)
  hapfreqs_1 <- c(hap_terms$f0_1,hap_terms$f1_1,hap_terms$f2_1,hap_terms$f3_1)
  rvals <- c(r10,r1x,rx0)
  Dvals <- c(D_10,D_1x,D_x0)
  
  frqs <- list("frq_0" = hapfreqs_0,
               "frq_1" = hapfreqs_1,
               "rvals" = rvals,
               "Dvals" = Dvals)
  
  return(frqs)
}

## Interpolation likelihood function

lik.function.n.log <- function(s, p, fem_counts, male_counts, freqs0, freqs1){
  
  #hap freqs
  f0_0 <- freqs0[1]
  f1_0 <- freqs0[2]
  f2_0 <- freqs0[3]
  f3_0 <- freqs0[4]
  f0_1 <- freqs1[1]
  f1_1 <- freqs1[2]
  f2_1 <- freqs1[3]
  f3_1 <- freqs1[4]
  
  #hap counts
  nF0 <- fem_counts[1]
  nF1 <- fem_counts[2]
  nF2 <- fem_counts[3]
  nF3 <- fem_counts[4]
  
  nM0 <- male_counts[1]
  nM1 <- male_counts[2]
  nM2 <- male_counts[3]
  nM3 <- male_counts[4]
  
  #alleles
  px <- p
  qx <- 1-px  
  
  
  #after selection (take logs)
  F0 <- ifelse(nF0 == 0, 0, log(((1-s*px)*f0_0)+((1+s*qx)*f0_1))*nF0)
  F1 <- ifelse(nF1 == 0, 0, log(((1-s*px)*f1_0)+((1+s*qx)*f1_1))*nF1)
  F2 <- ifelse(nF2 == 0, 0, log(((1-s*px)*f2_0)+((1+s*qx)*f2_1))*nF2)
  F3 <- ifelse(nF3 == 0, 0, log(((1-s*px)*f3_0)+((1+s*qx)*f3_1))*nF3)
  
  M0 <- ifelse(nM0 == 0, 0, log(((1+s*px)*f0_0)+((1-s*qx)*f0_1))*nM0)
  M1 <- ifelse(nM1 == 0, 0, log(((1+s*px)*f1_0)+((1-s*qx)*f1_1))*nM1)
  M2 <- ifelse(nM2 == 0, 0, log(((1+s*px)*f2_0)+((1-s*qx)*f2_1))*nM2)
  M3 <- ifelse(nM3 == 0, 0, log(((1+s*px)*f3_0)+((1-s*qx)*f3_1))*nM3)
  
  #log likelihood
  LL <- sum(F0,F1,F2,F3,M0,M1,M2,M3)
  
  return(LL)
}

#Optimization function for interpolation (viability)
likelihood.out.v <- function(data,window, pxhat, shrinkage=FALSE){
  
  #output DF
  output <- data.frame()
  
  #list of haplos from data
  haplist <- as.numeric(unique(unlist(data[data$Window == window,3])))
  expected_haps <- c(0,1,2,3)
  
  #variables 
  fem_counts <- NULL
  male_counts <- NULL
  freqs0 <- NULL
  freqs1 <- NULL
  px_hat <- pxhat
  
  #get counts
  window_counts <- data[data$Window == window,4:6]
  missing_haplotypes <- setdiff(expected_haps, window_counts$haplotype)
  missing_data <- data.frame(haplotype = missing_haplotypes,
                             female_counts = rep(0, length(missing_haplotypes)),
                             male_counts = rep(0, length(missing_haplotypes)))
  combined_data <- rbind(window_counts, missing_data)
  
  fem_counts <- as.numeric(unlist(combined_data[order(combined_data$haplotype), "female_counts"]))
  male_counts <- as.numeric(unlist(combined_data[order(combined_data$haplotype), "male_counts"]))
  
  #get freqs
  freqs <- hap_freqs(fem_counts, male_counts, p = px_hat, shrinkage = shrinkage)
  freqs0 <- freqs$frq_0
  freqs1 <- freqs$frq_1
  
  #Null ML
  ml_null<-lik.function.n.log(s=0,p=px_hat,fem_counts = fem_counts, male_counts = male_counts,
                              freqs0 = freqs0, freqs1 = freqs1)
  
  #optimize (dynamic intervals)
  ml <- suppressWarnings(optimize(lik.function.n.log, interval = c(-5,5), maximum = TRUE,
                                  p = px_hat, fem_counts = fem_counts, male_counts = male_counts,
                                  freqs0 = freqs0, freqs1 = freqs1))
  
  if (is.nan(ml$objective)) {
    max_s <- 0.5
    inp_s <- seq(-max_s,max_s,length=10000)
    vect.lik.func <- Vectorize(lik.function.n.log, vectorize.args = "s")
    int_t <- suppressWarnings(range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                                  fem_counts = fem_counts, male_counts = male_counts,
                                                                  freqs0 = freqs0, freqs1 = freqs1))]))
    if (any(is.nan(int_t)) || any(is.infinite(int_t)) || int_t[1] >= 0 || int_t[2] < 0) {
      inp_s <- seq(-max_s, max_s, length=100000) # Increase length to 100000
      int_t <- suppressWarnings(range(inp_s[is.finite(vect.lik.func(inp_s, p = px_hat,
                                                                    fem_counts = fem_counts, male_counts = male_counts,
                                                                    freqs0 = freqs0, freqs1 = freqs1))]))
    }
    
    ml <- optimize(lik.function.n.log, interval = c(int_t[1],int_t[2]), maximum = TRUE,
                   p = px_hat, fem_counts = fem_counts, male_counts = male_counts,
                   freqs0 = freqs0, freqs1 = freqs1)
  }

  ### get results
  s_int <- ml$maximum
  maxl_int <- ml$objective
  
  #calculate p-values (LRT)
  LRT_stat<- -2 * (as.numeric(ml$objective) - as.numeric(ml_null))
  p_val <- pchisq(abs(LRT_stat), df=1, lower.tail = FALSE)
  
  #output info
  chr <- unique(data$Chrom)
  pos <- unique(data[data$Window == window, ]$Leading_Position)
  snpid <- unique(data[data$Window == window, ]$Leading_SNP)
  lmaf <- unique(data[data$Window == window, ]$Leading_MAF)
  al0 <- unique(data[data$Window == window, ]$Leading_Allele0)
  al1 <- unique(data[data$Window == window, ]$Leading_Allele1)
  r10 <- freqs$rvals[1]
  r1x <- freqs$rvals[2]
  rx0 <- freqs$rvals[3]
  D10 <- freqs$Dvals[1]
  D1x <- freqs$Dvals[2]
  Dx0 <- freqs$Dvals[3]
  
  #Output results 
  out <- data.frame(chr,window,pos,snpid,lmaf,al0,al1,r10,r1x,rx0,D10,D1x,Dx0,   
                    s_int, maxl_int, p_val)
  
  colnames(out) <- (c("Chrom","Window","Leading_Position","Leading_SNP","Leading_MAF","Leading_Allele0","Leading_Allele1",
                      "r_10", "r_1x", "r_x0", "D_10","D_1x","D_x0",
                      "s_v_int","ML_v_int","LRT_p_value"))
  
  output <- rbind(output, out)
  
  return(output)
}
