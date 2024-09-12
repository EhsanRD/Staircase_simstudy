# Calculate empirical power for staircase designs
# Ehsan Rezaei (Ehsan.rezaeidarzi@monash.edu)

library(MASS)
library(lme4)
library(glmmTMB)
library(pbkrtest)
library(tidyr)
library(parameters)
library(lmerTest)
library(dplyr)
library(TeachingDemos)

setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Staircase_simstudy\\")

source('1. functions_sim.R')

gen_dat <- function(S, K, m, ICC, CAC, theta){
  # Generates a single dataset from a staircase design
  # and trial configuration,for exchangeable and block-exchangeable within-cluster  
  # correlation structure with given correlation parameters and treatment effect
  #
  # Inputs:
  #   S - number of unique treatment sequences
  #   K - number of clusters assigned to each sequence in SCdes
  #   m - cluster-period size
  #   ICC - within-period intracluster correlation
  #   CAC - cluster autocorrelation (CAC=1 returns exchangeable correlation)
  # Output:
  #   One dataset, as a vector with length equal to the total number of
  #   measurements
  
  # Determine covariance parameter values based on given correlation parameters
  # Assuming total variance of 1: sig2T = 1 = sig2C + sig2CP+ sig2e
  sig2C <- ICC*CAC
  sig2e <- 1 - ICC
  sig2CP<- ICC*(1-CAC)
  
  # Determine remaining trial configuration parameters
  # matrix
  
  nclust <- S*K# number of clusters
  
  # Get vector of treatment effect indicators from given design matrix
  SCdes <- SCdesmat(S,K,1,1)
  XvecSC <- rep(as.vector(t(SCdes)), each=m)
  Xvec <- XvecSC[!is.na(XvecSC)]
  #need to be updated
  
  # Set up indices for cluster, cluster-period, and period
  clusterind <- rep(1:nclust, each=2*m)
  clustperind <- rep(1:(nclust*2), each=m)
  # Get the column indices of non-NA values
  non_na_indices <- apply(!is.na(SCdes), 1, function(x) which(x))
  # Repeat the non-NA indices m times
  perind <- rep(non_na_indices, each = m)
  
  # Time period effects (Linear- Andrew adviced)
  betavec <- perind
  
  
  ## Block-exchangeable correlation
  C = rnorm(nclust,mean=0,sqrt(sig2C)) # cluster-level random effects (normalized)
  Cvec = rep(C,each=2*m)
  CP = rnorm(nclust*2,mean=0,sqrt(sig2CP))
  CPvec = rep(CP,each=m)
  
  e <- rnorm(nclust*2*m, mean=0, sd=sqrt(sig2e))
  # Combine terms to get outcome data for SC design
  
  Y  <-  betavec + Xvec*theta + Cvec + CPvec + e
  
  # Create data frame with everything needed for fitting the model
  dat <- data.frame(Y=Y, cluster=as.factor(clusterind), time=perind, clustper=as.factor(clustperind), treat=Xvec)
  
  return(dat)
}

#Fit both models

fitHHmodelSC <- function(dat,typ) {
  tryCatch({
    withCallingHandlers({
      w_convergence_HH  <- w_other_HH  <- e_all_HH <- 0
      if (typ=='cat'){
        lmerfit = lmer(Y ~ treat + as.factor(time) + (1|cluster), data=dat, REML=TRUE)
      }else if (typ=='lin'){
        lmerfit = lmer(Y ~ treat + time + (1|cluster), data=dat, REML=TRUE)
      }
      list(remlfit=lmerfit,
           w_convergence_HH = w_convergence_HH,
           w_other_HH = w_other_HH,
           e_all_HH = e_all_HH)
    },warning = function(w) {
      if(grepl("Model failed to converge", w$message, ignore.case = TRUE)) {
        w_convergence_HH <- 1
      }
      else {
        w_other_HH <- 1  
      }
    },error = function(e) {
      e_all_HH <- 1
    }
    
    )
  }) 
}

fitBEmodelSC <- function(dat,typ) {
  tryCatch({
    withCallingHandlers({
      w_convergence_BE  <- w_other_BE <- e_all_BE <- 0
      if (typ=='cat'){
        lmerfit =lmer(Y ~ treat + as.factor(time) + (1|cluster) + (1|clustper), data=dat, REML=TRUE)
      }else if (typ=='lin') {
        lmerfit =lmer(Y ~ treat + time + (1|cluster) + (1|clustper), data=dat, REML=TRUE)
      }
      list(remlfit=lmerfit,
           w_convergence_BE = w_convergence_BE,
           w_other_BE = w_other_BE,
           e_all_BE = e_all_BE)
    },warning = function(w) {
      if(grepl("Model failed to converge", w$message, ignore.case = TRUE)) {
        w_convergence_BE <- 1
      }
      else {
        w_other_BE <- 1  
      }
    },error = function(e) {
      e_all_BE <- 1
    }
    
    )
  }) 
}


fitmodels <- function(S, K, m, ICC, CAC, theta,typ){
  
  # Generates a single simulated trial dataset, fits corresponding model and
  # outputs the treatment effect estimate and standard error
  
  # Generate dataset
  dat <- gen_dat(S, K, m, ICC, CAC, theta)
  # Fit both models
  fitHHmodelSC_dat <- list()
  fitHHmodelSC_dat <- fitHHmodelSC(dat,typ)
  remlfit_HH    <- fitHHmodelSC_dat[[1]]
  est_trt_HH    <- ifelse(is.null(remlfit_HH), NA, fixef(remlfit_HH)['treat'])
  se_trt_HH     <- ifelse(is.null(remlfit_HH), NA, sqrt(vcov(remlfit_HH)['treat','treat']))
  IsSing_HH     <- ifelse(isSingular(remlfit_HH), 1, 0) 
  
  estvarclustr_HH <-  ifelse(is.null(remlfit_HH), NA, VarCorr(remlfit_HH)$cluster[1])
  estvarres_HH     <- ifelse(is.null(remlfit_HH), NA, sigma(remlfit_HH)^2)
  est_ICC_HH  <- (estvarclustr_HH)/(estvarclustr_HH+estvarres_HH)
  
  fitBEmodelSC_dat <- fitBEmodelSC(dat,typ)
  remlfit_BE    <- fitBEmodelSC_dat[[1]]
  est_trt_BE    <- ifelse(is.null(remlfit_BE), NA, fixef(remlfit_BE)['treat'])
  se_trt_BE     <- ifelse(is.null(remlfit_BE), NA, sqrt(vcov(remlfit_BE)['treat','treat']))
  IsSing_BE     <- ifelse(isSingular(remlfit_BE), 1, 0) 
  
  estvarclustr_BE   <- ifelse(is.null(remlfit_BE), NA, VarCorr(remlfit_BE)$cluster[1])
  estvarclustper_BE <- ifelse(is.null(remlfit_BE), NA, VarCorr(remlfit_BE)$clustper[1]) 
  estvarres_BE      <- ifelse(is.null(remlfit_BE), NA, sigma(remlfit_BE)^2)
  est_ICC_BE <- (estvarclustr_BE+estvarclustper_BE)/(estvarclustr_BE+estvarclustper_BE+estvarres_BE)
  est_CAC_BE <- (estvarclustr_BE)/(estvarclustr_BE+estvarclustper_BE)

  
  # Get adjusted SE & adjusted degree of freedom with Kenward-Roger correction 
  #also for non-singular fit
  adj_se_KR_HH <- adj_ddf_KR_HH <- adj_se_KR_NSing_HH <-adj_ddf_KR_NSing_HH <- NA
  adj_se_Sat_HH <- adj_ddf_Sat_HH <- NA
  if (!is.null(remlfit_HH) & (S<10 | K<10)){
    vcov0_HH <- vcov(remlfit_HH)
    vcovKR_HH <- vcovAdj(remlfit_HH)
    adj_se_KR_HH <- sqrt(vcovKR_HH['treat','treat'])
    adj_se_KR_NSing_HH <- ifelse(IsSing_HH==0, adj_se_KR_HH, NA)
    L1 <- rep(0, length(fixef(remlfit_HH)))
    L1[which(names(fixef(remlfit_HH))=='treat')] <- 1
    adj_ddf_KR_HH <- Lb_ddf(L1, vcov0_HH, vcovKR_HH)
    adj_ddf_KR_NSing_HH <- ifelse(IsSing_HH==0, adj_ddf_KR_HH, NA)
    
    all_Sat_ses_HH <- se_satterthwaite(remlfit_HH)
    adj_se_Sat_HH <- all_Sat_ses_HH[which(all_Sat_ses_HH$Parameter=='treat'),'SE']
    adj_ddf_Sat_HH <- dof_satterthwaite(remlfit_HH)['treat']
  } 
  
  adj_se_KR_BE <- adj_ddf_KR_BE <- adj_se_KR_NSing_BE <- adj_ddf_KR_NSing_BE <-NA 
  adj_se_Sat_BE <- adj_ddf_Sat_BE <- NA
  
  if (!is.null(remlfit_BE) & (S<10 | K<10)){
    
    vcov0_BE<- vcov(remlfit_BE)
    vcovKR_BE <- vcovAdj(remlfit_BE)
    adj_se_KR_BE <- sqrt(vcovKR_BE['treat','treat'])
    adj_se_KR_NSing_BE <- ifelse(IsSing_BE==0, adj_se_KR_BE, NA)
    L2 <- rep(0, length(fixef(remlfit_BE)))
    L2[which(names(fixef(remlfit_BE))=='treat')] <- 1
    adj_ddf_KR_BE <- Lb_ddf(L2, vcov0_BE, vcovKR_BE)
    adj_ddf_KR_NSing_BE <- ifelse(IsSing_BE==0, adj_ddf_KR_BE, NA)
    
    all_Sat_ses_BE <- se_satterthwaite(remlfit_BE)
    adj_se_Sat_BE <- all_Sat_ses_BE[which(all_Sat_ses_BE$Parameter=='treat'),'SE']
    adj_ddf_Sat_BE <- dof_satterthwaite(remlfit_BE)['treat']
  } 
  
  w_conv_HH <- fitHHmodelSC_dat$w_convergence_HH
  w_other_HH <-   fitHHmodelSC_dat$w_other_HH
  err_HH  <-   fitHHmodelSC_dat$e_all_HH
  
  w_conv_BE <- fitBEmodelSC_dat$w_convergence_BE
  w_other_BE <-   fitBEmodelSC_dat$w_other_BE
  err_BE  <-   fitBEmodelSC_dat$e_all_BE
  
  return(c(est_trt_HH, se_trt_HH,est_trt_BE, se_trt_BE,
           est_ICC_HH, est_ICC_BE, est_CAC_BE,
           adj_se_KR_HH,adj_ddf_KR_HH,adj_se_KR_NSing_HH, adj_ddf_KR_NSing_HH,
           adj_se_KR_BE,adj_ddf_KR_BE,adj_se_KR_NSing_BE, adj_ddf_KR_NSing_BE,
           adj_se_Sat_HH,adj_ddf_Sat_HH,adj_se_Sat_BE, adj_ddf_Sat_BE,
           w_conv_HH,w_conv_BE,w_other_HH,w_other_BE,
           IsSing_HH,IsSing_BE,err_HH,err_BE))
}

# Get estimates for nsim repeates
sim_res_fit <- function(nsim, S, K, m, ICC, CAC, theta,typ){
  # Calculates empirical power based on nsim simulated trial datasets
  # Generate trial dataset, fit both models, calculate rejection probability
  res_fit_mat <- replicate(nsim, fitmodels(S, K, m, ICC, CAC, theta,typ))
  res_fit_matx_t <-  t(res_fit_mat)
  #create a data frame convert res matrix to a data frame
  res_fit <- as.data.frame(res_fit_matx_t)
  colnames(res_fit) <- c("est_trt_HH", "se_trt_HH", "est_trt_BE", "se_trt_BE", 
                         "est_ICC_HH", "est_ICC_BE", "est_CAC_BE", 
                         "adj_se_KR_HH","adj_ddf_KR_HH","adj_se_KR_NSing_HH", "adj_ddf_KR_NSing_HH",
                         "adj_se_KR_BE","adj_ddf_KR_BE","adj_se_KR_NSing_BE", "adj_ddf_KR_NSing_BE",
                         "adj_se_Sat_HH","adj_ddf_Sat_HH","adj_se_Sat_BE", "adj_ddf_Sat_BE",
                         "w_conv_HH","w_conv_BE","w_other_HH","w_other_BE",
                         "IsSing_HH","IsSing_BE","err_HH","err_BE")
  
  res_est_fit <- res_fit %>%
    mutate(rep=1:nsim, nsim, S, K, m, ICC, CAC, theta,type=typ)
  # Reorder columns
  res_est_fit <- res_est_fit %>% 
    relocate (rep:type)
  # Save outputs
  save(res_est_fit, file = paste0("est_files\\estimates",
                                  "_nsim", nsim, "_S", S, "_K", K, "_m", m, "_ICC", ICC, "_CAC", CAC,
                                    "_theta", theta, "_", typ, ".RData"))
  
  return(res_est_fit)
  
}





