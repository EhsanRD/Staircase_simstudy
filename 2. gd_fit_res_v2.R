# Calculate empirical power for staircase designs
# Ehsan Rezaei (Ehsan.rezaeidarzi@monash.edu)
# S=6
# K=1
# m=4
# ICC=0.2
# CAC=0.5
# theta=0.15
# typ='cat'

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
  
  # Time period effects (Linear- Andrew advised)
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

fitHHmodelSC_c <- function(dat) {
  tryCatch({
    withCallingHandlers({
      w_convergence_HH_c  <- w_other_HH_c  <- e_all_HH_c <- 0
        lmerfit_c = lmer(Y ~ treat + as.factor(time) + (1|cluster), data=dat, REML=TRUE)
      list(remlfit_c=lmerfit_c,
           w_convergence_HH_c = w_convergence_HH_c,
           w_other_HH_c = w_other_HH_c,
           e_all_HH_c = e_all_HH_c
           )
    },warning = function(w) {
      if(grepl("Model failed to converge", w$message, ignore.case = TRUE)) {
        w_convergence_HH_c <- 1
      }
      else {
        w_other_HH_c <- 1  
      }
    },error = function(e) {
      e_all_HH_c <- 1
    }
    
    )
  }) 
}
fitHHmodelSC_l <- function(dat) {
  tryCatch({
    withCallingHandlers({
      w_convergence_HH_l  <- w_other_HH_l  <- e_all_HH_l <- 0
      lmerfit_l = lmer(Y ~ treat + time + (1|cluster), data=dat, REML=TRUE)
      list(remlfit_l=lmerfit_l,
           w_convergence_HH_l = w_convergence_HH_l,
           w_other_HH_l = w_other_HH_l,
           e_all_HH_l = e_all_HH_l)
    },warning = function(w) {
      if(grepl("Model failed to converge", w$message, ignore.case = TRUE)) {
        w_convergence_HH_l <- 1
      }
      else {
        w_other_HH_l <- 1  
      }
    },error = function(e) {
      e_all_HH_l <- 1
    }
    
    )
  }) 
}
fitBEmodelSC_c <- function(dat) {
  tryCatch({
    withCallingHandlers({
      w_convergence_BE_c  <- w_other_BE_c <- e_all_BE_c <- 0
        lmerfit_c =lmer(Y ~ treat + as.factor(time) + (1|cluster) + (1|clustper), data=dat, REML=TRUE)
      list(remlfit_c=lmerfit_c,
           w_convergence_BE_c = w_convergence_BE_c,
           w_other_BE_c = w_other_BE_c,
           e_all_BE_c = e_all_BE_c)
    },warning = function(w) {
      if(grepl("Model failed to converge", w$message, ignore.case = TRUE)) {
        w_convergence_BE_c <- 1
      }
      else {
        w_other_BE_c <- 1  
      }
    },error = function(e) {
      e_all_BE_c <- 1
    }
    
    )
  }) 
}

fitBEmodelSC_l <- function(dat) {
  tryCatch({
    withCallingHandlers({
      w_convergence_BE_l <- w_other_BE_l <- e_all_BE_l <- 0
      lmerfit_l =lmer(Y ~ treat + time + (1|cluster) + (1|clustper), data=dat, REML=TRUE)
      list(remlfit_l=lmerfit_l,
           w_convergence_BE_l = w_convergence_BE_l,
           w_other_BE_l = w_other_BE_l,
           e_all_BE_l = e_all_BE_l)
    },warning = function(w) {
      if(grepl("Model failed to converge", w$message, ignore.case = TRUE)) {
        w_convergence_BE_l <- 1
      }
      else {
        w_other_BE_l <- 1  
      }
    },error = function(e) {
      e_all_BE_l <- 1
    }
    
    )
  }) 
}



fitmodels <- function(S, K, m, ICC, CAC, theta){
  
  # Generates a single simulated trial dataset, fits corresponding model and
  # outputs the treatment effect estimate and standard error
  
  # Generate dataset
  dat <- gen_dat(S, K, m, ICC, CAC, theta)
  # Fit both models
  fitHHmodelSC_cdat <- list()
  fitHHmodelSC_cdat <- fitHHmodelSC_c(dat)
  remlfit_HH_c    <- fitHHmodelSC_cdat[[1]]
  est_trt_HH_c    <- ifelse(is.null(remlfit_HH_c), NA, fixef(remlfit_HH_c)['treat'])
  se_trt_HH_c     <- ifelse(is.null(remlfit_HH_c), NA, sqrt(vcov(remlfit_HH_c)['treat','treat']))
  IsSing_HH_c     <- ifelse(isSingular(remlfit_HH_c), 1, 0) 
  
  estvarclustr_HH_c <-  ifelse(is.null(remlfit_HH_c), NA, VarCorr(remlfit_HH_c)$cluster[1])
  estvarres_HH_c     <- ifelse(is.null(remlfit_HH_c), NA, sigma(remlfit_HH_c)^2)
  est_ICC_HH_c  <- (estvarclustr_HH_c)/(estvarclustr_HH_c+estvarres_HH_c)
  
  fitHHmodelSC_ldat <- list()
  fitHHmodelSC_ldat <- fitHHmodelSC_l(dat)
  remlfit_HH_l    <- fitHHmodelSC_ldat[[1]]
  est_trt_HH_l    <- ifelse(is.null(remlfit_HH_l), NA, fixef(remlfit_HH_l)['treat'])
  se_trt_HH_l     <- ifelse(is.null(remlfit_HH_l), NA, sqrt(vcov(remlfit_HH_l)['treat','treat']))
  IsSing_HH_l     <- ifelse(isSingular(remlfit_HH_l), 1, 0) 
  
  estvarclustr_HH_l <-  ifelse(is.null(remlfit_HH_l), NA, VarCorr(remlfit_HH_l)$cluster[1])
  estvarres_HH_l     <- ifelse(is.null(remlfit_HH_l), NA, sigma(remlfit_HH_l)^2)
  est_ICC_HH_l  <- (estvarclustr_HH_l)/(estvarclustr_HH_l+estvarres_HH_l)
  
  fitBEmodelSC_cdat <- list()
  fitBEmodelSC_cdat <- fitBEmodelSC_c(dat)
  remlfit_BE_c    <- fitBEmodelSC_cdat[[1]]
  est_trt_BE_c    <- ifelse(is.null(remlfit_BE_c), NA, fixef(remlfit_BE_c)['treat'])
  se_trt_BE_c     <- ifelse(is.null(remlfit_BE_c), NA, sqrt(vcov(remlfit_BE_c)['treat','treat']))
  IsSing_BE_c     <- ifelse(isSingular(remlfit_BE_c), 1, 0) 
  
  estvarclustr_BE_c   <- ifelse(is.null(remlfit_BE_c), NA, VarCorr(remlfit_BE_c)$cluster[1])
  estvarclustper_BE_c <- ifelse(is.null(remlfit_BE_c), NA, VarCorr(remlfit_BE_c)$clustper[1]) 
  estvarres_BE_c      <- ifelse(is.null(remlfit_BE_c), NA, sigma(remlfit_BE_c)^2)
  est_ICC_BE_c <- (estvarclustr_BE_c+estvarclustper_BE_c)/(estvarclustr_BE_c+estvarclustper_BE_c+estvarres_BE_c)
  est_CAC_BE_c <- (estvarclustr_BE_c)/(estvarclustr_BE_c+estvarclustper_BE_c)

  fitBEmodelSC_ldat <- list()
  fitBEmodelSC_ldat <- fitBEmodelSC_l(dat)
  remlfit_BE_l    <- fitBEmodelSC_ldat[[1]]
  est_trt_BE_l    <- ifelse(is.null(remlfit_BE_l), NA, fixef(remlfit_BE_l)['treat'])
  se_trt_BE_l     <- ifelse(is.null(remlfit_BE_l), NA, sqrt(vcov(remlfit_BE_l)['treat','treat']))
  IsSing_BE_l     <- ifelse(isSingular(remlfit_BE_l), 1, 0) 
  
  estvarclustr_BE_l   <- ifelse(is.null(remlfit_BE_l), NA, VarCorr(remlfit_BE_l)$cluster[1])
  estvarclustper_BE_l <- ifelse(is.null(remlfit_BE_l), NA, VarCorr(remlfit_BE_l)$clustper[1]) 
  estvarres_BE_l      <- ifelse(is.null(remlfit_BE_l), NA, sigma(remlfit_BE_l)^2)
  est_ICC_BE_l <- (estvarclustr_BE_l+estvarclustper_BE_l)/(estvarclustr_BE_l+estvarclustper_BE_l+estvarres_BE_l)
  est_CAC_BE_l <- (estvarclustr_BE_l)/(estvarclustr_BE_l+estvarclustper_BE_l)
  
  # Get adjusted SE & adjusted degree of freedom with Kenward-Roger correction 
  #also for non-singular fit
  adj_se_KR_HH_c <- adj_ddf_KR_HH_c <- adj_se_KR_NSing_HH_c <-adj_ddf_KR_NSing_HH_c <- NA
  adj_se_Sat_HH_c <- adj_ddf_Sat_HH_c <- NA

  
  if (!is.null(remlfit_HH_c) & (S<10 | K<10)){
    ###pbkrtest (KR)
    vcov0_HH_c <- vcov(remlfit_HH_c)
    vcovKR_HH_c <- vcovAdj(remlfit_HH_c)
    adj_se_KR_HH_c <- sqrt(vcovKR_HH_c['treat','treat'])
    adj_se_KR_NSing_HH_c <- ifelse(IsSing_HH_c==0, adj_se_KR_HH_c, NA)
    L1 <- rep(0, length(fixef(remlfit_HH_c)))
    L1[which(names(fixef(remlfit_HH_c))=='treat')] <- 1
    adj_ddf_KR_HH_c <- Lb_ddf(L1, vcov0_HH_c, vcovKR_HH_c)
    adj_ddf_KR_NSing_HH_c <- ifelse(IsSing_HH_c==0, adj_ddf_KR_HH_c, NA)
    ###parameters (Sat)
    all_Sat_ses_HH_c <- se_satterthwaite(remlfit_HH_c)
    adj_se_Sat_HH_c <- all_Sat_ses_HH_c[which(all_Sat_ses_HH_c$Parameter=='treat'),'SE']
    adj_ddf_Sat_HH_c <- dof_satterthwaite(remlfit_HH_c)['treat']
    # microbenchmark(adj_se_KR_HH_c <- sqrt(vcovKR_HH_c['treat','treat']),
    #                all_KR_ses_HH_c <- se_kenward(remlfit_HH_c) )
  } 
  
  # if (!is.null(remlfit_HH_c) & ((S==4 & K==1) | (S==4 & K==5) | (S==10 & K==1))){
  # ###parameters (KR)
  # all_KR_ses_HH_c <- se_kenward(remlfit_HH_c)
  # adj_se_KR_HH_c <- all_KR_ses_HH_c[which(all_KR_ses_HH_c$Parameter=='treat'),'SE']
  # adj_se_KR_NSing_HH_c <- ifelse(IsSing_HH_c==0, adj_se_KR_HH_c, NA)
  # adj_ddf_KR_HH_c <- dof_kenward(remlfit_HH_c)['treat']
  # adj_ddf_KR_NSing_HH_c <- ifelse(IsSing_HH_c==0, adj_ddf_KR_HH_c, NA)
  # }
  
  adj_se_KR_HH_l <- adj_ddf_KR_HH_l <- adj_se_KR_NSing_HH_l <-adj_ddf_KR_NSing_HH_l <- NA
  adj_se_Sat_HH_l <- adj_ddf_Sat_HH_l <- NA
  
  if (!is.null(remlfit_HH_l) & (S<10 | K<10)){
    ###pbkrtest (KR)
    vcov0_HH_l <- vcov(remlfit_HH_l)
    vcovKR_HH_l <- vcovAdj(remlfit_HH_l)
    adj_se_KR_HH_l <- sqrt(vcovKR_HH_l['treat','treat'])
    adj_se_KR_NSing_HH_l <- ifelse(IsSing_HH_l==0, adj_se_KR_HH_l, NA)
    L1 <- rep(0, length(fixef(remlfit_HH_l)))
    L1[which(names(fixef(remlfit_HH_l))=='treat')] <- 1
    adj_ddf_KR_HH_l <- Lb_ddf(L1, vcov0_HH_l, vcovKR_HH_l)
    adj_ddf_KR_NSing_HH_l <- ifelse(IsSing_HH_l==0, adj_ddf_KR_HH_l, NA)
    ###parameters (Sat)
    all_Sat_ses_HH_l <- se_satterthwaite(remlfit_HH_l)
    adj_se_Sat_HH_l <- all_Sat_ses_HH_l[which(all_Sat_ses_HH_l$Parameter=='treat'),'SE']
    adj_ddf_Sat_HH_l <- dof_satterthwaite(remlfit_HH_l)['treat']
    
  } 
  
  # if (!is.null(remlfit_HH_l) & ((S==4 & K==1) | (S==4 & K==5) | (S==10 & K==1))){
  #   ###parameters (KR)
  #   all_KR_ses_HH_l <- se_kenward(remlfit_HH_l)
  #   adj_se_KR_HH_l <- all_KR_ses_HH_l[which(all_KR_ses_HH_l$Parameter=='treat'),'SE']
  #   adj_se_KR_NSing_HH_l <- ifelse(IsSing_HH_l==0, adj_se_KR_HH_l, NA)
  #   adj_ddf_KR_HH_l <- dof_kenward(remlfit_HH_l)['treat']
  #   adj_ddf_KR_NSing_HH_l <- ifelse(IsSing_HH_l==0, adj_ddf_KR_HH_l, NA)  
  # }
  
  adj_se_KR_BE_c <- adj_ddf_KR_BE_c <- adj_se_KR_NSing_BE_c <- adj_ddf_KR_NSing_BE_c <-NA 
  adj_se_Sat_BE_c <- adj_ddf_Sat_BE_c <- NA
  
  if (!is.null(remlfit_BE_c) & (S<10 | K<10)){
    ###pbkrtest (KR)
    vcov0_BE_c<- vcov(remlfit_BE_c)
    vcovKR_BE_c <- vcovAdj(remlfit_BE_c)
    adj_se_KR_BE_c <- sqrt(vcovKR_BE_c['treat','treat'])
    adj_se_KR_NSing_BE_c <- ifelse(IsSing_BE_c==0, adj_se_KR_BE_c, NA)
    L2 <- rep(0, length(fixef(remlfit_BE_c)))
    L2[which(names(fixef(remlfit_BE_c))=='treat')] <- 1
    adj_ddf_KR_BE_c <- Lb_ddf(L2, vcov0_BE_c, vcovKR_BE_c)
    adj_ddf_KR_NSing_BE_c <- ifelse(IsSing_BE_c==0, adj_ddf_KR_BE_c, NA)
    ###parameters (Sat)
    all_Sat_ses_BE_c <- se_satterthwaite(remlfit_BE_c)
    adj_se_Sat_BE_c <- all_Sat_ses_BE_c[which(all_Sat_ses_BE_c$Parameter=='treat'),'SE']
    adj_ddf_Sat_BE_c <- dof_satterthwaite(remlfit_BE_c)['treat']
    
  } 
  
  # if (!is.null(remlfit_BE_c) & ((S==4 & K==1) | (S==4 & K==5) | (S==10 & K==1))){
  #   ###parameters (KR)
  #   all_KR_ses_BE_c <- se_kenward(remlfit_BE_c)
  #   adj_se_KR_BE_c <- all_KR_ses_BE_c[which(all_KR_ses_BE_c$Parameter=='treat'),'SE']
  #   adj_se_KR_NSing_BE_c <- ifelse(IsSing_BE_c==0, adj_se_KR_BE_c, NA)
  #   adj_ddf_KR_BE_c <- dof_kenward(remlfit_BE_c)['treat']
  #   adj_ddf_KR_NSing_BE_c <- ifelse(IsSing_BE_c==0, adj_ddf_KR_BE_c, NA)
  # }
  adj_se_KR_BE_l <- adj_ddf_KR_BE_l <- adj_se_KR_NSing_BE_l <- adj_ddf_KR_NSing_BE_l <-NA 
  adj_se_Sat_BE_l <- adj_ddf_Sat_BE_l <- NA
  
  if (!is.null(remlfit_BE_l) & (S<10 | K<10)){
    ###pbkrtest (KR)
    vcov0_BE_l<- vcov(remlfit_BE_l)
    vcovKR_BE_l <- vcovAdj(remlfit_BE_l)
    adj_se_KR_BE_l <- sqrt(vcovKR_BE_l['treat','treat'])
    adj_se_KR_NSing_BE_l <- ifelse(IsSing_BE_l==0, adj_se_KR_BE_l, NA)
    L2 <- rep(0, length(fixef(remlfit_BE_l)))
    L2[which(names(fixef(remlfit_BE_l))=='treat')] <- 1
    adj_ddf_KR_BE_l <- Lb_ddf(L2, vcov0_BE_l, vcovKR_BE_l)
    adj_ddf_KR_NSing_BE_l <- ifelse(IsSing_BE_l==0, adj_ddf_KR_BE_l, NA)

    ###parameters (Sat)
    all_Sat_ses_BE_l <- se_satterthwaite(remlfit_BE_l)
    adj_se_Sat_BE_l <- all_Sat_ses_BE_l[which(all_Sat_ses_BE_l$Parameter=='treat'),'SE']
    adj_ddf_Sat_BE_l <- dof_satterthwaite(remlfit_BE_l)['treat']
    
  } 
  
  # if (!is.null(remlfit_BE_l) & ((S==4 & K==1) | (S==4 & K==5) | (S==10 & K==1))){
  #   ###parameters (KR)
  #   all_KR_ses_BE_l <- se_kenward(remlfit_BE_l)
  #   adj_se_KR_BE_l <- all_KR_ses_BE_l[which(all_KR_ses_BE_l$Parameter=='treat'),'SE']
  #   adj_se_KR_NSing_BE_l <- ifelse(IsSing_BE_l==0, adj_se_KR_BE_l, NA)
  #   adj_ddf_KR_BE_l <- dof_kenward(remlfit_BE_l)['treat']
  #   adj_ddf_KR_NSing_BE_l <- ifelse(IsSing_BE_l==0, adj_ddf_KR_BE_l, NA)
  # }
  
  # str(dat)
  # 
  # dat$time <- as.integer(dat$time)
  # dat$cluster <- as.integer(dat$cluster)
  # 
  # mod <- Model$new(formula= Y ~ treat + time + (1 | gr(cluster)) + (1 | gr(cluster,time)),
  #                  data = dat,
  #                  family = gaussian())
  # 
  # mod$small_sample_correction("KR")
  # #improved correction given in Kenward & Roger (2009)
  # mod$small_sample_correction("KR2")
  # mod$small_sample_correction("sat")
  # kr <- mod$small_sample_correction(type="KR")
  # sqrt(diag(kr$vcov_beta))
  # kr$dof #degrees of freedom
  # 
  # mod <- Model$new(formula= Y ~ treat + time + (1 | gr(cluster)),
  #                  data = dat,
  #                  family = gaussian())
  # 
  # mod$small_sample_correction("KR")
  # #improved correction given in Kenward & Roger (2009)
  # mod$small_sample_correction("KR2")
  # mod$small_sample_correction("sat")
  # kr <- mod$small_sample_correction(type="KR")
  # sqrt(diag(kr$vcov_beta))
  # kr$dof #degrees of freedom
  
  #sqrt(mod$small_sample_correction("KR")$vcov_theta[2,2])
  
  w_conv_HH_c <- fitHHmodelSC_cdat$w_convergence_HH_c
  w_other_HH_c <-   fitHHmodelSC_cdat$w_other_HH_c
  err_HH_c  <-   fitHHmodelSC_cdat$e_all_HH_c
  
  w_conv_HH_l <- fitHHmodelSC_ldat$w_convergence_HH_l
  w_other_HH_l <-   fitHHmodelSC_ldat$w_other_HH_l
  err_HH_l  <-   fitHHmodelSC_ldat$e_all_HH_l
  
  w_conv_BE_c <- fitBEmodelSC_cdat$w_convergence_BE_c
  w_other_BE_c <-   fitBEmodelSC_cdat$w_other_BE_c
  err_BE_c  <-   fitBEmodelSC_cdat$e_all_BE_c
  
  w_conv_BE_l <- fitBEmodelSC_ldat$w_convergence_BE_l
  w_other_BE_l <-   fitBEmodelSC_ldat$w_other_BE_l
  err_BE_l  <-   fitBEmodelSC_ldat$e_all_BE_l
  
  return(c(est_trt_HH_c, se_trt_HH_c,est_trt_BE_c, se_trt_BE_c,
           est_ICC_HH_c, est_ICC_BE_c, est_CAC_BE_c,
           adj_se_KR_HH_c,adj_ddf_KR_HH_c,adj_se_KR_NSing_HH_c, adj_ddf_KR_NSing_HH_c,
           adj_se_KR_BE_c,adj_ddf_KR_BE_c,adj_se_KR_NSing_BE_c, adj_ddf_KR_NSing_BE_c,
           adj_se_Sat_HH_c,adj_ddf_Sat_HH_c,adj_se_Sat_BE_c, adj_ddf_Sat_BE_c,
           w_conv_HH_c,w_conv_BE_c,w_other_HH_c,w_other_BE_c,
           IsSing_HH_c,IsSing_BE_c,err_HH_c,err_BE_c,
           
           est_trt_HH_l, se_trt_HH_l,est_trt_BE_l, se_trt_BE_l,
           est_ICC_HH_l, est_ICC_BE_l, est_CAC_BE_l,
           adj_se_KR_HH_l,adj_ddf_KR_HH_l,adj_se_KR_NSing_HH_l, adj_ddf_KR_NSing_HH_l,
           #adj_se_KR2_HH_l,adj_ddf_KR2_HH_l,adj_se_KR2_NSing_HH_l, adj_ddf_KR2_NSing_HH_l,
           adj_se_KR_BE_l,adj_ddf_KR_BE_l,adj_se_KR_NSing_BE_l, adj_ddf_KR_NSing_BE_l,
           #adj_se_KR2_BE_l,adj_ddf_KR2_BE_l,adj_se_KR2_NSing_BE_l, adj_ddf_KR2_NSing_BE_l,
           adj_se_Sat_HH_l,adj_ddf_Sat_HH_l,adj_se_Sat_BE_l, adj_ddf_Sat_BE_l,
           w_conv_HH_l,w_conv_BE_l,w_other_HH_l,w_other_BE_l,
           IsSing_HH_l,IsSing_BE_l,err_HH_l,err_BE_l))
}

# Get estimates for nsim repeates
sim_res_fit <- function(nsim, S, K, m, ICC, CAC, theta){
  # Calculates empirical power based on nsim simulated trial datasets
  # Generate trial dataset, fit both models, calculate rejection probability
  res_fit_mat <- replicate(nsim, fitmodels(S, K, m, ICC, CAC, theta))
  res_fit_matx_t <-  t(res_fit_mat)
  #create a data frame convert res matrix to a data frame
  res_fit <- as.data.frame(res_fit_matx_t)
  colnames(res_fit) <- c("est_trt_HH_c", "se_trt_HH_c", "est_trt_BE_c", "se_trt_BE_c", 
                         "est_ICC_HH_c", "est_ICC_BE_c", "est_CAC_BE_c", 
                         "adj_se_KR_HH_c","adj_ddf_KR_HH_c","adj_se_KR_NSing_HH_c", "adj_ddf_KR_NSing_HH_c",
                         "adj_se_KR_BE_c","adj_ddf_KR_BE_c","adj_se_KR_NSing_BE_c", "adj_ddf_KR_NSing_BE_c",
                         "adj_se_Sat_HH_c","adj_ddf_Sat_HH_c","adj_se_Sat_BE_c", "adj_ddf_Sat_BE_c",
                         "w_conv_HH_c","w_conv_BE_c","w_other_HH_c","w_other_BE_c",
                         "IsSing_HH_c","IsSing_BE_c","err_HH_c","err_BE_c",
                         
                         "est_trt_HH_l", "se_trt_HH_l", "est_trt_BE_l", "se_trt_BE_l", 
                         "est_ICC_HH_l", "est_ICC_BE_l", "est_CAC_BE_l", 
                         "adj_se_KR_HH_l","adj_ddf_KR_HH_l","adj_se_KR_NSing_HH_l", "adj_ddf_KR_NSing_HH_l",
                         "adj_se_KR_BE_l","adj_ddf_KR_BE_l","adj_se_KR_NSing_BE_l", "adj_ddf_KR_NSing_BE_l",
                         "adj_se_Sat_HH_l","adj_ddf_Sat_HH_l","adj_se_Sat_BE_l", "adj_ddf_Sat_BE_l",
                         "w_conv_HH_l","w_conv_BE_l","w_other_HH_l","w_other_BE_l",
                         "IsSing_HH_l","IsSing_BE_l","err_HH_l","err_BE_l")
  
  res_est_fit <- res_fit %>%
    mutate(rep=1:nsim, nsim, S, K, m, ICC, CAC, theta)
  # Reorder columns
  res_est_fit <- res_est_fit %>% 
    relocate (rep:theta)
  # Save outputs
  save(res_est_fit, file = paste0("est_files_v2\\estimates",
                                  "_nsim", nsim, "_S", S, "_K", K, "_m", m, "_ICC", ICC, "_CAC", CAC,
                                    "_theta", theta, ".RData"))
  
  return(res_est_fit)
  
}





