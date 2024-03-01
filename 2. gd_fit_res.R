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
  # dat <- gen_dat(4, 1, 10, 0.01,1,0.15)
  # typ='lin'
  # typ='cat'
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
  #se_ICC_HH <- 
  
  fitBEmodelSC_dat <- fitBEmodelSC(dat,typ)
  remlfit_BE    <- fitBEmodelSC_dat[[1]]
  est_trt_BE    <- ifelse(is.null(remlfit_BE), NA, fixef(remlfit_BE)['treat'])
  se_trt_BE     <- ifelse(is.null(remlfit_BE), NA, sqrt(vcov(remlfit_BE)['treat','treat']))
  IsSing_BE     <- ifelse(isSingular(remlfit_BE), 1, 0) 
  
  estvarclustr_BE   <- ifelse(is.null(remlfit_BE), NA, VarCorr(remlfit_BE)$cluster[1])
  estvarclustper_BE <- ifelse(is.null(remlfit_BE), NA, VarCorr(remlfit_BE)$clustper[1]) 
  estvarres_BE      <- ifelse(is.null(remlfit_BE), NA, sigma(remlfit_BE)^2)
  est_ICC_BE <- (estvarclustr_BE+estvarclustper_BE)/(estvarclustr_BE+estvarclustper_BE+estvarres_BE)
  se_ICC_BE <-
  est_CAC_BE <- (estvarclustr_BE)/(estvarclustr_BE+estvarclustper_BE)
  #se_CAC_BE <- 
  
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
#a single simulation replicates fitmodels nsim times

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
  
  #Calculate rejection proportion for empirical power for HH
  #the proportion of times the test rejects the null hypothesis in the study.
  nsim_HH <- NA
  nsim_HH  <- sum(!is.na(res_fit$est_trt_HH) & !is.na(res_fit$se_trt_HH))
  pwr_HH   <- sum(abs(res_fit$est_trt_HH)/res_fit$se_trt_HH > 1.96, na.rm=TRUE)/nsim_HH
  varest_trt_HH<- var(res_fit$est_trt_HH, na.rm=TRUE)
  mest_trt_HH  <- mean(res_fit$est_trt_HH, na.rm=TRUE)
  mest_ICC_HH <- mean(res_fit$est_ICC_HH, na.rm=TRUE)
  varICC_HH  <- var(res_fit$est_ICC_HH, na.rm=TRUE)
  
  #Coverage for HH model without any adjustment
  num_cov_HH <- NA
  alpha <- (1 + 0.95)/2
  zstat_HH <- qnorm(alpha)
  #count the number of intervals that contain the true values.
  num_cov_HH <- sum(res_fit$est_trt_HH-(zstat_HH*res_fit$se_trt_HH) <= theta 
                    & res_fit$est_trt_HH+(zstat_HH*res_fit$se_trt_HH) >= theta,na.rm=T) 
  
  p_cov_HH <- ifelse(is.na(num_cov_HH/nsim_HH), NA, num_cov_HH/nsim_HH)
  
  #Calculate empirical power for HH with Kenward-Roger correction
  nsim_KR_HH<- NA
  tstat_KR_HH <- qt(alpha, res_fit$adj_ddf_KR_HH)
  nsim_KR_HH <- sum(!is.na(res_fit$est_trt_HH) & !is.na(res_fit$adj_se_KR_HH))
  pwr_KR_HH  <- sum(abs(res_fit$est_trt_HH)/res_fit$adj_se_KR_HH > tstat_KR_HH, na.rm=TRUE)/nsim_KR_HH
  
  #Calculate Coverage for HH with Kenward-Roger correction 
  num_KR_cov_HH<- NA
  #count the number of intervals that contain the true values.
  num_KR_cov_HH <- sum(res_fit$est_trt_HH-zstat_HH*(res_fit$adj_se_KR_HH) <= theta
                       & res_fit$est_trt_HH+zstat_HH*(res_fit$adj_se_KR_HH) >= theta,na.rm=T)
  p_KR_cov_HH <- ifelse(is.na(num_KR_cov_HH/nsim_KR_HH), NA, num_KR_cov_HH/nsim_KR_HH)
  
  #Calculate empirical power for HH with Kenward-Roger correction & non-singular fit
  nsim_KR_NSing_HH<- NA
  tstat_KR_NSing_HH <- qt(alpha, res_fit$adj_ddf_KR_NSing_HH)
  nsim_KR_NSing_HH <- sum(!is.na(res_fit$est_trt_HH) & !is.na(res_fit$adj_se_KR_NSing_HH))
  pwr_KR_NSing_HH  <- sum(abs(res_fit$est_trt_HH)/res_fit$adj_se_KR_NSing_HH > tstat_KR_NSing_HH, na.rm=TRUE)/nsim_KR_NSing_HH
  
  #Calculate Coverage for HH with Kenward-Roger correction & Non singular fit
  num_KR_NSing_cov_HH<- NA
  #count the number of intervals that contain the true values.
  num_KR_NSing_cov_HH <- sum(res_fit$est_trt_HH-zstat_HH*(res_fit$adj_se_KR_NSing_HH) <= theta
                          & res_fit$est_trt_HH+zstat_HH*(res_fit$adj_se_KR_NSing_HH) >= theta,na.rm=T)
  p_cov_KR_NSing_HH <- ifelse(is.na(num_KR_NSing_cov_HH/nsim_KR_NSing_HH), NA, num_KR_NSing_cov_HH/nsim_KR_NSing_HH)
  
  #Calculate empirical power for HH with Satterthwaite  correction
  nsim_Sat_HH<- NA
  tstat_Sat_HH <- qt(alpha, res_fit$adj_ddf_Sat_HH)
  nsim_Sat_HH  <- sum(!is.na(res_fit$est_trt_HH) & !is.na(res_fit$adj_se_Sat_HH))
  pwr_Sat_HH <- sum(abs(res_fit$est_trt_HH)/res_fit$adj_se_Sat_HH > tstat_Sat_HH,na.rm=TRUE)/nsim_Sat_HH
  
  #Calculate Coverage for HH with Satterthwaite  correction
  num_Sat_cov_HH<- NA
  #count the number of intervals that contain the true values.
  num_Sat_cov_HH <- sum(res_fit$est_trt_HH-zstat_HH*(res_fit$adj_se_Sat_HH) <= theta
                        & res_fit$est_trt_HH+zstat_HH*(res_fit$adj_se_Sat_HH) >= theta,na.rm=T)
  p_cov_Sat_HH <- ifelse(is.na(num_Sat_cov_HH/nsim_Sat_HH), NA, num_Sat_cov_HH/nsim_Sat_HH)
  
  #Calculate rejection proportion for empirical power for BE
  #the proportion of times the test rejects the null hypothesis in the study.
  nsim_BE<- NA
  nsim_BE  <-  sum(!is.na(res_fit$est_trt_BE) & !is.na(res_fit$se_trt_BE))
  pwr_BE   <- sum(abs(res_fit$est_trt_BE)/res_fit$se_trt_BE > 1.96, na.rm=TRUE)/nsim_BE
  varest_trt_BE<- var(res_fit$est_trt_BE, na.rm=TRUE)
  mest_trt_BE  <- mean(res_fit$est_trt_BE, na.rm=TRUE)
  mest_ICC_BE <- mean(res_fit$est_ICC_BE, na.rm=TRUE)
  varICC_BE  <- var(res_fit$est_ICC_BE, na.rm=TRUE)
  mest_CAC_BE <-  mean(res_fit$est_CAC_BE, na.rm=TRUE)
  var_CAC_BE  <- var(res_fit$est_CAC_BE, na.rm=TRUE)
  
  #Coverage for BE model without any adjustment
  num_cov_BE<- NA
  zstat_BE <- qnorm(alpha)
  #count the number of intervals that contain the true values.
  num_cov_BE <- sum(res_fit$est_trt_BE-zstat_BE*(res_fit$se_trt_BE) <= theta
                    & res_fit$est_trt_BE+zstat_BE*(res_fit$se_trt_BE) >= theta,na.rm=T)
  p_cov_BE <- ifelse(is.na(num_cov_BE/nsim_BE), NA, num_cov_BE/nsim_BE)
  
  
  #Calculate empirical power for BE with Kenward-Roger correction
  nsim_KR_BE<- NA
  tstat_KR_BE <- qt(alpha, res_fit$adj_ddf_KR_BE)
  nsim_KR_BE  <- sum(!is.na(res_fit$est_trt_BE) & !is.na(res_fit$adj_se_KR_BE))
  pwr_KR_BE <- sum(abs(res_fit$est_trt_BE)/res_fit$adj_se_KR_BE> tstat_KR_BE, na.rm=TRUE)/nsim_KR_BE
  
  
  #Calculate Coverage for BE with Kenward-Roger correction 
  num_KR_cov_BE<- NA
  #count the number of intervals that contain the true values.
  num_KR_cov_BE <- sum(res_fit$est_trt_BE-zstat_BE*(res_fit$adj_se_KR_BE) <= theta
                         & res_fit$est_trt_BE+zstat_BE*(res_fit$adj_se_KR_BE) >= theta,na.rm=T)
  p_KR_cov_BE <- ifelse(is.na(num_KR_cov_BE/nsim_KR_BE), NA, num_KR_cov_BE/nsim_KR_BE)
  
  #Calculate empirical power for BE with Kenward-Roger correction & Non singular fit 
  nsim_KR_NSing_BE<- NA
  tstat_KR_NSing_BE <- qt(alpha, res_fit$adj_ddf_KR_NSing_BE)
  nsim_KR_NSing_BE  <- sum(!is.na(res_fit$est_trt_BE) & !is.na(res_fit$adj_se_KR_NSing_BE))
  pwr_KR_NSing_BE <- sum(abs(res_fit$est_trt_BE)/res_fit$adj_se_KR_NSing_BE> tstat_KR_NSing_BE, na.rm=TRUE)/nsim_KR_NSing_BE 
  
  #Calculate Coverage for BE with Kenward-Roger correction & Non singular fit
  num_KR_NSing_cov_BE<- NA
  #count the number of intervals that contain the true values.
  num_KR_NSing_cov_BE <- sum(res_fit$est_trt_BE-zstat_BE*(res_fit$adj_se_KR_NSing_BE) <= theta
                    & res_fit$est_trt_BE+zstat_BE*(res_fit$adj_se_KR_NSing_BE) >= theta,na.rm=T)
  p_cov_KR_NSing_BE <- ifelse(is.na(num_KR_NSing_cov_BE/nsim_KR_NSing_BE), NA, num_KR_NSing_cov_BE/nsim_KR_NSing_BE)
  
  
  #Calculate empirical power for BE with Satterthwaite  correction
  nsim_Sat_BE<- NA
  tstat_Sat_BE <- qt(alpha, res_fit$adj_ddf_Sat_BE)
  nsim_Sat_BE  <- sum(!is.na(res_fit$est_trt_BE) & !is.na(res_fit$adj_se_Sat_BE))
  pwr_Sat_BE <- sum(abs(res_fit$est_trt_BE)/res_fit$adj_se_Sat_BE > tstat_Sat_BE,na.rm=TRUE)/nsim_Sat_BE
  
  
  #Calculate Coverage for BE with Satterthwaite  correction
  num_Sat_cov_BE<- NA
  #count the number of intervals that contain the true values.
  num_Sat_cov_BE <- sum(res_fit$est_trt_BE-zstat_BE*(res_fit$adj_se_Sat_BE) <= theta
                         & res_fit$est_trt_BE+zstat_BE*(res_fit$adj_se_Sat_BE) >= theta,na.rm=T)
  p_cov_Sat_BE <- ifelse(is.na(num_Sat_cov_BE/nsim_Sat_BE), NA, num_Sat_cov_BE/nsim_Sat_BE)
  
  
  #Type I error rate
  #error of concluding that there is a significant effect or difference
  #when there is, in fact, no such effect or difference.
  
  
  if (typ=='cat'){
    power  <- pow(VarSCcat(S,K,m, ICC, CAC),theta)
  } else if (typ=='lin') {
    power  <- pow(VarSClin(S,K,m, ICC, CAC),theta)
  }
  
  s_IsSing_HH <- sum(res_fit$IsSing_HH,na.rm=T)
  s_IsSing_BE <- sum(res_fit$IsSing_BE,na.rm=T)
  s_w_conv_HH <- sum(res_fit$w_conv_HH,na.rm=T)
  s_w_conv_BE <- sum(res_fit$w_conv_BE,na.rm=T)
  s_w_other_HH <- sum(res_fit$w_other_HH,na.rm=T)
  s_w_other_BE <- sum(res_fit$w_other_BE,na.rm=T)
  s_err_HH <- sum(res_fit$err_HH,na.rm=T)
  s_err_BE <- sum(res_fit$err_BE,na.rm=T)
  
  par_emp_vals <- data.frame(nsim, S, K, m, ICC, CAC, theta,type=typ,power,
           nsim_HH=nsim_HH, 
           mesttrt_HH=mest_trt_HH, varest_trt_HH=varest_trt_HH,
           nsim_BE=nsim_BE, 
           mesttrt_BE=mest_trt_BE,varest_trt_BE=varest_trt_BE,
           mestICC_HH=mest_ICC_HH,varICC_HH=varICC_HH, 
           mestICC_BE=mest_ICC_BE,varICC_BE=varICC_BE,
           mestCAC_BE=mest_CAC_BE,varCAC_BE=var_CAC_BE,
           sIsSing_HH=s_IsSing_HH,pow_HH=pwr_HH,
           powKR_HH=pwr_KR_HH, powKRnSing_HH=pwr_KR_NSing_HH,
           powSat_HH=pwr_Sat_HH,
           pcov_HH=p_cov_HH,pKRcov_HH=p_KR_cov_HH,
           pKRnSingcov_HH=p_cov_KR_NSing_HH,pcovSat_HH=p_cov_Sat_HH,
           sIsSing_BE=s_IsSing_BE,pow_BE=pwr_BE,
           powKR_BE=pwr_KR_BE,powKRnSing_BE=pwr_KR_NSing_BE,
           powSat_BE=pwr_Sat_BE,
           pcov_BE=p_cov_BE,pKRcov_BE=p_KR_cov_BE,
           pKRnSingcov_BE=p_cov_KR_NSing_BE,pSatcov_BE=p_cov_Sat_BE,
           swarconv_HH=s_w_conv_HH,swarconv_BE=s_w_conv_BE,
           swarother_HH=s_w_other_HH,swarother_BE=s_w_other_BE,
           serr_HH=s_err_HH,serr_BE=s_err_BE)
  

  
  res_est_fit <- res_fit %>%
    mutate(rep=1:nsim, nsim, S, K, m, ICC, CAC, theta,type=typ)
  
  # Reorder columns
  res_est_fit <- res_est_fit %>% 
    relocate (rep:type)
  
  return(list(par_emp_vals,res_est_fit))
}





