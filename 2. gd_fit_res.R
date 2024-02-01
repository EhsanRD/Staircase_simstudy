# Calculate empirical power for staircase designs
# Empirical power was measured as the proportion of trials for which a test would have identified
# a significant result (p value < .05) doi: 10.1177/1740774520940256.
# Ehsan Rezaei (Ehsan.rezaeidarzi@monash.edu)

library(MASS)
library(lme4)
library(glmmTMB)
library(pbkrtest)
library(tidyr)


setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Staircase_simstudy\\")

source('1. functions_sim.R')

gen_data <- function(S, K, m, ICC, CAC, theta){
  # Generates a single dataset from a staircase design
  # and trial configuration,for a block-exchangeable within-cluster  
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
  #inactive zero treatment effect statements and delete them later after Jess's thoughts on this
  #Y0 <-  betavec + Xvec*0 + Cvec + CPvec + e
  
  # Create data frame with everything needed for fitting the model
  dat <- data.frame(Y=Y, cluster=as.factor(clusterind), time=perind, clustper=as.factor(clustperind), treat=Xvec)
  #dat <- data.frame(Y=Y, Y_0=Y0, cluster=as.factor(clusterind), time=perind, clustper=as.factor(clustperind), treat=Xvec)
  # Get subset of data corresponding to embedded basic staircase design

  return(dat)
}

fitBEmodelSC <- function(dat){
  tryCatch(
    expr = {
      remlfit <- lmer(Y ~ treat + as.factor(time) + (1|cluster) + (1|clustper), data=dat, REML=TRUE)
    },
    error = function(e){
      message('Caught an error!')
      message(e)
      
      return(NULL)
    }
  )    
}


# fitBEmodelSC0 <- function(dat){
#   tryCatch(
#     expr = {
#       remlfit0 <- lmer(Y0 ~ treat + as.factor(time) + (1|cluster) + (1|clustper), data=dat, REML=TRUE)
#     },
#     error = function(e){
#       message('Caught an error!')
#       message(e)
#       
#       return(NULL)
#     }
#   )    
# }

# getSE <- function(remlfit_BE){
#   tryCatch(
#     expr = {
#       sqrt(vcov(remlfit_BE)[[1]]['treat','treat'])
#     },
#     error = function(e){
#       message('Caught an error!')
#       message(e)
#       
#       return(NA)
#     }
#   )
# }

fitmodels <- function(S, K, m, ICC, CAC, theta){
  
  #lmer for fitting linear mixed-effects models (LMMs)
  #glmmTMB for fitting generalized linear mixed-effects models (GLMMs)
  
  # Generates a single simulated trial dataset, fits corresponding model and
  # outputs the treatment effect estimate and standard error
  
  # Generate dataset
  dat <- gen_data(S, K, m, ICC, CAC, theta)
  #dat <- gen_data(10, 1, 100, 0.01, 0.95,0.15)
  
  # Fit both models
  remlfit_HH   <- lmer(Y ~ treat + as.factor(time) + (1|cluster), data=dat, REML=TRUE)
  est_trt_HH   <- fixef(remlfit_HH)['treat']
  se_trt_HH    <- sqrt(vcov(remlfit_HH)['treat','treat'])
  # Extract the number of parameters
  #num_par_HH <- length(unlist(lme4::fixef(remlfit_HH))) + length(unlist(lme4::ranef(remlfit_HH)))
  
  estvarclustr_HH <- VarCorr(remlfit_HH)$cluster[1]
  estvarres_HH     <- sigma(remlfit_HH)^2
  est_ICC_HH  <- (estvarclustr_HH)/(estvarclustr_HH+estvarres_HH)

  # remlfit0_HH  <- lmer(Y0 ~ treat +  as.factor(time)  + (1|cluster), data=dat, REML=TRUE)
  # est0_HH  <- fixef(remlfit0_HH)['treat']
  # se0_HH   <- sqrt(vcov(remlfit0_HH)['treat','treat'])
  # 
  # estvarclustr0_HH   <- VarCorr(remlfit0_HH)$cluster[1]
  # estvarres0_HH      <-  sigma(remlfit0_HH)^2
  # estICC0_HH <- (estvarclustr0_HH)/(estvarclustr0_HH+estvarres0_HH)

  remlfit_BE    <- fitBEmodelSC(dat)
  est_trt_BE    <- ifelse(is.null(remlfit_BE), NA, fixef(remlfit_BE)['treat'])
  se_trt_BE     <- ifelse(is.null(remlfit_BE), NA, sqrt(vcov(remlfit_BE)['treat','treat']))

  estvarclustr_BE   <- ifelse(is.null(remlfit_BE), NA, VarCorr(remlfit_BE)$cluster[1])
  estvarclustper_BE <- ifelse(is.null(remlfit_BE), NA, VarCorr(remlfit_BE)$clustper[1]) 
  estvarres_BE      <- ifelse(is.null(remlfit_BE), NA, sigma(remlfit_BE)^2)
  est_ICC_BE <- (estvarclustr_BE+estvarclustper_BE)/(estvarclustr_BE+estvarclustper_BE+estvarres_BE)
  est_CAC_BE <- (estvarclustr_BE)/(estvarclustr_BE+estvarclustper_BE)

  # #glmfit0      <- glmmTMB(Y0 ~ treat + time + (1|cluster) + (1|clustper), data=dat, REML=TRUE)
  # remlfit0_BE    <- fitBEmodelSC0(dat)
  # #estSC_BE   <- ifelse(is.null(remlfit_BE), NA, fixef(remlfit_BE)[[1]]['treat'])
  # est0_BE    <- ifelse(is.null(remlfit0_BE), NA, fixef(remlfit0_BE)['treat'])
  # #seSC_BE    <- ifelse(is.null(remlfit_BE), NA, getSE(remlfit_BE))
  # se0_BE     <- ifelse(is.null(remlfit0_BE), NA, sqrt(vcov(remlfit0_BE)['treat','treat']))
  # 
  # #estvarclustr_BE   <- ifelse(is.null(remlfit_BE), NA, VarCorr(remlfit_BE)[[c("cond","cluster")]][1]) 
  # estvarclustr0_BE   <- ifelse(is.null(remlfit0_BE), NA, VarCorr(remlfit0_BE)$cluster[1])
  # #estvarclustper_BE <- ifelse(is.null(remlfit_BE), NA, VarCorr(remlfit_BE)[[c("cond","clustper")]][1]) 
  # estvarclustper0_BE <- ifelse(is.null(remlfit0_BE), NA, VarCorr(remlfit0_BE)$clustper[1]) 
  # #estvarres_BE      <- ifelse(is.null(remlfit_BE), NA, attr(VarCorr(remlfit_BE)$cond, "sc")^2)
  # estvarres0_BE      <- ifelse(is.null(remlfit0_BE), NA, sigma(remlfit0_BE)^2)
  # estICC0_BE <- (estvarclustr0_BE+estvarclustper0_BE)/(estvarclustr0_BE+estvarclustper0_BE+estvarres0_BE)
  # estCAC0_BE <- (estvarclustr0_BE)/(estvarclustr0_BE+estvarclustper0_BE)

  # Get adjusted SE & adjusted degree of freedom with Kenward-Roger correction
  # Note: This only works for exchangeable models (lmer model fit objects).
   vcov0_HH <- vcov(remlfit_HH)
   vcovKR_HH <- vcovAdj(remlfit_HH)
   adj_se_HH <- sqrt(vcovKR_HH['treat','treat'])
   L1 <- rep(0, length(fixef(remlfit_HH)))
   L1[which(names(fixef(remlfit_HH))=='treat')] <- 1
   adj_ddf_HH <- Lb_ddf(L1, vcov0_HH, vcovKR_HH)
   
   vcov0_BE<- vcov(remlfit_BE)
   vcovKR_BE <- vcovAdj(remlfit_BE)
   adj_se_BE <- sqrt(vcovKR_BE['treat','treat'])
   L2 <- rep(0, length(fixef(remlfit_BE)))
   L2[which(names(fixef(remlfit_BE))=='treat')] <- 1
   adj_ddf_BE <- Lb_ddf(L2, vcov0_BE, vcovKR_BE)
   
  return(c(est_trt_HH, se_trt_HH,est_trt_BE, se_trt_BE, 
           est_ICC_HH, est_ICC_BE, est_CAC_BE, 
           adj_se_HH,adj_ddf_HH,adj_se_BE,adj_ddf_BE))
}

sim_res_fit <- function(nsim, S, K, m, ICC, CAC, theta){
  # Calculates empirical power based on nsim simulated trial datasets
  
  # Generate trial dataset, fit both models, calculate rejection probability
  res_fit_matx <- replicate(nsim, fitmodels(S, K, m, ICC, CAC, theta))
  #res_fit_matx <- replicate(5, fitmodels(10, 1, 100, 0.01, 0.95, 0.15))
  res_fit_matx_t <-  t(res_fit_matx)
  
  #create a data frame convert res matrix to a data frame
  res_fit <- as.data.frame(res_fit_matx_t)
  colnames(res_fit) <- c("est_trt_HH", "se_trt_HH", "est_trt_BE", "se_trt_BE", 
                        "est_ICC_HH", "est_ICC_BE", "est_CAC_BE", 
                        "adj_se_HH", "adj_ddf_HH", "adj_se_BE", "adj_ddf_BE")

  
  #Calculate rejection proportion for empirical power
  #the proportion of times the test rejects the null hypothesis in the study.
  
  nsim_HH  <- sum(!is.na(abs(res_fit$est_trt_HH)/res_fit$se_trt_HH > 1.96))
  pwr_HH   <- sum(abs(res_fit$est_trt_HH)/res_fit$se_trt_HH > 1.96, na.rm=TRUE)/nsim_HH
  varest_trt_HH<- var(res_fit$est_trt_HH, na.rm=TRUE)
  mest_trt_HH  <- mean(res_fit$est_trt_HH, na.rm=TRUE)
  mest_ICC_HH <- mean(res_fit$est_ICC_HH, na.rm=TRUE)
  varICC_HH  <- var(res_fit$est_ICC_HH, na.rm=TRUE)
  
  #Coverage 
  alpha <- (1 + 0.95)/2
  #dof: number of clusters minus the number of regression parameters (Fan Li 2020)? K*S-(T+1) (Fan Li 2018) ? or m*S*K*T-S*K-T (Grayling et al 2017)?
  #tstat_BE <- qt(alpha, dof)
  zstat_HH <- qnorm(alpha)
  #count the number of intervals that contain the true values.
  num_cov_HH <- sum(res_fit$est_trt_HH-(zstat_HH*res_fit$se_trt_HH) <= theta 
                  & res_fit$est_trt_HH+(zstat_HH*res_fit$se_trt_HH) >= theta) 
  p_cov_HH <- num_cov_HH/nsim_HH
  
  #Get t test statistic with KR-adjusted ddf
  tstat_KR_HH <- qt(alpha, res_fit$adj_ddf_HH)
  nsim_KR_HH  <- sum(!is.na(abs(res_fit$est_trt_HH)/res_fit$adj_se_HH > tstat_KR_HH))
  pwr_KR_HH <- sum(abs(res_fit$est_trt_HH)/res_fit$adj_se_HH > tstat_KR_HH)/nsim_KR_HH
  
  # # Construct 95% KR confidence intervals
  # theta_CI_KR_low  <- REML_theta_est - tstatSC * adj_SE
  # theta_CI_KR_high <- REML_theta_est + tstatSC * adj_SE
  
  #Type I error rate
  #error of concluding that there is a significant effect or difference
  #when there is, in fact, no such effect or difference.
  
  # nsim0_HH  <- sum(!is.na(abs(res[3,])/res[4,] > 1.96))
  # typeI_HH  <- sum(abs(res[3,])/res[4,] > 1.96, na.rm=TRUE)/nsim0_HH
  # varest0_HH<- var(res[3,], na.rm=TRUE)
  # mest0_HH  <- mean(res[3,], na.rm=TRUE)
  # mestICC0_HH <- mean(res[10,], na.rm=TRUE)
  # varICC0_HH  <- var(res[10,], na.rm=TRUE)
  
  nsim_BE  <- sum(!is.na(abs(res_fit$est_trt_BE)/res_fit$se_trt_BE > 1.96))
  pwr_BE   <- sum(abs(res_fit$est_trt_BE)/res_fit$se_trt_BE > 1.96, na.rm=TRUE)/nsim_BE
  varest_trt_BE<- var(res_fit$est_trt_BE, na.rm=TRUE)
  mest_trt_BE  <- mean(res_fit$est_trt_BE, na.rm=TRUE)
  mest_ICC_BE <- mean(res_fit$est_ICC_BE, na.rm=TRUE)
  varICC_BE  <- var(res_fit$est_ICC_BE, na.rm=TRUE)
  mest_CAC_BE <-  mean(res_fit$est_CAC_BE, na.rm=TRUE)
  varCAC_BE  <- var(res_fit$est_CAC_BE, na.rm=TRUE)
  
  #Coverage 
  #tstat_BE <- qt(alpha, dof)
  zstat_BE <- qnorm(alpha)
  #count the number of intervals that contain the true values.
  num_cov_BE <- sum(res_fit$est_trt_HH-zstat_BE*(res_fit$se_trt_BE) <= theta
                  & res_fit$est_trt_HH+zstat_BE*(res_fit$se_trt_BE) >= theta) 
  p_cov_BE <- num_cov_BE/nsim_BE
  
  # Calculate empirical power with Kenward-Roger correction
  tstat_KR_BE <- qt(alpha, res_fit$adj_ddf_BE)
  nsim_KR_BE  <- sum(!is.na(abs(res_fit$est_trt_BE)/res_fit$adj_se_BE > tstat_KR_BE))
  pwr_KR_BE <- sum(abs(res_fit$est_trt_BE)/res_fit$adj_se_BE> tstat_KR_BE)/nsim_KR_BE 
  
  #Type I error rate
  #error of concluding that there is a significant effect or difference
  #when there is, in fact, no such effect or difference.
  
  # nsim0_BE <- sum(!is.na(abs(res[7,])/res[8,] > 1.96))
  # typeI_BE <- sum(abs(res[7,])/res[8,] > 1.96, na.rm=TRUE)/nsim0_BE
  # est0_BE  <- var(res[7,], na.rm=TRUE)
  # mest0_BE <- mean(res[7,], na.rm=TRUE)
  # mestICC0_BE<- mean(res[12,], na.rm=TRUE)
  # mestCAC0_BE<-  mean(res[14,], na.rm=TRUE)
  # varICC0_BE <- var(res[12,], na.rm=TRUE)
  # mestCAC0_BE<-  mean(res[14,], na.rm=TRUE)
  # varCAC0_BE <- var(res[14,], na.rm=TRUE)
  
  power  <- pow(VarSCcat(S,K,m, ICC, CAC),theta)
  #power  <- pow(VarSClin(S,K,m, ICC, CAC),theta)
  
  par_emp_vals <- data.frame(nsim, S, K, m, ICC, CAC, theta,power,
                        nsim_HH=nsim_HH, power_HH=pwr_HH,
                        var_est_trt_HH=varest_trt_HH, m_est_trt_HH=mest_trt_HH,
                        nsim_BE=nsim_BE, power_BE=pwr_BE,
                        var_est_trt_BE=varest_trt_BE,m_est_trt_BE=mest_trt_BE,
                        m_est_ICC_HH=mest_ICC_HH, 
                        var_ICC_HH=varICC_HH, m_est_ICC_BE=mest_ICC_BE,
                        var_ICC_BE=varICC_BE,m_est_CAC_BE=mest_CAC_BE,
                        var_CAC_BE=varCAC_BE,power_KR_HH=pwr_KR_HH, 
                        power_KR_BE=pwr_KR_BE,p_cov_HH=p_cov_HH,
                        p_cov_BE=p_cov_BE)
  
  #rename empirical power to type I error if there is no treatment effect!
  # if (theta == 0) {
  #   par_emp_vals <- par_emp_vals %>%
  #     rename(typeI_BE = power_BE, typeI_HH = power_HH,
  #            typeI_KR_BE =power_KR_BE, typeI_KR_HH = power_KR_HH)
  # }
  # 
  #if theta==0 the empirical power calculation corresponds to type I error.
  

  return(par_emp_vals)
}




