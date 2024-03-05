#generate summary results
summary_res_fit <- function(nsim, S, K, m, ICC, CAC, theta,typ){
  
  # load estimates
  config <- paste0(
    '_nsim', nsim,
    '_S', S,
    '_K', K,
    '_m', m,
    '_ICC', ICC,
    '_CAC', CAC,
    '_theta', theta,
    '_', typ
  )
  
  infile=paste0(
    'est_files//',
    'estimates',
    config,
    '.RData'
  )
  load(file=infile)
  
  # Calculate rejection proportion for empirical power for HH
  # the proportion of times the test rejects the null hypothesis in the study.
  nsim_HH <- NA
  nsim_HH  <- sum(!is.na(res_est_fit$est_trt_HH) & !is.na(res_est_fit$se_trt_HH))
  pwr_HH   <- sum(abs(res_est_fit$est_trt_HH)/res_est_fit$se_trt_HH > 1.96, na.rm=TRUE)/nsim_HH
  varest_trt_HH<- var(res_est_fit$est_trt_HH, na.rm=TRUE)
  mest_trt_HH  <- mean(res_est_fit$est_trt_HH, na.rm=TRUE)
  mest_ICC_HH <- mean(res_est_fit$est_ICC_HH, na.rm=TRUE)
  varICC_HH  <- var(res_est_fit$est_ICC_HH, na.rm=TRUE)
  
  #Coverage for HH model without any adjustment
  num_cov_HH <- NA
  alpha <- (1 + 0.95)/2
  zstat_HH <- qnorm(alpha)
  #count the number of intervals that contain the true values.
  num_cov_HH <- sum(res_est_fit$est_trt_HH-(zstat_HH*res_est_fit$se_trt_HH) <= theta 
                    & res_est_fit$est_trt_HH+(zstat_HH*res_est_fit$se_trt_HH) >= theta,na.rm=T) 
  
  p_cov_HH <- ifelse(is.na(num_cov_HH/nsim_HH), NA, num_cov_HH/nsim_HH)
  
  #Calculate empirical power for HH with Kenward-Roger correction
  nsim_KR_HH<- NA
  tstat_KR_HH <- qt(alpha, res_est_fit$adj_ddf_KR_HH)
  nsim_KR_HH <- sum(!is.na(res_est_fit$est_trt_HH) & !is.na(res_est_fit$adj_se_KR_HH))
  pwr_KR_HH  <- sum(abs(res_est_fit$est_trt_HH)/res_est_fit$adj_se_KR_HH > tstat_KR_HH, na.rm=TRUE)/nsim_KR_HH
  
  #Calculate Coverage for HH with Kenward-Roger correction 
  num_KR_cov_HH<- NA
  #count the number of intervals that contain the true values.
  num_KR_cov_HH <- sum(res_est_fit$est_trt_HH-zstat_HH*(res_est_fit$adj_se_KR_HH) <= theta
                       & res_est_fit$est_trt_HH+zstat_HH*(res_est_fit$adj_se_KR_HH) >= theta,na.rm=T)
  p_KR_cov_HH <- ifelse(is.na(num_KR_cov_HH/nsim_KR_HH), NA, num_KR_cov_HH/nsim_KR_HH)
  
  #Calculate empirical power for HH with Kenward-Roger correction & non-singular fit
  nsim_KR_NSing_HH<- NA
  tstat_KR_NSing_HH <- qt(alpha, res_est_fit$adj_ddf_KR_NSing_HH)
  nsim_KR_NSing_HH <- sum(!is.na(res_est_fit$est_trt_HH) & !is.na(res_est_fit$adj_se_KR_NSing_HH))
  pwr_KR_NSing_HH  <- sum(abs(res_est_fit$est_trt_HH)/res_est_fit$adj_se_KR_NSing_HH > tstat_KR_NSing_HH, na.rm=TRUE)/nsim_KR_NSing_HH
  
  #Calculate Coverage for HH with Kenward-Roger correction & Non singular fit
  num_KR_NSing_cov_HH<- NA
  #count the number of intervals that contain the true values.
  num_KR_NSing_cov_HH <- sum(res_est_fit$est_trt_HH-zstat_HH*(res_est_fit$adj_se_KR_NSing_HH) <= theta
                             & res_est_fit$est_trt_HH+zstat_HH*(res_est_fit$adj_se_KR_NSing_HH) >= theta,na.rm=T)
  p_cov_KR_NSing_HH <- ifelse(is.na(num_KR_NSing_cov_HH/nsim_KR_NSing_HH), NA, num_KR_NSing_cov_HH/nsim_KR_NSing_HH)
  
  #Calculate empirical power for HH with Satterthwaite  correction
  nsim_Sat_HH<- NA
  tstat_Sat_HH <- qt(alpha, res_est_fit$adj_ddf_Sat_HH)
  nsim_Sat_HH  <- sum(!is.na(res_est_fit$est_trt_HH) & !is.na(res_est_fit$adj_se_Sat_HH))
  pwr_Sat_HH <- sum(abs(res_est_fit$est_trt_HH)/res_est_fit$adj_se_Sat_HH > tstat_Sat_HH,na.rm=TRUE)/nsim_Sat_HH
  
  #Calculate Coverage for HH with Satterthwaite  correction
  num_Sat_cov_HH<- NA
  #count the number of intervals that contain the true values.
  num_Sat_cov_HH <- sum(res_est_fit$est_trt_HH-zstat_HH*(res_est_fit$adj_se_Sat_HH) <= theta
                        & res_est_fit$est_trt_HH+zstat_HH*(res_est_fit$adj_se_Sat_HH) >= theta,na.rm=T)
  p_cov_Sat_HH <- ifelse(is.na(num_Sat_cov_HH/nsim_Sat_HH), NA, num_Sat_cov_HH/nsim_Sat_HH)
  
  #Calculate rejection proportion for empirical power for BE
  #the proportion of times the test rejects the null hypothesis in the study.
  nsim_BE<- NA
  nsim_BE  <-  sum(!is.na(res_est_fit$est_trt_BE) & !is.na(res_est_fit$se_trt_BE))
  pwr_BE   <- sum(abs(res_est_fit$est_trt_BE)/res_est_fit$se_trt_BE > 1.96, na.rm=TRUE)/nsim_BE
  varest_trt_BE<- var(res_est_fit$est_trt_BE, na.rm=TRUE)
  mest_trt_BE  <- mean(res_est_fit$est_trt_BE, na.rm=TRUE)
  mest_ICC_BE <- mean(res_est_fit$est_ICC_BE, na.rm=TRUE)
  varICC_BE  <- var(res_est_fit$est_ICC_BE, na.rm=TRUE)
  mest_CAC_BE <-  mean(res_est_fit$est_CAC_BE, na.rm=TRUE)
  var_CAC_BE  <- var(res_est_fit$est_CAC_BE, na.rm=TRUE)
  
  #Coverage for BE model without any adjustment
  num_cov_BE<- NA
  zstat_BE <- qnorm(alpha)
  #count the number of intervals that contain the true values.
  num_cov_BE <- sum(res_est_fit$est_trt_BE-zstat_BE*(res_est_fit$se_trt_BE) <= theta
                    & res_est_fit$est_trt_BE+zstat_BE*(res_est_fit$se_trt_BE) >= theta,na.rm=T)
  p_cov_BE <- ifelse(is.na(num_cov_BE/nsim_BE), NA, num_cov_BE/nsim_BE)
  
  
  #Calculate empirical power for BE with Kenward-Roger correction
  nsim_KR_BE<- NA
  tstat_KR_BE <- qt(alpha, res_est_fit$adj_ddf_KR_BE)
  nsim_KR_BE  <- sum(!is.na(res_est_fit$est_trt_BE) & !is.na(res_est_fit$adj_se_KR_BE))
  pwr_KR_BE <- sum(abs(res_est_fit$est_trt_BE)/res_est_fit$adj_se_KR_BE> tstat_KR_BE, na.rm=TRUE)/nsim_KR_BE
  
  
  #Calculate Coverage for BE with Kenward-Roger correction 
  num_KR_cov_BE<- NA
  #count the number of intervals that contain the true values.
  num_KR_cov_BE <- sum(res_est_fit$est_trt_BE-zstat_BE*(res_est_fit$adj_se_KR_BE) <= theta
                       & res_est_fit$est_trt_BE+zstat_BE*(res_est_fit$adj_se_KR_BE) >= theta,na.rm=T)
  p_KR_cov_BE <- ifelse(is.na(num_KR_cov_BE/nsim_KR_BE), NA, num_KR_cov_BE/nsim_KR_BE)
  
  #Calculate empirical power for BE with Kenward-Roger correction & Non singular fit 
  nsim_KR_NSing_BE<- NA
  tstat_KR_NSing_BE <- qt(alpha, res_est_fit$adj_ddf_KR_NSing_BE)
  nsim_KR_NSing_BE  <- sum(!is.na(res_est_fit$est_trt_BE) & !is.na(res_est_fit$adj_se_KR_NSing_BE))
  pwr_KR_NSing_BE <- sum(abs(res_est_fit$est_trt_BE)/res_est_fit$adj_se_KR_NSing_BE> tstat_KR_NSing_BE, na.rm=TRUE)/nsim_KR_NSing_BE 
  
  #Calculate Coverage for BE with Kenward-Roger correction & Non singular fit
  num_KR_NSing_cov_BE<- NA
  #count the number of intervals that contain the true values.
  num_KR_NSing_cov_BE <- sum(res_est_fit$est_trt_BE-zstat_BE*(res_est_fit$adj_se_KR_NSing_BE) <= theta
                             & res_est_fit$est_trt_BE+zstat_BE*(res_est_fit$adj_se_KR_NSing_BE) >= theta,na.rm=T)
  p_cov_KR_NSing_BE <- ifelse(is.na(num_KR_NSing_cov_BE/nsim_KR_NSing_BE), NA, num_KR_NSing_cov_BE/nsim_KR_NSing_BE)
  
  
  #Calculate empirical power for BE with Satterthwaite  correction
  nsim_Sat_BE<- NA
  tstat_Sat_BE <- qt(alpha, res_est_fit$adj_ddf_Sat_BE)
  nsim_Sat_BE  <- sum(!is.na(res_est_fit$est_trt_BE) & !is.na(res_est_fit$adj_se_Sat_BE))
  pwr_Sat_BE <- sum(abs(res_est_fit$est_trt_BE)/res_est_fit$adj_se_Sat_BE > tstat_Sat_BE,na.rm=TRUE)/nsim_Sat_BE
  
  
  #Calculate Coverage for BE with Satterthwaite  correction
  num_Sat_cov_BE<- NA
  #count the number of intervals that contain the true values.
  num_Sat_cov_BE <- sum(res_est_fit$est_trt_BE-zstat_BE*(res_est_fit$adj_se_Sat_BE) <= theta
                        & res_est_fit$est_trt_BE+zstat_BE*(res_est_fit$adj_se_Sat_BE) >= theta,na.rm=T)
  p_cov_Sat_BE <- ifelse(is.na(num_Sat_cov_BE/nsim_Sat_BE), NA, num_Sat_cov_BE/nsim_Sat_BE)
  
  
  #Type I error rate
  #error of concluding that there is a significant effect or difference
  #when there is, in fact, no such effect or difference.
  
  
  if (typ=='cat'){
    power  <- pow(VarSCcat(S,K,m, ICC, CAC),theta)
  } else if (typ=='lin') {
    power  <- pow(VarSClin(S,K,m, ICC, CAC),theta)
  }
  
  s_IsSing_HH <- sum(res_est_fit$IsSing_HH,na.rm=T)
  s_IsSing_BE <- sum(res_est_fit$IsSing_BE,na.rm=T)
  s_w_conv_HH <- sum(res_est_fit$w_conv_HH,na.rm=T)
  s_w_conv_BE <- sum(res_est_fit$w_conv_BE,na.rm=T)
  s_w_other_HH <- sum(res_est_fit$w_other_HH,na.rm=T)
  s_w_other_BE <- sum(res_est_fit$w_other_BE,na.rm=T)
  s_err_HH <- sum(res_est_fit$err_HH,na.rm=T)
  s_err_BE <- sum(res_est_fit$err_BE,na.rm=T)
  
  summary_emp_vals <- data.frame(nsim, S, K, m, ICC, CAC, theta,type=typ,power,
                                 nsim_HH=nsim_HH, 
                                 mesttrt_HH=mest_trt_HH, varest_trt_HH=varest_trt_HH,
                                 mestICC_HH=mest_ICC_HH,varICC_HH=varICC_HH, 
                                 nsim_BE=nsim_BE, 
                                 mesttrt_BE=mest_trt_BE,varest_trt_BE=varest_trt_BE,
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
                                 swarconv_HH=s_w_conv_HH, swarother_HH=s_w_other_HH,
                                 serr_HH=s_err_HH,swarconv_BE=s_w_conv_BE,
                                 swarother_BE=s_w_other_BE, serr_BE=s_err_BE)
  
  
  # Add the remaining performance measures
  summary_emp_vals  <- summary_emp_vals %>%
    mutate(Biastrt_HH = mesttrt_HH - theta, Biastrt_BE = mesttrt_BE - theta,
           BiasICC_HH = mestICC_HH - ICC, BiasICC_BE = mestICC_BE - ICC,
           BiasCAC_BE = mestCAC_BE - CAC)
  
  return(summary_emp_vals)
}

# List all files in the directory
file_name=as.list(dir(path = "est_files/", pattern="estimates_*"))

# Initialize an empty list to store results
smry_totl <- data.frame()
# Loop through each RData file
for (i in seq_along(file_name)) {
  # Load the RData file
  est_list <- file_name[[i]]
  load(file=paste0("est_files/",est_list))
  unique_val <- sapply(res_est_fit, unique)
  smry <- summary_res_fit(unique_val$nsim, unique_val$S, unique_val$K, unique_val$m, 
                          unique_val$ICC, unique_val$CAC, unique_val$theta,unique_val$typ)
  # Add the result to the data frame
  smry_totl <- rbind(smry_totl, smry)
}

#Save the data frames to a csv file
write.csv(smry_totl, "summary.csv", row.names = FALSE)
