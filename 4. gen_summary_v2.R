#generate summary results
summary_res_fit <- function(nsim, S, K, m, ICC, CAC, theta){
  ############################
  ##########Categorical (HH)##
  ############################
  # Calculate rejection proportion for empirical power for HH
  # the proportion of times the test rejects the null hypothesis in the study.
  nsim_HH_c <- NA
  nsim_HH_c  <- sum(!is.na(res_est_fit$est_trt_HH_c) & !is.na(res_est_fit$se_trt_HH_c))
  pow_HH_c   <- sum(abs(res_est_fit$est_trt_HH_c)/res_est_fit$se_trt_HH_c > 1.96, na.rm=TRUE)/nsim_HH_c
  
  # calculate mean & var for treatment effect, ICC, ddf
  mest_trt_HH_c  <- mean(res_est_fit$est_trt_HH_c, na.rm=TRUE)
  varest_trt_HH_c<- var(res_est_fit$est_trt_HH_c, na.rm=TRUE)
  mest_ICC_HH_c <- mean(res_est_fit$est_ICC_HH_c, na.rm=TRUE)
  varICC_HH_c  <- var(res_est_fit$est_ICC_HH_c, na.rm=TRUE)
  madj_ddf_KR_HH_c  <- mean(res_est_fit$adj_ddf_KR_HH_c, na.rm=TRUE)
  madj_ddf_KR_NSing_HH_c  <- mean(res_est_fit$adj_ddf_KR_NSing_HH_c, na.rm=TRUE)
  madj_ddf_Sat_HH_c  <- mean(res_est_fit$adj_ddf_Sat_HH_c, na.rm=TRUE)
  
  # Coverage for HH model without any adjustment
  num_cov_HH_c <- NA
  alpha <- (1 + 0.95)/2
  zstat_HH_c <- qnorm(alpha)
  # count the number of intervals that contain the true values.
  num_cov_HH_c <- sum(res_est_fit$est_trt_HH_c-(zstat_HH_c*res_est_fit$se_trt_HH_c) <= theta 
                    & res_est_fit$est_trt_HH_c+(zstat_HH_c*res_est_fit$se_trt_HH_c) >= theta,na.rm=T) 
  p_cov_HH_c <- ifelse(is.na(num_cov_HH_c/nsim_HH_c), NA, num_cov_HH_c/nsim_HH_c)
  # Calculate empirical power for HH with Kenward-Roger correction
  nsim_KR_HH_c<- NA
  tstat_KR_HH_c <- qt(alpha, res_est_fit$adj_ddf_KR_HH_c)
  nsim_KR_HH_c <- sum(!is.na(res_est_fit$est_trt_HH_c) & !is.na(res_est_fit$adj_se_KR_HH_c))
  pow_KR_HH_c  <- sum(abs(res_est_fit$est_trt_HH_c)/res_est_fit$adj_se_KR_HH_c > tstat_KR_HH_c, na.rm=TRUE)/nsim_KR_HH_c
  
  # Calculate Coverage for HH with Kenward-Roger correction 
  num_KR_cov_HH_c<- NA
  # count the number of intervals that contain the true values.
  num_KR_cov_HH_c <- sum(res_est_fit$est_trt_HH_c-tstat_KR_HH_c*(res_est_fit$adj_se_KR_HH_c) <= theta
                       & res_est_fit$est_trt_HH_c+tstat_KR_HH_c*(res_est_fit$adj_se_KR_HH_c) >= theta,na.rm=T)
  p_cov_KR_HH_c <- ifelse(is.na(num_KR_cov_HH_c/nsim_KR_HH_c), NA, num_KR_cov_HH_c/nsim_KR_HH_c)
  
  # Calculate empirical power for HH with Kenward-Roger correction & non-singular fit
  nsim_KR_NSing_HH_c<- NA
  tstat_KR_NSing_HH_c <- qt(alpha, res_est_fit$adj_ddf_KR_NSing_HH_c)
  nsim_KR_NSing_HH_c <- sum(!is.na(res_est_fit$est_trt_HH_c) & !is.na(res_est_fit$adj_se_KR_NSing_HH_c))
  pow_KR_NSing_HH_c  <- sum(abs(res_est_fit$est_trt_HH_c)/res_est_fit$adj_se_KR_NSing_HH_c > tstat_KR_NSing_HH_c, na.rm=TRUE)/nsim_KR_NSing_HH_c
  
  # Calculate Coverage for HH with Kenward-Roger correction & Non singular fit
  num_KR_NSing_cov_HH_c<- NA
  # count the number of intervals that contain the true values.
  num_KR_NSing_cov_HH_c <- sum(res_est_fit$est_trt_HH_c-tstat_KR_NSing_HH_c*(res_est_fit$adj_se_KR_NSing_HH_c) <= theta
                             & res_est_fit$est_trt_HH_c+tstat_KR_NSing_HH_c*(res_est_fit$adj_se_KR_NSing_HH_c) >= theta,na.rm=T)
  p_cov_KR_NSing_HH_c <- ifelse(is.na(num_KR_NSing_cov_HH_c/nsim_KR_NSing_HH_c), NA, num_KR_NSing_cov_HH_c/nsim_KR_NSing_HH_c)
  
  # Calculate empirical power for HH with Satterthwaite  correction
  nsim_Sat_HH_c<- NA
  tstat_Sat_HH_c <- qt(alpha, res_est_fit$adj_ddf_Sat_HH_c)
  nsim_Sat_HH_c  <- sum(!is.na(res_est_fit$est_trt_HH_c) & !is.na(res_est_fit$adj_se_Sat_HH_c))
  pow_Sat_HH_c <- sum(abs(res_est_fit$est_trt_HH_c)/res_est_fit$adj_se_Sat_HH_c > tstat_Sat_HH_c,na.rm=TRUE)/nsim_Sat_HH_c
  
  # Calculate Coverage for HH with Satterthwaite  correction
  num_Sat_cov_HH_c<- NA
  # count the number of intervals that contain the true values.
  num_Sat_cov_HH_c <- sum(res_est_fit$est_trt_HH_c-tstat_Sat_HH_c*(res_est_fit$adj_se_Sat_HH_c) <= theta
                        & res_est_fit$est_trt_HH_c+tstat_Sat_HH_c*(res_est_fit$adj_se_Sat_HH_c) >= theta,na.rm=T)
  p_cov_Sat_HH_c <- ifelse(is.na(num_Sat_cov_HH_c/nsim_Sat_HH_c), NA, num_Sat_cov_HH_c/nsim_Sat_HH_c)
  
  ############################
  ##########linear (HH)#######
  ############################
  
  nsim_HH_l <- NA
  nsim_HH_l  <- sum(!is.na(res_est_fit$est_trt_HH_l) & !is.na(res_est_fit$se_trt_HH_l))
  pow_HH_l   <- sum(abs(res_est_fit$est_trt_HH_l)/res_est_fit$se_trt_HH_l > 1.96, na.rm=TRUE)/nsim_HH_l
  
  # calculate mean & var for treatment effect, ICC, ddf
  mest_trt_HH_l  <- mean(res_est_fit$est_trt_HH_l, na.rm=TRUE)
  varest_trt_HH_l<- var(res_est_fit$est_trt_HH_l, na.rm=TRUE)
  mest_ICC_HH_l <- mean(res_est_fit$est_ICC_HH_l, na.rm=TRUE)
  varICC_HH_l  <- var(res_est_fit$est_ICC_HH_l, na.rm=TRUE)
  madj_ddf_KR_HH_l  <- mean(res_est_fit$adj_ddf_KR_HH_l, na.rm=TRUE)
  madj_ddf_KR_NSing_HH_l  <- mean(res_est_fit$adj_ddf_KR_NSing_HH_l, na.rm=TRUE)
  madj_ddf_Sat_HH_l  <- mean(res_est_fit$adj_ddf_Sat_HH_l, na.rm=TRUE)
  
  
  # Coverage for HH model without any adjustment
  num_cov_HH_l <- NA
  alpha <- (1 + 0.95)/2
  zstat_HH_l <- qnorm(alpha)
  # count the number of intervals that contain the true values.
  num_cov_HH_l <- sum(res_est_fit$est_trt_HH_l-(zstat_HH_l*res_est_fit$se_trt_HH_l) <= theta 
                    & res_est_fit$est_trt_HH_l+(zstat_HH_l*res_est_fit$se_trt_HH_l) >= theta,na.rm=T) 
  
  p_cov_HH_l <- ifelse(is.na(num_cov_HH_l/nsim_HH_l), NA, num_cov_HH_l/nsim_HH_l)
  
  # Calculate empirical power for HH with Kenward-Roger correction
  nsim_KR_HH_l<- NA
  tstat_KR_HH_l <- qt(alpha, res_est_fit$adj_ddf_KR_HH_l)
  nsim_KR_HH_l <- sum(!is.na(res_est_fit$est_trt_HH_l) & !is.na(res_est_fit$adj_se_KR_HH_l))
  pow_KR_HH_l  <- sum(abs(res_est_fit$est_trt_HH_l)/res_est_fit$adj_se_KR_HH_l > tstat_KR_HH_l, na.rm=TRUE)/nsim_KR_HH_l
  
  # Calculate Coverage for HH with Kenward-Roger correction 
  num_KR_cov_HH_l<- NA
  # count the number of intervals that contain the true values.
  num_KR_cov_HH_l <- sum(res_est_fit$est_trt_HH_l-tstat_KR_HH_l*(res_est_fit$adj_se_KR_HH_l) <= theta
                       & res_est_fit$est_trt_HH_l+tstat_KR_HH_l*(res_est_fit$adj_se_KR_HH_l) >= theta,na.rm=T)
  p_cov_KR_HH_l <- ifelse(is.na(num_KR_cov_HH_l/nsim_KR_HH_l), NA, num_KR_cov_HH_l/nsim_KR_HH_l)
  
  # Calculate empirical power for HH with Kenward-Roger correction & non-singular fit
  nsim_KR_NSing_HH_l<- NA
  tstat_KR_NSing_HH_l <- qt(alpha, res_est_fit$adj_ddf_KR_NSing_HH_l)
  nsim_KR_NSing_HH_l <- sum(!is.na(res_est_fit$est_trt_HH_l) & !is.na(res_est_fit$adj_se_KR_NSing_HH_l))
  pow_KR_NSing_HH_l  <- sum(abs(res_est_fit$est_trt_HH_l)/res_est_fit$adj_se_KR_NSing_HH_l > tstat_KR_NSing_HH_l, na.rm=TRUE)/nsim_KR_NSing_HH_l
  
  # Calculate Coverage for HH with Kenward-Roger correction & Non singular fit
  num_KR_NSing_cov_HH_l<- NA
  # count the number of intervals that contain the true values.
  num_KR_NSing_cov_HH_l <- sum(res_est_fit$est_trt_HH_l-tstat_KR_NSing_HH_l*(res_est_fit$adj_se_KR_NSing_HH_l) <= theta
                             & res_est_fit$est_trt_HH_l+tstat_KR_NSing_HH_l*(res_est_fit$adj_se_KR_NSing_HH_l) >= theta,na.rm=T)
  p_cov_KR_NSing_HH_l <- ifelse(is.na(num_KR_NSing_cov_HH_l/nsim_KR_NSing_HH_l), NA, num_KR_NSing_cov_HH_l/nsim_KR_NSing_HH_l)
  
  # Calculate empirical power for HH with Satterthwaite  correction
  nsim_Sat_HH_l<- NA
  tstat_Sat_HH_l <- qt(alpha, res_est_fit$adj_ddf_Sat_HH_l)
  nsim_Sat_HH_l  <- sum(!is.na(res_est_fit$est_trt_HH_l) & !is.na(res_est_fit$adj_se_Sat_HH_l))
  pow_Sat_HH_l <- sum(abs(res_est_fit$est_trt_HH_l)/res_est_fit$adj_se_Sat_HH_l > tstat_Sat_HH_l,na.rm=TRUE)/nsim_Sat_HH_l
  
  # Calculate Coverage for HH with Satterthwaite  correction
  num_Sat_cov_HH_l<- NA
  # count the number of intervals that contain the true values.
  num_Sat_cov_HH_l <- sum(res_est_fit$est_trt_HH_l-tstat_Sat_HH_l*(res_est_fit$adj_se_Sat_HH_l) <= theta
                        & res_est_fit$est_trt_HH_l+tstat_Sat_HH_l*(res_est_fit$adj_se_Sat_HH_l) >= theta,na.rm=T)
  p_cov_Sat_HH_l <- ifelse(is.na(num_Sat_cov_HH_l/nsim_Sat_HH_l), NA, num_Sat_cov_HH_l/nsim_Sat_HH_l)
 
  ############################
  ##########Categorical (BE)##
  ############################
  
  # Calculate rejection proportion for empirical power for BE
  # the proportion of times the test rejects the null hypothesis in the study.
  nsim_BE_c<- NA
  nsim_BE_c  <-  sum(!is.na(res_est_fit$est_trt_BE_c) & !is.na(res_est_fit$se_trt_BE_c))
  pow_BE_c   <- sum(abs(res_est_fit$est_trt_BE_c)/res_est_fit$se_trt_BE_c > 1.96, na.rm=TRUE)/nsim_BE_c
  
  # calculate mean & var for treatment effect, ICC, CAC, and ddf
  mest_trt_BE_c  <- mean(res_est_fit$est_trt_BE_c, na.rm=TRUE)
  varest_trt_BE_c<- var(res_est_fit$est_trt_BE_c, na.rm=TRUE)
  mest_ICC_BE_c <- mean(res_est_fit$est_ICC_BE_c, na.rm=TRUE)
  varICC_BE_c  <- var(res_est_fit$est_ICC_BE_c, na.rm=TRUE)
  mest_CAC_BE_c <-  mean(res_est_fit$est_CAC_BE_c, na.rm=TRUE)
  var_CAC_BE_c  <- var(res_est_fit$est_CAC_BE_c, na.rm=TRUE)
  madj_ddf_KR_BE_c <- mean(res_est_fit$adj_ddf_KR_BE_c, na.rm=TRUE)
  madj_ddf_KR_NSing_BE_c  <- mean(res_est_fit$adj_ddf_KR_NSing_BE_c, na.rm=TRUE)
  madj_ddf_Sat_BE_c <- mean(res_est_fit$adj_ddf_Sat_BE_c, na.rm=TRUE)
  
  # Coverage for BE model without any adjustment
  num_cov_BE_c<- NA
  zstat_BE_c <- qnorm(alpha)
  # count the number of intervals that contain the true values.
  num_cov_BE_c <- sum(res_est_fit$est_trt_BE_c-zstat_BE_c*(res_est_fit$se_trt_BE_c) <= theta
                    & res_est_fit$est_trt_BE_c+zstat_BE_c*(res_est_fit$se_trt_BE_c) >= theta,na.rm=T)
  p_cov_BE_c <- ifelse(is.na(num_cov_BE_c/nsim_BE_c), NA, num_cov_BE_c/nsim_BE_c)
  
  
  # Calculate empirical power for BE with Kenward-Roger correction
  nsim_KR_BE_c<- NA
  tstat_KR_BE_c <- qt(alpha, res_est_fit$adj_ddf_KR_BE_c)
  nsim_KR_BE_c  <- sum(!is.na(res_est_fit$est_trt_BE_c) & !is.na(res_est_fit$adj_se_KR_BE_c))
  pow_KR_BE_c <- sum(abs(res_est_fit$est_trt_BE_c)/res_est_fit$adj_se_KR_BE_c> tstat_KR_BE_c, na.rm=TRUE)/nsim_KR_BE_c
  
  
  # Calculate Coverage for BE with Kenward-Roger correction 
  num_KR_cov_BE_c<- NA
  # count the number of intervals that contain the true values.
  num_KR_cov_BE_c <- sum(res_est_fit$est_trt_BE_c-tstat_KR_BE_c*(res_est_fit$adj_se_KR_BE_c) <= theta
                       & res_est_fit$est_trt_BE_c+tstat_KR_BE_c*(res_est_fit$adj_se_KR_BE_c) >= theta,na.rm=T)
  p_cov_KR_BE_c <- ifelse(is.na(num_KR_cov_BE_c/nsim_KR_BE_c), NA, num_KR_cov_BE_c/nsim_KR_BE_c)
  
  # Calculate empirical power for BE with Kenward-Roger correction & Non singular fit 
  nsim_KR_NSing_BE_c<- NA
  tstat_KR_NSing_BE_c <- qt(alpha, res_est_fit$adj_ddf_KR_NSing_BE_c)
  nsim_KR_NSing_BE_c  <- sum(!is.na(res_est_fit$est_trt_BE_c) & !is.na(res_est_fit$adj_se_KR_NSing_BE_c))
  pow_KR_NSing_BE_c <- sum(abs(res_est_fit$est_trt_BE_c)/res_est_fit$adj_se_KR_NSing_BE_c> tstat_KR_NSing_BE_c, na.rm=TRUE)/nsim_KR_NSing_BE_c 
  
  # Calculate Coverage for BE with Kenward-Roger correction & Non singular fit
  num_KR_NSing_cov_BE_c<- NA
  # count the number of intervals that contain the true values.
  num_KR_NSing_cov_BE_c <- sum(res_est_fit$est_trt_BE_c-tstat_KR_NSing_BE_c*(res_est_fit$adj_se_KR_NSing_BE_c) <= theta
                             & res_est_fit$est_trt_BE_c+tstat_KR_NSing_BE_c*(res_est_fit$adj_se_KR_NSing_BE_c) >= theta,na.rm=T)
  p_cov_KR_NSing_BE_c <- ifelse(is.na(num_KR_NSing_cov_BE_c/nsim_KR_NSing_BE_c), NA, num_KR_NSing_cov_BE_c/nsim_KR_NSing_BE_c)
  
  
  # Calculate empirical power for BE with Satterthwaite  correction
  nsim_Sat_BE_c<- NA
  tstat_Sat_BE_c <- qt(alpha, res_est_fit$adj_ddf_Sat_BE_c)
  nsim_Sat_BE_c  <- sum(!is.na(res_est_fit$est_trt_BE_c) & !is.na(res_est_fit$adj_se_Sat_BE_c))
  pow_Sat_BE_c <- sum(abs(res_est_fit$est_trt_BE_c)/res_est_fit$adj_se_Sat_BE_c > tstat_Sat_BE_c,na.rm=TRUE)/nsim_Sat_BE_c
  
  
  # Calculate Coverage for BE with Satterthwaite  correction
  num_Sat_cov_BE_c<- NA
  # count the number of intervals that contain the true values.
  num_Sat_cov_BE_c <- sum(res_est_fit$est_trt_BE_c-tstat_Sat_BE_c*(res_est_fit$adj_se_Sat_BE_c) <= theta
                        & res_est_fit$est_trt_BE_c+tstat_Sat_BE_c*(res_est_fit$adj_se_Sat_BE_c) >= theta,na.rm=T)
  p_cov_Sat_BE_c <- ifelse(is.na(num_Sat_cov_BE_c/nsim_Sat_BE_c), NA, num_Sat_cov_BE_c/nsim_Sat_BE_c)
  
  
  #Type I error rate
  # error of concluding that there is a significant effect or difference
  # when there is, in fact, no such effect or difference.
  
  ############################
  ##########linear (BE)#######
  ############################
  
  # Calculate rejection proportion for empirical power for BE
  # the proportion of times the test rejects the null hypothesis in the study.
  nsim_BE_l<- NA
  nsim_BE_l  <-  sum(!is.na(res_est_fit$est_trt_BE_l) & !is.na(res_est_fit$se_trt_BE_l))
  pow_BE_l   <- sum(abs(res_est_fit$est_trt_BE_l)/res_est_fit$se_trt_BE_l > 1.96, na.rm=TRUE)/nsim_BE_l
  
  # calculate mean & var for treatment effect, ICC, CAC, and ddf
  mest_trt_BE_l  <- mean(res_est_fit$est_trt_BE_l, na.rm=TRUE)
  varest_trt_BE_l<- var(res_est_fit$est_trt_BE_l, na.rm=TRUE)
  mest_ICC_BE_l <- mean(res_est_fit$est_ICC_BE_l, na.rm=TRUE)
  varICC_BE_l  <- var(res_est_fit$est_ICC_BE_l, na.rm=TRUE)
  mest_CAC_BE_l <-  mean(res_est_fit$est_CAC_BE_l, na.rm=TRUE)
  var_CAC_BE_l  <- var(res_est_fit$est_CAC_BE_l, na.rm=TRUE)
  madj_ddf_KR_BE_l <- mean(res_est_fit$adj_ddf_KR_BE_l, na.rm=TRUE)
  madj_ddf_KR_NSing_BE_l  <- mean(res_est_fit$adj_ddf_KR_NSing_BE_l, na.rm=TRUE)
  madj_ddf_Sat_BE_l <- mean(res_est_fit$adj_ddf_Sat_BE_l, na.rm=TRUE)
  
  # Coverage for BE model without any adjustment
  num_cov_BE_l<- NA
  zstat_BE_l <- qnorm(alpha)
  # count the number of intervals that contain the true values.
  num_cov_BE_l <- sum(res_est_fit$est_trt_BE_l-zstat_BE_l*(res_est_fit$se_trt_BE_l) <= theta
                    & res_est_fit$est_trt_BE_l+zstat_BE_l*(res_est_fit$se_trt_BE_l) >= theta,na.rm=T)
  p_cov_BE_l <- ifelse(is.na(num_cov_BE_l/nsim_BE_l), NA, num_cov_BE_l/nsim_BE_l)
  
  
  # Calculate empirical power for BE with Kenward-Roger correction
  nsim_KR_BE_l<- NA
  tstat_KR_BE_l <- qt(alpha, res_est_fit$adj_ddf_KR_BE_l)
  nsim_KR_BE_l  <- sum(!is.na(res_est_fit$est_trt_BE_l) & !is.na(res_est_fit$adj_se_KR_BE_l))
  pow_KR_BE_l <- sum(abs(res_est_fit$est_trt_BE_l)/res_est_fit$adj_se_KR_BE_l> tstat_KR_BE_l, na.rm=TRUE)/nsim_KR_BE_l
  
  
  # Calculate Coverage for BE with Kenward-Roger correction 
  num_KR_cov_BE_l<- NA
  # count the number of intervals that contain the true values.
  num_KR_cov_BE_l <- sum(res_est_fit$est_trt_BE_l-tstat_KR_BE_l*(res_est_fit$adj_se_KR_BE_l) <= theta
                       & res_est_fit$est_trt_BE_l+tstat_KR_BE_l*(res_est_fit$adj_se_KR_BE_l) >= theta,na.rm=T)
  p_cov_KR_BE_l <- ifelse(is.na(num_KR_cov_BE_l/nsim_KR_BE_l), NA, num_KR_cov_BE_l/nsim_KR_BE_l)
  
  # Calculate empirical power for BE with Kenward-Roger correction & Non singular fit 
  nsim_KR_NSing_BE_l<- NA
  tstat_KR_NSing_BE_l <- qt(alpha, res_est_fit$adj_ddf_KR_NSing_BE_l)
  nsim_KR_NSing_BE_l  <- sum(!is.na(res_est_fit$est_trt_BE_l) & !is.na(res_est_fit$adj_se_KR_NSing_BE_l))
  pow_KR_NSing_BE_l <- sum(abs(res_est_fit$est_trt_BE_l)/res_est_fit$adj_se_KR_NSing_BE_l> tstat_KR_NSing_BE_l, na.rm=TRUE)/nsim_KR_NSing_BE_l 
  
  # Calculate Coverage for BE with Kenward-Roger correction & Non singular fit
  num_KR_NSing_cov_BE_l<- NA
  # count the number of intervals that contain the true values.
  num_KR_NSing_cov_BE_l <- sum(res_est_fit$est_trt_BE_l-tstat_KR_NSing_BE_l*(res_est_fit$adj_se_KR_NSing_BE_l) <= theta
                             & res_est_fit$est_trt_BE_l+tstat_KR_NSing_BE_l*(res_est_fit$adj_se_KR_NSing_BE_l) >= theta,na.rm=T)
  p_cov_KR_NSing_BE_l <- ifelse(is.na(num_KR_NSing_cov_BE_l/nsim_KR_NSing_BE_l), NA, num_KR_NSing_cov_BE_l/nsim_KR_NSing_BE_l)
  
  
  # Calculate empirical power for BE with Satterthwaite  correction
  nsim_Sat_BE_l<- NA
  tstat_Sat_BE_l <- qt(alpha, res_est_fit$adj_ddf_Sat_BE_l)
  nsim_Sat_BE_l  <- sum(!is.na(res_est_fit$est_trt_BE_l) & !is.na(res_est_fit$adj_se_Sat_BE_l))
  pow_Sat_BE_l <- sum(abs(res_est_fit$est_trt_BE_l)/res_est_fit$adj_se_Sat_BE_l > tstat_Sat_BE_l,na.rm=TRUE)/nsim_Sat_BE_l
  
  
  # Calculate Coverage for BE with Satterthwaite  correction
  num_Sat_cov_BE_l<- NA
  # count the number of intervals that contain the true values.
  num_Sat_cov_BE_l <- sum(res_est_fit$est_trt_BE_l-tstat_Sat_BE_l*(res_est_fit$adj_se_Sat_BE_l) <= theta
                        & res_est_fit$est_trt_BE_l+tstat_Sat_BE_l*(res_est_fit$adj_se_Sat_BE_l) >= theta,na.rm=T)
  p_cov_Sat_BE_l <- ifelse(is.na(num_Sat_cov_BE_l/nsim_Sat_BE_l), NA, num_Sat_cov_BE_l/nsim_Sat_BE_l)
  
  
  #Type I error rate
  # error of concluding that there is a significant effect or difference
  # when there is, in fact, no such effect or difference.
  
    power_c  <- pow(VarSCcat(S,K,m, ICC, CAC),theta)
    power_l  <- pow(VarSClin(S,K,m, ICC, CAC),theta)
  
  s_IsSing_HH_c <- sum(res_est_fit$IsSing_HH_c,na.rm=T)
  s_IsSing_BE_c <- sum(res_est_fit$IsSing_BE_c,na.rm=T)
  s_w_conv_HH_c <- sum(res_est_fit$w_conv_HH_c,na.rm=T)
  s_w_conv_BE_c <- sum(res_est_fit$w_conv_BE_c,na.rm=T)
  s_w_other_HH_c <- sum(res_est_fit$w_other_HH_c,na.rm=T)
  s_w_other_BE_c <- sum(res_est_fit$w_other_BE_c,na.rm=T)
  s_err_HH_c <- sum(res_est_fit$err_HH_c,na.rm=T)
  s_err_BE_c <- sum(res_est_fit$err_BE_c,na.rm=T)
  
  s_IsSing_HH_l <- sum(res_est_fit$IsSing_HH_l,na.rm=T)
  s_IsSing_BE_l <- sum(res_est_fit$IsSing_BE_l,na.rm=T)
  s_w_conv_HH_l <- sum(res_est_fit$w_conv_HH_l,na.rm=T)
  s_w_conv_BE_l <- sum(res_est_fit$w_conv_BE_l,na.rm=T)
  s_w_other_HH_l <- sum(res_est_fit$w_other_HH_l,na.rm=T)
  s_w_other_BE_l <- sum(res_est_fit$w_other_BE_l,na.rm=T)
  s_err_HH_l <- sum(res_est_fit$err_HH_l,na.rm=T)
  s_err_BE_l <- sum(res_est_fit$err_BE_l,na.rm=T)
  
  summary_emp_vals <- data.frame(nsim, S, K, m, ICC, CAC, theta,power_c,power_l,
                                 nsim_HH_c=nsim_HH_c, 
                                 nsimKRNSing_HH_c=nsim_KR_NSing_HH_c,
                                 mesttrt_HH_c=mest_trt_HH_c, varesttrt_HH_c=varest_trt_HH_c,
                                 mestICC_HH_c=mest_ICC_HH_c,varICC_HH_c=varICC_HH_c, 
                                 madjddfKR_HH_c=madj_ddf_KR_HH_c,madjddfKRNSing_HH_c=madj_ddf_KR_NSing_HH_c,
                                 madjddfSat_HH_c=madj_ddf_Sat_HH_c,
                                 nsim_BE_c=nsim_BE_c, 
                                 nsimKRNSing_BE_c=nsim_KR_NSing_BE_c,
                                 mesttrt_BE_c=mest_trt_BE_c,varesttrt_BE_c=varest_trt_BE_c,
                                 mestICC_BE_c=mest_ICC_BE_c,varICC_BE_c=varICC_BE_c,
                                 mestCAC_BE_c=mest_CAC_BE_c,varCAC_BE_c=var_CAC_BE_c,
                                 madjddfKR_BE_c=madj_ddf_KR_BE_c,madjddfKRNSing_BE_c=madj_ddf_KR_NSing_BE_c,
                                 madjddfSat_BE_c=madj_ddf_Sat_BE_c,
                                 sIsSing_HH_c=s_IsSing_HH_c,pow_HH_c=pow_HH_c,
                                 powKR_HH_c=pow_KR_HH_c, powKRnSing_HH_c=pow_KR_NSing_HH_c,
                                 powSat_HH_c=pow_Sat_HH_c,
                                 cov_HH_c=p_cov_HH_c,covKR_HH_c=p_cov_KR_HH_c,
                                 covKRnSing_HH_c=p_cov_KR_NSing_HH_c,covSat_HH_c=p_cov_Sat_HH_c,
                                 sIsSing_BE_c=s_IsSing_BE_c,pow_BE_c=pow_BE_c,
                                 powKR_BE_c=pow_KR_BE_c,powKRnSing_BE_c=pow_KR_NSing_BE_c,
                                 powSat_BE_c=pow_Sat_BE_c,
                                 cov_BE_c=p_cov_BE_c,covKR_BE_c=p_cov_KR_BE_c,
                                 covKRnSing_BE_c=p_cov_KR_NSing_BE_c,covSat_BE_c=p_cov_Sat_BE_c,
                                 swarconv_HH_c=s_w_conv_HH_c, swarother_HH_c=s_w_other_HH_c,
                                 serr_HH_c=s_err_HH_c,swarconv_BE_c=s_w_conv_BE_c,
                                 swarother_BE_c=s_w_other_BE_c, serr_BE_c=s_err_BE_c,
                                 nsim_HH_l=nsim_HH_l, 
                                 nsimKRNSing_HH_l=nsim_KR_NSing_HH_l,
                                 mesttrt_HH_l=mest_trt_HH_l, varesttrt_HH_l=varest_trt_HH_l,
                                 mestICC_HH_l=mest_ICC_HH_l,varICC_HH_l=varICC_HH_l, 
                                 madjddfKR_HH_l=madj_ddf_KR_HH_l,madjddfKRNSing_HH_l=madj_ddf_KR_NSing_HH_l,
                                 madjddfSat_HH_l=madj_ddf_Sat_HH_l,
                                 nsim_BE_l=nsim_BE_l, 
                                 nsimKRNSing_BE_l=nsim_KR_NSing_BE_l,
                                 mesttrt_BE_l=mest_trt_BE_l,varesttrt_BE_l=varest_trt_BE_l,
                                 mestICC_BE_l=mest_ICC_BE_l,varICC_BE_l=varICC_BE_l,
                                 mestCAC_BE_l=mest_CAC_BE_l,varCAC_BE_l=var_CAC_BE_l,
                                 madjddfKR_BE_l=madj_ddf_KR_BE_l,madjddfKRNSing_BE_l=madj_ddf_KR_NSing_BE_l,
                                 madjddfSat_BE_l=madj_ddf_Sat_BE_l,
                                 sIsSing_HH_l=s_IsSing_HH_l,pow_HH_l=pow_HH_l,
                                 powKR_HH_l=pow_KR_HH_l, powKRnSing_HH_l=pow_KR_NSing_HH_l,
                                 powSat_HH_l=pow_Sat_HH_l,
                                 cov_HH_l=p_cov_HH_l,covKR_HH_l=p_cov_KR_HH_l,
                                 covKRnSing_HH_l=p_cov_KR_NSing_HH_l,covSat_HH_l=p_cov_Sat_HH_l,
                                 sIsSing_BE_l=s_IsSing_BE_l,pow_BE_l=pow_BE_l,
                                 powKR_BE_l=pow_KR_BE_l,powKRnSing_BE_l=pow_KR_NSing_BE_l,
                                 powSat_BE_l=pow_Sat_BE_l,
                                 cov_BE_l=p_cov_BE_l,covKR_BE_l=p_cov_KR_BE_l,
                                 covKRnSing_BE_l=p_cov_KR_NSing_BE_l,covSat_BE_l=p_cov_Sat_BE_l,
                                 swarconv_HH_l=s_w_conv_HH_l, swarother_HH_l=s_w_other_HH_l,
                                 serr_HH_l=s_err_HH_l,swarconv_BE_l=s_w_conv_BE_l,
                                 swarother_BE_l=s_w_other_BE_l, serr_BE_l=s_err_BE_l)
  
  # Add the remaining performance measures
  summary_emp_vals  <- summary_emp_vals %>%
    mutate(Biastrt_HH_c = mesttrt_HH_c - theta, Biastrt_BE_c = mesttrt_BE_c - theta,
           BiasICC_HH_c = mestICC_HH_c - ICC, BiasICC_BE_c = mestICC_BE_c - ICC,
           BiasCAC_BE_c = mestCAC_BE_c - CAC,
           MCSE_HH_c=sqrt(varesttrt_HH_c/nsim_HH_c), MCSE_BE_c=sqrt(varesttrt_BE_c/nsim_BE_c),
           Biastrt_HH_l = mesttrt_HH_l - theta, Biastrt_BE_l = mesttrt_BE_l - theta,
           BiasICC_HH_l = mestICC_HH_l - ICC, BiasICC_BE_l = mestICC_BE_l - ICC,
           BiasCAC_BE_l = mestCAC_BE_l - CAC,
           MCSE_HH_l=sqrt(varesttrt_HH_l/nsim_HH_l), MCSE_BE_l=sqrt(varesttrt_BE_l/nsim_BE_l),
           MCSEpow_HH_c=sqrt(pow_HH_c*(1-pow_HH_c)/nsim_HH_c),
           MCSEcov_HH_c=sqrt(cov_HH_c*(1-cov_HH_c)/nsim_HH_c) ,
           MCSEpowKR_HH_c=sqrt(powKR_HH_c*(1-powKR_HH_c)/nsim_HH_c) ,
           MCSEcovKR_HH_c=sqrt(covKR_HH_c*(1-covKR_HH_c)/nsim_HH_c),
           MCSEpowSat_HH_c=sqrt(powSat_HH_c*(1-powSat_HH_c)/nsim_HH_c),
           MCSEcovSat_HH_c=sqrt(covSat_HH_c*(1-covSat_HH_c)/nsim_HH_c) ,
           MCSEpow_HH_l=sqrt(pow_HH_l*(1-pow_HH_l)/nsim_HH_l) ,
           MCSEcov_HH_l=sqrt(cov_HH_l*(1-cov_HH_l)/nsim_HH_l) ,
           MCSEpowKR_HH_l=sqrt(powKR_HH_l*(1-powKR_HH_l)/nsim_HH_l) ,
           MCSEcovKR_HH_l=sqrt(covKR_HH_l*(1-covKR_HH_l)/nsim_HH_l) ,
           MCSEpowSat_HH_l=sqrt(powSat_HH_l*(1-powSat_HH_l)/nsim_HH_l) ,
           MCSEcovSat_HH_l=sqrt(covSat_HH_l*(1-covSat_HH_l)/nsim_HH_l) ,
           MCSEpow_BE_c= sqrt(pow_BE_c*(1-pow_BE_c)/nsim_BE_c),
           MCSEcov_BE_c= sqrt(cov_BE_c*(1-cov_BE_c)/nsim_BE_c),
           MCSEpowKR_BE_c=sqrt(powKR_BE_c*(1-powKR_BE_c)/nsim_BE_c) ,
           MCSEcovKR_BE_c= sqrt(covKR_BE_c*(1-covKR_BE_c)/nsim_BE_c),
           MCSEpowSat_BE_c= sqrt(powSat_BE_c*(1-powSat_BE_c)/nsim_BE_c),
           MCSEcovSat_BE_c= sqrt(covSat_BE_c*(1-covSat_BE_c)/nsim_BE_c),
           MCSEpow_BE_l= sqrt(pow_BE_l*(1-pow_BE_l)/nsim_BE_l),
           MCSEcov_BE_l= sqrt(cov_BE_l*(1-cov_BE_l)/nsim_BE_l),
           MCSEpowKR_BE_l= sqrt(powKR_BE_l*(1-powKR_BE_l)/nsim_BE_l),
           MCSEcovKR_BE_l= sqrt(covKR_BE_l*(1-covKR_BE_l)/nsim_BE_l),
           MCSEpowSat_BE_l= sqrt(powSat_BE_l*(1-powSat_BE_l)/nsim_BE_l),
           MCSEcovSat_BE_l= sqrt(covSat_BE_l*(1-covSat_BE_l)/nsim_BE_l),
           )
  
  return(summary_emp_vals)
}
#################################################################
#######get summary for all data files of estimations#############
#################################################################
# List all files in the directory
file_name=as.list(dir(path = "est_files_v2/", pattern="estimates_*"))

# Initialize an empty list to store results
smry_totl <- data.frame()

#Loop through each RData file
for (i in seq_along(file_name)) {
  # Load the RData file
  est_list <- file_name[[i]]
  load(file=paste0("est_files_v2/",est_list))
  unique_val <- sapply(res_est_fit, unique)
  smry <- summary_res_fit(unique_val$nsim, unique_val$S, unique_val$K, unique_val$m,
                          unique_val$ICC, unique_val$CAC, unique_val$theta)
  # Add the result to the data frame
  smry_totl <- rbind(smry_totl, smry)
}

#reshape 
#smry_totl to smry_totl_long

smry_totl_long <- smry_totl %>%
  pivot_longer(
    cols = ends_with("_c") | ends_with("_l"),
    names_to = c(".value", "type"),
    names_pattern = "(.*)_(c|l)"
  )

#Save the data frames to a csv file
write.csv(smry_totl_long, "summary_v2_MCSE.csv", row.names = FALSE)
