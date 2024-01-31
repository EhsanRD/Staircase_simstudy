library(dplyr)
########################################
#run the whole simulation & get results#
########################################

#run the simulation for one single setting

#sim_res(nsim, S, K, m, ICC, CAC, theta)

reslt <- sim_res(50, 10,5,20, 0.05, 0.8, 0.15)
saveRDS(reslt, file = paste0('results/', 'par_emp_vals', '.RDS'))

vals <- readRDS(file = 'results/par_emp_vals.RDS')

perfs <- vals %>%
  mutate(Bias_trt_HH=m_est_HH-theta,Bias_trt_BE=m_est_BE-theta,
                    Bias_ICC_HH=m_estICC_HH-ICC,
                    Bias_ICC_BE=m_estICC_BE-ICC, 
                    Bias_CAC_BE=m_estCAC_BE-CAC)


