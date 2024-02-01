library(dplyr)
########################################
#run the whole simulation & get results#
########################################

#run the simulation for one single setting

#sim_res(nsim, S, K, m, ICC, CAC, theta)

pow(VarSCcat(10, 1, 100, 0.01, 0.8),0)
reslt <- sim_res_fit(1000, 10, 1, 100, 0.01, 0.8,0)
#reslt <- sim_res_fit(1000, 10, 1, 100, 0.01, 0.95,0)

saveRDS(reslt, file = paste0('results/', 'par_emp_vals', '.RDS'))

vals <- readRDS(file = 'results/par_emp_vals.RDS')

perfs <- vals %>%
  mutate(Bias_trt_HH=m_est_trt_HH-theta,Bias_trt_BE=m_est_trt_BE-theta,
                    Bias_ICC_HH=m_est_ICC_HH-ICC,
                    Bias_ICC_BE=m_est_ICC_BE-ICC, 
                    Bias_CAC_BE=m_est_CAC_BE-CAC)



vals <- readRDS(file = 'results/par_emp_vals.RDS')
vals2 <- readRDS(file = 'results/par_emp_vals2.RDS')
