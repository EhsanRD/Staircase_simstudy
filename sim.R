
library(dplyr)
library(parallel)
library(purrr)
library(foreach)
library(doParallel)
########################################
#run the whole simulation & get results#
#MU00212779
########################################
source('1. functions_sim.R')
source('2. gd_fit_res.R')

# Your parameter vectors
nsim_values <- 5
S_values <- c(4,10)
K_values <-c(1, 5, 10)  
m_values <- c(10, 50, 100)  
ICC_values <-c(0.01, 0.05, 0.1, 0.2) 
CAC_values <- 0.5#c(1, 0.95, 0.8, 0.5) 
theta_values <-c(0, 0.15) 
type_values <- c("cat","lin")

# Create a list of all combinations of parameter values
param_combinations <- expand.grid(
  nsim = nsim_values,
  S = S_values,
  K = K_values,
  m = m_values,
  ICC = ICC_values,
  CAC = CAC_values,
  theta = theta_values,
  typ=type_values
)

#Function to run sim_res_fit for a specific combination of parameters
#The duration time for each configuration might include the overhead of setting up and managing the parallel processes.
#Note: the durt is not accurate!
run_all_sim <- function(params) {
  start_time <- format(Sys.time(), "%a %b %d %X %Y %Z")
  timeStart<-Sys.time()
  reslt <- sim_res_fit(params$nsim, params$S, params$K, params$m, params$ICC, params$CAC, params$theta,params$typ)
  end_time <- format(Sys.time(), "%a %b %d %X %Y %Z")
  timeEnd<-Sys.time()
  durt <- as.numeric(difftime(timeEnd, timeStart, units = "mins"))
  reslt_main<- reslt[[1]]
  reslt_est <- reslt[[2]]
  return(list(reslt_main,start_time,end_time,durt,reslt_est))
}
#split the data frame param_combinations into a list of data frames
#based on the values in the sequence from 1 to the number of rows in param_combinations
sim_lst <- split(param_combinations, seq(nrow(param_combinations)))

#Using sockets with parLapply: The socket approach launches a new version of R on each core. dept.stat.lsa.umich.edu
numCores <- detectCores()
#system.time({

clust <- makeCluster(numCores)
# Load necessary packages in the worker processes
clusterEvalQ(clust, { # Load any other necessary packages here
  library(MASS)
  library(lme4)
  library(glmmTMB)
  library(pbkrtest)
  library(tidyr)
  library(parameters)
  library(lmerTest)
  library(dplyr)
})

clusterExport(clust, c("sim_res_fit",
                       "fitmodels","gen_dat","SCdesmat","fitHHmodelSC",
                       "fitBEmodelSC","pow","VarSCcat")) # Export necessary functions to the cluster
set.seed(958624)
#Perform parallel computation
res_lst <- parLapply(clust, sim_lst, run_all_sim)


#avoid running R sessions in different clusters in the background.
stopCluster(clust)
#})
# set.seed(958624)
# res_lst <- lapply(sim_lst, run_all_sim)
res_plst<- lapply(res_lst, function(x) x[[1]])

# Combine the results into a single data frame
res_df <- do.call(rbind, res_plst)
resall_df<- res_df %>%
  mutate(Biastrt_HH=mesttrt_HH-theta,Biastrt_BE=mesttrt_BE-theta,
         BiasICC_HH=mestICC_HH-ICC,
         BiasICC_BE=mestICC_BE-ICC,
         BiasCAC_BE=mestCAC_BE-CAC)

resall_df$startmin <- unlist(lapply(res_lst, function(x) x[[2]]))
resall_df$endmin <- unlist(lapply(res_lst, function(x) x[[3]]))
resall_df$durtmin <- unlist(lapply(res_lst, function(x) x[[4]]))
resest_df <- do.call(rbind, lapply(res_lst, function(x) x[[5]]))


#Save the data frames to a csv file
write.csv(resall_df, "result_overall.csv", row.names = FALSE)
write.csv(resest_df, "result_indivdiual_estimates.csv", row.names = FALSE)



