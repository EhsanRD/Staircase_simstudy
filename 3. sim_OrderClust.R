rm(list=ls())
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
#run the simulation for one single setting

# Your parameter vectors
nsim_values <- 2
S_values <-c(4,10)
K_values <- c(1, 5, 10)  
m_values <- 10#c(10, 50, 100)  
ICC_values <-c(0.01, 0.05, 0.1, 0.2) 
CAC_values <- 0.5#c(1, 0.95, 0.8, 0.5) 
theta_values <- 0#c(0, 0.15) 
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
#Function to filter combinations based on the specified order
filter_combinations <- function(order) {
  #filter the complete set of parameter combinations based on a specific order of S, K, and m values
  filter(param_combinations, S == order[1] & K == order[2] & m == order[3])
}

# Specify the desired order
desired_order <- list(
  # c(10, 5, 100),
  #  c(4, 10, 100),
  # c(4, 5, 100),
  # c(10, 5, 50),
  # c(10, 1, 100),
  # c(4, 10, 50),
  # c(10, 10, 100),
  # c(4, 5, 50),
  # c(10, 10, 50),
  # c(4, 1, 100),
  # c(10, 1, 50),
  c(10, 5, 10),
  # c(4, 1, 50),
  c(4, 10, 10),
  c(4, 5, 10),
  c(4, 1, 10),
  c(10, 10, 10),
  c(10, 1, 10)
)
run_all_sim <- function(params) {
  start_time <- format(Sys.time(), "%a %b %d %X %Y %Z")
  timeStart<-Sys.time()
  reslt <- sim_res_fit(params$nsim, params$S, params$K, params$m, params$ICC, params$CAC, params$theta,params$typ)
  end_time <- format(Sys.time(), "%a %b %d %X %Y %Z")
  timeEnd<-Sys.time()
  durt <- as.numeric(difftime(timeEnd, timeStart, units = "mins"))
  tim_df <- cbind(params$nsim, params$S, params$K, params$m, params$ICC, params$CAC, params$theta,as.character(params$typ),start_time,end_time,durt)
  return(list(reslt,tim_df))
}
# Initialize an empty list to store results
all_est_lst <- list()
numCores <- detectCores()
# Loop through each desired order
for (gpparams in desired_order) {
  
  #Filter combinations based on the specific order
  filtered_combinations <- filter_combinations(gpparams)
  
  #Split the filtered combinations into a list of data frames
  sim_lst <- split(filtered_combinations, seq_len(nrow(filtered_combinations)))
  
  # Perform parallel computation for the current set of combinations
  clust <- makeCluster(numCores)
  clusterEvalQ(clust, {
    library(MASS)
    library(lme4)
    library(glmmTMB)
    library(pbkrtest)
    library(tidyr)
    library(parameters)
    library(lmerTest)
    library(dplyr)
  })
  clusterExport(clust, c("sim_res_fit", "fitmodels", "gen_dat", "SCdesmat",
                         "fitHHmodelSC", "fitBEmodelSC", "pow" ,"VarSCcat","VarSClin"))
  
  est_lst <- parLapply(clust, sim_lst, run_all_sim)
  stopCluster(clust)
  
  # Append the results to the list
  all_est_lst <- c(all_est_lst, est_lst)
}


# Combine the results into a single data frame
est_df <- do.call(rbind, lapply(all_est_lst, function(x) x[[1]]))
#write.csv(est_df, "all_estimates.csv", row.names = FALSE)

sve_t <- do.call(rbind, lapply(all_est_lst, function(x) x[[2]]))
adj_names <- c("nsim", "S", "K", "m", "ICC", "CAC", "theta", "type","start_time","end_time","durt")
colnames(sve_t) <- adj_names
#write.csv(sve_t, "time_records.csv", row.names = FALSE)


