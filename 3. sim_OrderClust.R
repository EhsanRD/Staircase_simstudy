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
nsim_values <- 10
S_values <- c(4,10)
K_values <-c(1, 5, 10)  
m_values <- c(10, 50, 100)  
ICC_values <-c(0.01,0.05)#c(0.01, 0.05, 0.1, 0.2) 
CAC_values <- 1#c(1, 0.95, 0.8, 0.5) 
theta_values <- c(0, 0.15) 
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
  c(10, 5, 100),
  c(4, 10, 100),
  c(4, 5, 100),
  c(10, 5, 50),
  c(10, 1, 100),
  c(4, 10, 50),
  c(10, 10, 100),
  c(4, 5, 50),
  c(10, 10, 50),
  c(4, 1, 100),
  c(10, 1, 50),
  c(10, 5, 10),
  c(4, 1, 50),
  c(4, 10, 10),
  c(4, 5, 10),
  c(4, 1, 10),
  c(10, 10, 10),
  c(10, 1, 10)
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
# Initialize an empty list to store results
all_res_lst <- list()
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
                         "fitHHmodelSC", "fitBEmodelSC", "pow","VarSClin" ,"VarSCcat"))
  
  res_lst <- parLapply(clust, sim_lst, run_all_sim)
  stopCluster(clust)
  
  # Append the results to the list
  all_res_lst <- c(all_res_lst, res_lst)
}

# Combine the results into a single data frame
res_df <- do.call(rbind, lapply(all_res_lst, function(x) x[[1]]))
resall_df <- res_df %>%
  mutate(Biastrt_HH = mesttrt_HH - theta, Biastrt_BE = mesttrt_BE - theta,
         BiasICC_HH = mestICC_HH - ICC,
         BiasICC_BE = mestICC_BE - ICC,
         BiasCAC_BE = mestCAC_BE - CAC)

resall_df$startmin <- unlist(lapply(all_res_lst, function(x) x[[2]]))
resall_df$endmin <- unlist(lapply(all_res_lst, function(x) x[[3]]))
resall_df$durtmin <- unlist(lapply(all_res_lst, function(x) x[[4]]))
resest_df <- do.call(rbind, lapply(all_res_lst, function(x) x[[5]]))

# Save the data frame to a CSV file
write.csv(resall_df, "testing1_parLapply_sim_res_18mainconfigs.csv", row.names = FALSE)
write.csv(resest_df, "testing2_parLapply_sim_res_18mainconfigs.csv", row.names= FALSE)
          