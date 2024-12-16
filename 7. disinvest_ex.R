setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Staircase_simstudy\\")
dis_dat <- read.csv("Disinvestment example/pmed.1002412.s002_Patient-level data.csv")
subdis_dat <- subset(dis_dat, site ==1 & study1==0, select = c(dhvicavglos, sw_step,no_we_exposure,
                        new_we_exposure, study1, site, index_ward,los))
##colnames(subdis_dat)[colnames(subdis_dat) == "no_we_exposure"] <- "we_exposure_0"
colnames(subdis_dat)[colnames(subdis_dat) == "new_we_exposure"] <- "treat"
colnames(subdis_dat)[colnames(subdis_dat) == "dhvicavglos"] <- "Y"
subdis_dat$time=subdis_dat$sw_step-7

original_values <- c(5, 1, 6, 2, 4, 3)
mapped_values <- c(1, 2, 3, 4, 5, 6)

subdis_dat$cluster <- mapped_values[match(subdis_dat$index_ward, original_values)]
#limit cluster-period cells to create a staircase design 
subdis_dat2 <- subset(subdis_dat, (cluster ==1 & (time==1 | time==2)) | (cluster ==2 & (time==2 | time==3)) 
                      | (cluster ==3 & (time==3 | time==4)) | (cluster ==4 & (time==4 | time==5))  
                      | (cluster ==5 & (time==5 | time==6)) | (cluster ==6 & (time==6 | time==7)) )

subdis_dat2$clustertime= subdis_dat2$cluster*subdis_dat2$time
sort(unique(subdis_dat2$clustertime))

original_values2 <- c( 1 , 2 , 4 , 6 , 9 ,12 ,16 ,20, 25 ,30, 36, 42)
mapped_values2 <- c(1  ,2  ,3  ,4,5  ,6  ,7 , 8  ,9 ,10 ,11 ,12)
subdis_dat2$clustper <- mapped_values2[match(subdis_dat2$clustertime, original_values2)]
library(dplyr)
count_summary <- subdis_dat2 %>%
  group_by(cluster, time) %>%  # Group by 'time' and 'cluster'
  summarise(count = n())   
count_summary
#generate results
dat<-subdis_dat2 
#95% CI HH Cat
est_trt_HH_c-1.96*(se_trt_HH_c)
est_trt_HH_c+1.96*(se_trt_HH_c)
#95% CI BE Cat
est_trt_BE_c-1.96*(se_trt_BE_c)
est_trt_BE_c+1.96*(se_trt_BE_c)
#95% CI HH Lin
est_trt_HH_l-1.96*(se_trt_HH_l)
est_trt_HH_l+1.96*(se_trt_HH_l)
#95% CI BE Lin
est_trt_BE_l-1.96*(se_trt_BE_l)
est_trt_BE_l+1.96*(se_trt_BE_l)

#95% CI HH Cat SAT
tstat_Sat_HH_c <- qt(alpha, adj_ddf_Sat_HH_c)
est_trt_HH_c-tstat_Sat_HH_c*adj_se_Sat_HH_c
est_trt_HH_c+tstat_Sat_HH_c*adj_se_Sat_HH_c
#95% CI BE Cat SAT
tstat_Sat_BE_c <- qt(alpha, adj_ddf_Sat_BE_c)
est_trt_BE_c-tstat_Sat_BE_c*adj_se_Sat_BE_c
est_trt_BE_c+tstat_Sat_BE_c*adj_se_Sat_BE_c
#95% CI HH Lin SAT
tstat_Sat_HH_l <- qt(alpha, adj_ddf_Sat_HH_l)
est_trt_HH_l-tstat_Sat_HH_l*adj_se_Sat_HH_l
est_trt_HH_l+tstat_Sat_HH_l*adj_se_Sat_HH_l
#95% CI BE Lin SAT
tstat_Sat_BE_l <- qt(alpha, adj_ddf_Sat_BE_l)
est_trt_BE_l-tstat_Sat_BE_l*adj_se_Sat_BE_l
est_trt_BE_l+tstat_Sat_BE_l*adj_se_Sat_BE_l

#95% CI HH Cat KR
tstat_KR_HH_c <- qt(alpha, adj_ddf_KR_HH_c)
est_trt_HH_c-tstat_KR_HH_c*adj_se_KR_HH_c
est_trt_HH_c+tstat_KR_HH_c*adj_se_KR_HH_c
#95% CI BE Cat KR
tstat_KR_BE_c <- qt(alpha, adj_ddf_KR_BE_c)
est_trt_BE_c-tstat_KR_BE_c*adj_se_KR_BE_c
est_trt_BE_c+tstat_KR_BE_c*adj_se_KR_BE_c
#95% CI HH Lin KR
tstat_KR_HH_l <- qt(alpha, adj_ddf_KR_HH_l)
est_trt_HH_l-tstat_KR_HH_l*adj_se_KR_HH_l
est_trt_HH_l+tstat_KR_HH_l*adj_se_KR_HH_l
#95% CI BE Lin KR
tstat_KR_BE_l <- qt(alpha, adj_ddf_KR_BE_l)
est_trt_BE_l-tstat_KR_BE_l*adj_se_KR_BE_l
est_trt_BE_l+tstat_KR_BE_l*adj_se_KR_BE_l
