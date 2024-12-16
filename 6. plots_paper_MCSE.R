setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Staircase_simstudy\\")
source('5. dat_plots.R')
setwd("G:\\Shared drives\\Ehsan PhD work\\paper_3\\Figures_v2_MCSE\\")
############POWER
ylow_1 <- 0
yup_1  <- 100
fval_1 <- cal_val(ylow_1, yup_1)
ylow_2 <- 0
yup_2  <- 50 #20 CS
fval_2 <- cal_val(ylow_2, yup_2)


# Adjust power values
smry_df$pow = smry_df$pow * 100

smry_df$pow_HH_l=(smry_df$pow_HH-2*smry_df$MCSEpow_HH)*100
smry_df$pow_HH_u=(smry_df$pow_HH+2*smry_df$MCSEpow_HH)*100
smry_df$pow_HH = smry_df$pow_HH * 100

smry_df$pow_BE_l=(smry_df$pow_BE-2*smry_df$MCSEpow_BE)*100
smry_df$pow_BE_u=(smry_df$pow_BE+2*smry_df$MCSEpow_BE)*100
smry_df$pow_BE = smry_df$pow_BE * 100

smry_df$powKR_l_HH=(smry_df$powKR_HH-2*smry_df$MCSEpowKR_HH)*100
smry_df$powKR_u_HH=(smry_df$powKR_HH+2*smry_df$MCSEpowKR_HH)*100
smry_df$powKR_HH = smry_df$powKR_HH * 100

smry_df$powKR_l_BE=(smry_df$powKR_BE-2*smry_df$MCSEpowKR_BE)*100
smry_df$powKR_u_BE=(smry_df$powKR_BE+2*smry_df$MCSEpowKR_BE)*100
smry_df$powKR_BE = smry_df$powKR_BE * 100

smry_df$powSat_l_HH=(smry_df$powSat_HH-2*smry_df$MCSEpowSat_HH)*100
smry_df$powSat_u_HH=(smry_df$powSat_HH+2*smry_df$MCSEpowSat_HH)*100
smry_df$powSat_HH = smry_df$powSat_HH * 100

smry_df$powSat_l_BE=(smry_df$powSat_BE-2*smry_df$MCSEpowSat_BE)*100
smry_df$powSat_u_BE=(smry_df$powSat_BE+2*smry_df$MCSEpowSat_BE)*100
smry_df$powSat_BE = smry_df$powSat_BE * 100

# replace KR (or Sat) correction with uncorrected power estimation for S<10 or K<10

smry_df$pow_HH_2_l<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$powKR_l_HH,smry_df$pow_HH_l)
smry_df$pow_HH_2_u<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$powKR_u_HH,smry_df$pow_HH_u)
smry_df$pow_HH_2<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$powKR_HH,smry_df$pow_HH)

smry_df$pow_BE_2_l<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$powKR_l_BE,smry_df$pow_BE_l)
smry_df$pow_BE_2_u<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$powKR_u_BE,smry_df$pow_BE_u)
smry_df$pow_BE_2<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$powKR_BE,smry_df$pow_BE)

#CS correctly specified
smry_df$pow_HH_2_l_CS<-ifelse(smry_df$CAC == 1, smry_df$pow_HH_2_l,  -999) 
smry_df$pow_HH_2_u_CS<-ifelse(smry_df$CAC == 1, smry_df$pow_HH_2_u,  -999) 
smry_df$pow_HH_2_CS <- ifelse(smry_df$CAC == 1, smry_df$pow_HH_2,  -999) 

smry_df$pow_BE_2_CS_l <- ifelse(smry_df$CAC != 1,  smry_df$pow_BE_2_l, -999) 
smry_df$pow_BE_2_CS_u <- ifelse(smry_df$CAC != 1,  smry_df$pow_BE_2_u, -999) 
smry_df$pow_BE_2_CS <- ifelse(smry_df$CAC != 1,  smry_df$pow_BE_2, -999) 

#ICS correctly specified
smry_df$pow_HH_2_l_ICS<-ifelse(smry_df$CAC != 1, smry_df$pow_HH_2_l,  -999) 
smry_df$pow_HH_2_u_ICS<-ifelse(smry_df$CAC != 1, smry_df$pow_HH_2_u,  -999) 
smry_df$pow_HH_2_ICS <- ifelse(smry_df$CAC != 1, smry_df$pow_HH_2,  -999) 

smry_df$pow_BE_2_ICS_l <- ifelse(smry_df$CAC == 1,  smry_df$pow_BE_2_l, -999) 
smry_df$pow_BE_2_ICS_u <- ifelse(smry_df$CAC == 1,  smry_df$pow_BE_2_u, -999) 
smry_df$pow_BE_2_ICS <- ifelse(smry_df$CAC == 1,  smry_df$pow_BE_2, -999) 

smry_df$pow_HH_3_l<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$powSat_l_HH,smry_df$pow_HH_l)
smry_df$pow_HH_3_u<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$powSat_u_HH,smry_df$pow_HH_u)
smry_df$pow_HH_3<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$powSat_HH,smry_df$pow_HH)

smry_df$pow_BE_3_l<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$powSat_l_BE,smry_df$pow_BE_l)
smry_df$pow_BE_3_u<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$powSat_u_BE,smry_df$pow_BE_u)
smry_df$pow_BE_3<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$powSat_BE,smry_df$pow_BE)

#CS correctly specified
smry_df$pow_HH_3_l_CS<-ifelse(smry_df$CAC == 1, smry_df$pow_HH_3_l,  -999) 
smry_df$pow_HH_3_u_CS<-ifelse(smry_df$CAC == 1, smry_df$pow_HH_3_u,  -999) 
smry_df$pow_HH_3_CS <- ifelse(smry_df$CAC == 1, smry_df$pow_HH_3,  -999) 

smry_df$pow_BE_3_CS_l <- ifelse(smry_df$CAC != 1,  smry_df$pow_BE_3_l, -999) 
smry_df$pow_BE_3_CS_u <- ifelse(smry_df$CAC != 1,  smry_df$pow_BE_3_u, -999) 
smry_df$pow_BE_3_CS <- ifelse(smry_df$CAC != 1,  smry_df$pow_BE_3, -999) 

#ICS correctly specified
smry_df$pow_HH_3_l_ICS<-ifelse(smry_df$CAC != 1, smry_df$pow_HH_3_l,  -999) 
smry_df$pow_HH_3_u_ICS<-ifelse(smry_df$CAC != 1, smry_df$pow_HH_3_u,  -999) 
smry_df$pow_HH_3_ICS <- ifelse(smry_df$CAC != 1, smry_df$pow_HH_3,  -999) 

smry_df$pow_BE_3_ICS_l <- ifelse(smry_df$CAC == 1,  smry_df$pow_BE_3_l, -999) 
smry_df$pow_BE_3_ICS_u <- ifelse(smry_df$CAC == 1,  smry_df$pow_BE_3_u, -999) 
smry_df$pow_BE_3_ICS <- ifelse(smry_df$CAC == 1,  smry_df$pow_BE_3, -999) 

# Reorder the CAC factor levels
smry_df$CAC <- factor(smry_df$CAC, levels = c(1, 0.95, 0.8, 0.5))
#smry_df <- smry_df[order(smry_df$CAC), ]

# Plotting function
plot_nested_loop <- function(output_filename, data_subset, y_variable, y_variable_name,legend_labels) {
  png(file = output_filename, width = 1200, height = 900) 
  
  nested_loop_plot(resdf = data_subset, 
                   x = "m",
                   grid_rows = "CAC", steps = c("ICC","K","S"),
                   steps_y_base = ifelse(y_variable == "Power (%)",fval_1$syb,fval_2$syb), 
                   steps_y_height =ifelse(y_variable == "Power (%)",fval_1$syh,fval_2$syh),  
                   steps_y_shift = ifelse(y_variable == "Power (%)",fval_1$sysh,fval_2$sysh), 
                   x_name = "m", y_name = y_variable_name,
                   spu_x_shift = 60,line_size = 1, line_alpha = 1,
                   point_shapes =c(20, 19, 20,20,19,20,15),point_size =2,
                   colors = c("#5D3FD3","#F08080","#FF6B6B","#F08080","#A7E9DF","#4ECDC4","#A7E9DF"), #CS
                   #colors = c("#A7E9DF","#4ECDC4","#A7E9DF","#F08080","#FF6B6B","#F08080","#5D3FD3"), #ICS
                    line_linetypes = c(1,1,1,1,1,1),
                   steps_color = "dimgrey",
                   steps_values_annotate = TRUE, steps_annotation_size = 3,
                   steps_annotation_nudge =0.4,
                   steps_annotation_color = "dimgrey", ylim = c(ylow_1, ifelse(y_variable == "Power (%)", yup_1, yup_2)),
                   hline_intercept = c(ylow_1,ifelse(y_variable == "Power (%)", NA, 5) ,
                                       ifelse(y_variable == "Power (%)", yup_1, 0)),
                   
                   y_expand_add = c(ifelse(y_variable == "Power (%)",fval_1$yexp1,fval_2$yexp1), 
                                    ifelse(y_variable == "Power (%)",fval_1$yexp2,fval_2$yexp2)),
                   legend_labels = legend_labels,base_size = 16,
                   steps_names = c("ICC","K", "S"),
                   post_processing = list(
                     add_custom_theme = list(
                       axis.text.x = element_text(angle = -90, vjust = 0.5, size = 10,face = "bold"), 
                       axis.text.y = element_text(face = "bold", size = 14), 
                       strip.text = element_text(face = "bold", size = 14),
                       axis.title.y = element_text(face = "bold", size = 16)
                     ))
      )
  
}
# Plotting for Type I error (%) using Sat correction  
#categorical #TYPE I ERROR #SAT
smry_df_sub <- smry_df[(smry_df$theta ==0 & 
                          smry_df$`Time effect`=='categorical'), c("CAC","ICC","S", "K", "m", "pow_HH_3_l_CS","pow_HH_3_CS","pow_HH_3_u_CS",
                                                                   "pow_BE_3_CS_l","pow_BE_3_CS","pow_BE_3_CS_u")]

plot_nested_loop(output_filename = "Fig2_estTypeISat_CatbyCACs.png", 
                 data_subset = smry_df_sub, y_variable = "Type I error (%)", y_variable_name="Type I error (%)",
                 legend_labels = c("+ 2 * MCSE","Exch model W/ Sat","- 2 * MCSE",
                                   "+ 2 * MCSE", "BE model W/ Sat", "- 2 * MCSE",
                                   "Theoretical value"))
dev.off()
write.csv(smry_df_sub, "smry_df_sub_Type1CS.csv", row.names = FALSE)


#categorical #POWER #SAT
smry_df_sub <- smry_df[(smry_df$theta > 0 & 
                          smry_df$`Time effect`=='categorical'), c("CAC","ICC","S", "K", "m","pow","pow_HH_3_l_CS","pow_HH_3_CS","pow_HH_3_u_CS",
                                                                   "pow_BE_3_CS_l","pow_BE_3_CS","pow_BE_3_CS_u")]
plot_nested_loop(output_filename = "Fig3_estpowSat_CatbyCACs.png", 
                 data_subset = smry_df_sub, y_variable = "Power (%)", y_variable_name="Power (%)",
                 legend_labels =c("Theoretical value","+ 2 * MCSE","Exch model W/ Sat","- 2 * MCSE",
                                  "+ 2 * MCSE", "BE model W/ Sat", "- 2 * MCSE"))
dev.off()
write.csv(smry_df_sub, "smry_df_sub_powCS.csv", row.names = FALSE)
#categorical #POWER #SAT #ICS
smry_df_sub <- smry_df[(smry_df$theta ==0 &
                          smry_df$`Time effect`=='categorical'), c("CAC","ICC","S", "K", "m","pow_BE_3_ICS_l","pow_BE_3_ICS","pow_BE_3_ICS_u",
                                                                   "pow_HH_3_l_ICS","pow_HH_3_ICS","pow_HH_3_u_ICS")]
plot_nested_loop(output_filename = "Fig4_estTypeISat_CatbyCACs.png", 
                 data_subset = smry_df_sub, y_variable = "Type I error (%)", y_variable_name="Type I error (%)",
                 legend_labels = c("+ 2 * MCSE", "BE model W/ Sat", "- 2 * MCSE",
                                   "+ 2 * MCSE","Exch model W/ Sat","- 2 * MCSE", "Theoretical value"))
dev.off()
write.csv(smry_df_sub, "smry_df_sub_Type1ICS.csv", row.names = FALSE)


#Plotting for TypeI error (%) using KR correction  
#categorical #TYPE I ERROR # KR 
smry_df_sub <- smry_df[(smry_df$theta ==0 & 
                          smry_df$`Time effect`=='categorical'), c("CAC","ICC","S", "K", "m", "pow_HH_2_l_CS","pow_HH_2_CS","pow_HH_2_u_CS",
                                                                   "pow_BE_2_CS_l","pow_BE_2_CS","pow_BE_2_CS_u")]
plot_nested_loop(output_filename = "FigS4_estTypeIKR_CatbyCACs.png", 
                 data_subset = smry_df_sub, y_variable = "Type I error (%)", y_variable_name="Type I error (%)",
                 legend_labels = c("+ 2 * MCSE","Exch model W/ KR","- 2 * MCSE",
                                   "+ 2 * MCSE", "BE model W/ KR", "- 2 * MCSE",
                                   "Theoretical value"))
dev.off()
#categorical #POWER # KR 
smry_df_sub <- smry_df[(smry_df$theta > 0 & 
                          smry_df$`Time effect`=='categorical'), c("CAC","ICC","S", "K", "m","pow","pow_HH_2_l_CS","pow_HH_2_CS","pow_HH_2_u_CS",
                                                                   "pow_BE_2_CS_l","pow_BE_2_CS","pow_BE_2_CS_u")]
plot_nested_loop(output_filename = "FigS5_estpowKR_CatbyCACs.png", 
                 data_subset = smry_df_sub, y_variable = "Power (%)", y_variable_name="Power (%)",
                 legend_labels =c("Theoretical value","+ 2 * MCSE","Exch model W/ KR","- 2 * MCSE",
                                  "+ 2 * MCSE", "BE model W/ KR", "- 2 * MCSE"))
dev.off()
#categorical #TYPE I ERROR # KR #ICS
smry_df_sub <- smry_df[(smry_df$theta ==0 & 
                          smry_df$`Time effect`=='categorical'), c("CAC","ICC","S", "K", "m","pow_BE_2_ICS_l","pow_BE_2_ICS","pow_BE_2_ICS_u",
                                                                   "pow_HH_2_l_ICS","pow_HH_2_ICS","pow_HH_2_u_ICS")]
plot_nested_loop(output_filename = "FigS6_estTypeIKR_CatbyCACs.png", 
                 data_subset = smry_df_sub, y_variable = "Type I error (%)", y_variable_name="Type I error (%)",
                 legend_labels = c("+ 2 * MCSE", "BE model W/ KR", "- 2 * MCSE",
                                   "+ 2 * MCSE","Exch model W/ KR","- 2 * MCSE", "Theoretical value"))
dev.off()
#linear #TYPE I ERROR # SAT
smry_df_sub <- smry_df[(smry_df$theta ==0 & 
                          smry_df$`Time effect`=='linear'), c("CAC","ICC","S", "K", "m", "pow_HH_3_l_CS","pow_HH_3_CS","pow_HH_3_u_CS",
                                                              "pow_BE_3_CS_l","pow_BE_3_CS","pow_BE_3_CS_u")]
plot_nested_loop(output_filename = "FigS12_estTypeISat_LinbyCACs.png", 
                 data_subset = smry_df_sub, y_variable = "Type I error (%)", y_variable_name="Type I error (%)",
                 legend_labels = c("+ 2 * MCSE","Exch model W/ Sat","- 2 * MCSE",
                                   "+ 2 * MCSE", "BE model W/ Sat", "- 2 * MCSE",
                                   "Theoretical value"))
dev.off()
#linear #POWER # SAT
smry_df_sub <- smry_df[(smry_df$theta > 0 & 
                          smry_df$`Time effect`=='linear'), c("CAC","ICC","S", "K", "m","pow","pow_HH_3_l_CS","pow_HH_3_CS","pow_HH_3_u_CS",
                                                              "pow_BE_3_CS_l","pow_BE_3_CS","pow_BE_3_CS_u")]
plot_nested_loop(output_filename = "FigS13_estpowSat_LinbyCACs.png", 
                 data_subset = smry_df_sub, y_variable = "Power (%)", y_variable_name="Power (%)",
                 legend_labels =c("Theoretical value","+ 2 * MCSE","Exch model W/ Sat","- 2 * MCSE",
                                  "+ 2 * MCSE", "BE model W/ Sat", "- 2 * MCSE"))
dev.off()
#linear #TYPE I ERROR # SAT  #ICS
smry_df_sub <- smry_df[(smry_df$theta ==0 & 
                          smry_df$`Time effect`=='linear'), c("CAC","ICC","S", "K", "m", "pow_BE_3_ICS_l","pow_BE_3_ICS","pow_BE_3_ICS_u",
                                                              "pow_HH_3_l_ICS","pow_HH_3_ICS","pow_HH_3_u_ICS")]
plot_nested_loop(output_filename = "FigS14_estTypeISat_LinbyCACs.png", 
                 data_subset = smry_df_sub, y_variable = "Type I error (%)", y_variable_name="Type I error (%)",
                 legend_labels = c("+ 2 * MCSE", "BE model W/ Sat", "- 2 * MCSE",
                                   "+ 2 * MCSE","Exch model W/ Sat","- 2 * MCSE", "Theoretical value"))
dev.off()

#linear #TYPE I ERROR # KR
smry_df_sub <- smry_df[(smry_df$theta ==0 & 
                          smry_df$`Time effect`=='linear'), c("CAC","ICC","S", "K", "m",  "pow_HH_2_l_CS","pow_HH_2_CS","pow_HH_2_u_CS",
                                                              "pow_BE_2_CS_l","pow_BE_2_CS","pow_BE_2_CS_u")]
plot_nested_loop(output_filename = "FigS16_estTypeIKR_LinbyCACs.png", 
                 data_subset = smry_df_sub, y_variable = "Type I error (%)", y_variable_name="Type I error (%)",
                 legend_labels =c("+ 2 * MCSE","Exch model W/ KR","- 2 * MCSE",
                                  "+ 2 * MCSE", "BE model W/ KR", "- 2 * MCSE",
                                  "Theoretical value"))
dev.off()
#linear #POWER # KR
smry_df_sub <- smry_df[(smry_df$theta > 0 & 
                          smry_df$`Time effect`=='linear'), c("CAC","ICC","S", "K", "m" ,"pow","pow_HH_2_l_CS","pow_HH_2_CS","pow_HH_2_u_CS",
                                                              "pow_BE_2_CS_l","pow_BE_2_CS","pow_BE_2_CS_u")]
plot_nested_loop(output_filename = "FigS17_estpowKR_LinbyCACs.png", 
                 data_subset = smry_df_sub, y_variable = "Power (%)", y_variable_name="Power (%)",
                 legend_labels =c("Theoretical value","+ 2 * MCSE","Exch model W/ KR","- 2 * MCSE",
                                  "+ 2 * MCSE", "BE model W/ KR", "- 2 * MCSE" ))
dev.off()

#linear #TYPE I ERROR # KR  #ICS
smry_df_sub <- smry_df[(smry_df$theta ==0 & 
                          smry_df$`Time effect`=='linear'), c("CAC","ICC","S", "K", "m", "pow_BE_2_ICS_l","pow_BE_2_ICS","pow_BE_2_ICS_u",
                                                              "pow_HH_2_l_ICS","pow_HH_2_ICS","pow_HH_2_u_ICS")]
plot_nested_loop(output_filename = "FigS18_estTypeIKR_LinbyCACs.png", 
                 data_subset = smry_df_sub, y_variable = "Type I error (%)", y_variable_name="Type I error (%)",
                 legend_labels = c("+ 2 * MCSE", "BE model W/ KR", "- 2 * MCSE",
                                   "+ 2 * MCSE","Exch model W/ KR","- 2 * MCSE", "Theoretical value"))
dev.off()
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
# treatment effect bias (mean) ICC bias (mean)
setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Staircase_simstudy\\")
source('5. dat_plots.R')
setwd("G:\\Shared drives\\Ehsan PhD work\\paper_3\\Figures_v2_MCSE\\")
ylow_1 <- -0.05 #appendix
yup_1  <- 0.05 #appendix
# ylow_1 <- -0.025 
# yup_1  <- 0.025 
fval_1 <- cal_val(ylow_1, yup_1)
ylow_2 <- -0.1
yup_2  <-0.05
fval_2 <- cal_val(ylow_2, yup_2)

smry_df$Biastrt_HH_l = smry_df$Biastrt_HH - 2*smry_df$MCSE_HH     
smry_df$Biastrt_HH_u = smry_df$Biastrt_HH + 2*smry_df$MCSE_HH 
smry_df$Biastrt_BE_l = smry_df$Biastrt_BE - 2*smry_df$MCSE_BE     
smry_df$Biastrt_BE_u = smry_df$Biastrt_BE + 2*smry_df$MCSE_BE

#CS correctly specified
smry_df$Biastrt_HH_l_CS<-ifelse(smry_df$CAC == 1, smry_df$Biastrt_HH_l,  -999) 
smry_df$Biastrt_HH_u_CS<-ifelse(smry_df$CAC == 1, smry_df$Biastrt_HH_u,  -999) 
smry_df$Biastrt_HH_CS <- ifelse(smry_df$CAC == 1, smry_df$Biastrt_HH,  -999) 

smry_df$Biastrt_BE_l_CS <- ifelse(smry_df$CAC != 1,  smry_df$Biastrt_BE_l, -999) 
smry_df$Biastrt_BE_u_CS <- ifelse(smry_df$CAC != 1,  smry_df$Biastrt_BE_u, -999) 
smry_df$Biastrt_BE_CS <- ifelse(smry_df$CAC != 1,  smry_df$Biastrt_BE, -999) 

#ICS correctly specified
smry_df$Biastrt_HH_l_ICS<-ifelse(smry_df$CAC != 1, smry_df$Biastrt_HH_l,  -999) 
smry_df$Biastrt_HH_u_ICS<-ifelse(smry_df$CAC != 1, smry_df$Biastrt_HH_u,  -999) 
smry_df$Biastrt_HH_ICS <- ifelse(smry_df$CAC != 1, smry_df$Biastrt_HH,  -999) 

smry_df$Biastrt_BE_l_ICS <- ifelse(smry_df$CAC == 1,  smry_df$Biastrt_BE_l, -999) 
smry_df$Biastrt_BE_u_ICS <- ifelse(smry_df$CAC == 1,  smry_df$Biastrt_BE_u, -999) 
smry_df$Biastrt_BE_ICS <- ifelse(smry_df$CAC == 1,  smry_df$Biastrt_BE, -999) 


# Reorder the CAC factor levels
smry_df$CAC <- factor(smry_df$CAC, levels = c(1, 0.95, 0.8, 0.5))

# Plotting function
plot_nested_loop <- function(output_filename, data_subset, y_variable,theta_val, y_variable_name,legend_labels) {
  png(file = output_filename, width = 1200, height = 900) 
  
  nested_loop_plot(resdf = data_subset, 
                   x = "m",
                   grid_rows = "CAC", steps = c("ICC","K","S"),
                   steps_y_base = ifelse(y_variable == "trt",fval_1$syb,fval_2$syb), 
                   steps_y_height =ifelse(y_variable == "trt",fval_1$syh,fval_2$syh),  
                   steps_y_shift = ifelse(y_variable == "trt",fval_1$sysh,fval_2$sysh), 
                   x_name = "m", y_name = y_variable_name,
                   spu_x_shift = 60,line_size = 1, line_alpha = 1,
                   point_shapes =c(20, 19, 20,20,19,20,15),point_size =2,
                   #colors = scales::brewer_pal(palette = "Set2"), #appendix
                   colors = c("#F08080","#FF6B6B","#F08080","#A7E9DF","#4ECDC4","#A7E9DF"), #CS
                   #colors = c("#A7E9DF","#4ECDC4","#A7E9DF","#F08080","#FF6B6B","#F08080","#5D3FD3"), #ICS
                   line_linetypes = c(1,1,1,1,1,1),
                   steps_color = "dimgrey",
                   steps_values_annotate = TRUE, steps_annotation_size = 3,
                   steps_annotation_nudge = 0.4,
                   steps_annotation_color = "dimgrey", ylim = c(ifelse(y_variable == "trt",ylow_1,ylow_2),
                                                                ifelse(y_variable == "trt",yup_1,yup_2)),
                   hline_intercept =0, #c(ifelse(y_variable == "trt",ylow_1,ylow_2), 
                   # ifelse(y_variable == "trt",yup_1,yup_2)),
                   y_expand_add = c(ifelse(y_variable == "trt",fval_1$yexp1,fval_2$yexp1), 
                                    ifelse(y_variable == "trt",fval_1$yexp2,fval_2$yexp2)), 
                   legend_name = ifelse(theta_val == "zero","Method","Method"), 
                   steps_names = c("ICC","K", "S"),
                   legend_labels = legend_labels,base_size = 16,
                   post_processing = list(
                     add_custom_theme = list(
                       axis.text.x = element_text(angle = -90, vjust = 0.5, size = 10,face = "bold"), 
                       axis.text.y = element_text(face = "bold", size = 14), 
                       strip.text = element_text(face = "bold", size = 14),
                       axis.title.y = element_text(face = "bold", size = 16)
                     ))
  )
  
}
# Plotting treatment effect bias trt & icc for each theta  
#categorical #CS
smry_df_sub <- smry_df[smry_df$theta >0  &
                         smry_df$`Time effect`=='categorical', c("CAC","ICC","S", "K", "m",
                                    "Biastrt_HH_l_CS","Biastrt_HH_CS","Biastrt_HH_u_CS",
                                    "Biastrt_BE_l_CS","Biastrt_BE_CS","Biastrt_BE_u_CS")]
plot_nested_loop(output_filename = "FigS2_esttrtbias_theta15_CatbyCACs.png",
                 data_subset = smry_df_sub, y_variable = "trt",theta_val="non-zero", y_variable_name= "Treatment effect bias (mean)",
                 legend_labels = c("+ 2 * MCSE","Exch model","- 2 * MCSE",
                                   "+ 2 * MCSE", "BE model", "- 2 * MCSE"))
dev.off()
#categorical #ICS
smry_df_sub <- smry_df[smry_df$theta >0  &
                         smry_df$`Time effect`=='categorical', c("CAC","ICC","S", "K", "m",
                                                                 "Biastrt_BE_l_ICS","Biastrt_BE_ICS","Biastrt_BE_u_ICS",
                                                                 "Biastrt_HH_l_ICS","Biastrt_HH_ICS","Biastrt_HH_u_ICS")]
plot_nested_loop(output_filename = "FigS3_esttrtbias_theta15_CatbyCACs.png",
                 data_subset = smry_df_sub, y_variable = "trt",theta_val="non-zero", y_variable_name= "Treatment effect bias (mean)",
                 legend_labels = c("+ 2 * MCSE","BE model","- 2 * MCSE",
                                   "+ 2 * MCSE", "Exch model", "- 2 * MCSE"))
dev.off()

#linear
smry_df_sub <- smry_df[smry_df$theta >0  &
                         smry_df$`Time effect`=='linear', c("CAC","ICC","S", "K", "m",
                                                            "Biastrt_HH_l_CS","Biastrt_HH_CS","Biastrt_HH_u_CS",
                                                            "Biastrt_BE_l_CS","Biastrt_BE_CS","Biastrt_BE_u_CS")]
plot_nested_loop(output_filename = "FigS20_esttrtbias_theta15_LinbyCACs.png",
                 data_subset = smry_df_sub, y_variable = "trt",theta_val="non-zero", y_variable_name= "Treatment effect bias (mean)",
                 legend_labels = c("+ 2 * MCSE","Exch model","- 2 * MCSE",
                                   "+ 2 * MCSE", "BE model", "- 2 * MCSE"))
dev.off()

#linear
smry_df_sub <- smry_df[smry_df$theta >0  &
                         smry_df$`Time effect`=='linear', c("CAC","ICC","S", "K", "m",
                                                            "Biastrt_BE_l_ICS","Biastrt_BE_ICS","Biastrt_BE_u_ICS",
                                                            "Biastrt_HH_l_ICS","Biastrt_HH_ICS","Biastrt_HH_u_ICS")]
plot_nested_loop(output_filename = "FigS21_esttrtbias_theta15_LinbyCACs.png",
                 data_subset = smry_df_sub, y_variable = "trt",theta_val="non-zero", y_variable_name= "Treatment effect bias (mean)",
                 legend_labels = c("+ 2 * MCSE","BE model","- 2 * MCSE",
                                    "+ 2 * MCSE", "Exch model", "- 2 * MCSE"))
dev.off()
# smry_df_sub <- smry_df[smry_df$theta>0 & smry_df$m!=50 & smry_df$K!=5 & (smry_df$ICC==0.05 | smry_df$ICC==0.2) & 
#                          smry_df$`Time effect`=='categorical' & smry_df$CAC==0.8, c("CAC","ICC","S", "K", "m", "Biastrt_HH","Biastrt_BE")]
# plot_nested_loop(output_filename = "Fig3_esttrtbias_theta15_CAC0.8.png", 
#                  data_subset = smry_df_sub, y_variable = "trt", theta_val="non-zero",y_variable_name= "Treatment effect bias (mean)",
#                  legend_labels = c("Exchangeable model",  "Block-exchangeable model"))
# dev.off()
# 
# smry_df_sub <- smry_df[smry_df$theta>0 & smry_df$m!=50 & smry_df$K!=5 & (smry_df$ICC==0.05 | smry_df$ICC==0.2) & 
#                          smry_df$`Time effect`=='categorical' & smry_df$CAC==1, c("CAC","ICC","S", "K", "m", "Biastrt_HH","Biastrt_BE")]
# plot_nested_loop(output_filename = "Fig3_esttrtbias_theta15_CAC1.png", 
#                  data_subset = smry_df_sub, y_variable = "trt", theta_val="non-zero",y_variable_name= "Treatment effect bias (mean)",
#                  legend_labels = c("Exchangeable model",  "Block-exchangeable model"))
# dev.off()

#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#coverage
setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Staircase_simstudy\\")
source('5. dat_plots.R')
setwd("G:\\Shared drives\\Ehsan PhD work\\paper_3\\Figures_v2_MCSE\\")
ylow <- 50
yup  <- 100
fval <- cal_val(ylow, yup)


# Adjust coverage values
smry_df$cov_HH_l=(smry_df$cov_HH-2*smry_df$MCSEcov_HH)*100
smry_df$cov_HH_u=(smry_df$cov_HH+2*smry_df$MCSEcov_HH)*100
smry_df$cov_HH = smry_df$cov_HH * 100

smry_df$cov_BE_l=(smry_df$cov_BE-2*smry_df$MCSEcov_BE)*100
smry_df$cov_BE_u=(smry_df$cov_BE+2*smry_df$MCSEcov_BE)*100
smry_df$cov_BE = smry_df$cov_BE * 100

smry_df$covKR_l_HH=(smry_df$covKR_HH-2*smry_df$MCSEcovKR_HH)*100
smry_df$covKR_u_HH=(smry_df$covKR_HH+2*smry_df$MCSEcovKR_HH)*100
smry_df$covKR_HH = smry_df$covKR_HH * 100

smry_df$covKR_l_BE=(smry_df$covKR_BE-2*smry_df$MCSEcovKR_BE)*100
smry_df$covKR_u_BE=(smry_df$covKR_BE+2*smry_df$MCSEcovKR_BE)*100
smry_df$covKR_BE = smry_df$covKR_BE * 100

smry_df$covSat_l_HH=(smry_df$covSat_HH-2*smry_df$MCSEcovSat_HH)*100
smry_df$covSat_u_HH=(smry_df$covSat_HH+2*smry_df$MCSEcovSat_HH)*100
smry_df$covSat_HH = smry_df$covSat_HH * 100

smry_df$covSat_l_BE=(smry_df$covSat_BE-2*smry_df$MCSEcovSat_BE)*100
smry_df$covSat_u_BE=(smry_df$covSat_BE+2*smry_df$MCSEcovSat_BE)*100
smry_df$covSat_BE = smry_df$covSat_BE * 100


# replace KR (or Sat) correction with uncorrected cover estimation for S<10 or K<10
smry_df$cov_HH_2_l<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$covKR_l_HH,smry_df$cov_HH_l)
smry_df$cov_HH_2_u<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$covKR_u_HH,smry_df$cov_HH_u)
smry_df$cov_HH_2<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$covKR_HH,smry_df$cov_HH)

smry_df$cov_BE_2_l<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$covKR_l_BE,smry_df$cov_BE_l)
smry_df$cov_BE_2_u<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$covKR_u_BE,smry_df$cov_BE_u)
smry_df$cov_BE_2<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$covKR_BE,smry_df$cov_BE)

#ICS correctly specified
smry_df$cov_HH_2_l_ICS<-ifelse(smry_df$CAC != 1, smry_df$cov_HH_2_l,  -999) 
smry_df$cov_HH_2_u_ICS<-ifelse(smry_df$CAC != 1, smry_df$cov_HH_2_u,  -999) 
smry_df$cov_HH_2_ICS <- ifelse(smry_df$CAC != 1, smry_df$cov_HH_2,  -999) 

smry_df$cov_BE_2_ICS_l <- ifelse(smry_df$CAC == 1,  smry_df$cov_BE_2_l, -999) 
smry_df$cov_BE_2_ICS_u <- ifelse(smry_df$CAC == 1,  smry_df$cov_BE_2_u, -999) 
smry_df$cov_BE_2_ICS <- ifelse(smry_df$CAC == 1,  smry_df$cov_BE_2, -999) 

smry_df$cov_HH_3_l<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$covSat_l_HH,smry_df$cov_HH_l)
smry_df$cov_HH_3_u<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$covSat_u_HH,smry_df$cov_HH_u)
smry_df$cov_HH_3<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$covSat_HH,smry_df$cov_HH)

smry_df$cov_BE_3_l<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$covSat_l_BE,smry_df$cov_BE_l)
smry_df$cov_BE_3_u<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$covSat_u_BE,smry_df$cov_BE_u)
smry_df$cov_BE_3<-ifelse(smry_df$S<10 | smry_df$K<10,smry_df$covSat_BE,smry_df$cov_BE)

#ICS correctly specified
smry_df$cov_HH_3_l_ICS<-ifelse(smry_df$CAC != 1, smry_df$cov_HH_3_l,  -999) 
smry_df$cov_HH_3_u_ICS<-ifelse(smry_df$CAC != 1, smry_df$cov_HH_3_u,  -999) 
smry_df$cov_HH_3_ICS <- ifelse(smry_df$CAC != 1, smry_df$cov_HH_3,  -999) 

smry_df$cov_BE_3_ICS_l <- ifelse(smry_df$CAC == 1,  smry_df$cov_BE_3_l, -999) 
smry_df$cov_BE_3_ICS_u <- ifelse(smry_df$CAC == 1,  smry_df$cov_BE_3_u, -999) 
smry_df$cov_BE_3_ICS <- ifelse(smry_df$CAC == 1,  smry_df$cov_BE_3, -999) 

# Reorder the CAC factor levels
smry_df$CAC <- factor(smry_df$CAC, levels = c(1, 0.95, 0.8, 0.5))

# Plotting function
plot_nested_loop <- function(output_filename, data_subset, theta_val, y_variable_name,legend_labels) {
  png(file = output_filename, width = 1200, height = 900) 
  
  nested_loop_plot(resdf = data_subset, 
                   x = "m",
                   grid_rows = "CAC", steps = c("ICC","K","S"),
                   steps_y_base = fval$syb, 
                   steps_y_height =fval$syh,  
                   steps_y_shift = fval$sysh, 
                   x_name = "m", y_name = y_variable_name,
                   spu_x_shift = 60,line_size = 1, line_alpha = 1,
                   #point_shapes =c(15, 19, 17),point_size =3,
                   #colors = c("#000080","#006400","#DC143C"),#scales::brewer_pal(palette = "Set1"), for #1
                   #colors = c("#FF6B6B","dimgrey"), #for #2
                   #colors = c("#4ECDC4","dimgrey"), #for #2
                   # line_linetypes = c(1,1,3),
                   point_shapes =c(20, 19, 20,20,19,20),point_size =2,
                   #colors = c("#F08080","#FF6B6B","#F08080","#A7E9DF","#4ECDC4","#A7E9DF","#5D3FD3"), #CS
                   colors = c("#A7E9DF","#4ECDC4","#A7E9DF","#F08080","#FF6B6B","#F08080","#5D3FD3"), #ICS
                   line_linetypes = c(1,1,1,1,1,1),
                   steps_color = "dimgrey",
                   steps_values_annotate = TRUE, steps_annotation_size = 3,
                   steps_annotation_nudge = 0.4,
                   steps_annotation_color = "dimgrey", ylim = c(ylow,yup),
                   hline_intercept = 95,
                   y_expand_add = c(fval$yexp1,fval$yexp2),
                   legend_name = ifelse(theta_val == "zero","Method","Method"),
                   legend_labels = legend_labels,base_size = 16,
                   steps_names = c("ICC","K", "S"),
                   post_processing = list(
                     add_custom_theme = list(
                       axis.text.x = element_text(angle = -90, vjust = 0.5, size = 10,face = "bold"), 
                       axis.text.y = element_text(face = "bold", size = 14), 
                       strip.text = element_text(face = "bold", size = 14),
                       axis.title.y = element_text(face = "bold", size = 16)
                     ))
  )
  
}

# Plotting for Coverage (%) using KR correction 
# Plotting for Coverage (%) using Sat correction  
smry_df_sub <- smry_df[(smry_df$theta > 0 & 
                          smry_df$`Time effect`=='categorical'), c("CAC","ICC","S", "K", "m", "cov_BE_3_ICS_l","cov_BE_3_ICS","cov_BE_3_ICS_u",
                                                                   "cov_HH_3_l_ICS","cov_HH_3_ICS","cov_HH_3_u_ICS")]
plot_nested_loop(output_filename = "Fig5_estcovSat_theta15_CatbyCACs.png", 
                 data_subset = smry_df_sub,theta_val = "non-zero",y_variable_name="Coverage (%)",
                 legend_labels = c("+ 2 * MCSE", "BE model W/ Sat", "- 2 * MCSE",
                                   "+ 2 * MCSE","Exch model W/ Sat","- 2 * MCSE"))
dev.off()
write.csv(smry_df_sub, "smry_df_sub_covICS.csv", row.names = FALSE)

#1
#categorical
smry_df_sub <- smry_df[(smry_df$theta > 0 & 
                          smry_df$`Time effect`=='categorical'), c("CAC","ICC","S", "K", "m","cov_BE_2_ICS_l","cov_BE_2_ICS","cov_BE_2_ICS_u",
                                                                   "cov_HH_2_l_ICS","cov_HH_2_ICS","cov_HH_2_u_ICS")]
plot_nested_loop(output_filename = "FigS7_estcovKR_theta15_CatbyCACs.png", 
                 data_subset = smry_df_sub,theta_val = "non-zero",y_variable_name="Coverage (%)",
                 legend_labels = c("+ 2 * MCSE", "BE model W/ KR", "- 2 * MCSE",
                                   "+ 2 * MCSE","Exch model W/ KR","- 2 * MCSE"))
dev.off()
#linear
smry_df_sub <- smry_df[(smry_df$theta > 0 & 
                          smry_df$`Time effect`=='linear'), c("CAC","ICC","S", "K", "m", "cov_BE_3_ICS_l","cov_BE_3_ICS","cov_BE_3_ICS_u",
                                                              "cov_HH_3_l_ICS","cov_HH_3_ICS","cov_HH_3_u_ICS")]
plot_nested_loop(output_filename = "FigS15_estcovSat_theta15_LinbyCACs.png", 
                 data_subset = smry_df_sub,theta_val = "non-zero",y_variable_name="Coverage (%)",
                 legend_labels = c("+ 2 * MCSE", "BE model W/ Sat", "- 2 * MCSE",
                                   "+ 2 * MCSE","Exch model W/ Sat","- 2 * MCSE"))
dev.off()

smry_df_sub <- smry_df[(smry_df$theta > 0 & 
                          smry_df$`Time effect`=='linear'), c("CAC","ICC","S", "K", "m", "cov_BE_2_ICS_l","cov_BE_2_ICS","cov_BE_2_ICS_u",
                                                              "cov_HH_2_l_ICS","cov_HH_2_ICS","cov_HH_2_u_ICS")]
plot_nested_loop(output_filename = "FigS19_estcovKR_theta15_LinbyCACs.png", 
                 data_subset = smry_df_sub,theta_val = "non-zero",y_variable_name="Coverage (%)",
                 legend_labels = c("+ 2 * MCSE", "BE model W/ KR", "- 2 * MCSE",
                                   "+ 2 * MCSE","Exch model W/ KR","- 2 * MCSE"))
dev.off()
# here we need to change the colors
#2 Plotting for HH using Sat vs uncorrected method for coverage (%)
smry_df_sub <- smry_df[(smry_df$theta > 0 &
                          smry_df$`Time effect`=='categorical'), c("CAC","ICC","S", "K", "m", "covSat_HH","cov_HH")]
plot_nested_loop(output_filename = "FigS8_estcovSat_vs_cov_theta15_HH_CatbyCACs.png",
                 data_subset = smry_df_sub,theta_val = "non-zero", y_variable_name= "Coverage (%)",
                 legend_labels = c("Exch model W/ Sat",  "Exch model W/O Sat"))
dev.off()
# Plotting for BE using Sat vs uncorrected method for coverage (%)
#3 here we need to change the colors
smry_df_sub <- smry_df[(smry_df$theta > 0 &
                          smry_df$`Time effect`=='categorical'), c("CAC","ICC","S", "K", "m", "covSat_BE","cov_BE")]
plot_nested_loop(output_filename = "FigS9_estcovSat_vs_cov_theta15_BE_CatbyCACs.png",
                 data_subset = smry_df_sub,theta_val = "non-zero", y_variable_name= "Coverage (%)",
                 legend_labels = c("BE model W/ Sat","BE model W/O Sat"))
dev.off()
#2 Plotting for HH using KR  vs uncorrected method for coverage (%)
smry_df_sub <- smry_df[(smry_df$theta > 0 &
                          smry_df$`Time effect`=='categorical'), c("CAC","ICC","S", "K", "m", "covKR_HH","cov_HH")]
plot_nested_loop(output_filename = "FigS10_estcovKR_vs_cov_theta15_HH_CatbyCACs.png",
                 data_subset = smry_df_sub,theta_val = "non-zero", y_variable_name= "Coverage (%)",
                 legend_labels = c("Exch model W/ KR",  "Exch model W/O KR"))
dev.off()
# Plotting for BE using KR vs uncorrected method for coverage (%)
#3 here we need to change the colors
smry_df_sub <- smry_df[(smry_df$theta > 0 &
                          smry_df$`Time effect`=='categorical'), c("CAC","ICC","S", "K", "m", "covKR_BE","cov_BE")]
plot_nested_loop(output_filename = "FigS11_estcovKR_vs_cov_theta15_BE_CatbyCACs.png",
                 data_subset = smry_df_sub,theta_val = "non-zero", y_variable_name= "Coverage (%)",
                 legend_labels = c("BE model W/ KR","BE model W/O KR"))
dev.off()




#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#Singular
setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Staircase_simstudy\\")
source('5. dat_plots.R')
setwd("G:\\Shared drives\\Ehsan PhD work\\paper_3\\Figures_v2_MCSE\\")

ylow <- 0
yup <- 1000
fval <- cal_val(ylow, yup)

# Plotting function
plot_nested_loop <- function(output_filename, data_subset, theta_val, y_variable_name,legend_labels) {
  png(file = output_filename, width = 1200, height = 900) 
  
  nested_loop_plot(resdf = data_subset, 
                   x = "m",
                   grid_rows = "CAC", grid_cols = "Time effect", steps = c("ICC","K","S"),
                   steps_y_base = fval$syb, 
                   steps_y_height =fval$syh,  
                   steps_y_shift = fval$sysh, 
                   x_name = "m", y_name = y_variable_name,
                   spu_x_shift = 60,,line_size = 1, line_alpha = 1,
                   point_shapes =c(19,19),point_size =2,
                   colors = c("#FF6B6B","#4ECDC4"),
                   line_linetypes = c(1, 1),
                   steps_color = "dimgrey",
                   steps_values_annotate = TRUE, steps_annotation_size = 3,
                   steps_annotation_nudge = 0.4,
                   steps_annotation_color = "dimgrey", ylim = c(ylow,yup),
                   y_expand_add = c(fval$yexp1,fval$yexp2),
                   legend_name = ifelse(theta_val == "zero","Method","Method"),
                   legend_labels = legend_labels,base_size = 16,
                   steps_names = c("ICC","K", "S"),
                   post_processing = list(
                     add_custom_theme = list(
                       axis.text.x = element_text(angle = -90, vjust = 0.5, size = 7,face = "bold"), 
                       axis.text.y = element_text(face = "bold", size = 14), 
                       strip.text = element_text(face = "bold", size = 14),
                       axis.title.y = element_text(face = "bold", size = 16) 
                     ))
  )
  
}

# Reorder the CAC factor levels
smry_df$CAC <- factor(smry_df$CAC, levels = c(1, 0.95, 0.8, 0.5))
# Plotting for number of non-Singular fit using KR correction
smry_df_sub <- smry_df[(smry_df$theta > 0 ), c("Time effect","CAC","ICC","S", "K", "m", "sIsSing_HH", "sIsSing_BE")]
plot_nested_loop(output_filename = "FigS01_sSingKR_theta15.png",
                 data_subset = smry_df_sub,theta_val="non-zero", y_variable_name= "Number of Singular fit",
                 legend_labels = c("Exch model", "BE model"))
dev.off()

