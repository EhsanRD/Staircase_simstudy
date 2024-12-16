rm(list=ls())
# make a structure for the loop plot values
cal_val <- function(ylow, yup) {
  yexp1 <- 0.75 * (yup - ylow)
  yexp2 <- 0.05 * (yup - ylow)
  syb <- (-0.12) * (yup - ylow) + ylow
  syh <- 0.0533 * yexp1
  sysh <- 0.1733 * yexp1
  
  return(list(yexp1 = yexp1, yexp2 = yexp2, syb = syb, syh = syh, sysh = sysh))
}

library(looplot)
library(gridExtra)

setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Staircase_simstudy\\")

# Read data
smry_df <- read.csv(file='summary_v2_MCSE.csv', header= T)

# Adjust type column
smry_df$type <- ifelse(smry_df$type %in% 'c', "categorical", "linear")

# Rename columns
names(smry_df)[names(smry_df) == "type"]  <- "Time effect"
names(smry_df)[names(smry_df) == "power"]  <- "pow"
