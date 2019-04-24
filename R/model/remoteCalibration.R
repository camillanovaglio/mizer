

# CN run calibration in pearcey 
# thsi file will be uploaded in the server (you need to keep a copy as the file will be deleted in the server)
# calibration will be run on the server
# outputs will be store on the server and saved on mac
# all directory need to be saved accordingly 

##### load inputs as in starting community, effort etc. 
rm(list=ls())
setwd("/home/nov017/flush1/R_model")
load("dataForRemoteCalibration.RData")

##### load all files and libraries to run model
library(tidyverse)
library(devtools)
library(plyr) 
library(Rcpp) # this will allow you to run the inner_project_loop 
library(reshape2)
library(inline)
library(ggplot2)
library(vegan)
library(mizer) # if inner_loop error, you need to upload mizer
# install.packages("mizer")

source("MizerParams-class.R") 
source("MizerSim-class.R")
source("plots.R")
source("project_methods.R")
source("project.R")
source("RcppExports.R") 
source("selectivity_funcs.R")
source("summary_methods.R")
source("wrapper_functions.R")
source("Calibration_function.R")

##### set optimisation criteria
Initialcomm = sim2
kappa = kappa # 0.001
df_param = df_param

# for mizerParams()
theta = theta2
fleetDynamics = FALSE
selectivity_params = NA
catchability = NA
target = NA
dt =  0.5

# for project() 
kappa_ben = kappa  
kappa_alg = 0 
w_pp_cutoff = w_pp_cutoff 
min_w_bb = min_w_bb 
w_bb_cutoff = w_bb_cutoff
management = FALSE
effort = matrix_effort[301:400,] # constant effort (for 100 times) as per effort used in sim2
rownames(effort)<-seq(1:100)
price = NA
cost = NA

# for error functions
meantsteps = 10 
extinction_test = TRUE 
extinct_threshold = 0.01
diet_steps = 0


##### calibration options:

### option A - calibrate kappa, r_max and Q
logR <- log10(Initialcomm@params@species_params$r_max)
logQ <- log10(df_param$catchability)
logK <- log10(kappa)
logParams <- c(logR, logQ, logK)
# lower and upper are needed for calibration with Q as this value should not be above 1 
lower = c(log10(Initialcomm@params@species_params$r_max*0.01),log10(rep(0.001, nrow(df_param))), log10(kappa*0.01))
upper = c(log10(Initialcomm@params@species_params$r_max/0.01),log10(rep(1, nrow(df_param))), log10(kappa/0.01))

# run calibration
optimizer_count=0 # Initialize count of function evaluations
optim_SEA   <- optim(par = logParams,
                     lower=lower,
                     upper=upper,
                     method ="L-BFGS-B",
                     fn = calibrate,
                     SSBcali = TRUE, 
                     rankCali = TRUE, 
                     Q = TRUE, 
                     K = TRUE,
                     R = TRUE)


# save calibrated param
save(optim_SEA, file = "optim_SEA_kappaQrmax_yieldSsbRank.Rdata")

optimizer_count=0 
optim_SEA   <- optim(par = logParams,
                     lower=lower,
                     upper=upper,
                     method ="L-BFGS-B",
                     fn = calibrate,
                     SSBcali = FALSE, 
                     rankCali = TRUE,
                     Q = TRUE, 
                     K = TRUE,
                     R = TRUE)

save(optim_SEA, file = "optim_SEA_kappaQrmax_yieldRank.Rdata")

optimizer_count=0 
optim_SEA   <- optim(par = logParams,
                     lower=lower,
                     upper=upper,
                     method ="L-BFGS-B",
                     fn = calibrate,
                     SSBcali = FALSE, 
                     rankCali = FALSE,
                     Q = TRUE, 
                     K = TRUE,
                     R = TRUE)

save(optim_SEA, file = "optim_SEA_kappaQrmax_yield.Rdata")


### option B- calibrate kappa and r_max
logParams <- c(log10(Initialcomm@params@species_params$r_max),log10(kappa))

optimizer_count=0 
optim_SEA   <- optim(par = logParams,
                     method ="L-BFGS-B",
                     fn = calibrate,
                     SSBcali = TRUE, 
                     rankCali = TRUE,
                     Q = FALSE, 
                     K = TRUE,
                     R = TRUE)
save(optim_SEA, file = "optim_SEA_kappaRmax_yieldSsbRank.Rdata")

optim_SEA   <- optim(par = logParams,
                     method ="L-BFGS-B",
                     fn = calibrate,
                     SSBcali = FALSE, 
                     rankCali = TRUE,
                     Q = FALSE, 
                     K = TRUE,
                     R = TRUE)
save(optim_SEA, file = "optim_SEA_kappaRmax_yieldRank.Rdata")

optim_SEA   <- optim(par = logParams,
                     method ="L-BFGS-B",
                     fn = calibrate,
                     SSBcali = TRUE, 
                     rankCali = FALSE,
                     Q = FALSE, 
                     K = TRUE,
                     R = TRUE)
save(optim_SEA, file = "optim_SEA_kappaRmax_yieldSsb.Rdata")

optim_SEA   <- optim(par = logParams,
                     method ="L-BFGS-B",
                     fn = calibrate,
                     SSBcali = FALSE, 
                     rankCali = FALSE,
                     Q = FALSE, 
                     K = TRUE,
                     R = TRUE)
save(optim_SEA, file = "optim_SEA_kappaRmax_yield.Rdata")
