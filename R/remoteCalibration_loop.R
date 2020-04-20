

# CN run calibration in pearcey 
# thsi file will be uploaded in the server (you need to keep a copy as the file will be deleted in the server)
# calibration will be run on the server
# outputs will be store on the server and saved on mac
# all directory need to be saved accordingly 

##### load inputs as in starting community, effort etc. 
rm(list=ls())
#setwd("/home/nov017/flush1/R_model")
load("dataForRemoteCalibration.RData")
.libPaths("~/local/R_libs/")

# only parse the command line arguments following the script name
argv <- commandArgs(TRUE)
# print(argv)

# CN not working
# argv <- c("SSBcali_value", "rankCali_valu", "Q_value", "K_value", "R_value", "changeSSBCali_value", "output_file_name") 
# str(argv)

# check for required num.rows argument
if (length(argv) < 6)
    stop("Missing argument: SSBcali, rankCali, Q, K, R, output_file_name")
  
# Arguments have type character, so coerce to logical - should cast a 1 to TRUE, 1 to FALSE
SSBcali_value <- as.logical(as.numeric(argv[1]))
rankCali_value <- as.logical(as.numeric(argv[2]))
Q_value <- as.logical(as.numeric(argv[3]))
K_value <- as.logical(as.numeric(argv[4]))
R_value <- as.logical(as.numeric(argv[5]))
# changeSSBCali_value <- as.logical(as.numeric(argv[6]))

# # try only the one you are interested in 
# SSBcali_value <- TRUE
# rankCali_value <- FALSE
# Q_value <- FALSE
# K_value <- TRUE
# R_value <- TRUE

output_file_name <- argv[6] # will be a string by default.

print(K_value)
##### load all files and libraries to run model
#library(tidyverse)
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
meantsteps = 30 
extinction_test = TRUE 
extinct_threshold = 0.9
diet_steps = 0


##### calibration options:

### option A - calibrate the kappas (kappa & kappa_ben), r_max and Q

# kappa, r_max and Q
logR <- log10(Initialcomm@params@species_params$r_max)
logQ <- log10(df_param$catchability)
logK <- log10(kappa)
logK_ben <- log10(kappa_ben)
logParams <- c(logR, logQ, logK, logK_ben)

# # lower and upper are needed for calibration with Q as this value should not be above 1
lower = c(log10(Initialcomm@params@species_params$r_max*0.01),log10(rep(0.001, nrow(df_param))), log10(kappa*0.01), log10(kappa_ben*0.01))
upper = c(log10(Initialcomm@params@species_params$r_max/0.01),log10(rep(1, nrow(df_param))), log10(kappa/0.01), log10(kappa_ben/0.01))

# trial here initial community is sim2, but kappa, rmax nd q values are those per calibration with extinction... see if all tests work
# use calibrated values from a colibration where spp got extinct 
# load("optim_SEA_kappa_SSB_K.RData")# ... EXTINCTION!
# kappa<-10^optim_SEA$par[length(optim_SEA$par)]
# r_max<-10^optim_SEA$par[1:nrow(df_param)]
# Q<-10^optim_SEA$par[(nrow(df_param)+1):(length(optim_SEA$par)-1)]
# df_param<-df_param
# kappa_ben<-kappa
# df_param$r_max<-r_max
# df_param$catchability<-Q
# 
# SSBcali_value <- TRUE
# rankCali_value <- TRUE
# Q_value <- TRUE
# K_value <- TRUE
# R_value <- TRUE

# run calibration
optimizer_count=0 # Initialize count of function evaluations
optim_SEA   <- optim(par = logParams,
                     lower=lower,
                     upper=upper,
                     method ="L-BFGS-B",
                     control=list(maxit=10000), # REPORT=1, trace=6),
                     fn = calibrate,
                     SSBcali = SSBcali_value, 
                     rankCali = rankCali_value, 
                     Q = Q_value, 
                     K = K_value,
                     R = R_value,
                     # changeSSBCali = changeSSBCali_value,
                     extinction_test.f = extinction_test) # make sure extinction test is running?!


# save calibrated param
save(optim_SEA, file = output_file_name)


