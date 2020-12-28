
# CN run calibration in pearcey 
# thsi file will be uploaded in the server (you need to keep a copy as the file will be deleted in the server)
# calibration will be run on the server
# outputs will be store on the server and saved on mac
# all directory need to be saved accordingly 

##### load inputs as in starting community, effort etc. 
rm(list=ls())
load("dataForRemoteCalibration_FD.RData")
.libPaths("~/local/R_libs/")

# only parse the command line arguments following the script name
argv <- commandArgs(TRUE)

# check for required num.rows argument
if (length(argv) < 5)
    stop("Missing argument: yieldCali, effortCali, P, F, output_file_name")
  
# Arguments have type character, so coerce to logical - should cast a 1 to TRUE, 1 to FALSE
yieldCali_value <- as.logical(as.numeric(argv[1]))
effortCali_value <- as.logical(as.numeric(argv[2]))
P_value <- as.logical(as.numeric(argv[3]))
F_value <- as.logical(as.numeric(argv[4]))

output_file_name <- argv[5] # will be a string by default.

print(yieldCali_value)
print(effortCali_value)
print(P_value)
print(F_value)


##### load all files and libraries to run model
library(devtools)
library(plyr) 
library(dplyr) 
library(Rcpp) 
library(reshape2)
library(inline)
library(ggplot2)
library(vegan)
library(mizer) 

source("MizerParams-class.R") 
source("MizerSim-class.R")
source("plots.R")
source("project_methods.R")
source("project.R")
source("RcppExports.R") 
source("selectivity_funcs.R")
source("summary_methods.R")
source("wrapper_functions.R")
source("Calibration_function_FD.R")

##### set optimisation criteria
Initialcomm = sim_calibrated # sim_FD # see notes inside calibration_function_FD on why 
kappa = kappa3
kappa_ben = kappa_ben
initial_effort = initial_effort

# turn fishing on 
df_param = df_param3
theta = theta
selectivity_params = df_selParam_new
catchability = df_target
target = df_target
fleetDynamics = TRUE
dt =  0.5

# for project function 
kappa_alg = 0 
w_pp_cutoff = w_pp_cutoff 
min_w_bb = min_w_bb 
w_bb_cutoff = w_bb_cutoff
management = TRUE
effort = 0
price = df_price_mean2
cost = df_cost_mean2
Blevel_management = Blevel_management # could be changed to a combination

# for error functions 
# meantsteps = 10 
# extinction_test = TRUE
# extinct_threshold = 0.01
meantsteps = 30 
extinction_test = TRUE 
extinct_threshold = 0.9

diet_steps = 0

scaling_effort = 10000

##### calibration options:

### option A - calibrate the kappas (kappa & kappa_ben), r_max and Q
logParams <- c(log10(ke_fleet$ke), log10(scaling_price)) 

# lower and upper are needed for calibration with Q as this value should not be above 1
lower = c(log10(ke_fleet$ke*0.01),log10(scaling_price*0.001))
upper = c(log10(ke_fleet$ke/0.01),log10(scaling_price/0.001))

# run calibration
optimizer_count=0 
optim_SEA   <- optim(par = logParams,
                     lower=lower,
                     upper=upper,
                     method ="L-BFGS-B",
                     control=list(maxit=10000), # REPORT=1, trace=6),
                     fn = calibrate,
                     P = P_value, 
                     F = F_value,
                     yieldCali = yieldCali_value, 
                     effortCali = effortCali_value, 
                     extinction_test.f = extinction_test)


# save calibrated param
save(optim_SEA, file = output_file_name)


