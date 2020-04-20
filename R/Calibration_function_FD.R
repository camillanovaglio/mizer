


#----------------------------------------------
# calibration function for FD
#----------------------------------------------

calibrate <- function(logParams,
                      Initialcomm.f = Initialcomm,
                      meantsteps.f = meantsteps, 
                      extinction_test.f = extinction_test, 
                      df_param.f = df_param, 
                      ke_fleet.f = ke_fleet, 
                      scaling_price.f = scaling_price, 
                      scaling_effort.f = scaling_effort,
                      extinct_threshold.f = extinct_threshold, 
                      P = P, 
                      F = F, 
                      yieldCali = yieldCali, 
                      effortCali = effortCali){
  
  # logParams = logParams
  # Initialcomm.f = Initialcomm
  # meantsteps.f = meantsteps
  # extinction_test.f = extinction_test
  # df_param.f = df_param
  # ke_fleet.f = ke_fleet
  # scaling_price.f = scaling_price
  # scaling_effort.f = scaling_effort
  # extinct_threshold.f = extinct_threshold
  # P = TRUE
  # F = TRUE
  # yieldCali = TRUE
  # effortCali = TRUE
  
  optimizer_count <- optimizer_count + 1                                       
  assign("optimizer_count", optimizer_count, pos = .GlobalEnv)
  
  # mizer params 
  params <- MizerParams(df_param.f, interaction = theta, kappa = kappa, kappa_ben = kappa_ben, kappa_alg = kappa_alg, w_pp_cutoff = w_pp_cutoff, min_w_bb = min_w_bb, w_bb_cutoff = w_bb_cutoff, fleetDynamics = fleetDynamics, selectivity_params = df_selParam_new, catchability = df_target, target = df_target)
  
  # chenge param to calibrate
  if(P == TRUE){
    scaling_price.f <- 10^(logParams[length(logParams)]) 
  }
  if(F == TRUE){
    ke_fleet.f$ke <- 10^(logParams[-length(logParams)])  
  }
  
  # update initial abundance to those at equlibrium - abundance of initial community
  initial_n = Initialcomm.f@n[dim(Initialcomm.f@n)[1],,]
  initial_n_pp = Initialcomm.f@n_pp[dim(Initialcomm.f@n_pp)[1],]
  initial_n_bb = Initialcomm.f@n_bb[dim(Initialcomm.f@n_bb)[1],]
  
  # update initial effort 
  # initial_effort = Initialcomm.f@effortOut[dim(Initialcomm.f@effortOut)[1],]
  
  # projection 
  sim <- project(params, effort = effort, dt = dt, fleetDynamics = fleetDynamics, management = management, price =price, cost = cost, diet_steps = diet_steps, ke = ke_fleet.f, initial_n = initial_n, initial_n_pp = initial_n_pp, initial_n_bb = initial_n_bb, initial_effort = initial_effort, scaling_price = scaling_price.f, Blevel_management = Blevel_management)
  
  # plot_CN(sim)
  # plotFleet(sim)
  # 
  # # if Initialcomm.f== sim_FD and initial effort is per sim_FD, the community does not change when I change fishing parameters (scaling_price and ke) because managemment keeps effort from increasing. however, effort can decrease if revenue decreases - i.e. scaling_price.f2 = scaling_price.f/1000)
  # scaling_price.f2 = scaling_price.f*1000
  # ke_fleet.f2<-ke_fleet.f
  # ke_fleet.f2$ke <- ke_fleet.f$ke*100
  # 
  # # if I don't statr from sim_FD, but I use sim_calibrated as starting point instead and I keep effort as the initial starting value? this is the current trial  
  # 
  # initial_effort = rep(0.1,5)
  # initial_n_new = Initialcomm.f@n[1,,]
  # initial_n_pp_new = Initialcomm.f@n_pp[1,]
  # initial_n_bb_new = Initialcomm.f@n_bb[1,]
  # 
  # # projection 
  # sim2 <- project(params, effort = effort, dt = dt, fleetDynamics = fleetDynamics, management = management, price =price, cost = cost, diet_steps = diet_steps, ke = ke_fleet.f2, initial_n = initial_n_new, initial_n_pp = initial_n_pp_new, initial_n_bb = initial_n_bb_new, initial_effort = initial_effort, scaling_price = scaling_price.f2, Blevel_management = Blevel_management)
  # 
  # plot_CN(sim2)
  # plotFleet(sim2)
  
  # error basaed on yields 
  # NOTE: if yield is FALSE, then the below cannot be done as error is not specified as argument
  
  if(yieldCali == TRUE){
    ye <- Yielderror_FD(sim, meantsteps = meantsteps.f) 
    error <- ye
  }

  # add effort error
  if(effortCali == TRUE){
    ee <- Efforterror_FD(sim, meantsteps = meantsteps.f) 
    error <- ye + ee
  }
  
  # extinction
  if(extinction_test.f==TRUE){
    extinct <- Extinct_test(sim, extinct_threshold = extinct_threshold.f, meantsteps = meantsteps.f)
    error <- error + extinct
  }
  
  # fleet inactivity - I need to penalise combinations that make fleets inactive - though it should be included in yield and effort error... increase error if any yield or any effort == 0? 
  
  if(0 %in% sim@effortOut[dim(sim@effortOut)[1],] == TRUE){
    error <- error + 500
  }

  print(paste(optimizer_count, "Error:", error))
  print(logParams)
  return(error)
  
}

#----------------------------------------------
# yelderror and extinction function for calibration
#----------------------------------------------

Yielderror_FD <- function(Newcomm, 
                       meantsteps){ 
  
  # Newcomm = sim
  # meantsteps = meantsteps.f
  
  # calibrated to yield fleet
  # Yieldobs <- compare$yobs 
  Yieldobs<- compare[which(compare$theme=="yield"), "obs"]
  
  yield<-Newcomm@yield
  yield<-rowSums(aperm(yield,c(4,1,2,3)),dims=3)
  yield<-rowSums(yield,dims=2)

  if (is.na(meantsteps)){
    Yieldhat <- yield[,dim(yield)[2]] 
  }else{
    Yieldhat <- apply(yield[,(dim(yield)[2]-meantsteps+1):dim(yield)[2]],1,mean)
  }
  
  # yield difference and error (form Reum)
  ydif<- log10(Yieldhat)-log10(Yieldobs)
  # cat("\nthis yield= ",ydif)
  
  # relplace all inf/NaN/NA with NA before sum() - problems otherwise 
  # BUT inf/NaN/NA are given if Yieldhat = 0 (fleet not active) and replacing this with NA means that the sum() below either does not consider 0 yields and calcualtes errors based on active fleets alone or, if returns 0 (perfect fit) if all fleets report yield = 0 - i.e. not active. Instead set the difference to a number so that the error increase for any non-active fleet (e.g. 2 - so taht if each fleet is inactive you get an error of 20). Same is done in Efforterror_FD and extinction_test
  
  if(length(ydif[!is.finite(ydif)])>0){
    ydif[!is.finite(ydif)]<-2
  }
  
  Yielderr<-sum(ydif^2,na.rm=TRUE) 
  
  return(Yielderr)
}

Efforterror_FD <- function(Newcomm, 
                       meantsteps){ 
  
  # calibrated to effort fleet
  # EffortObs <- compare$eobs
  EffortObs <- compare[which(compare$theme=="effort"), "obs"]
  effort<-Newcomm@effortOut
  
  if (is.na(meantsteps)){
    efforthat <- effort[dim(effort)[1],] 
  }else{
    efforthat <- apply(effort[(dim(effort)[1]-meantsteps+1):dim(effort)[1],],2,mean) #
  }
  
  # scaling for effort 
  efforthat<-efforthat/scaling_effort
  
  # effort difference and error (form Reum)
  edif<- log10(efforthat)-log10(EffortObs)
  
  # relplace all inf/NaN/NA with NA before sum() - problems otherwise 
  if(length(edif[!is.finite(edif)])>0){
    edif[!is.finite(edif)]<-2
  }
  
  Eerr<-sum(edif^2,na.rm=TRUE) 
  
  return(Eerr)
}

# extinction - see Calibration_function.R for explanation

Extinct_test <- function(sim, 
                         extinct_threshold, 
                         meantsteps){
  
  Biomass <- getBiomass(sim)
  BiomassInit<- Biomass[1,] # option 1 - if you consider extinction 
  BiomassInit<- Biomass[dim(Biomass)[1]-meantsteps,] 
  BiomassEnd<-  Biomass[dim(Biomass)[1],]
  
  relB<-(BiomassEnd/BiomassInit) 
  relB<-relB[relB<extinct_threshold] 
  extinct<- sum(1/relB)
  
  # if some spp are extinct or declining, this error needs to be big
  if(is.finite(extinct) & extinct>0){
    extinct = 500
  }
  
  # if this returs: 
  if(is.infinite(extinct)){
    extinct = 500  
  }
  
  return(extinct)
}






















