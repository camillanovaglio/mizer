


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
                      scaling_effort.f = scaling_effort){
  
  # # # trial
  # logParams <- c(log10(df_param3$r_max),log10(kappa3))
  # logParams <- c(log10(ke_fleet$ke), log10(scaling_price))
  # Initialcomm.f = Initialcomm
  # meantsteps.f = meantsteps
  # extinction_test.f = extinction_test
  # extinct_threshold.f = extinct_threshold
  # df_param.f = df_param
  # ke_fleet.f = ke_fleet
  # scaling_price.f = scaling_price
  # scaling_effort.f = scaling_effort
  
  optimizer_count <- optimizer_count + 1                                       
  assign("optimizer_count", optimizer_count, pos = .GlobalEnv)
  
  # kappa
  # kappa <- 10^(logParams[length(logParams)])
  
  # mizer params 
  params <- MizerParams(df_param.f, interaction = theta, kappa = kappa, kappa_ben = kappa_ben, kappa_alg = kappa_alg, w_pp_cutoff = w_pp_cutoff, min_w_bb = min_w_bb, w_bb_cutoff = w_bb_cutoff, fleetDynamics = fleetDynamics, selectivity_params = df_selParam, catchability = df_target, target = df_target)
  
  # r_max or ke and scaling
  # params@species_params$r_max <- 10^(logParams[-length(logParams)])
  scaling_price.f <- 10^(logParams[length(logParams)])
  ke_fleet.f$ke <- 10^(logParams[-length(logParams)])
  
  
  # # trial errors ------
  # error0<-c(-0.5093663, -1.6893985, -1.1628286, 182582.5647805, -0.8170459, 182601.8888313)
  # erroeinf<-c( -0.8150214, -1.4521220, -1.1633706, 81.6695910, -0.9667953, 89.2898399)
  # erroBig<-c(-0.8147440,-1.4523374,  -1.1633701, 247.3012384,  -0.9666585, 254.9311093)
  # 
  # scaling_price.f <- 10^(erroeinf[length(erroeinf)])
  # ke_fleet.f$ke <- 10^(erroeinf[-length(erroeinf)])
  # # end trials -------
  
  
  # update initial abundance to those at equlibrium - abundance of initial community
  initial_n = Initialcomm.f@n[dim(Initialcomm.f@n)[1],,]
  initial_n_pp = Initialcomm.f@n_pp[dim(Initialcomm.f@n_pp)[1],]
  initial_n_bb = Initialcomm.f@n_bb[dim(Initialcomm.f@n_bb)[1],]
  # update initial effort 
  initial_effort = Initialcomm.f@effortOut[dim(Initialcomm.f@effortOut)[1],]
  
  # projection 
  sim <- project(params, t_max = t_max, effort = effort, dt = dt, fleetDynamics = fleetDynamics, management = management, price =price, cost = cost, diet_steps = diet_steps, ke = ke_fleet.f, initial_n = initial_n, initial_n_pp = initial_n_pp, initial_n_bb = initial_n_bb, initial_effort = initial_effort, scaling_price = scaling_price.f)
  
  # plot_CN(sim)
  # plotFleet(sim)
  # sim@profit
  
  
  # error basaed on yields 
  ye <- Yielderror_FD(sim, meantsteps = meantsteps.f) 
  error <- ye
  
  # add effort error
  ee <- Efforterror_FD(sim, meantsteps = meantsteps.f) 
  error <- ye + ee 
  
  # extinction
  if(extinction_test.f==TRUE){
    extinct <- Extinct_test(sim, extinct_threshold = extinct_threshold.f)
    error <- error + extinct
  }
  
  # fleet inactivity - I need to penalise combinations that make fleets inactive - though it should be included in yield and effort error... maybe increase error if any yield or any effort ==0? 
  
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
  
  # calibrate to spp
  # Yieldobs <- Newcomm@params@species_params$catchComm_t
  # yield<-Newcomm@yield
  # yield<-rowSums(yield,dims=3) # sum over fleet
  # yield<-rowSums(yield,dims=2) # sum over size

  # calibrated to yield fleet
  Yieldobs <- compare$yobs # these are mean yield over 1995-2005
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
    # ydif[!is.finite(ydif)]<-NA # instead of 
    ydif[!is.finite(ydif)]<-2
  }
  
  Yielderr<-sum(ydif^2,na.rm=TRUE) # if things are not working - check this too sum(is.finite(ydif^2))
  # cat("\nthe new value= ", yiel_try, "\t the old value= ", Yielderr)
  
  return(Yielderr)
}

Efforterror_FD <- function(Newcomm, 
                       meantsteps){ 
  
  # calibrated to effort fleet
  EffortObs <- compare$eobs
  effort<-Newcomm@effortOut
  
  if (is.na(meantsteps)){
    efforthat <- effort[dim(effort)[1],] #*scaling # to use if kappa =0.001
  }else{
    efforthat <- apply(effort[(dim(effort)[1]-meantsteps+1):dim(effort)[1],],2,mean) #*scaling
  }
  
  # scaling for effort 
  efforthat<-efforthat/scaling_effort
  
  # effort difference and error (form Reum)
  edif<- log10(efforthat)-log10(EffortObs)
  
  # relplace all inf/NaN/NA with NA before sum() - problems otherwise 
  # see explanation abive for changes 
  if(length(edif[!is.finite(edif)])>0){
    # ydif[!is.finite(ydif)]<-NA # instead of 
    edif[!is.finite(edif)]<-2
  }
  
  Eerr<-sum(edif^2,na.rm=TRUE) 
  
  return(Eerr)
}

# extinction

Extinct_test <- function(sim, 
                         extinct_threshold){
  
  Biomass <- getBiomass(sim)
  BiomassInit<- Biomass[1,]
  BiomassEnd<-  Biomass[dim(Biomass)[1],]
  
  relB<-(BiomassEnd/BiomassInit)
  relB<-relB[relB<extinct_threshold]
  extinct<- sum(1/relB)
  
  # if this returs: 
  # Inf -> the final biomass of all spp below extinct_threshold = 0 - this needs to change. see below solution
  # big number -> the final biomass of all spp below extinct_threshold = very small number compared to initial - this is OK 
  if(is.infinite(extinct)){
    extinct = 10 # this is simply a high number compared to the expected error ~ 3.5 given the conditions used for this model 
  }
  
  return(extinct)
}




