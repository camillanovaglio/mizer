


#----------------------------------------------
# calibration function
#----------------------------------------------

calibrate <- function(logParams,
                      Initialcomm.f = Initialcomm,
                      meantsteps.f = meantsteps, 
                      extinction_test.f = extinction_test, 
                      df_param.f = df_param,
                      extinct_threshold.f = extinct_threshold, 
                      
                      # according to param to calibrate
                      K = K, # not sure how to specify param when kappa dn r_max are not present... 
                      Q = Q,
                      R = R,
                      
                      
                      # according to data use: 
                      SSBcali = SSBcali, 
                      rankCali = rankCali){
  
  # # trial
  # logParams = logParams
  # Initialcomm.f = Initialcomm
  # meantsteps.f = meantsteps
  # extinction_test.f = extinction_test
  # extinct_threshold.f = extinct_threshold
  # df_param.f = df_param
  # Q=TRUE
  # K=TRUE
  # R=TRUE
  # SSBcali = TRUE
  # rankCali = TRUE
  
  optimizer_count <- optimizer_count + 1                                       
  assign("optimizer_count", optimizer_count, pos = .GlobalEnv)
  
  # kappa and Q
  if(K == TRUE){
    kappa <- 10^(logParams[length(logParams)]) 
  }
  if(Q == TRUE){
    # catchability - it should be here as it's transformed into mizer object in MizerParams and used as object in project
    df_param.f$catchability<-10^(logParams[(nrow(df_param)+1):(length(logParams)-1)])
  }
  
  # mizer params 
  params <- MizerParams(df_param.f, interaction = theta, kappa = kappa, kappa_ben = kappa_ben, kappa_alg = kappa_alg, w_pp_cutoff = w_pp_cutoff, min_w_bb = min_w_bb, w_bb_cutoff = w_bb_cutoff, fleetDynamics = fleetDynamics, selectivity_params = df_selParam, catchability = df_Q, target = df_target)
  
  # update initial abundance to those at equlibrium - abundance of initial community
  # if(!is.na(Initialcomm.f)){
    initial_n = Initialcomm.f@n[dim(Initialcomm.f@n)[1],,]
    initial_n_pp = Initialcomm.f@n_pp[dim(Initialcomm.f@n_pp)[1],]
    initial_n_bb = Initialcomm.f@n_bb[dim(Initialcomm.f@n_bb)[1],]
  # } else{
  #   initial_n = initial_n
  #   initial_n_pp = initial_n_pp
  #   initial_n_bb = initial_n_bb
  # }
  
  # r_max (also when you are also calibrating Q) 
    if(R == TRUE){
      params@species_params$r_max <- 10^(logParams[1:nrow(df_param.f)]) 
    }
  
  # projection 
  # if you specified an initial community, then project starting from those abundances
  # if(!is.na(Initialcomm.f)){
    sim <- project(params, t_max = t_max, effort = effort, dt = dt, fleetDynamics = fleetDynamics, management = management, price =price, cost = cost, diet_steps = diet_steps, initial_n = initial_n, initial_n_pp = initial_n_pp)
  # }else{
    # sim <- project(params, t_max = t_max, effort = effort, dt = dt, fleetDynamics = fleetDynamics, management = management, price =price, cost = cost, diet_steps = diet_steps)
  # }
  
  # plot(sim)
  
  # error basaed on yields (always done) 
  error <- Yielderror(sim, meantsteps = meantsteps.f) 
  
  # error basaed on changes in SSB before/after fishing ?
  # cSSBe <- changeSSBerror(sim, meantsteps = meantsteps.f)
  
  # error basaed on SSB
  if(SSBcali==TRUE){
    se<-SSBerror(sim, meantsteps = meantsteps.f)
    error <- error + se
  }
  
  # error basaed on rank abundance 
  if(rankCali==TRUE){
    re<-rankerror(sim, meantsteps = meantsteps.f)
    error <- error + re
  }
  
  # # error basaed on rank abundance 
  # re<-rankerror(sim, meantsteps = meantsteps.f)
  # 
  # # total error 
  # error <- ye + se + re # + cSSBe
  # 
  # # error based on ssb if you specify the argument
  # 
  # if(SSBcali==TRUE){
  #   se<-SSBerror(sim, meantsteps = meantsteps.f)
  #   error <- ye + se
  # }
  
  # extinction
  if(extinction_test.f==TRUE){
    extinct <- Extinct_test(sim, extinct_threshold = extinct_threshold.f)
    error <- error + extinct
  }
  
  print(paste(optimizer_count, "Error:", error))
  print(logParams)
  return(error)
  
}

#----------------------------------------------
# yelderror and extinction function for calibration
#----------------------------------------------

Yielderror <- function(Newcomm, 
                       meantsteps){ 
  
  # # # trial 
  # Newcomm = sim
  # meantsteps = meantsteps
  
  # do not compare yield from mycto - not there - and perch - there but not to be considered 
  toDelete<-c("myctophids","helicolenus barathri")
  temp<-Newcomm@params@species_params
  temp<-temp[-which(temp$species %in% toDelete),]
  Yieldobs <- temp$catchComm_t

  # modelled yield
  yield<-getYield(Newcomm)
  # do not compare yield from mycto - not there - and perch - there but not to be considered as above
  yield<-yield[, !colnames(yield) %in% toDelete]
    
  if (is.na(meantsteps)){
    Yieldhat <- yield[dim(yield)[1],] #*scaling # to use if kappa =0.001
  }else{
    Yieldhat <- apply(yield[(dim(yield)[1]-meantsteps+1):dim(yield)[1],],2,mean) #*scaling
  }
  
  # yield difference and error (form Reum)
  ydif<- log10(Yieldhat)-log10(Yieldobs)
  # cat("\nthis yield= ",ydif)
  
  # relplace all inf/NaN/NA with NA before sum() - problems otherwise 
  if(length(ydif[!is.finite(ydif)])>0){
    ydif[!is.finite(ydif)]<-NA
  }
  
  Yielderr<-sum(ydif^2,na.rm=TRUE) # if things are not working - check this too sum(is.finite(ydif^2))
  # cat("\nthe new value= ", yiel_try, "\t the old value= ", Yielderr)
  
  # note from Reum
  # Yielderr<-sum(is.finite(ydif^2), na.rm=TRUE)/var(log10(na.omit(Yieldobs[Yieldobs>0]))) # (from Reum) Standardize the RSS by the variance of the yield observations. We're doing this because we are combining RSS from SSB (or biomass) and they are one a slightly larger scale, but we whish to weight the data types the same in the objective function (see Hilborn and Walters 1992).
  
  return(Yielderr)
}

## changes in SSB

changeSSBerror <- function(Newcomm, 
                       meantsteps){ 
  
  # # # trial 
  # Newcomm = sim
  # meantsteps = meantsteps
  
  # do not compare yield from mycto - not there - and perch - there but not to be considered 
  toDelete<-c("myctophids","helicolenus barathri","nototodarus gouldi","trachurus declivis")
  temp<-Newcomm@params@species_params
  temp<-temp[-which(temp$species %in% toDelete),]
  changeObs <- temp$changesSSB
  
  # you need to run it with an effort matrix of 1/2 no fishing and 1/2 fishing 
  changeMod<-round((getSSB(Newcomm)[dim(effort)[1],]/getSSB(Newcomm)[(dim(effort)[1]/2)-1,])*100)
  changeMod<-changeMod[- which(names(changeMod) %in% toDelete)]
  
  # yield difference and error (form Reum)
  changedif<- log10(changeMod)-log10(changeObs)

  # relplace all inf/NaN/NA with NA before sum() - problems otherwise 
  if(length(changedif[!is.finite(changedif)])>0){
    changedif[!is.finite(changedif)]<-NA
  }
  
  changeSSBerr<-sum(changedif^2,na.rm=TRUE) 
  
  return(changeSSBerr)
}

## SSB 

SSBerror <- function(Newcomm, 
                       meantsteps){ 
  
  # trial 
  # Newcomm = sim2
  # meantsteps = meantsteps
  
  # compare only spp for which you have values 
  toKeep<-Newcomm@params@species_params[which(!is.na(Newcomm@params@species_params$ssbObs)),]
  
  # observed SSB 
  SSBobs <- toKeep$ssbObs
  
  # modelled SSB
  ssb<-getSSB(Newcomm)
  
  if (is.na(meantsteps)){
    ssbhat <- ssb[dim(ssb)[1],] 
  }else{
    ssbhat <- apply(ssb[(dim(ssb)[1]-meantsteps+1):dim(ssb)[1],],2,mean) 
  }
  
  # compare only ssb of spp you have info for 
  ssbhat<-ssbhat[names(ssbhat) %in% rownames(toKeep)]
  
  # yield difference and error (form Reum)
  ssbdif<- log10(ssbhat)-log10(SSBobs)
  
  # relplace all inf/NaN/NA with NA before sum() - problems otherwise 
  if(length(ssbdif[!is.finite(ssbdif)])>0){
    ssbdif[!is.finite(ssbdif)]<-NA
  }
  
  ssberr<-sum(ssbdif^2,na.rm=TRUE) 
  
  return(ssberr)
}

## rank abundance  

rankerror <- function(Newcomm, 
                     meantsteps){ 
  
  # # trial 
  # Newcomm = sim
  # meantsteps = meantsteps
  
  # observed rank 
  rankObs <- Newcomm@params@species_params$rankAb
  
  # modelled rank
  rank<-getBiomass(Newcomm)
  
  if (is.na(meantsteps)){
    rankhat <- rank[dim(rank)[1],] 
  }else{
    rankhat <- apply(rank[(dim(rank)[1]-meantsteps+1):dim(rank)[1],],2,mean) 
  }
  
  rank<-rank(-rankhat)
  
  # yield difference and error (form Reum)
  rankdif<- log10(rank)-log10(rankObs)
  
  # relplace all inf/NaN/NA with NA before sum() - problems otherwise 
  if(length(rankdif[!is.finite(rankdif)])>0){
    rankdif[!is.finite(rankdif)]<-NA
  }
  
  rankerr<-sum(rankdif^2,na.rm=TRUE) 
  
  return(rankerr)
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







