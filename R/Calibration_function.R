


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
                      # changeSSBCali = changeSSBCali, 
                      # cpueCali = cpueCali){
  
  # trial - testing calibration outputs
  # logParams = optim_SEA$par
  # logParams = logParams
  # Initialcomm.f = sim2
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
    kappa <- 10^(logParams[length(logParams)-1]) 
    kappa_ben <- 10^(logParams[length(logParams)])
  }
  if(Q == TRUE){
    df_param.f$catchability<-10^(logParams[(nrow(df_param)+1):(length(logParams)-2)])
  }
  
  # mizer params 
  params <- MizerParams(df_param.f, interaction = theta, kappa = kappa, kappa_ben = kappa_ben, kappa_alg = kappa_alg, w_pp_cutoff = w_pp_cutoff, min_w_bb = min_w_bb, w_bb_cutoff = w_bb_cutoff, fleetDynamics = fleetDynamics, selectivity_params = df_selParam, catchability = df_Q, target = df_target)
  
  initial_n = Initialcomm.f@n[dim(Initialcomm.f@n)[1],,]
  initial_n_pp = Initialcomm.f@n_pp[dim(Initialcomm.f@n_pp)[1],]
  initial_n_bb = Initialcomm.f@n_bb[dim(Initialcomm.f@n_bb)[1],]
 
  # r_max 
    if(R == TRUE){
      params@species_params$r_max <- 10^(logParams[1:nrow(df_param.f)]) 
    }
  
  sim <- project(params, t_max = t_max, effort = effort, dt = dt, fleetDynamics = fleetDynamics, management = management, price =price, cost = cost, diet_steps = diet_steps, initial_n = initial_n, initial_n_pp = initial_n_pp)
  # plot(sim)
  
  # unfished community
  effort_unfished<-effort
  effort_unfished[]<-0
  sim_unfished<-project(params, t_max = t_max, effort = effort_unfished, dt = dt, fleetDynamics = fleetDynamics, management = management, price =price, cost = cost, diet_steps = diet_steps, initial_n = initial_n, initial_n_pp = initial_n_pp)
  # plot(sim_unfished)
  
  # error basaed on yields (always done) 
  ye <- Yielderror(sim, meantsteps = meantsteps.f) 
  error <- ye
  
  # error basaed on SSB
  if(SSBcali==TRUE){
    se <- SSBerror(sim, sim_unfished, effort, effort_unfished, meantsteps = meantsteps.f)
    error <- error + se
  }
  
  # error basaed on rank abundance 
  if(rankCali==TRUE){
    re <-rankerror(sim, meantsteps = meantsteps.f)
    error <- error + re
  }
  
  # # error basaed on changes in SSB  
  # if(changeSSBCali==TRUE){
  #   ch_re <-changeSSBerror(sim, sim_unfished, effort, effort_unfished, meantsteps = meantsteps.f)
  #   error <- error + ch_re
  # }
  
  # extinction
  if(extinction_test.f==TRUE){
    extinct <- Extinct_test(sim, extinct_threshold = extinct_threshold.f,meantsteps = meantsteps.f)
    error <- error + extinct
  }
  
  # extinction for unfished community
  if(extinction_test.f==TRUE){
    extinct_unfished <- Extinct_test(sim_unfished, extinct_threshold = extinct_threshold.f,meantsteps = meantsteps.f)
    error <- error + extinct_unfished
  }
  
  
  # total error and weights 
  # you can't udd them up here because you have lots of if() statments 
  # error <- ye + se + re + extinct + extinct_unfished
  
  # how to weight the error !?? 
  # extinction error is big compared to all the others, and that's OK as we don't want this to happen
  
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
  # Yielderr<-sum(is.finite(ydif^2), na.rm=TRUE)/var(log10(na.omit(Yieldobs[Yieldobs>0]))) # Standardize the RSS by the variance of the yield observations. We're doing this because we are combining RSS from SSB (or biomass) and they are one a slightly larger scale, but we whish to weight the data types the same in the objective function (see Hilborn and Walters 1992).
  
  return(Yielderr)
}

## changes in SSB

# changeSSBerror <- function(Newcomm,
#                            Newcomm2,
#                            effort, 
#                            effort_unfished,
#                            meantsteps){ 
#   
#   # # trial
#   # Newcomm = sim
#   # Newcomm2 = sim_unfished
#   # meantsteps = meantsteps
#   
#   # do not compare yield from mycto - not there - and perch - there but not to be considered 
#   toDelete<-c("myctophids","helicolenus barathri","nototodarus gouldi","trachurus declivis")
#   temp<-Newcomm@params@species_params
#   temp<-temp[-which(temp$species %in% toDelete),]
#   changeObs <- temp$changesSSB
#   
#   # you need to run it with an effort matrix of 1/2 no fishing and 1/2 fishing 
#   # changeMod<-round((getSSB(Newcomm)[dim(effort)[1],]/getSSB(Newcomm)[(dim(effort)[1]/2)-1,])*100)
#   # changeMod<-changeMod[- which(names(changeMod) %in% toDelete)]
#   
#   # OR 
#   changeMod<-round((getBiomass(Newcomm)[dim(effort)[1],]/getBiomass(Newcomm2)[(dim(effort_unfished)[1]),])*100)
#   changeMod<-changeMod[- which(names(changeMod) %in% toDelete)]
#   
#   # yield difference and error (form Reum)
#   changedif<- log10(changeMod)-log10(changeObs)
# 
#   # relplace all inf/NaN/NA with NA before sum() - problems otherwise 
#   if(length(changedif[!is.finite(changedif)])>0){
#     changedif[!is.finite(changedif)]<-NA
#   }
#   
#   changeSSBerr<-sum(changedif^2,na.rm=TRUE) 
#   
#   return(changeSSBerr)
# }

## SSB 

SSBerror <- function(Newcomm, 
                     meantsteps, 
                     Newcomm2,
                     effort, 
                     effort_unfished){ 
  
  # # trial 
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
  
  # error changesSSB - added inside this error becaseu outside would need to add orgument and the scripts stop working in pearcey 
  
  # do not compare yield from mycto - not there - and perch - there but not to be considered 
  toDelete<-c("myctophids","helicolenus barathri","nototodarus gouldi","trachurus declivis")
  temp<-Newcomm@params@species_params
  temp<-temp[-which(temp$species %in% toDelete),]
  changeObs <- temp$changesSSB
    
  changeMod<-round((getBiomass(Newcomm)[dim(effort)[1],]/getBiomass(Newcomm2)[(dim(effort_unfished)[1]),])*100)
  changeMod<-changeMod[- which(names(changeMod) %in% toDelete)]
    
  changedif<- log10(changeMod)-log10(changeObs)
    
  # relplace all inf/NaN/NA with NA before sum() - problems otherwise 
  if(length(changedif[!is.finite(changedif)])>0){
    changedif[!is.finite(changedif)]<-NA
  }
    
  changeSSBerr<-sum(changedif^2,na.rm=TRUE) 
    
  ssberr<-ssberr+changeSSBerr
  
  return(ssberr)
  
  }

## rank abundance  

rankerror <- function(Newcomm, 
                     meantsteps){ 
  
  # # # trial 
  # Newcomm = sim2
  # meantsteps = 10
  
  #### instead used for CPUE comparison - given taht rank abundance does not work and is not to trust much 
  
  # compare only spp for which you have values 
  # toKeep<-Newcomm@params@species_params[which(!is.na(Newcomm@params@species_params$rankAb)),]
  toKeep<-Newcomm@params@species_params[which(!is.na(Newcomm@params@species_params$cpue_gm3)),]
  
  # observed rank 
  # rankObs <- toKeep$rankAb
  rankObs <- toKeep$cpue_gm3
  
  # modelled rank
  rank<-getBiomass(Newcomm)
  
  if (is.na(meantsteps)){
    rankhat <- rank[dim(rank)[1],] 
  }else{
    rankhat <- apply(rank[(dim(rank)[1]-meantsteps+1):dim(rank)[1],],2,mean) 
  }
  
  # only for rank abundance 
  # rank<-rank(-rankhat)
  
  # compare only rank of spp you have info for 
  # rank<-rank[names(rank) %in% rownames(toKeep)]
  rankhat<-rankhat[names(rankhat) %in% rownames(toKeep)]
  
  # yield difference and error (form Reum)
  # rankdif<- log10(rank)-log10(rankObs)
  rankdif<- log10(rankhat)-log10(rankObs)
  
  # relplace all inf/NaN/NA with NA before sum() - problems otherwise 
  if(length(rankdif[!is.finite(rankdif)])>0){
    rankdif[!is.finite(rankdif)]<-NA
  }
  
  rankerr<-sum(rankdif^2,na.rm=TRUE) 
  
  return(rankerr)
}

# extinction

Extinct_test <- function(sim, 
                         extinct_threshold, 
                         meantsteps){
  
  # when considering declines in biomass instead of extinction - what happens if a species is extinct before you caonsider declines? possibly not happeining - otherwise may need to increase meantsteps
  
  # # trial 
  # sim = sim_unfished
  # extinct_threshold =  extinct_threshold
  # meantsteps = meantsteps
  
  Biomass <- getBiomass(sim)
  BiomassInit<- Biomass[1,] # option 1 - if you consider extinction 
  BiomassInit<- Biomass[dim(Biomass)[1]-meantsteps,] # option 2 # if you consider declines - biomass (30?) time steps before last time step - this should depend on t_max (set to 100 now) and on when communitiy stabilise (after the first few runs) 
  BiomassEnd<-  Biomass[dim(Biomass)[1],]
  
  relB<-(BiomassEnd/BiomassInit) # option 1 and 2
  relB<-relB[relB<extinct_threshold] # option 1 and 2 but different extinct_threshold values: 0.01 for extinction and e.g. 0.9 (biomass of last time step 90% biomass of n time steps before)
  extinct<- sum(1/relB)
  
  # if some spp are extinct or declining, this error needs to be big
  if(is.finite(extinct) & extinct>0){
    extinct = 500
  }

  # if this returs: 
  if(is.infinite(extinct)){
    extinct = 500 # this is simply a high number compared to the expected error ~ 3.5 given the conditions used for this model 
  }
  
  return(extinct)
}







