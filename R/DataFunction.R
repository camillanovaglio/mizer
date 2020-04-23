
# function for SEAmodel_run.r

########### select speceis ----

SpSelect<-function(df_param){
  
  df_param <- as.data.frame(df_param) # colnames are as the one reqiored by mizer
  
  # delete some spp
  all<-c("lanternfish","royal red prawn","blacktip cucumberfish","school whiting","sardine","arrow squid","bigeye ocean perch","mackerel","gurnards and latchet", "redfish","deepwater sharks","jackass morwong", "tiger flathead","dories and oreos","ocean jacket", "sea mullet","bight redfish","ribaldo", "blue warehou","elephantfish","orange roughy", "eastern rock lobster","blue grenadier","gulper sharks","silver warehou","black cardinalfish","southern chimaera","gemfish","silver trevally","sharpnose sevengill shark","eastern angelshark","snapper and breams","pink ling", "common sawshark","gummy shark","school shark","blueye trevalla","bight skate")
  
  toKeep<-c("lanternfish","arrow squid","bigeye ocean perch","mackerel","school whiting","redfish","deepwater sharks","jackass morwong","tiger flathead","dories and oreos", "blue warehou","orange roughy","blue grenadier","silver warehou","gemfish","pink ling","common sawshark", "gummy shark","school shark") #!!! bight redfish was included
  
  df_param<-df_param[which(df_param$spCommon %in% toKeep),]
  
  # add catcah data for squid, perch - there is no catch, and mackerel - don't have the fisheries so use AFMA info
  datCalibration_add<-df_log_spp %>%
    filter(YEAR>1994 & YEAR <2006, SPC_NAME %in% c("nototodarus gouldi","helicolenus barathri","trachurus declivis")) %>%
    group_by(SPC_NAME, YEAR) %>%
    dplyr::summarise(catchComm = sum(TOT_CATCH_KG, na.rm=TRUE)) %>%
    group_by(SPC_NAME) %>%
    dplyr::summarise(catchComm_t = round(mean(catchComm, na.rm=TRUE)/1e3)) 
  
  df_param[which(df_param$spCommon == "arrow squid"),"catchComm_t"]<-1307 
  
  return(df_param = df_param)
}

########### calculate Fmort ----

CreateFmort<-function(df_param, all_mortality, t_max, areaEco){

  names<-df_param[,which(colnames(df_param) %in% c("species","spCommon"))]
  
  # add OR mortality from SA as this one was from coco and included only cascade OR 
  or_fmort<-read.csv("/Users/nov017/Multispecies_sizebasedmodels/SouthEastAustralia/input/originalData/StockAssessment/OR_Fmort.csv")
  or_fmort$spe<-as.factor("hoplostethus atlanticus")
  all_mortality<-all_mortality[-which(all_mortality$spe=="hoplostethus atlanticus"),]
  all_mortality<-rbind(or_fmort, as.data.frame(all_mortality))
  
  # create an Fmort 
  Fmort <- all_mortality %>% 
    `colnames<-`(c("species","year","Fmort")) %>% 
    mutate(species = as.character(species)) %>% 
    right_join(names) %>% 
    ungroup() %>% 
    select(-spCommon)
  Fmort<-acast(Fmort, year ~ species)
  Fmort<-Fmort[,df_param$species]
  Fmort<-Fmort[as.character(1995:2017),]
  
  # add Fmort for ssp with no values - these will become reference values for the time variant bit below. mean values will change too  
  Fmort[Fmort==0]<-NA
  Fmort[,"mustelus antarcticus"]<-rep(0.13,dim(Fmort)[1]) # assessment to come, here ERA 
  Fmort[,"galeorhinus galeus"]<-rep(0.17,dim(Fmort)[1]) # assessment to come, here ERA
  FmortSqualus<-mean(c(0.47,0.4,0.6,0.43,0.21,0.26,0.3,0.5,0.35)) #,0.14,0.13))
  Fmort[,"squalus spp."]<-rep(FmortSqualus,dim(Fmort)[1])
  FmortSaw<-0.2
  Fmort[,"pristiophorus cirratus"]<-rep(FmortSaw,dim(Fmort)[1])
  Fmort_dories<-mean(c(0.45,0.31,0.1,0.35,0.37))
  Fmort[,"zeus faber"]<-rep(Fmort_dories,dim(Fmort)[1])
  
  # add background mortality to mesopelagic
  Fmort[,"nototodarus gouldi"]<- rep(0.01, dim(Fmort)[1])
  Fmort[,"trachurus declivis"]<- rep(0.01, dim(Fmort)[1])
  Fmort[,"helicolenus barathri"]<- rep(0.05, dim(Fmort)[1]) # state fishery? no catches
  
  # build trends in Fmort for thoses spp that don't have them - use initial values values as given above and relative changes in observed catches to build relative changes in Fmort from initial values  
  
  # get OBSERVED YIELD 
  yieldObs_timeVariant<-datValidationYieldSpp %>% 
    filter(species %in% df_param$species) %>% 
    mutate(catchComm_t = (catchComm_t*1000000)/areaEco) # transfor in g m3
  
  # add catches for tracurus - these are from reports - 2010 2018 
  tr<-data.frame(species = "trachurus declivis", year = seq(1995, 2017), catchComm_t = c(100, 300, 500, 9800, 5000, 1000, 100, 300, 2300, 3800,2000, 800, 700, 200, 300, 200, 50,50,50,300,6000,4000,3000))
  
  # update mean cathes in df_param for calibration 
  df_param[which(df_param$spCommon == "mackerel"),"catchComm_t"]<-mean(tr[1:11,"catchComm_t"])
  
  # change units adn add to observed trends in yields
  tr$catchComm_t<-(tr$catchComm_t*1000000)/areaEco
  tr$species<-as.character(tr$species)
  yieldObs_timeVariant<-yieldObs_timeVariant %>% 
    filter(species != "trachurus declivis") %>% 
    bind_rows(tr)
  
  yieldObs_timeVariant<-acast(yieldObs_timeVariant, year ~ species)
  yieldObs_timeVariant<-yieldObs_timeVariant[-nrow(yieldObs_timeVariant),]
  yieldObs_timeVariant<-yieldObs_timeVariant[,df_param$species]
  
  # get relative changes in observed yield to be used as Fmort
  relative_Fmort<-yieldObs_timeVariant
  relative_Fmort[is.na(relative_Fmort)]<-0
  ref<-relative_Fmort[1,]
  
  # OPTION of calibrationna dbefore adjusting r_max etc. 
  # # ref points a bit low for these spp resulting in unrealistically high Fmort 
  # ref["pristiophorus cirratus"]<-relative_Fmort[5,"pristiophorus cirratus"]
  # ref["mustelus antarcticus"]<-relative_Fmort[6,"mustelus antarcticus"] 
  # ref["galeorhinus galeus"]<-relative_Fmort[3,"galeorhinus galeus"] 
  ref["squalus spp."]<-relative_Fmort[6,"squalus spp."]
  ref["trachurus declivis"]<-relative_Fmort[6,"trachurus declivis"]
  
  # to increase fishing mortality
  ref["pristiophorus cirratus"]<-relative_Fmort[18,"pristiophorus cirratus"]
  # need to change reference level... 
  ref["mustelus antarcticus"]<-relative_Fmort[3,"mustelus antarcticus"]
  ref["galeorhinus galeus"]<-relative_Fmort[18,"galeorhinus galeus"]
  # ref["squalus spp."]<-relative_Fmort[8,"squalus spp."] 
  # ref["trachurus declivis"]<-relative_Fmort[1,"trachurus declivis"]
  
  relative_Fmort<-sweep(relative_Fmort,2,ref,"/")
  
  Fmort[is.na(Fmort)]<-0
  relative_Fmort[,c("sillago flindersi","centroberyx affinis" ,"nemadactylus macropterus","platycephalus richardsoni","seriolella brama","hoplostethus atlanticus","macruronus novaezelandiae", "seriolella punctata","rexea solandri","genypterus blacodes","helicolenus barathri")]<-1
  
  Fmort<-relative_Fmort*Fmort
  
  # write.csv(Fmort, file = "/Users/nov017/Desktop/Fmort.csv",row.names=FALSE)
  
  
  
  # CAN DELETE depending on trials in SEAmodel_run.Rmd

  # add missing Fmort values for recent years based on catches as above
  trial<-acast(datValidationYieldSpp, year ~ species)
  trial<-trial[-nrow(trial),]
  trial<-trial[,df_param$species]
  trial[is.na(trial)]<-0
  # dim(trial)

  Fmort2<-Fmort

  # seriorella
  ref<-trial[14,] # per seriorella b
  trial_s<-sweep(trial,2,ref,"/")
  Fmort2[14:23,"seriolella brama"]<-Fmort2[14,"seriolella brama"]*trial_s[14:23,"seriolella brama"]

  # sillago
  ref<-trial[22,]
  trial_s<-sweep(trial,2,ref,"/")
  Fmort2[22:23,"sillago flindersi"]<-Fmort2[22,"sillago flindersi"]*trial_s[22:23,"sillago flindersi"]

  # centroberyx affinis
  ref<-trial[22,]
  trial_s<-sweep(trial,2,ref,"/")
  Fmort2[22:23,"centroberyx affinis"]<-Fmort2[22,"centroberyx affinis"]*trial_s[22:23,"centroberyx affinis"]

  # nemadactylus macropterus
  ref<-trial[20,]
  trial_s<-sweep(trial,2,ref,"/")
  Fmort2[20:23,"nemadactylus macropterus"]<-Fmort2[20,"nemadactylus macropterus"]*trial_s[20:23,"nemadactylus macropterus"]

  # platycephalus richardsoni
  ref<-trial[21,]
  trial_s<-sweep(trial,2,ref,"/")
  Fmort2[21:23,"platycephalus richardsoni"]<-Fmort2[21,"platycephalus richardsoni"]*trial_s[21:23,"platycephalus richardsoni"]

  # macruronus novaezelandiae
  ref<-trial[18,]
  trial_s<-sweep(trial,2,ref,"/")
  Fmort2[18:23,"macruronus novaezelandiae"]<-Fmort2[18,"macruronus novaezelandiae"]*trial_s[18:23,"macruronus novaezelandiae"]

  # seriolella punctata
  ref<-trial[20,]
  trial_s<-sweep(trial,2,ref,"/")
  Fmort2[20:23,"seriolella punctata"]<-Fmort2[20,"seriolella punctata"]*trial_s[20:23,"seriolella punctata"]

  # rexea solandri
  ref<-trial[21,]
  trial_s<-sweep(trial,2,ref,"/")
  Fmort2[21:23,"rexea solandri"]<-Fmort2[21,"rexea solandri"]*trial_s[21:23,"rexea solandri"]

  # genypterus blacodes
  ref<-trial[19,]
  trial_s<-sweep(trial,2,ref,"/")
  Fmort2[19:23,"genypterus blacodes"]<-Fmort2[19,"genypterus blacodes"]*trial_s[19:23,"genypterus blacodes"]

  Fmort<-Fmort2

  # end CAN DELETE
  
  
  
  
  
  # calculate mean Fmort constant across years 
  Fmort_mean<-Fmort[as.character(1995:2005),] 
  Fmort_mean<-colMeans(Fmort_mean, na.rm=TRUE)
  Fmort_mean[is.nan(Fmort_mean)]<-0
  Fmort_mean<-unname(Fmort_mean)
  
  matrix_effort<-t(matrix(c(rep(0,300*19), rep(Fmort_mean, 300*19)), nrow = length(df_param$species), ncol = t_max))
  rownames(matrix_effort)<-seq(1:t_max)
  colnames(matrix_effort)<-df_param$species
  
  return(list(df_param = df_param, Fmort = Fmort, Fmort_mean = Fmort_mean, matrix_effort = matrix_effort, yieldObs_timeVariant = yieldObs_timeVariant))
  
}

########### calculate theta ----

CalcTheta<-function(df_param, df_ismp_spp){
  
  # Method from old version:
  # source("/Users/nov017/multispecies_sizebasedmodels/SouthEastAustralia/R/GETtheta.R")
  # res = GETtheta(df_kapCatch, df_param)
  
  load("/Users/nov017/Multispecies_sizebasedmodels/SouthEastAustralia/input/SESSF_ISMP_Cleaned_noCanatus.RData")

  sp_name<-df_param %>%
    select("species","spCommon") %>%
    `colnames<-`(c("SPC_NAME","spCommon"))

  occ<-df_ismp_spp %>%
    right_join(sp_name) %>%
    filter(!is.na(ACTIVITY_ID)) %>%
    select(ACTIVITY_ID,SPC_NAME,SPC_WT) %>%
    group_by(ACTIVITY_ID, SPC_NAME) %>%
    dplyr::summarise(SPC_WT = sum(SPC_WT))

  n_tows<-length(unique(occ$ACTIVITY_ID))
  
  occ<-as.data.frame(occ)
  occ<-acast(occ, SPC_NAME~ACTIVITY_ID)
  occ[is.na(occ)]<-0
  occ[occ>0]<-1 # presence-absence data instead of abundance (as used before)

  # Probability of cooccurrence
  library("cooccur")
  cooccur<-cooccur(mat=occ, type = "spp_site", spp_names = TRUE)
  cooccur<-cooccur$results

  # probability of cooccurrence = expected cooccurrence/n_tows
  cooccur<-cooccur[, which(colnames(cooccur) %in% c("sp1_name","sp2_name", "prob_cooccur"))]

  # need to increase values otherwise there is no interaction
  cooccur$prob_cooccur = cooccur$prob_cooccur*10
  cooccur$prob_cooccur = ifelse(cooccur$prob_cooccur >1, 1, cooccur$prob_cooccur)

  cooccur$sp1_name<-as.character(cooccur$sp1_name)
  cooccur$sp2_name<-as.character(cooccur$sp2_name)

  # add sp not in ismp
  add<-df_param[-which(df_param$species %in% unique(c(cooccur$sp1_name,cooccur$sp2_name))),"species"]
  add<-expand.grid(sp1_name = add, sp2_name= unique(c(add, cooccur$sp1_name, cooccur$sp2_name)))
  add$prob_cooccur<-1
  cooccur<-rbind(cooccur, add)

  # add the diagonal - cannibalism
  add2<-data.frame(sp1_name=unique(rownames(occ)), sp2_name=unique(rownames(occ)), prob_cooccur=0.6)
  cooccur<-rbind(cooccur, add2)

  # rearrange occurr to include all combinations
  cooccur1<-cooccur
  colnames(cooccur)<-c("prob_cooccur","sp2_name", "sp1_name")
  cooccur<-rbind(cooccur, cooccur1)
  cooccur<-unique(cooccur)

  # order cooccur as df_param sp list
  df_param<-df_param[order(df_param$w_inf),]
  cooccur$sp2_name<-factor(cooccur$sp2_name, c(paste(df_param$species)))
  cooccur$sp1_name<-factor(cooccur$sp1_name, c(paste(df_param$species)))

  # trnsform in matrix
  library(reshape2)
  theta<-acast(cooccur, sp1_name~sp2_name, value.var="prob_cooccur")

  # taking into account food ecologies:

  # theta[,c("myctophids","nototodarus gouldi","trachurus declivis")]<- 1 # Always preys
  # theta[c("myctophids","nototodarus gouldi","trachurus declivis"),]<- 0 # never predators - no predation on small fish of all spp (larvae)
  # Planktivory, bottom feeders and sp not feeding on other fishes. (this should relay on PP and BB)
  # theta[c("nemadactylus macropterus","seriolella punctata","seriolella brama"),]<-0 # bigone starv after calibration - should adjust beta

  # 2 using diet data to inform some spp

  # from GETpesciData.R in GETbetaWrap.R in SEA_EcolParam
  source("/Users/nov017/multispecies_sizebasedmodels/SouthEastAustralia/R/GETpesciData.R")
  res = GETpesciData(df_param)
  df_pesci <- res$df_pesci
  df_pesci2<-df_pesci %>%
    mutate(pred = tolower(pred)) %>%
    mutate(prey = tolower(prey)) %>%
    filter(pred %in% df_param$species) %>%
    mutate(n = 1) %>%
    select("pred","prey","n") %>%
    group_by(pred,prey) %>%
    dplyr::summarise(n = sum(n)) # %>%
    # spread(prey,n) # better 1 by 1

  toMatch<-c("myctophids","sillago","flindersi","nototodarus","gouldi","helicolenus","barathri","trachurus","declivis","centroberyx","affinis","squalus","nemadactylus","macropterus","platycephalus","richardsoni","zeus","faber","seriolella","brama","hoplostethus","atlanticus","macruronus","novaezelandiae","punctata","rexea","solandri","genypterus","blacodes","pristiophorus","cirratus","mustelus","antarcticus","galeorhinus","galeus")

  matches <- unique (grep(paste(toMatch,collapse="|"),
                          df_pesci2$prey, value=TRUE))

  df_pesci2<-df_pesci2 %>%
    filter(prey %in% matches)

  # increase these interactions
  theta["centroberyx affinis","helicolenus barathri"]<-1
  theta["hoplostethus atlanticus","macruronus novaezelandiae"]<-1
  theta["macruronus novaezelandiae","macruronus novaezelandiae"]<-1
  theta["platycephalus richardsoni","genypterus blacodes"]<-1
  theta["platycephalus richardsoni","sillago flindersi"]<-1
  theta["platycephalus richardsoni","platycephalus richardsoni"]<-1
  theta["squalus spp.","platycephalus richardsoni"]<-1
  theta["zeus faber","centroberyx affinis"]<-1
  theta["zeus faber","helicolenus barathri"]<-1
  theta["zeus faber","sillago flindersi"]<-1
  theta["zeus faber","macruronus novaezelandiae"]<-1
  
  # write.csv(theta,"/Users/nov017/Desktop/feedingTrial/theta_trial.csv")
  
  # explore PPMR
  # i=6  
  # bell<-function(wp) exp((-(log(df_param2[i,"w_inf"]/(wp*df_param2[i,"beta"])))^2)/2*(df_param2[i,"sigma"])^2)
  # 
  # ggplot(data.frame(wp=c(0:10)), aes(wp)) +
  #   stat_function(fun = bell, geom="line")+
  #   labs(x = "Prey weight", y="Predator preference", caption = 'Prey-weight selection function')
  
  return(theta = theta)
  
}

########### modify df_param ----

modDf_Param<-function(df_param, kappa){
  
  df_param <- df_param %>%
    
    # OPTION 1-3
    # H and KS (and GAMMA)
    # default h (constant for maximum food intake/consumption) = [3*k_vb/(alfa*f0)]*W_inf^(1/3)
    # default ks (constant for standard metabolism and activity) = 20% of h
    # default gamma (search parameter) = f0*h*beta^(2-lambda)/((1-f0)*2*3.14*kappa*sigma) # vignette pp 10 vignette
    mutate(h = 50 * (w_inf/1000)^0.15) %>% # h = constant for maximum food intake
    mutate(ks = 20 * w_inf^(-0.25)) %>% # ks = constant for standard metabolism and activity
    mutate(ks = ifelse(species == "myctophids", 5, ks)) %>% # however this relationship gives Trachinops (in my case myctophids) very high ks, so I decrease it a bit
    # OPTION 5 
    # mutate(ks = 0.15*h) %>% # closer to mizer h*0.2
    # h drives everithing - try keeping h as per Asta - problems - see option 5 documentation  
    # select(-c("h","ks")) %>%
    
    # R MAX and Q
    mutate(r_max = kappa * w_inf^(-1.5)) %>% # if you are using r_max instead, here is Julia's formula for deep sea
    
    # BETA and SIGMA
    # OPTION 4 - use your values
    # mutate(beta = beta/exp((2-2*n+q)/sigma^2)) %>% # using correnction factor in Bluanchard paper but not sure as she divides PPMR by this nnumber instead of beta - nothing changes anyway...
    # OPTION 1-3
    mutate(beta = ifelse(species %in% c("nemadactylus macropterus","seriolella brama","seriolella punctata"), 1000,100)) %>%
    mutate(sigma = ifelse(sigma> 2, 2, sigma)) %>%
    
    # PP and BB
    mutate(avail_PP = ifelse(species %in% c("myctophids","trachurus declivis"),1,0.3)) %>%
    mutate(avail_BB = ifelse(species %in% c("myctophids","trachurus declivis"),0,0.7))
  
  # according to Beth's
  df_param <- df_param %>%
    
    # R MAX and Q to fix abundances
    # to macth catches observed vs modelled catches and to have a Q param that can be calibrated
    # mutate(scalar_r_max = c(30, 1, 150, 4, 160, 10, 6, 3, 40, 10, 1, 1, 5, 1, 1, 1, 0.1, 0.2, 1)) %>% # OPTION 1
    mutate(scalar_r_max = c(30, 1, 150, 4, 160, 10, 0.5, 3, 40, 20, 1, 10, 15, 10, 1, 10, 0.1, 0.2, 1)) %>% # OPTION 2 # decrease squalus, increase dories, OR, blue granadier and seriorella p. and genypt
    # mutate(scalar_r_max = c(30, 1, 150, 4, 160, 10, 6, 3, 40, 10, 1, 10, 20, 1, 1, 1, 0.1, 0.2, 1)) %>% # OPTION 3 # increase blue gran and OR
    mutate(catchability = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)) %>%
    # mutate(catchability = c(1, 0.5, 0.5, 1, 1, 0.5, 0.5, 0.5, 0.8, 1, 1, 1, 1, 0.5, 1, 1, 1, 1.5, 0.5)) %>%
    mutate(r_max = r_max * scalar_r_max) %>%
    
    # W_MIN
    # to distribute food over size classes - i.e. incresing w_min increases feeding on other spp rather than squids, because abundances at this size, which are the ones required by other spp, increase
    mutate(w_min = ifelse(species == "myctophids", 0.001, ifelse(w_min < 1, 0.005, w_min))) %>%
    mutate(avail_PP = ifelse(species %in% c("myctophids","trachurus declivis","nototodarus gouldi"),0.6,0.2)) %>%
    mutate(avail_PP = ifelse(species %in% c("myctophids"),1,avail_PP)) %>%
    mutate(avail_BB = ifelse(species %in% c("myctophids","trachurus declivis","nototodarus gouldi"),0,0.4))
  
  # erepro? what's the meaning of this param? and what's the def value? mizer guide says 1,
  # df_param <- df_param %>%
  #   mutate(erepro = ifelse(species == "nototodarus gouldi", 0.2, ifelse(species %in% c("pristiophorus cirratus","mustelus antarcticus","galeorhinus galeus"), 0.05, 0.1))) # use opposite specification - "squalus spp." very problematic... this changes things... 
  # df_param$erepro<-1 # thsi does not change things a lot, but for squalus: higher decreases when fished if rerpro = 0.5 instead of 1 - i.e. lower reproduction? no changes for the other sharks  
  
  return(df_param = df_param)
  
}

########### modify r_max ----

modQ<-function(df_param){
  
  # # use the same catcahbility that you've explore using shiny and when the FD was on... 
  # df_param[2,"catchability"]<-0.4 # whiting # OK
  # df_param[3,"catchability"]<-0.5 # squid # great
  # df_param[6,"catchability"]<-0.3 # redfish # OK
  # df_param[7,"catchability"]<-0.1 # deepw shark # change effort ?
  # df_param[8,"catchability"]<-0.3 # morwong # OK
  # df_param[9,"catchability"]<-0.05 # platy # OK
  # df_param[11,"catchability"]<-0.02 # b war. # OK
  # df_param[13,"catchability"]<-0.8 # blue g. # great
  # df_param[14,"catchability"]<-0.2 # silver w. # great
  # df_param[18,"catchability"]<-1 # gummy
  # df_param[19,"catchability"]<-0.1 # school # great
  
  # use the r_max that you've explored using shiny
  # df_param<-df_param %>%
  #   mutate(r_max = r_max*2)
 
  df_param<-df_param %>% 
    mutate(r_max = case_when(spCommon == "lanternfish" ~ r_max,
                             spCommon == "school whiting" ~ r_max,
                             spCommon == "arrow squid" ~ r_max,
                             spCommon == "bigeye ocean perch" ~ r_max,
                             spCommon == "mackerel" ~ r_max*1.5,
                             spCommon == "redfish" ~ r_max*1.5,
                             spCommon == "deepwater sharks" ~ r_max,
                             spCommon == "jackass morwong" ~ r_max,
                             spCommon == "tiger flathead" ~ r_max*0.7,
                             spCommon == "dories and oreos" ~ r_max*1.5,
                             spCommon == "blue warehou" ~ r_max,
                             spCommon == "orange roughy" ~ r_max*1.5,
                             spCommon == "blue grenadier" ~ r_max,
                             spCommon == "silver warehou" ~ r_max,
                             spCommon == "gemfish" ~ r_max*1.2,
                             spCommon == "pink ling" ~ r_max*1.3,
                             spCommon == "common sawshark" ~ r_max*2,
                             spCommon == "gummy shark" ~ r_max*7,
                             spCommon == "school shark" ~ r_max*0.4))
  
  return(df_param = df_param)
}

########### modify Fmort ----

modFmort<-function(matrix_effort){
  
  # option 1
  # cahnges to effort matrix according to Beth - see comments below. not sure the 2nd group had this combindation
  
  # # essential to crush sharks a bit, also high abundance and low yields after fishing
  # matrix_effort[,"squalus spp."]<-matrix_effort[,"squalus spp."]*3.5
  # matrix_effort[,"mustelus antarcticus"]<-matrix_effort[,"mustelus antarcticus"]*2 # historical depletion
  # matrix_effort[,"galeorhinus galeus"]<-matrix_effort[,"galeorhinus galeus"]*1.5 # historical depletion
  # 
  # # too hign abundance after fishing and too low yields
  # matrix_effort[,"hoplostethus atlanticus"]<-matrix_effort[,"hoplostethus atlanticus"]*6
  # matrix_effort[,"macruronus novaezelandiae"]<-matrix_effort[,"macruronus novaezelandiae"]*3
  # matrix_effort[,"centroberyx affinis"]<-matrix_effort[,"centroberyx affinis"]*2
  # 
  # # too low yields but alos low abundance after fishing...
  # matrix_effort[,"genypterus blacodes"]<-matrix_effort[,"genypterus blacodes"]*1.5
  # matrix_effort[,"seriolella punctata"]<-matrix_effort[,"seriolella punctata"]*1.5
  
  # # option 2-4
  matrix_effort[,"squalus spp."]<-matrix_effort[,"squalus spp."]*2 # 3.5 
  # matrix_effort[,"mustelus antarcticus"]<-matrix_effort[,"mustelus antarcticus"]*2 # historical depletion
  # matrix_effort[,"galeorhinus galeus"]<-matrix_effort[,"galeorhinus galeus"]*1.5 # historical depletion
  matrix_effort[,"macruronus novaezelandiae"]<-matrix_effort[,"macruronus novaezelandiae"]*3
  
  return(matrix_effort = matrix_effort)
}

########### modify theta ----

modTheta<-function(theta){
  
  # decrease predation on flathead that is declining a lot after fishing
  theta["platycephalus richardsoni",]<-theta["platycephalus richardsoni",]*0.3
  
  # improve interaction - increase feeding on forage
  theta[,"helicolenus barathri"]<-theta[,"nototodarus gouldi"]
  theta[,"trachurus declivis"]<-theta[,"nototodarus gouldi"]
  
  # OPTION 3 and 7 - add this 
  # macro-squalus interaction: the higher macrurus abundance, the lower squalus biomass as predation here is high 
  theta["macruronus novaezelandiae","squalus spp."]<-theta["macruronus novaezelandiae","squalus spp."]*0.5
  
  # # OPTION ? to consider as these spp are starving and according to literature - but more problems if I do this
  # theta["macruronus novaezelandiae","myctophids"]<-1
  # theta["macruronus novaezelandiae","nototodarus gouldi"]<-1
  # theta["macruronus novaezelandiae","trachurus declivis"]<-1
  # 
  # theta["hoplostethus atlanticus","myctophids"]<-1
  # theta["hoplostethus atlanticus","nototodarus gouldi"]<-1
  # theta["hoplostethus atlanticus","trachurus declivis"]<-1
  #
  # theta["galeorhinus galeus","nototodarus gouldi"]<-1
  # # AFMA Fish, squid, octopus and benthic invertebrates.
  # # also lower cannobalism from 0.6 to 0.2 for these spp?
  # theta["macruronus novaezelandiae","macruronus novaezelandiae"]<-0.2
  # theta["hoplostethus atlanticus","hoplostethus atlanticus"]<-0.2
  # theta["galeorhinus galeus","galeorhinus galeus"]<-0.2
  
  # ALSO! squalus is too dependent on small pelagice - i.e. when they increase due to unfish situation - squalus increases too. 
  
  # plot predatro-prey matrix
  theta_plot<-as.data.frame(theta)
  theta_plot$predator<-rownames(theta_plot)
  theta_plot<-melt(theta_plot)
  theta_plot$predator<-factor(theta_plot$predator, levels = levels(theta_plot$variable))
  
  p2 <- ggplot(theta_plot, aes(variable, predator)) + 
    geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = "white",high = "blue")+
    theme_bw()+
    theme(text = element_text(size=18),
          # axis.title.y = element_text(vjust=0.4),
          # axis.title.x = element_text(vjust=0.3),
          axis.text.x = element_text(angle=90, hjust=0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  return(list(theta = theta, p_theta = p2))
  
}

########### plot mod vs obs eql ----

compareEql<-function(areaEco, df_param, sim, matrix_effort, df_log_spp, calibration){  
  
  # sim = sim_calibrated
  # matrix_effort = effort
  
  # 1. prepare data 
  # this could be a separate step/function but you'd need to separate all calls using sim as argument. 
  
  ### comparison with ATLANTIS abundance 
  # see compareB0.xls for values from Beth 
  compare_ab<-data.frame(species = df_param$species, obs = round(c(2149027.8, 45485.12221, 1594790.423,0, 598937.9633, 249000.592, 953398.2439, 14866.67109, 57652.67201, 745870.2626, 189868.1142, 485537.1347, 102083.8022, 32721.00727, 10667.09547, 19279.65539, 0, 21350.21696, 207126.9421)), mod = round((getBiomass(sim)[dim(matrix_effort)[1],]*areaEco)/1000000))
  compare_ab$theme<-"compare_ab"
  compare_ab<-compare_ab[match(df_param$species, compare_ab$species),]
  
  ### CHANGES IN SSB - changes in biomass due to fishing and compare with SA 
  compare_changeSSB<-data.frame(species = df_param$species, obs = c(NA, 50, NA, NA, NA, 15, 30, 40, 40, 40, 20, 15, 100, 25, 15, 25, 30, 30, 25), mod=round((getSSB(sim)[dim(matrix_effort)[1],]/getSSB(sim)[(dim(matrix_effort)[1]/2)-1,])*100))
  compare_changeSSB$theme<-"compare_changeSSB"
  compare_changeSSB<-compare_changeSSB[match(df_param$species, compare_changeSSB$species),]
  
  ### RANK abundance using CPUE from ISMP data 
  load("/Users/nov017/Multispecies_sizebasedmodels/SouthEastAustralia/input/SESSF_ISMP_Cleaned.RData")
  
  sampling<-df_ismp_spp_LF %>% 
    filter(SPC_NAME %in% df_param$species,YEAR >1994 & YEAR <2006) %>%
    filter(!SPC_NAME %in% c("nototodarus gouldi","trachurus declivis")) %>% 
    group_by(SPC_NAME) %>% 
    dplyr::summarise(n = length(unique(ACTIVITY_ID)))
  
  rank<-df_ismp_spp_LF %>% 
    filter(SPC_NAME %in% df_param$species,YEAR >1994 & YEAR <2006) %>% 
    filter(!SPC_NAME %in% c("nototodarus gouldi","trachurus declivis")) %>% 
    mutate(ab = WT_CLASS*SPC_CNT) %>% 
    group_by(SPC_NAME) %>% 
    dplyr::summarise(ab = sum(ab)) %>% 
    right_join(sampling) %>% 
    mutate(ab = ab/n) %>% 
    select(-n) 
  
  colnames(rank)<-c("species","rankAb")
  rank<-rank %>% 
    left_join(df_param) %>% 
    select(c(species,rankAb)) %>% 
    mutate(rankAb = rank(-rankAb)) %>% 
    mutate(rankAb = rankAb+3) # you are not considering mesopelagics
  
  rankMod<-data.frame(species = df_param$species, mod = getBiomass(sim)[dim(matrix_effort)[1],]) %>% 
    mutate(mod = rank(-mod))
  
  compare_rank<-merge(rank, rankMod, all=TRUE)
  compare_rank$theme<-"rank"
  compare_rank<-compare_rank[match(df_param$species, compare_rank$species),]
  
  ### YIELD 
  yobs<-sim@params@species_params$catchComm_t # these are already rescaled in gm3 outside teh function
  ymod<- getYield(sim) # modelled yields are in g m3 
  ymod<-ymod[dim(matrix_effort)[1],] # last time step 
  
  compare_y<-data.frame(species = df_param$species,obs = yobs, mod = ymod) 
  compare_y$theme<-"compare_y"

  # SSB 
  # data from Punt files (coco's folder)
  ssbObs<-read.csv("/Users/nov017/Multispecies_sizebasedmodels/SouthEastAustralia/input/originalData/StockAssessment/corentin/SS_outputs/all_info.csv")
  
  # what about OR east stock?! need to take values from the recovery chapter as in here values are for cascade only
  ssbOR<-read.csv("/Users/nov017/Publications/Smith2018_recovery/camilla/data/20170724_SSB2.csv")
  ssbOR<-ssbOR %>%  filter(sp == "Orange Roughy East")  
  colnames(ssbOR)<-c("Year","Spawning.biomass","spCommon")
  ssbOR$spCommon<-as.factor("orange roughy")
  ssbOR$species<-as.factor("hoplostethus atlanticus") # values are for female only and doubled below 

  # delete cascade
  ssbObs<-ssbObs %>% filter(spCommon != "orange roughy cascade")

  # add OR east stock 
  ssbObs<-merge(ssbObs, ssbOR, all =TRUE)
  
  ssbObs<-ssbObs %>%
    select(Year,Spawning.biomass, species, spCommon) %>% 
    mutate(Spawning.biomass = ifelse(Spawning.biomass == -1, NA, Spawning.biomass)) %>%
    group_by(species, Year) %>%
    dplyr::summarize(Spawning.biomass = sum(Spawning.biomass, na.rm=TRUE)) %>% 
    `colnames<-`(c("species", "year", "SSB")) %>% 
    mutate(SSB = ifelse(SSB == 0, NA, SSB)) %>% 
    mutate(gender = ifelse(species %in% c("seriolella brama","centroberyx affinis","seriolella punctata"),1,2)) %>% # for some sp the SSB is that of female and male, wherease for others is that of female only. The one with modelGender == 2 are the one with only female and their values is doubled for consistency. see jemery table 
    mutate(SSB = ifelse(gender == 2, SSB*2, SSB)) %>% 
    select(-gender)
  
  ssbObsCalibration<-ssbObs %>% 
    filter(year > 1994 & year< 2006) %>%
    group_by(species) %>%
    dplyr::summarize(ssb = mean(SSB, na.rm=TRUE)) %>% 
    mutate(ssbObs = ssb*scaling_total) %>% # convert biomass in g m3
    select(-ssb) %>% 
    mutate(species = as.character(species))
  
  ssbm<-getSSB(sim)
  ssbm<-data.frame(species = df_param$species, ssbm=ssbm[dim(matrix_effort)[1],]) 
  ssbm$species<-as.character(ssbm$species)
  
  allSSB<- merge(ssbObsCalibration, ssbm, all=TRUE)
  allSSB$theme<-"SSB"
  colnames(allSSB)<-c("species","obs","mod","theme")

  # CPUE based on logbook data 
  catch <- df_log_spp %>%
    filter(SPC_NAME %in% df_param$species,
           YEAR > 1994 & YEAR< 2006) %>%
    group_by(SPC_NAME, YEAR) %>%
    dplyr::summarize(catch_g = sum(TOT_CATCH_KG, na.rm=TRUE)*1000) %>% # total g in system 
    group_by(SPC_NAME) %>%
    dplyr::summarize(catch_g = mean(catch_g, na.rm=TRUE)) # mean catch in g between 1995-2005 by spp
  
  effort_obs<- df_log_spp %>%
    filter(SPC_NAME %in% df_param$species,
           YEAR > 1994 & YEAR< 2006) %>%
    group_by(SPC_NAME, YEAR) %>%
    dplyr::summarize(opn = n()) %>%
    group_by(SPC_NAME) %>%
    dplyr::summarize(opn = mean(opn, na.rm=TRUE)*14000*10*20) # mean m3 trawled for each spp
  
  cpue<-catch %>%
    left_join(effort_obs) %>%
    mutate(cpue_gm3 = catch_g/opn) %>%  # these should be g m3
    select(SPC_NAME, cpue_gm3) %>% 
    `colnames<-`(c("species", "cpue_gm3")) %>% 
    right_join(select(df_param, c(species, spCommon))) %>% 
    select(-spCommon)
  
  compare_cpue <- data.frame(species = df_param$species, obs = cpue$cpue_gm3, mod = getBiomass(sim)[dim(matrix_effort)[1],]*10000)
  compare_cpue <- data.frame(species = df_param$species, obs = cpue$cpue_gm3/10000, mod = getBiomass(sim)[dim(matrix_effort)[1],])
  compare_cpue[c(3,5), "obs"]<-NA 
  compare_cpue$theme<-"compare_CPUE"
  
  #### merge and plot all data to compare
  compare<-merge(compare_y,allSSB, all=TRUE) # yield abd SSB
  # compare<-merge(compare,compare_ab, all=TRUE ) # Atlantis
  # compare<-merge(compare, compare_changeSSB, all=TRUE) # cahnges in SSB
  # compare<-merge(compare, compare_cpue, all=TRUE) # cpue from logbooks 
  # colnames(compare_rank)<-c("species","obs","mod","theme") 
  # compare<-merge(compare, compare_rank, all=TRUE)
  compare<-compare[-which(compare$species %in% c("myctophids","helicolenus barathri")),]
  compare<-droplevels(compare)

# 2. plot 

  p<-ggplot(compare, aes(x = log(obs), y = log(mod), label = species))+
    geom_point()+
    geom_text()+
    geom_smooth(method='lm',formula=y~x)+
    geom_abline(intercept = 0, slope = 1)+
    facet_wrap(~theme, scale= "free")  

# 3. add info in df_param and sim for calibration routine. you need this in many cases. some is before calibration, others are for chacking. so you add all of this info before calibration in SEAmodel_run.Rmd  

  if(calibration == TRUE){

    # # TO DO IF decide to use changes in SSB in calibration
    df_param$changesSSB<-compare_changeSSB[,"obs"]
    sim@params@species_params$changesSSB<-df_param$changesSSB
    # # TO DO IF decide to use rank abundance in calibration. also in sim as that will be the starting community
    df_param$rankAb<-compare_rank$rankAb
    sim@params@species_params$rankAb<-df_param$rankAb
    # TO DO IF decide for this option 
    temp<-df_param %>% 
      left_join(ssbObsCalibration)
    df_param$ssbObs<-temp$ssbObs
    sim@params@species_params$ssbObs<-df_param$ssbObs
    # # TO DO IF decide to use CPUE from logbook instead of catches 
    df_param$cpue_gm3<-cpue$cpue_gm3/10000
    sim@params@species_params$cpue_gm3<-df_param$cpue_gm3
    
  }
  
  return (list (df_param = df_param, sim = sim, plotYield_sp = p, compare = compare, ssbObs = ssbObs))
  
}

########### plot mod vs obs time ----

compareTrends<-function(sim,sim_FD,fleetDynamics,type,yieldObs_timeVariant,ssbObs,rescale,areaEco){
  
  # need to fix for the FD part: 
  # other inputs 
  # sim_FD = sim_fitted_FD
  # time of interest? e.g. 1995:2017 for the Fmort part and 2006:2017 for the FD
  # comparison option: scale = 'fleet_yield', "fleet_effort". Need "species"? this component needs both sim and sim_FD for the FD part as you want to compare them all 
  # areaEco, datValidationYield and datValidationEffort
  
  ######## comparison options
  # 1. compare yield by spp: sim vs sim_FD vs observed
  # 2. compare yield by fleet: sim_FD vs observed 
  # 3. compare effort by fleet: sim_FD vs observed
  
  # 4. also need to add option to do it for SSB instead of yield for plot option 2 and 3
  
  ######## scaling options: 
  # 1. rescale yield 0-1 for each species independently - most importnat for trends  
  # 2. rescale yield 0-1 across species - i.e. relative abundance - meaningless
  # 3. rescale modelled and observed yield 0-1 together for each speceis - most important for magnitude of agreements 
  # 4. raw values - rescale = 0 - most meaningful as overall picture
  
  ######## comparison option 1 

  # sim = sim_fitted # this might be NA if you don't care about
  # sim_FD = sim_FD_bmsy # this needs to be NA in
  # fleetDynamics = TRUE
  # type = "SSB"
  # yieldObs_timeVariant
  # ssbObs
  # rescale = 0
  # areaEco

  # sim 
  if(type == "yield"){
    yield_sim<-getYield(sim)
    yield_sim<-yield_sim[as.character(1995:2017),] 
  }else {
    yield_sim<-getSSB(sim)
    yield_sim<-yield_sim[as.character(1995:2017),] 
  }

  # function to rescale data according to rescale option 
  fn_rescale<-function(df){
    if(rescale == 1){
      library(matrixStats)
      max = colMaxs(df, na.rm =TRUE)
      min = colMins(df, na.rm =TRUE)
      ratio = max-min
      df <- sweep(df,2,min,"-")
      df <- sweep(df,2,ratio,"/")
    } else if(rescale == 2){
      max = max(df, na.rm =TRUE) 
      min = min(df, na.rm =TRUE)
      ratio = max-min
      df <- sweep(df,2,min,"-")
      df <- sweep(df,2,ratio,"/")
    } else if(rescale == 3){
      # owful code but working...
      trial<-df %>%
        spread(species,yield)
      trial2<-trial %>%
        select(-c("year","color")) %>%
        as.matrix()
      library(matrixStats)
      max = colMaxs(trial2, na.rm =TRUE) 
      min = colMins(trial2, na.rm =TRUE)
      ratio = max-min
      trial2 <- sweep(trial2,2,min,"-")
      trial2 <- sweep(trial2,2,ratio,"/")
      trial3<-cbind(trial[,c(1,2)], as.data.frame(trial2))
      df<-trial3 %>%
        gather(species, yield, -c(year, color)) 
    }
    return(df = df)
  }
  
  if (rescale == 1 | rescale == 2){
    yield_sim<-fn_rescale(yield_sim)
  }
  
  yield_sim<-as.data.frame(yield_sim)
  yield_sim$year<-rownames(yield_sim)
  plot_yield_sim<-yield_sim %>% 
    gather(species, yield, -year) %>% 
    mutate(year = as.numeric(year)) %>% 
    mutate(color = "modelled")
  
  if(fleetDynamics == TRUE){
    
    # sim_FD  
    if (type == "yield"){
      yield_sim_FD_spp<-sim_FD@yield
      yield_sim_FD_spp<-rowSums(aperm(yield_sim_FD_spp,c(1,2,4,3)),dims=3) # sum over size
      yield_sim_FD_spp<-rowSums(yield_sim_FD_spp,dims=2) 
      yield_sim_FD_spp<-yield_sim_FD_spp[as.character(2006:2017),] 
    }else{
      yield_sim_FD_spp<-getSSB(sim_FD)
      yield_sim_FD_spp<-yield_sim_FD_spp[as.character(2006:2017),] 
    }
    
    # rescale
    if (rescale == 1 | rescale == 2){
      yield_sim_FD_spp<-fn_rescale(yield_sim_FD_spp)
    }
    
    yield_sim_FD_spp<-as.data.frame(yield_sim_FD_spp)
    yield_sim_FD_spp$year<-rownames(yield_sim_FD_spp)
    plot_yield_sim_FD<-yield_sim_FD_spp %>% 
      gather(species, yield, -year) %>% 
      mutate(year = as.numeric(year)) %>% 
      mutate(color = "modelled_FD")
    
  } 
  
  # observed
  
  # yieldObs_timeVariant in g m-3, see CreateFmort
  # or SSB here trnsformed in g m3
  if(type == "yield"){
    observed<-yieldObs_timeVariant
  }else {
    observed<-ssbObs %>% 
      filter(species %in% df_param$species) %>% 
      mutate(SSB = (SSB*1000000)/areaEco)
    observed<-acast(observed, year ~ species) 
    observed<-observed[as.character(1995:2017),]  
  }

  # rescale
  if (rescale == 1 | rescale == 2){
    observed<-fn_rescale(observed)
  }
  
  observed[observed==0]<-NA # not sure about this ...
  observed<-as.data.frame(observed)
  observed$year<-rownames(observed)
  
  plot_observed<-observed %>% 
    gather(species, yieldObs, -year) %>% 
    mutate(year = as.numeric(year)) %>% 
    mutate(color = "observed") %>% 
    `colnames<-`(c("year", "species", "yield", "color")) 
  
  if(fleetDynamics==FALSE){
    plot_y<-rbind(plot_yield_sim, plot_observed)
  }else{
    plot_y<-rbind(plot_yield_sim, plot_yield_sim_FD, plot_observed) 
  }
  
  # rescale option 3
  if(rescale == 3){
    plot_y<-fn_rescale(plot_y)
  }
  
  # fit lm and calcualte R^2
  if (fleetDynamics==FALSE){
    fit<-plot_yield_sim %>% 
      select(-color) %>% 
      `colnames<-`(c("year","species", "modelled"))  
  }else{
    fit<-plot_yield_sim_FD %>% 
      select(-color) %>% 
      `colnames<-`(c("year","species", "modelled"))
  }
  
  fit2<-plot_observed %>% 
    select(-color) %>% 
    `colnames<-`(c("year","species", "observed"))
  r<-merge(fit, fit2, all =TRUE)
  
  library(broom)
  r <- r %>% 
    filter(!is.na(observed),!is.na(modelled)) %>% 
    group_by(species) %>%
    do(model = lm(modelled ~ observed, data = .)) 
  Coef = glance(r, model)
  
  plot_y<-plot_y %>% 
    right_join(df_param[,c("species","spCommon")]) %>% 
    filter(!species %in% c("myctophids","helicolenus barathri")) %>%
    mutate(yield = ifelse(species %in% c("nototodarus gouldi","trachurus declivis" ) & yield==0, NA, yield))
  
  # if(type == "SSB"){
  #   plot_y<-plot_y %>% 
  #   filter(species %in% c("sillago flindersi","centroberyx affinis","nemadactylus macropterus","platycephalus richardsoni", "centroberyx gerrardi","seriolella brama", "hoplostethus atlanticus","macruronus novaezelandiae", "seriolella punctata","rexea solandri", "genypterus blacodes")) # don't have observed
  # }
  
  # arrange spp in order of size
  plot_y$spCommon<-factor(plot_y$spCommon, levels = c(df_param$spCommon))
  plot_y$color<- as.factor(plot_y$color)
  
  if(fleetDynamics==TRUE){
    plot_y$color<- ordered(plot_y$color, levels = c("observed","modelled","modelled_FD") )
  }else{
    plot_y$color<- ordered(plot_y$color, levels = c("observed","modelled") )
  }
  
  # plot real numbers
  if(rescale == 0){
    plot_y$yield<-plot_y$yield/1000000 # to tonnes 
    plot_y$yield<-plot_y$yield*areaEco # in ecosystem modelled 
  }
  
  # plot trends in modelled and observed catches 
  plot_trends<-function(df){
    p <- ggplot(df) + 
      geom_line(aes(x=year, y = yield, colour=color), size=1) +
      scale_x_continuous(name = "Year") +
      theme_bw() +
      theme(legend.title = element_blank(),
            text = element_text(size=12),
            legend.text=element_text(size=12),
            axis.title.y = element_text(size=15, vjust=0.4),
            axis.title.x = element_text(size=15, vjust=0.3),
            axis.text.y = element_text(size=10, hjust=0.5),
            axis.text.x = element_text(size=10, angle = 90, hjust=0.5),
            panel.grid.major = element_blank(),
            strip.background = element_rect(colour="black", fill= NA)) 
    
    if(type == "yield"){
      p<-p+scale_y_continuous(name = "Yield [tonnes]") 
    }else {
      p<-p+scale_y_continuous(name = "SSB [tonnes]")
    }
  }
  
  # call the plotting function
  # plot_y$year<- as.factor(plot_y$year)
  plotYield<-plot_trends(plot_y)
  plotYield<-plotYield+facet_wrap(~spCommon)
  
  ######## comparison option 2
  
  if(fleetDynamics == TRUE){ 
    
    yield_sim_FD_fl<-sim_FD@yield
    yield_sim_FD_fl<-rowSums(aperm(yield_sim_FD_fl,c(1,2,4,3)),dims=3) # sum over size
    trial = yield_sim_FD_fl
    yield_sim_FD_fl<-rowSums(aperm(yield_sim_FD_fl,c(1,3,2)),dims=2) # sum over speceis
    yield_sim_FD_fl<-yield_sim_FD_fl[as.character(2006:2017),] 
    
    # rescale
    if (rescale == 1 | rescale == 2){
      yield_sim_FD_fl<-fn_rescale(yield_sim_FD_fl)
    }
    
    yield_sim_FD_fl<-as.data.frame(yield_sim_FD_fl)
    yield_sim_FD_fl$year<-rownames(yield_sim_FD_fl)
    yield_sim_FD_fl<-yield_sim_FD_fl %>% 
      gather(fleet, yield, -year) %>% 
      mutate(year = as.numeric(year)) %>% 
      mutate(color = "modelled_FD")
    
    # compare observed catches at fleet level
    c_temp<-datValidationYield %>% 
      filter(metier %in% df_main_metier) %>% 
      mutate(metier = case_when(metier == "GHAT - Southern Shark Gillnet" ~ "SSG",
                                metier == "South East Trawl Fishery - Danish Seine" ~ "SED",
                                metier == "South East Trawl Fishery - Otter trawl deepSlope" ~ "SET-DS",
                                metier == "South East Trawl Fishery - Otter trawl shelf" ~ "SET-SH",
                                metier == "South East Trawl Fishery - Otter trawl upperSlope" ~ "SET-US")) %>% 
      mutate(catch = (catchComm_t*1000000)/areaEco) %>% # from t vol eco to g m3 
      select(-catchComm_t)
    
    colnames(c_temp)<-c("fleet","year", "yield")
    c_temp<-c_temp[,c(2,1,3)]

    c_temp<-c_temp %>% 
      spread(fleet, yield) %>% 
      column_to_rownames("year")
    c_temp<-as.matrix(c_temp)
    
    # rescale
    if (rescale == 1 | rescale == 2){
      c_temp<-fn_rescale(c_temp)
    }
    
    c_temp<-as.data.frame(c_temp)
    c_temp$year<-rownames(c_temp)
    
    plot_observed<-c_temp %>% 
      gather(species, yield, -year) %>% 
      mutate(year = as.numeric(year)) %>% 
      mutate(color = "observed") %>% 
      `colnames<-`(c("year", "fleet", "yield", "color")) 
    
    plot_y_FD_fl<-rbind(yield_sim_FD_fl, plot_observed)
    plot_y_FD_fl$color<-as.factor(plot_y_FD_fl$color)
    plot_y_FD_fl$color<-ordered(plot_y_FD_fl$color, levels = c("observed","modelled_FD"))
    
    # rescale option 3
    if(rescale == 3){
      plot_y_FD_fl <- plot_y_FD_fl %>% dplyr::rename(species = fleet) 
      plot_y_FD_fl<-fn_rescale(plot_y_FD_fl)
      plot_y_FD_fl <- plot_y_FD_fl %>% dplyr::rename(fleet = species) 
    }
    
    # plot real numbers
    if(rescale == 0){
      plot_y_FD_fl$yield<-plot_y_FD_fl$yield/1000000 # to tonnes 
      plot_y_FD_fl$yield<-plot_y_FD_fl$yield*areaEco # in ecosystem modelled 
      trial[,,]<-trial[,,]/1000000
      trial[,,]<-trial[,,]*areaEco
    }
    
    # plot trends in modelled and observed catches 
    # plot_y_FD_fl$year<- as.factor(plot_y_FD_fl$year)
    plotYield_FD_fl<-plot_trends(plot_y_FD_fl)
    plotYield_FD_fl<-plotYield_FD_fl+facet_wrap(~fleet, nrow = 1)
  
  ######## comparison option 3
  
    effort_sim_FD_fl<-sim_FD@effortOut
    effort_sim_FD_fl<-effort_sim_FD_fl[as.character(2006:2016),] 
    
    # rescale
    if (rescale == 1 | rescale == 2){
      effort_sim_FD_fl<-fn_rescale(effort_sim_FD_fl)
    }
    
    effort_sim_FD_fl<-as.data.frame(effort_sim_FD_fl)
    effort_sim_FD_fl$year<-rownames(effort_sim_FD_fl)
    effort_sim_FD_fl<-effort_sim_FD_fl %>% 
      gather(fleet, effort, -year) %>% 
      mutate(year = as.numeric(year)) %>% 
      mutate(color = "modelled_FD")
    
    # compare observed vs modelled effort
    # units: observed effort is in opn y-1 for the whole system. here transformed in m-3 y-1. then scaling factor applied as modelled effort high compared to observed effort. Though high modelled effort results in the observed catches.  
    e_temp<-datValidationEffort %>% 
      filter(metier %in% df_main_metier) %>% 
      mutate(metier = case_when(metier == "GHAT - Southern Shark Gillnet" ~ "SSG",
                                metier == "South East Trawl Fishery - Danish Seine" ~ "SED",
                                metier == "South East Trawl Fishery - Otter trawl deepSlope" ~ "SET-DS",
                                metier == "South East Trawl Fishery - Otter trawl shelf" ~ "SET-SH",
                                metier == "South East Trawl Fishery - Otter trawl upperSlope" ~ "SET-US")) %>% 
      mutate(opn = opn*scaling_cost_area) %>% # from opn to area trawled in the ecosystem
      mutate(opn = opn/areaEco) %>% # divided by the m3 of the system -> becomes area trawled in 1m3  
      mutate(opn = ifelse(metier == "SED", opn*scaling_effort[1], ifelse(metier == "SET-DS", opn*scaling_effort[2], ifelse(metier == "SET-SH",opn*scaling_effort[3], ifelse(metier == "SET-US", opn*scaling_effort[4], opn*scaling_effort[5])))))
    
    
    colnames(e_temp)<-c("fleet","year", "effort")
    e_temp<-e_temp[,c(2,1,3)]
    
    e_temp<-e_temp %>% 
      spread(fleet, effort) %>% 
      column_to_rownames("year")
    e_temp<-as.matrix(e_temp)
    
    # rescale
    if (rescale == 1 | rescale == 2){
      e_temp<-fn_rescale(e_temp)
    }
    
    e_temp<-as.data.frame(e_temp)
    e_temp$year<-rownames(e_temp)
    
    plot_observed<-e_temp %>% 
      gather(fleet, effort, -year) %>% 
      mutate(year = as.numeric(year)) %>% 
      mutate(color = "observed") %>% 
      `colnames<-`(c("year", "fleet", "effort", "color")) 

    plot_e_FD_fl<-rbind(effort_sim_FD_fl, plot_observed)
    plot_e_FD_fl$color<-as.factor(plot_e_FD_fl$color)
    plot_e_FD_fl$color<-ordered(plot_e_FD_fl$color, levels = c("observed","modelled_FD"))
    
    plot_e_FD_fl<-plot_e_FD_fl %>% 
      `colnames<-`(c("year","fleet","yield","color"))
    
    # rescale option 3
    if(rescale == 3){
      plot_e_FD_fl <- plot_e_FD_fl %>% dplyr::rename(species = fleet) 
      plot_e_FD_fl<-fn_rescale(plot_e_FD_fl)
      plot_e_FD_fl <- plot_e_FD_fl %>% dplyr::rename(fleet = species) 
    }
    
    # plot real numbers - opn in ecosystem
    # note that this is observed values in operations in the ecosystem as per original data, but could be also be transformed in m or km fished more in line with what has been used in the model
    if(rescale == 0){
      plot_e_FD_fl$yield<-plot_e_FD_fl$yield/scaling_cost_area # from area trawled to opn  
      plot_e_FD_fl$yield<-plot_e_FD_fl$yield*areaEco # in ecosystem modelled 
    }
    
    # call the plotting function 
    # plot_e_FD_fl$year<- as.factor(plot_e_FD_fl$year)
    plotEffort_FD_fl<-plot_trends(plot_e_FD_fl)
    plotEffort_FD_fl<-plotEffort_FD_fl+facet_wrap(~fleet, nrow = 1)
    plotEffort_FD_fl<-plotEffort_FD_fl+scale_y_continuous(name = "Effort [opn]")
  }
  
  # return time series of catcah and SSB too? orbetter is these are outputs from the model  
  if(fleetDynamics == TRUE){
    return(list(plotYield_sp = plotYield, DataPlotYield = plot_y, plotYield_fl = plotYield_FD_fl, plotEffort_fl = plotEffort_FD_fl, DataPlotYield_fl = plot_y_FD_fl, DataPlotEffort_fl = plot_e_FD_fl, trial = trial))
  }else{
    return(list(plotYield_sp = plotYield, DataPlotYield = plot_y))
  }
  
}

########### plot biomass ---- 

plotBiomass_CN <- function(sim, sim_FD ,df_param, rescale){
 
  # # this is the same as plotBiomass() but it's your version and it's similar to the ones aboves
  # sim = sim_fitted
  # sim_FD = sim_FD_bmsy
  # df_param
  # rescale = 0
  
  b <- getBiomass(sim)
  b<-b[as.character(1995:2016),] 
  bm <- reshape2::melt(b)
  bm$color <- "modelled"
  b2 <- getBiomass(sim_FD)
  bm2 <- reshape2::melt(b2)
  bm2$color <- "modelled_FD"
  bm <- rbind(bm, bm2)
  colnames(bm)<-c("time", "Species","value","color")
  bm$Species <- as.factor(bm$Species)
  min_value <- 1e-20
  spec_bm <- bm[bm$value >= min_value, ]

  # common names and trsnform in tonnes per ecosystem 

  spec_bm$Species<-as.character(spec_bm$Species)
  df_param2<-df_param
  
  colnames(df_param2)<-c("Species","spCommon")
  spec_bm<-spec_bm %>% 
    left_join(df_param2[,c(1,2)])
  
  # arrange spp in order of size
  spec_bm$spCommon<-factor(spec_bm$spCommon, levels = c(df_param2$spCommon))
  spec_bm$color<- as.factor(spec_bm$color)
  spec_bm$color<- ordered(spec_bm$color, levels = c("modelled","modelled_FD") )
  
  if(rescale == 0 ){
    spec_bm$value<-spec_bm$value/1000000 # to tonnes 
    spec_bm$value<-spec_bm$value*areaEco # in ecosystem modelled 
  }
  
  p <- ggplot(spec_bm) + 
    geom_line(aes(x=time, y = value, group = color, color = color), size=1) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Biomass [tonnes]") +
    theme_bw() +
    theme(legend.title = element_blank(),
          text = element_text(size=12),
          legend.text=element_text(size=12),
          axis.title.y = element_text(size=15, vjust=0.4),
          axis.title.x = element_text(size=15, vjust=0.3),
          axis.text.y = element_text(size=10, hjust=0.5),
          axis.text.x = element_text(size=10, angle = 90, hjust=0.5),
          panel.grid.major = element_blank(),
          strip.background = element_rect(colour="black", fill= NA)) + 
    facet_wrap(~spCommon)
  # facet_wrap(~spCommon, scale = "free_y")
  
  return(list(p = p, BioData = spec_bm))
}














########### get msy ----

getBlevel<-function(sim_calibrated_unfished,matrix_effort_nofishing,constant_effort,sim_calibrated,df_param,sim_fitted){
  
  # *Bref* needs an unfished community at equilibrium (sim_calibrated_unfished)
  # *Bmsy* needs a fished community at equilibrium (sim_calibrated)
  # *Bhist* needs a fished commnity time variant (sim_fitted)
  
  ### calculate *Bref* at unfished community
  # calcaulte abundance at last time step - these are Bref 
  bioUnfished<-getBiomass(sim_calibrated_unfished)[dim(matrix_effort_nofishing)[1],]
  bioUnfished<-as.data.frame(bioUnfished)
  bioUnfished$species<-rownames(bioUnfished)
  
  # calcaulte abundance fully fished 
  BrefFished <- getBiomass(sim_calibrated_unfished)[dim(matrix_effort_nofishing)[1]/2,]
  BrefFished <- as.data.frame(BrefFished)
  BrefFished$species<-rownames(BrefFished) 
  
  # calculate *Bmsy*
  effort_y<-constant_effort[1:20,]
  out_y<-data.frame()
  out_e<-data.frame()
  out_b<-data.frame()
  
  # start from community at equilibrium (sim_calibrated) and fish each spp for 10 time steps effort (as in sim_calibrated). BUT change effort for each spp one at the time
  for(i in 1:ncol(effort_y)){
    
    prova <- effort_y
    # prova[,i]<-seq(from = 0.01, to = 0.5, length.out = dim(effort_y)[1])
    prova[,i]<-seq(from = 0.001, to = 1, length.out = dim(effort_y)[1])
    
    # NOTE here whether this is sim_calibrated@params or params changes a lot.... 
    sim_check <- project(sim_calibrated@params, effort = prova, dt = dt, fleetDynamics = FALSE, multiFleet = FALSE, management = FALSE, price = NA, cost = NA, diet_steps = 0, initial_n = sim_calibrated@n[dim(sim_calibrated@n)[1],,], initial_n_pp = sim_calibrated@n_pp[dim(sim_calibrated@n_pp)[1],], initial_n_bb = sim_calibrated@n_bb[dim(sim_calibrated@n_bb)[1],]) 
    
    check_y<-getYield(sim_check)
    check_b<-getBiomass(sim_check)
    
    out_y<-rbind(out_y, check_y[,i])
    out_b<-rbind(out_b, check_b[,i])
    out_e<-rbind(out_e, prova[,i])
  }
  
  # Yield by species and time steps
  out_y<-as.data.frame(t(out_y))
  colnames(out_y)<-df_param$species
  out_y <- out_y %>% 
    gather("species", "yield") %>% 
    mutate(time = seq(1:dim(effort_y)[1]))
  
  # Biomass by species and time steps
  out_b<-as.data.frame(t(out_b))
  colnames(out_b)<-df_param$species
  out_b <- out_b %>% 
    gather("species", "biomass") %>% 
    mutate(time = seq(1:dim(effort_y)[1]))
  
  # Effort by species and time steps
  out_e<-as.data.frame(t(out_e))
  colnames(out_e)<-df_param$species
  
  # combined dataset to calculate Bmsy and for plotting 
  check_yield <- out_e %>% 
    gather("species", "effort")%>% 
    mutate(time = seq(1:dim(effort_y)[1])) %>% 
    left_join(out_y) %>% 
    left_join(out_b) %>% 
    filter(species != "myctophids") %>% 
    mutate(species = factor(species, levels=c(df_param$species)))
  
  # caclaulte Bmsy 
  Bmsy<-split(check_yield, check_yield$species)
  msy<-do.call("rbind",lapply(Bmsy, function(x) x[which.max(x$yield),]))
  msy<-msy %>% 
    select(-time) %>% 
    `colnames<-`(c("species","Fmort","yiled","Bmsy")) 
  rownames(msy) <- NULL
  
  # for plotting 
  check_yield_plot<- check_yield %>% 
    gather(variable, value, -c(species,effort,time))
  
  # plot 
  p<-ggplot()+ 
    geom_smooth(data = filter(check_yield_plot, variable == "yield"), aes(y = value, x = effort, group = variable, color = variable))+
    geom_vline(data = msy, mapping = aes(xintercept = Fmort))+
    scale_y_continuous(name = "Yield") +
    scale_x_continuous(name = "Effort") +
    scale_color_manual(values = c("blue"))+
    theme_bw() +
    theme(text = element_text(size=12),
          legend.text=element_text(size=12),
          legend.position = "none", 
          axis.title.y = element_text(size=15, vjust=0.4),
          axis.title.x = element_text(size=15, vjust=0.3),
          axis.text.y = element_text(size=10, hjust=0.5),
          axis.text.x = element_text(size=10, angle = 90, hjust=0.5),
          panel.grid.major = element_blank(),
          strip.background = element_rect(colour="black", fill= NA))+
    facet_wrap(~species, scale="free")
  
  # caclualte *Bhist* 
  hist_y<-getYield(sim_fitted)
  hist_y<-as.data.frame(hist_y)
  hist_y$time<-rownames(hist_y)
  hist_y<-hist_y %>% 
    gather(species, yield, "myctophids":"galeorhinus galeus")
  
  hist_b<-getBiomass(sim_fitted)
  hist_b<-as.data.frame(hist_b)
  hist_b$time<-rownames(hist_b)
  hist_b<-hist_b %>% 
    gather(species, biomass, "myctophids":"galeorhinus galeus") %>% 
    left_join(hist_y)
  
  Bhist<-split(hist_b, hist_b$species)
  Bhist<-do.call("rbind",lapply(Bhist, function(x) x[which.max(x$yield),]))
  Bhist<-Bhist %>% 
    select(-time) %>% 
    `colnames<-`(c("species","Bhist","yiled"))
  rownames(Bhist) <- NULL
  
  # results 
  return(list(bioUnfished = bioUnfished, msy=msy, p = p, Bhist = Bhist))
}
