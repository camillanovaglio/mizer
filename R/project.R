# Project method for the size based modelling package mizer

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

# project can dispatch with effort being different classes (missing, numeric,
# array).

#' @useDynLib mizer
#' @importFrom Rcpp sourceCpp
NULL


#' project method for the size based modelling
#' 
#' Runs the size-based model simulation and projects the size based model
#' through time. \code{project} is called using an object of type
#' \linkS4class{MizerParams} and an object that contains the effort of the
#' fishing gears through time. The method returns an object of type
#' \linkS4class{MizerSim} which can then be explored with a range of summary and
#' plotting methods.
#' 
#' @param params A \linkS4class{MizerParams} object
#' @param effort The effort of each fishing gear through time. See notes below.
#' @param t_max The maximum time the projection runs for. The default value is
#'   100. However, this argument is not needed if an array is used for the
#'   \code{effort} argument, in which case this argument is ignored. See notes
#'   below.
#' @param dt Time step of the solver. The default value is 0.1.
#' @param t_save The frequency with which the output is stored. The default
#'   value is 1. Must be an integer multiple of dt.
#' @param initial_n The initial populations of the species. By default the 
#'   \code{initial_n} slot of the \linkS4class{MizerParams} argument is used.
#'   See the notes below.
#' @param initial_n_pp The initial population of the plankton spectrum. It
#'   should be a numeric vector of the same length as the \code{w_full} slot of
#'   the \code{MizerParams} argument. By default the \code{initial_n_pp} slot of the
#'   \linkS4class{MizerParams} argument is used.
#'   ##AAsp##
#' @param initial_n_bb The initial population of the benthic spectrum. At the moment it is
#'  just follows the same code as for the plankton spectrum 
#' @param initial_n_aa The initial population of the algal spectrum. At the moment it is
#'  just follows the same code as for the plankton spectrum 
#'  @param shiny_progress A shiny progress object used to update shiny progress bar.
#'   Default NULL.
#' @param ... Currently unused.
#' 
#' @note The \code{effort} argument specifies the level of fishing effort during
#' the simulation. It can be specified in three different ways: \itemize{ \item
#' A single numeric value. This specifies the effort of all fishing gears which
#' is constant through time (i.e. all the gears have the same constant effort). 
#' \item A numerical vector which has the same length as the number of fishing
#' gears. The vector must be named and the names must correspond to the gear
#' names in the \code{MizerParams} object. The values in the vector specify the
#' constant fishing effort of each of the fishing gears, i.e. the effort is
#' constant through time but each gear may have a different fishing effort. 
#' \item A numerical array with dimensions time step x gear. This specifies the
#' fishing effort of each gear at each time step.  The first dimension, time,
#' must be named numerically and contiguously. The second dimension of the array
#' must be named and the names must correspond to the gear names in the
#' \code{MizerParams} argument. The value for the effort for a particular time
#' is used during the interval from that time to the next time in the array.}
#' 
#' If effort is specified as an array then the smallest time in the array is 
#' used as the initial time for the simulation. Otherwise the initial time is
#' set to 0. Also, if the effort is an array then the \code{t_max} argument is 
#' ignored and the maximum simulation time is the largest time of the effort
#' array.
#' 
#' The \code{initial_n} argument is a matrix with dimensions species x size. 
#' It specifies the abundances of the species at the initial time. The
#' order of species must be the same as in the \code{MizerParams} argument. If
#' the initial population is not specified, the argument is set by default by
#' the \code{get_initial_n} function which is set up for a North Sea model.
#' @return An object of type \code{MizerSim}.
#' @export
#' @seealso \code{\link{MizerParams}}
#' @examples
#' \dontrun{
#' # Data set with different fishing gears
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # With constant fishing effort which is different for each gear
#' effort <- c(Industrial = 0, Pelagic = 1, Beam = 0.5, Otter = 0.5)
#' sim <- project(params, t_max = 20, effort = effort)
#' # With fishing effort that varies through time for each gear
#' gear_names <- c("Industrial","Pelagic","Beam","Otter")
#' times <- seq(from = 1, to = 10, by = 1)
#' effort_array <- array(NA, dim = c(length(times), length(gear_names)),
#'     dimnames = list(time = times, gear = gear_names))
#' effort_array[,"Industrial"] <- 0.5
#' effort_array[,"Pelagic"] <- seq(from = 1, to = 2, length = length(times))
#' effort_array[,"Beam"] <- seq(from = 1, to = 0, length = length(times))
#' effort_array[,"Otter"] <- seq(from = 1, to = 0.5, length = length(times))
#' sim <- project(params, effort = effort_array)
#' }
#' 
project <- function(params, effort = 0,  t_max = 100, dt = 0.1, t_save=1,
                    temperature = rep(params@t_ref, times = t_max), # what do we do with the t_ref?
                    initial_n = params@initial_n,

                    # AA
                    initial_n_pp = params@initial_n_pp, 
                    initial_n_bb = params@initial_n_bb,
                    initial_n_aa = params@initial_n_aa,
                    # CN
                    fleetDynamics, management, price = price, 
                    cost = cost, scaling_price = scaling_price, ke = ke,
                    
                    shiny_progress = NULL, 
                    # AA
                    diet_steps=10, ...) {  #default number of years (steps?) to calcualte diet for 


  # # # trial
  # params = params_FD
  # effort = 0
  # # effort = relative_effort
  # t_max = dim(df_price_mean)[1]
  # dt = 0.5
  # t_save=1
  # temperature = rep(params@t_ref, times = t_max)
  # initial_n = initial_n
  # initial_n_pp = initial_n_pp
  # initial_n_bb = initial_n_bb
  # initial_n_aa = params@initial_n_aa
  # fleetDynamics = TRUE
  # management = TRUE
  # price = df_price_mean2
  # scaling_price = scaling_price
  # cost = df_cost_mean2
  # shiny_progress = NULL
  # diet_steps=0
  # 
    validObject(params)

    # Do we need to create an effort array?
    
    if (is.vector(effort)) {
        no_gears <- dim(params@catchability)[1]
        if ((length(effort) > 1) & (length(effort) != no_gears)) {
            stop("Effort vector must be the same length as the number of fishing gears\n")
        }
        # If more than 1 gear need to check that gear names match
        gear_names <- dimnames(params@catchability)[[1]]
        effort_gear_names <- names(effort)
        if (length(effort) == 1 & is.null(effort_gear_names)) {
            effort_gear_names <- gear_names
        }
        if (!all(gear_names %in% effort_gear_names)) {
            gear_names_error_message <- paste("Gear names in the MizerParams object (", paste(gear_names, collapse=", "), ") do not match those in the effort vector.", sep="")
            stop(gear_names_error_message)
        }
        # Set up the effort array transposed so we can use the recycling rules
        time_dimnames <- signif(seq(from = 0, to = t_max, by = dt), 3)
        effort <- t(array(effort, dim = c(no_gears, length(time_dimnames)), 
                          dimnames = list(gear = effort_gear_names, time = time_dimnames)))
    }
    
    # Check that number and names of gears in effort array is same as in MizerParams object
    no_gears <- dim(params@catchability)[1]
    if(dim(effort)[2] != no_gears){
        no_gears_error_message <- paste("The number of gears in the effort array (length of the second dimension = ", dim(effort)[2], ") does not equal the number of gears in the MizerParams object (", no_gears, ").", sep="")
        stop(no_gears_error_message)
    }
    gear_names <- dimnames(params@catchability)[[1]]
    if(!all(gear_names %in% dimnames(effort)[[2]])){
        gear_names_error_message <- paste("Gear names in the MizerParams object (", paste(gear_names, collapse=", "), ") do not match those in the effort array.", sep="")
        stop(gear_names_error_message)
    }
    # Sort effort array to match order in MizerParams
    effort <- effort[,gear_names, drop=FALSE]
    
    # Blow up time dimension of effort array
    # i.e. effort might have been passed in using time steps of 1, but actual dt = 0.1, so need to blow up
    if (is.null(dimnames(effort)[[1]])){
        stop("The time dimname of the effort argument must be numeric.")
    }
    if (any(is.na(as.numeric(dimnames(effort)[[1]])))){
        stop("The time dimname of the effort argument must be numeric.")
    }
    time_effort <- as.numeric(dimnames(effort)[[1]])
    if (is.unsorted(time_effort)) {
        stop("The time dimname of the effort argument should be increasing.")
    }
    
    t_max <- time_effort[length(time_effort)]
    # Blow up effort so that rows are dt spaced
    time_effort_dt <- seq(from = time_effort[1], to = t_max, by = dt)
    effort_dt <- t(array(NA, dim = c(length(time_effort_dt), dim(effort)[2]), dimnames=list(time = time_effort_dt, dimnames(effort)[[2]])))
    for (i in 1:(length(time_effort)-1)){
        effort_dt[,time_effort_dt >= time_effort[i]] <- effort[i,]
    }
    effort_dt <- t(effort_dt)
    
    dim(effort_dt)
    
    ## AA
    ## expand temperature vector in the required number of timesteps 
    ###TODO?####
    ## at the moment we just repeat yearly values for the entire year, no smothing or interpolation is used
    if (length(temperature) != t_max) {
      stop("your temperature input vector is not the same length as t_max")
    }
    
    time_temperature_dt <- rep(temperature, length = t_max/dt, each = 1/dt) # works if t_max = length(temperature)
    x_axis <- seq(length.out=(t_max/dt),from =1)   # = time vector
    # need smoothing?
    # myData <- data.frame("y" = time_temperature_dt, "x" = x_axis) # create dataframe for smoothing (not sure if needed)
    # temperature_dt <- matrix(predict(loess(y~x, myData, span = 0.1)), dimnames = list(x_axis, "temperature")) # temperature vector following dt
    
    temperature_dt <- matrix(time_temperature_dt, dimnames = list(x_axis, "temperature")) # without smoothing

    dim(time_temperature_dt)
    
    #arrays with scalar values for all time, species and size
    metTempScalar <- array(NA, dim = c(dim(params@species_params)[1], length(params@w), length(temperature_dt)), dimnames = list(params@species_params$species,params@w,temperature_dt)) 
    matTempScalar <- array(NA, dim = c(dim(params@species_params)[1], length(params@w), length(temperature_dt)), dimnames = list(params@species_params$species,params@w,temperature_dt)) 
    morTempScalar <- array(NA, dim = c(dim(params@species_params)[1], length(params@w), length(temperature_dt)), dimnames = list(params@species_params$species,params@w,temperature_dt)) 
    intTempScalar <- array(NA, dim = c(dim(params@species_params)[1], length(params@w), length(temperature_dt)), dimnames = list(params@species_params$species,params@w,temperature_dt)) 

    for(iSpecies in 1:dim(params@species_params)[1])
    {
      metTempScalar[iSpecies,,] <-  tempFun(temperature = temperature_dt[,1], t_ref = params@t_ref, 
                                            Ea = params@species_params$ea_met[iSpecies], 
                                            c_a = params@species_params$ca_met[iSpecies], w = params@w)
      
      matTempScalar[iSpecies,,] <-  tempFun(temperature = temperature_dt[,1], t_ref = params@t_ref, 
                                            Ea = params@species_params$ea_mat[iSpecies], 
                                            c_a = params@species_params$ca_mat[iSpecies],  w = params@w)
      
      morTempScalar[iSpecies,,] <-  tempFun(temperature = temperature_dt[,1], t_ref = params@t_ref, 
                                            Ea = params@species_params$ea_mor[iSpecies], 
                                            c_a = params@species_params$ca_mor[iSpecies],  w = params@w)
      
      intTempScalar[iSpecies,,] <-  tempFun(temperature = temperature_dt[,1], t_ref = params@t_ref, 
                                            Ea = params@species_params$ea_int[iSpecies], 
                                            c_a = params@species_params$ca_int[iSpecies],  w = params@w)
    }


    # CN create effortOut, yield, revenue and profit matrices and specify price and costs 
    # All new arguments have the blown up time dimension as effort does
    # CN moved from below:
    # Handy things
    no_sp <- nrow(params@species_params) # number of species
    no_w <- length(params@w) # number of fish size bins
    idx <- 2:no_w
    # CN added: 
    fleet = params@fleet
    species = params@species_params$species
    w = params@w
 
    if (fleetDynamics == TRUE){ # if it's FALSE you use effort matrix or t_max as specified above - so no need of else{}
      
      # old version when only mean prices were considered 
      # effortOut_dt<-array(0,dim = c(length(time_effort_dt),length(fleet)),dimnames=list(time = time_effort_dt, fleet = fleet)) # effort output of the model 
      # yield_dt <- array(0,dim = c(length(time_effort_dt),no_sp,no_w,length(fleet)), dimnames = list(time = time_effort_dt, species, w, fleet)) # B*Q*S*E
      # revenue_dt<-array(0,dim = c(length(time_effort_dt),length(fleet)),dimnames=list(time = time_effort_dt, fleet = fleet)) # yield * price
      # profit_dt<-array(0,dim=c(length(time_effort_dt),length(fleet)),dimnames=list(time= time_effort_dt, fleet=fleet)) # revenu - costXeffort
      # F_dt <- array(0,dim = c(length(time_effort_dt),no_sp,no_w,length(fleet)), dimnames = list(time = time_effort_dt, species, w, fleet)) # f_mort_gear matrix. time X sp X w X gear 
      # # dim(F_dt)
      # 
      # # here prices are constant and you may need to change this bit when using the time variant option. See SizeBasedModel_FLEET for rer. also price shold be use to determine the runs es effort is right now (?)
      # # order price according to speceis. Note that the input price need to have an entry for each spp. 
      # # price = df_price_mean
      # price<-price[order(match(price$species, params@species_params$species)), ]
      # price_dt<-do.call("rbind",replicate(length(time_effort_dt),price,simplify = FALSE))
      # price_dt$time <-sort(rep(time_effort_dt,nrow(price))) 
      # 
      # # cost = df_cost_opn_mean
      # cost<-cost %>% 
      #   group_by(subfleet) %>% # calcualte total costs - this chould be changed in a later version 
      #   dplyr::summarise(cost = sum(cost))
      # cost_dt<-do.call("rbind",replicate(length(time_effort_dt),cost,simplify = FALSE))
      # cost_dt$time <-sort(rep(time_effort_dt,nrow(cost))) 
      # cost_dt<-acast(cost_dt, formula = time ~ subfleet, value.var = "cost") # matrix of time X fleet
      # cost_dt<-cost_dt[,params@fleet]
      
      # trial - time variant price and cost
      
      ###### matrices should be given as inputs (as for effort) - this bit should be outside project
      
      # # if df_price_mean
      # t_max = 5
      # price = df_price_mean
      # price <- t(matrix(rep(df_price_mean$price,t_max), ncol = t_max, nrow = length(df_price_mean$species), dimnames = list(df_price_mean$species,1:t_max)))
      # price<-price[,params@species_params$species]
      # 
      # # if df_price             
      # price = spread(df_price, species, price)
      # rownames(price)<-price$year
      # price<-price[,-1]
      # class(price)
      # price<-as.matrix(price)
      # price<-price[, params@species_params$species] # order spp
      
      ###### 
      # NOW price dictate all dimentions
      # price_dt
      # effort_dt (needs to be adjusted even if not used )
      # cost_dt
      # temperature_dt
      # all other matrices as profits_dt etc... do a serch of _dt
      
      # Blow up time dimension of price array
      # i.e. effort might have been passed in using time steps of 1, but actual dt = 0.1, so need to blow up
      if (is.null(dimnames(price)[[1]])){
        stop("The time dimname of the price argument must be numeric.")
      }
      if (any(is.na(as.numeric(dimnames(price)[[1]])))){
        stop("The time dimname of the price argument must be numeric.")
      }
      time_price <- as.numeric(dimnames(price)[[1]])
      if (is.unsorted(time_price)) {
        stop("The time dimname of the price argument should be increasing.")
      }
      
      t_max <- time_price[length(time_price)]
      # Blow up price so that rows are dt spaced
      time_price_dt <- seq(from = time_price[1], to = t_max, by = dt)
      price_dt <- t(array(NA, dim = c(length(time_price_dt), dim(price)[2]), dimnames=list(time = time_price_dt, dimnames(price)[[2]])))
      for (i in 1:(length(time_price)-1)){ # I m not sure why this is -1
        price_dt[,time_price_dt >= time_price[i]] <- price[i,]
      }
      price_dt <- t(price_dt)
      
      dim(price_dt)
      
      # effort_dt[,1]
      # price_dt[,1]
      ### NOTE that price_dt was a df not a matrix so will need to change things from here on (line 653)
    
 
      ##### trial 
      # need to adjust effort_dt so that it can have the same time dimention as price. effort dt is a vector of 0 if fleetDynamics is on and is blown up to a matrix using the function argument t_max (see above), which her is overwrittent by the price dimention 
      effort_dt <- t(array(NA, dim = c(length(time_price_dt), dim(effort)[2]), dimnames=list(time = time_price_dt, dimnames(effort)[[2]])))
      # effort_dt <- t(array(NA, dim = c(length(time_price_dt), length(species)), dimnames=list(time = time_price_dt, species))) # no
      effort_dt <- t(effort_dt)
      dim(effort_dt)
      
      ### end trial 
      
      
      
      # # then cost but with same size than price and no setting size for other things 
      # cost = df_cost_opn
      # df_cost_opn_mean # this is different.... 
      # str(cost)
      # cost<-cost %>% 
      #   group_by(subfleet, year) %>% # calcualte total costs - this chould be changed in a later version 
      #   dplyr::summarise(cost = sum(param_value, na.rm=TRUE))
      # cost<-acast(cost, formula = year ~ subfleet, value.var = "cost") # matrix of time X fleet
      # cost<-cost[,params@fleet]
      # # cost need to have the same years as price!!! just a trial here 
      # addCost<-cost[8:9,]
      # rownames(addCost)<-c(2015, 2016)
      # cost<-rbind(cost, addCost)
      
      ###### 
      
      # blow up costs with same dimension as price
      cost_dt <- t(array(NA, dim = c(length(time_price_dt), dim(cost)[2]), dimnames=list(time = time_price_dt, dimnames(cost)[[2]])))
      for (i in 1:(length(time_price)-1)){ # I m not sure why this is -1
        cost_dt[,time_price_dt >= time_price[i]] <- cost[i,]
      }
      cost_dt <- t(cost_dt)
      
      dim(cost_dt)
      
      # then create all other matrices but with price determining their dimention
      
      effortOut_dt<-array(0,dim = c(length(time_price_dt),length(fleet)),dimnames=list(time = time_price_dt, fleet = fleet)) # effort output of the model 
      yield_dt <- array(0,dim = c(length(time_price_dt),no_sp,no_w,length(fleet)), dimnames = list(time = time_price_dt, species, w, fleet)) # B*Q*S*E
      revenue_dt<-array(0,dim = c(length(time_price_dt),length(fleet)),dimnames=list(time = time_price_dt, fleet = fleet)) # yield * price
      profit_dt<-array(0,dim=c(length(time_price_dt),length(fleet)),dimnames=list(time= time_price_dt, fleet=fleet)) # revenu - costXeffort
      F_dt <- array(0,dim = c(length(time_price_dt),no_sp,no_w,length(fleet)), dimnames = list(time = time_price_dt, species, w, fleet)) # f_mort_gear matrix. time X sp X w X gear 
      dim(effortOut_dt)
      
      # I am not sure how temperature will merge with fleetDynamics as you also have time_temperature_dt etc. - should the temperature dimention be included here and rescalde to the price dimention instead? yes you need to do so otherwise temp had a different dimention
      # copy asta's stuff befeore effortOut_dt for consistency
      
      ## AA
      ## expand temperature vector in the required number of timesteps 
      ###TODO?####
      ## at the moment we just repeat yearly values for the entire year, no smothing or interpolation is used
      
      ## CN added this - not sure though 
      temperature = rep(params@t_ref, times = t_max)
      
      if (length(temperature) != t_max) {
        stop("your temperature input vector is not the same length as t_max")
      }
      
      time_temperature_dt <- rep(temperature, length = t_max/dt, each = 1/dt) # works if t_max = length(temperature)
      # length(time_temperature_dt)
      x_axis <- seq(length.out=(t_max/dt),from =1)   # = time vector
      # length(x_axis)
      # need smoothing?
      # myData <- data.frame("y" = time_temperature_dt, "x" = x_axis) # create dataframe for smoothing (not sure if needed)
      # temperature_dt <- matrix(predict(loess(y~x, myData, span = 0.1)), dimnames = list(x_axis, "temperature")) # temperature vector following dt
      
      temperature_dt <- matrix(time_temperature_dt, dimnames = list(x_axis, "temperature")) # without smoothing
      
      # CN should these be given once and moved outside the fleetdynamis if and at the end?  
      
      #arrays with scalar values for all time, species and size
      metTempScalar <- array(NA, dim = c(dim(params@species_params)[1], length(params@w), length(temperature_dt)), dimnames = list(params@species_params$species,params@w,temperature_dt)) 
      matTempScalar <- array(NA, dim = c(dim(params@species_params)[1], length(params@w), length(temperature_dt)), dimnames = list(params@species_params$species,params@w,temperature_dt)) 
      morTempScalar <- array(NA, dim = c(dim(params@species_params)[1], length(params@w), length(temperature_dt)), dimnames = list(params@species_params$species,params@w,temperature_dt)) 
      intTempScalar <- array(NA, dim = c(dim(params@species_params)[1], length(params@w), length(temperature_dt)), dimnames = list(params@species_params$species,params@w,temperature_dt)) 
      
      for(iSpecies in 1:dim(params@species_params)[1])
      {
        metTempScalar[iSpecies,,] <-  tempFun(temperature = temperature_dt[,1], t_ref = params@t_ref, 
                                              Ea = params@species_params$ea_met[iSpecies], 
                                              c_a = params@species_params$ca_met[iSpecies], w = params@w)
        
        matTempScalar[iSpecies,,] <-  tempFun(temperature = temperature_dt[,1], t_ref = params@t_ref, 
                                              Ea = params@species_params$ea_mat[iSpecies], 
                                              c_a = params@species_params$ca_mat[iSpecies],  w = params@w)
        
        morTempScalar[iSpecies,,] <-  tempFun(temperature = temperature_dt[,1], t_ref = params@t_ref, 
                                              Ea = params@species_params$ea_mor[iSpecies], 
                                              c_a = params@species_params$ca_mor[iSpecies],  w = params@w)
        
        intTempScalar[iSpecies,,] <-  tempFun(temperature = temperature_dt[,1], t_ref = params@t_ref, 
                                              Ea = params@species_params$ea_int[iSpecies], 
                                              c_a = params@species_params$ca_int[iSpecies],  w = params@w)
      }
       
    } # end of fleetdynamics == TRUE
    
    # Make the MizerSim object with the right size
    # We only save every t_save steps
    # Divisibility test needs to be careful about machine rounding errors,
    # see https://github.com/sizespectrum/mizer/pull/2
    if((t_save < dt) || !isTRUE(all.equal((t_save - round(t_save / dt) * dt), 0)))
        stop("t_save must be a positive multiple of dt")
    t_skip <- round(t_save/dt)
    
    # CN if fleetdynamics is on, effort is not given and the price matrix gives the dimention instead (it overwrites t_max) 
    if(fleetDynamics==TRUE){
      t_dimnames_index <- seq(1, to = length(time_price_dt), by = t_skip)
      t_dimnames <- time_price_dt[t_dimnames_index]
    }else{
      t_dimnames_index <- seq(1, to = length(time_effort_dt), by = t_skip)
      t_dimnames <- time_effort_dt[t_dimnames_index]
    }
    
    # CN mizer uses effort and effort_dt. The fleetDynamics version needs these parameters as well (did not want to touch them) but also uses effortOut and effortOut_dt. the first is used to save results (tmax times) and the second to update and store effort values at each time step in the loop (tmax/dt times)  
    # I also created yield_dt, revenue_dt and profit_dt and yield revenue and profit matrices, using the same concept, but I am not using the dt version of these matrices for much (beside the below if statement)
    # t_dimnames now is either from effort or price depending on whethere you are giving the model an effort matrix or a price matrix when fleetdynamics is on. 
    sim <- MizerSim(params, t_dimnames = t_dimnames) 

    # CN otherwise it messes with effort/price dimensions 
    # if(fleetDynamics==FALSE){
      # Fill up the effort array
      sim@effort[] <- effort_dt[t_dimnames_index,] # save only some but not all dt steps 
    # }
    
    # AA
    #attach temperature dimension and all the rate scalars to the sim object, so it is stored and can be viewed later
    sim@temperature <- temperature_dt
    sim@metTempScalar <- metTempScalar
    sim@matTempScalar <- matTempScalar
    sim@morTempScalar <- morTempScalar
    sim@intTempScalar <- intTempScalar


    # CN effortOut etc are created in MizerSim, now you change the time indexing using the previously created effortOut_dt etc. - again I am not sure I actually need the dt version of yields etc for much  
    if(fleetDynamics==TRUE){
      sim@effortOut[] <- effortOut_dt[t_dimnames_index,]
      sim@yield[] <- yield_dt[t_dimnames_index,,,]
      sim@revenue[] <- revenue_dt[t_dimnames_index,]
      sim@profit[] <- profit_dt[t_dimnames_index,]
      sim@F[] <- F_dt[t_dimnames_index,,,]
    }
    
    # dim(sim@F)

    # Set initial population
    sim@n[1,,] <- initial_n 
    sim@n_pp[1,] <- initial_n_pp
    sim@n_bb[1,] <- initial_n_bb
    sim@n_aa[1,] <- initial_n_aa

    # Handy things

    # AA
    no_w_full<- length(sim@params@w_full) # full size spectrum (including background spectra)
    
    # CN this 3 lines have been moved above
    # no_sp <- nrow(sim@params@species_params) # number of species
    # no_w <- length(sim@params@w) # number of fish size bins
    # idx <- 2:no_w

    # Hacky shortcut to access the correct element of a 2D array using 1D notation
    # This references the egg size bracket for all species, so for example
    # n[w_minidx_array_ref] = n[,w_min_idx]
    w_min_idx_array_ref <- (sim@params@w_min_idx - 1) * no_sp + (1:no_sp)
    
    # sex ratio - DO SOMETHING LATER WITH THIS
    sex_ratio <- 0.5
    
    #create a matrix for diet comparison. For prey it has the number of columns set at no_sp+3 because we have 3 background spectra
    sim@diet_comp<-array(0, c(no_sp, no_w, no_sp + 3, no_w_full), 
                         dimnames=list( predator=as.character(params@species_params$species), pred_size = params@w, 
                                        prey = c(as.character(params@species_params$species), "plankton", "benthos", "algae"),
                                        prey_size = params@w_full))
    
    
    # Matrices for solver
    A <- matrix(0, nrow = no_sp, ncol = no_w)
    B <- matrix(0, nrow = no_sp, ncol = no_w)
    S <- matrix(0, nrow = no_sp, ncol = no_w)
    
    # initialise n and nPP
    # We want the first time step only but cannot use drop as there may only be a single species
    n <- array(sim@n[1, , ], dim = dim(sim@n)[2:3])
    dimnames(n) <- dimnames(sim@n)[2:3]
    n_pp <- sim@n_pp[1, ]

    n_bb <- sim@n_bb[1, ]
    n_aa <- sim@n_aa[1, ]

    # CN 
    if(fleetDynamics == TRUE){
      t_steps <- dim(price_dt)[1] - 1
    }else{
      t_steps <- dim(effort_dt)[1] - 1
    }
    

    # Set up progress bar
    # CN turn off when running calibration 
    # pb <- progress::progress_bar$new(
    #     format = "[:bar] :percent ETA: :eta",
    #     total = length(t_dimnames_index), width = 60)
    # if (hasArg(shiny_progress)) {
    #     # We have been passed a shiny progress object
    #     shiny_progress$set(message = "Running simulation", value = 0)
    #     proginc <- 1/length(t_dimnames_index)
    # }
    
    # If storing diet compisiton, make diet_comp_all array ahead of main loop 
    if (diet_steps>0){
      diet_comp_all<- array(0, dim(sim@diet_comp))
    }

    for (i_time in 1:t_steps) {
      # print(i_time)

    # CN we could initialise effort for the fleetDynamics too but the first time step is 0 anyway 
    
    
      # trial   
      # i_time=1
      
        # Do it piece by piece to save repeatedly calling methods
        # Calculate amount E_{a,i}(w) of available food
        avail_energy <- getAvailEnergy(sim@params, n = n, n_pp = n_pp, n_bb = n_bb, n_aa = n_aa)
        
        # Calculate amount f_i(w) of food consumed
        feeding_level <- getFeedingLevel(sim@params, n = n, n_pp = n_pp, n_bb = n_bb, n_aa = n_aa, 
                                         avail_energy = avail_energy)
        
        # Calculate the predation rate
        pred_rate <- getPredRate(sim@params, n = n, n_pp = n_pp, n_bb = n_bb, n_aa = n_aa, 
                                 intakeScalar = sim@intTempScalar[,,i_time], feeding_level = feeding_level)
        
        # Calculate predation mortality on fish \mu_{p,i}(w)
        # AA
        m2 <- getPredMort(sim@params, pred_rate = pred_rate, intakeScalar = sim@intTempScalar[,,i_time])
        
# Calculate mortality on the plankton spectrum
        m2_background <- getPlanktonMort(sim@params, n = n, n_pp = n_pp, n_bb = n_bb, n_aa = n_aa, 
                                         intakeScalar = sim@intTempScalar[,,i_time], pred_rate = pred_rate)

                #Calculate mortality of the benthis spectrum 
        m2_benthos <- getBenthosMort(sim@params, n = n, n_pp = n_pp, n_bb = n_bb, n_aa = n_aa, 
                                         pred_rate = pred_rate)
        m2_algae <- getAlgalMort(sim@params, n = n, n_pp = n_pp, n_bb = n_bb, n_aa = n_aa, 
                                     pred_rate = pred_rate)
        
        # Calculate the resources available for reproduction and growth

        e <- getEReproAndGrowth(sim@params, n = n, n_pp = n_pp, n_bb = n_bb, n_aa = n_aa, 
                                intakeScalar = sim@intTempScalar[,,i_time], metScalar = sim@metTempScalar[,,i_time], 

                                feeding_level = feeding_level)
        # CN the fleetDynamic needs more outputs of fishing mortalities use to calcualte parameters for the fleetDynamics eqn, hence functions used below are the _CN versions (sligly different). Also the argument of these functions is effortOut_dt, which is a model output (instead of input) and is updated every time step  
        
        if(fleetDynamics==TRUE){

          # error checking - ask Javi
          # if(i_time==102)browser()
          # calcualte biomass at this time step and important to calucalte yield in fleetDynamics eqn
          B_itime<-sweep(n, 2, sim@params@w * sim@params@dw, "*") # is n updated every time step?
          
          # if you are using eqn with no effort at the denominator you need a starting value
          # effort<-ifelse(i_time ==1, 0.1, effortOut_dt[i_time,]) 
          
          # Calculate total mortality \mu_i(w)
          # dim(e)
          # e[1,]<-NA # Error in mu_S[mu_S < 0] <- x[x < 0] : 
          # NAs are not allowed in subscripted assignments
          # how can e be NA?  getEReproAndGrowth
          
          z <- getMort_CN(sim@params, n = n, n_pp = n_pp, n_bb = n_bb, n_aa = n_aa, 
                     intakeScalar = sim@intTempScalar[,,i_time], metScalar = sim@metTempScalar[,,i_time], 
                     morScalar = sim@morTempScalar[,,i_time], effort = effortOut_dt[i_time,], e = e, m2 = m2)

          # new outputs needed to calcualte yield in fleetDynamics eqn
          # Fmortality = effort * catchability * selectivity 
          F_itime <- getFMort_CN(sim@params, effort = effortOut_dt[i_time,])$f_mort_gear # this will also be saved as F (fishing mortality)
          F_itime <-aperm(F_itime, c(2,3,1)) # as for yield and all the other matrices 
          # Fmortality = catchability * selectivity
          F_itime_FL <- getFMort_CN(sim@params, effort = effortOut_dt[i_time,])$f_mort_gear_FL
          F_itime_FL <-aperm(F_itime_FL, c(2,3,1))

        }else{
          
         # AA
         # Moved total mortality calculation after the e calculation betcause we need e for stravation
         # Calculate total mortality \mu_i(w)
         z <- getMort(sim@params, n = n, n_pp = n_pp, n_bb = n_bb, n_aa = n_aa, 
                     intakeScalar = sim@intTempScalar[,,i_time], metScalar = sim@metTempScalar[,,i_time], 
                     morScalar = sim@morTempScalar[,,i_time], effort = effort_dt[i_time,], e = e, m2 = m2)
        }
        
        # Calculate the resources for reproduction
        e_repro <- getERepro(sim@params, n = n, n_pp = n_pp, n_bb = n_bb, n_aa = n_aa, e = e, 
                             intakeScalar = sim@intTempScalar[,,i_time], metScalar = sim@metTempScalar[,,i_time])
        
        # Calculate the growth rate g_i(w)
        # AA
        e_growth <- getEGrowth(sim@params, n = n, n_pp = n_pp, n_bb = n_bb, n_aa = n_aa, 
                               intakeScalar = sim@intTempScalar[,,i_time], metScalar = sim@metTempScalar[,,i_time], e_repro = e_repro, e = e) 
  
        # R_{p,i}
        # AA
        rdi <- getRDI(sim@params, n = n, n_pp = n_pp, n_bb = n_bb, n_aa = n_aa, 
                      intakeScalar = sim@intTempScalar[,,i_time], metScalar = sim@metTempScalar[,,i_time],
                      e_repro = e_repro, sex_ratio = sex_ratio)
        # R_i
        rdd <- getRDD(sim@params, n = n, n_pp = n_pp, n_bb = n_bb, n_aa = n_aa, rdi = rdi, 
                      intakeScalar = sim@intTempScalar[,,i_time],metScalar = sim@metTempScalar[,,i_time])

        # check if diet output is needed and if so call the diet comparison function 
        diet_store <- tail(t_dimnames_index, diet_steps) %in% (i_time + 1)
        
        if (any(diet_store)) {
          diet_comp_all[]<- getDietComp(sim@params, n=n,  n_pp=n_pp, n_bb = n_bb, n_aa = n_aa, 
                                        diet_comp_all=diet_comp_all, diet_steps=diet_steps, 
                                        intakeScalar = sim@intTempScalar[,,i_time])
          sim@diet_comp[]<-(diet_comp_all/diet_steps) + sim@diet_comp
        }

        # CN this is still a mistery to me... need to ask Asta

        # Iterate species one time step forward:
        # See Ken's PDF
        # A_{ij} = - g_i(w_{j-1}) / dw_j dt
        A[,idx] <- sweep(-e_growth[,idx-1,drop=FALSE]*dt, 2, sim@params@dw[idx], "/")
        # B_{ij} = 1 + g_i(w_j) / dw_j dt + \mu_i(w_j) dt
        B[,idx] <- 1 + sweep(e_growth[,idx,drop=FALSE]*dt,2,sim@params@dw[idx],"/") + z[,idx,drop=FALSE]*dt
        # S_{ij} <- N_i(w_j)
        S[,idx] <- n[,idx,drop=FALSE]
        # Boundary condition upstream end (recruitment)
        B[w_min_idx_array_ref] <- 1+e_growth[w_min_idx_array_ref]*dt/sim@params@dw[sim@params@w_min_idx]+z[w_min_idx_array_ref]*dt
        # Update first size group of n
        n[w_min_idx_array_ref] <- (n[w_min_idx_array_ref] + rdd*dt/sim@params@dw[sim@params@w_min_idx]) / B[w_min_idx_array_ref]
        # Update n
        # for (i in 1:no_sp) # number of species assumed small, so no need to vectorize this loop over species
        #     for (j in (sim@params@w_min_idx[i]+1):no_w)
        #         n[i,j] <- (S[i,j] - A[i,j]*n[i,j-1]) / B[i,j]

        # NOTE CN: Error in .Call("_mizer_inner_project_loop", PACKAGE = "mizer", no_sp,  :
        # "_mizer_inner_project_loop" not available for .Call() for package "mizer"
        n <- inner_project_loop(no_sp = no_sp, no_w = no_w, n = n,
                                A = A, B = B, S = S,
                                w_min_idx = sim@params@w_min_idx)

        # Dynamics of plankton spectrum uses a semi-chemostat model (de Roos - ask Ken)
        # We use the exact solution under the assumption of constant mortality during timestep
        ###TODO?#### 
        # now we need to scale r_pp (not rr_pp!) with the temperature response. 
        # this could be done at the start again where we create rr_pp object. this means it will have an extra time dimension 
        # or we calcualte rr_pp scalar here at every time step, which might be simpler 
        tmp <- (sim@params@rr_pp * sim@params@cc_pp / (sim@params@rr_pp + m2_background))
        n_pp <- tmp - (tmp - n_pp) * exp(-(sim@params@rr_pp + m2_background) * dt)

        # CN add effortDynamics and management
        # should this be getEffort in project_methods? 
        
        ##AAsp####
        # Dynamics of benthic spectrum uses a semi-chemostat model 
        # currently it follows exactly the same rules as plankton but has it's own parameters
        tmp <- (sim@params@rr_bb * sim@params@cc_bb / (sim@params@rr_bb + m2_benthos))
        n_bb <- tmp - (tmp - n_bb) * exp(-(sim@params@rr_bb + m2_benthos) * dt)
        n_bb[sim@params@initial_n_bb == 0] <- 0 # destroy what's below (and above) threshold sizes
        
        # Dynamics of the algal spectrum uses a semi-chemostat model 
        # currently it follows exactly the same rules as plankton but has it's own parameters
        tmp <- (sim@params@rr_aa * sim@params@cc_aa / (sim@params@rr_aa + m2_algae))
        n_aa <- tmp - (tmp - n_aa) * exp(-(sim@params@rr_aa + m2_algae) * dt)
        n_aa[sim@params@initial_n_aa == 0] <- 0 # destroy what's below (and above) threshold sizes
        ##AAsp##
        
        # CN                       
        if(fleetDynamics==TRUE){

          # calcualte yield, revenue and profit - all based on biomass at current time step (B_itime) as per EULER method

          # yield (spp X w X fleet) yield by fleet and total
          yield_itime <-sweep(F_itime,c(1,2), B_itime, "*")*scaling_price

          # yield for fleetdyncamic eqn (no effort and considering target matrix)
          yield_itime_FL <- sweep(F_itime_FL,c(1,2), B_itime, "*")*scaling_price
          # do not consider target is target == catchability and catchability was used to calculate F_itime
          # yield_itime_FL <- sweep(yield_itime_FL,c(1,3), t(sim@params@target), "*")

          # price
          # price_itime<-split(price_dt, price_dt$time)
          # price_itime<-price_itime[[i_time]]
          price_itime<-price_dt[i_time,]

          # revenue ()
          # revenue_itime<-sweep(yield_itime, c(1), price_itime$price,"*")
          revenue_itime<-sweep(yield_itime, c(1), price_itime,"*")
          revenue_itime_tot<-apply(revenue_itime, 3, sum)

          # revenue for fleetdyncamic eqn
          # revenue_itime_FL<-sweep(yield_itime_FL, c(1), price_itime$price,"*")
          revenue_itime_FL<-sweep(yield_itime_FL, c(1), price_itime,"*")
          revenue_itime_FL_tot<-apply(revenue_itime_FL, 3, sum)

          # cost
          cost_itime<-cost_dt[i_time,]

          # profit
          profit_itime<-revenue_itime_tot - cost_itime*effortOut_dt[i_time,] # NOTE: it was effort_dt ?!
          profit_itime_FL<-revenue_itime_FL_tot - cost_itime

          # fleet adjusment constant
          # ke= 0.00000001
          ke= ke

          # effort will be calcuated in loop and depends on the management component
          Effort_itime_change <-list()
          Effort_itime <-list()

          if(management == TRUE){

            # biomass reference points for management at start of projections (i_time = 1) or as given 
            
            if("bioUnfished" %in% colnames(params@species_params) == TRUE){
              Blim<-params@species_params$Bref
            }else{
              Blim<-rowSums(sweep(sim@n[1,,],2,sim@params@w * sim@params@dw,"*"))
            }
            
            Blim<-data.frame(species = species, bioUnfished = Blim, bio20 = Blim*0.2, bio40 = Blim*0.4, bio48 = Blim*0.48)
            
            # biomass level
            Blevel<-data.frame(species = species, bioLevel = rowSums(B_itime))

            # biomass level vs biomass ref points
            Bio<-merge(Blim, Blevel)
            Bio$bio48check<-ifelse(Bio$bioLevel > Bio$bio48, "Above","Below")
            Bio$bio40check<-ifelse(Bio$bioLevel > Bio$bio40, "Above","Below")
            Bio$bio20check<-ifelse(Bio$bioLevel > Bio$bio20, "Above","Below")
            Bio<-melt(Bio[,c("species","bio48check","bio40check","bio20check")], id.vars = "species")
            colnames(Bio)<-c("species","bioLim","bioStatus")

            # which spp is under management for each fishery?
            # you can use sim@params@target or sim@params@catchability
            target<-as.data.frame(t(sim@params@target))
            target$species<-rownames(target)
            target<-melt(target,id.vars = "species")
            colnames(target)<-c("species","fleet","target")

            # merge information on biomass level and biomass ref points with target/bycatch spp
            Bio<-merge(Bio,target)
            Bio<-split(Bio, Bio$fleet)
            Bio<-lapply(Bio, function(x) x[-which(x$target == 0 | x$bioStatus =="Above"),])

            # calcualte effort of each fleet for the next time step and based on profits
            # whether effort is allowed to increase or forced to decrease is based on whether the biomass of speceis under management have reached the treshold biomass (e.g. 20, 40 or 48%). How much effort changes is given by profits, which are calcualted considering only the target spp of each fleet as we are considering fishers expectation on target

            for(i in 1:length(Bio)){

              # i=1
              if(is.na(match('bio48check',Bio[[i]]$bioLim)) & is.na(match('bio40check',Bio[[i]]$bioLim)) & is.na(match('bio20check',Bio[[i]]$bioLim)))
              {
                # if no spp below any bioref limts
                # print("all good")
                Effort_itime_change[i]<-ke[i,"ke"]*(profit_itime_FL[[i]])
                # instead use eqn with no effort at the denominator (and no target matrix)
                # Effort_itime_change[i]<-ke[i,"ke"]*(profit_itime[[i]])
                Effort_itime[i]<-effortOut_dt[i_time,][i] + (Effort_itime_change[[i]] * dt)

              }else if(!is.na(match('bio48check',Bio[[i]]$bioLim)) & is.na(match('bio20check',Bio[[i]]$bioLim)) & is.na(match('bio40check',Bio[[i]]$bioLim))){
                # if some specied below 48%
                # print("some below 48%")
                if(profit_itime[[i]]>0) { # if fishing is profitable, do not allow any increase in effort
                  Effort_itime_change[i]<-0
                  Effort_itime[i]<-effortOut_dt[i_time,][i] + (Effort_itime_change[[i]] * dt)
                }else{ # if fishing is not profitable and it is decreasing, let it decrease
                  Effort_itime_change[i]<-ke[i,"ke"]*(profit_itime_FL[[i]])
                  # instead use eqn with no effort at the denominator (and no target matrix)
                  # Effort_itime_change[i]<-ke[i,"ke"]*(profit_itime[[i]])
                  Effort_itime[i]<-effortOut_dt[i_time,][i] + (Effort_itime_change[[i]] * dt)
                }

              }else if(!is.na(match('bio48check',Bio[[i]]$bioLim)) & !is.na(match('bio20check',Bio[[i]]$bioLim)) & is.na(match('bio40check',Bio[[i]]$bioLim))){
                # if some species below 40%
                # print("some below 40%")
                Effort_itime_change[i]<-NA
                Effort_itime[i]<-effortOut_dt[i_time,][i] * (1 - 0.2*dt) # force a descrease in effort

              }else{
                # if some speceis below 20%
                # print("some below 20%")
                Effort_itime_change[i]<-0
                Effort_itime[i]<-0
              }

              names(Effort_itime_change)[i]<-names(Bio)[i]
              names(Effort_itime)[i]<-names(Bio)[i]
              
              # do not allow for negative effort
              Effort_itime[[i]]<-ifelse(is.na(Effort_itime[[i]]) | Effort_itime[[i]]<0 | is.infinite(Effort_itime[[i]]),0, Effort_itime[[i]])

            } # end of for loop

          }else{ # management = FALSE

            for(i in 1:length(fleet)){

              Effort_itime_change[i]<-ke[i,"ke"]*(profit_itime_FL[[i]])
              # instead use eqn with no effort at the denominator (and no target matrix)
              # Effort_itime_change[i]<-ke[i,"ke"]*(profit_itime[[i]])
              Effort_itime[i]<-effortOut_dt[i_time,][i] + (Effort_itime_change[[i]] * dt)

              names(Effort_itime_change)[i]<-fleet[i]
              names(Effort_itime)[i]<-fleet[i]
              Effort_itime[[i]]<-ifelse(is.na(Effort_itime[[i]]) | Effort_itime[[i]]<0 | is.infinite(Effort_itime[[i]]),0, Effort_itime[[i]])

            } # end of fo loop

          } # end of managemetn = FALSE

          # change format of effort
          Effort_itime_change<-t(do.call("rbind", Effort_itime_change))
          Effort_itime<-t(do.call("rbind", Effort_itime))

          # update effort for next time step
          effortOut_dt[i_time+1, ]<-Effort_itime

        } # end of fleetDynamics = TRUE
       
        # Store results only every t_step steps.
        store <- t_dimnames_index %in% (i_time + 1)
        if (any(store)) {

          # Advance progress bar
          # CN you can turn off bar when running calibration 
          # pb$tick()
          # if (hasArg(shiny_progress)) {
          #   shiny_progress$inc(amount = proginc)
          # }
          
          # # Store result
          # sim@n[which(store), , ] <- n
          # sim@n_pp[which(store), ] <- n_pp


            # Advance progress bar
          #  pb$tick() #Asta commented this out, because this produces huge amount of bars at every time step
            # if (hasArg(shiny_progress)) {
            #     shiny_progress$inc(amount = proginc)
            # }
            # Store result
            sim@n[which(store), , ] <- n 
            sim@n_pp[which(store), ] <- n_pp
            ##AAsp######
            sim@n_bb[which(store), ] <- n_bb
            sim@n_aa[which(store), ] <- n_aa
            ##AAsp#  
          
          # CN store fleetDynamics results
          if(fleetDynamics == TRUE){

            # sim@effort
            sim@effortOut[which(store), ] <- Effort_itime
            sim@yield[which(store),,,] <- yield_itime
            sim@revenue[which(store), ] <- revenue_itime_tot
            sim@profit[which(store), ] <- profit_itime
            sim@F[which(store),,,] <- F_itime
            
          }

        }

    } # end of for loop

    return(sim)

 } # end of function


#' Calculate initial population abundances for the community populations
#' 
#' This function uses the model parameters and other parameters to calculate 
#' initial population abundances for the community populations. These initial 
#' abundances should be reasonable guesses at the equilibrium values. The 
#' returned population can be passed to the \code{project} method.
#' 
#' @param params The model parameters. An object of type \linkS4class{MizerParams}.
#' @param a A parameter with a default value of 0.35.
#' @param n0_mult Multiplier for the abundance at size 0. Default value is
#'   kappa/1000.
#' @export
#' @return A matrix (species x size) of population abundances.
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' params <- MizerParams(NS_species_params_gears)
#' init_n <- get_initial_n(params)
#' }
get_initial_n <- function(params, n0_mult = NULL, a = 0.35) {
    if (!is(params,"MizerParams"))
        stop("params argument must of type MizerParams")
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    initial_n <- array(NA, dim = c(no_sp, no_w))
    dimnames(initial_n) <- dimnames(params@intake_max)
    # N = N0 * Winf^(2*n-q-2+a) * w^(-n-a)
    # Reverse calc n and q from intake_max and search_vol slots (could add get_n as method)
    n <- (log(params@intake_max[,1] / params@species_params$h) / log(params@w[1]))[1]
    q <- (log(params@search_vol[,1] / params@species_params$gamma) / log(params@w[1]))[1]
    # Guessing at a suitable n0 value based on kappa - this was figured out using trial and error and should be updated
    if (is.null(n0_mult)) {
        lambda <- 2 + q - n
        kappa <- params@cc_pp[1] / (params@w_full[1]^(-lambda))
        n0_mult <- kappa / 1000
    }
    initial_n[] <- unlist(tapply(params@w, 1:no_w, function(wx,n0_mult,w_inf,a,n,q)
        n0_mult * w_inf^(2 * n - q - 2 + a) * wx^(-n - a),
        n0_mult = n0_mult, w_inf = params@species_params$w_inf, a=a, n=n, q=q))
    #set densities at w > w_inf to 0
    initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_inf) w_inf<wx, w_inf=params@species_params$w_inf))] <- 0
    # Also any densities at w < w_min set to 0
    initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_min)w_min>wx, w_min=params@species_params$w_min))] <- 0    
    return(initial_n)
}
