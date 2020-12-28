
## Shiny for Fmort model ----

library(shiny)

fn_shiny<-function(df_param, theta, kappa, kappa_ben, kappa_alg, w_pp_cutoff, min_w_bb, w_bb_cutoff, fleetDynamics, selectivity_params, catchability, target, effort, dt, management, price, cost, diet_steps, rescale){
  
  min = 0
  max = 10
  setp = 0.1
  
  ui <- fluidPage(
    
    titlePanel("Checking catches"),
    
    sidebarLayout(
      sidebarPanel(
        
        # cacthability
        sliderInput(inputId = "Q_lanternfish", label = "lanternfish", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_whiting", label = "whiting", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_squid", label = "squid", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_perch", label = "perch", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_mackerel", label = "mackerel", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_redfish", label = "redfish", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_deepShark", label = "deepShark", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_morwong", label = "morwong", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_flathead", label = "flathead", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_dory", label = "dory", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_bluew", label = "bluew", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_or", label = "or", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_blueg", label = "blueg", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_silverw", label = "silverw", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_gemfish", label = "gemfish", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_ling", label = "ling", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_sawshark", label = "sawshark", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_gummy", label = "gummy", min = min, max = max, value = 1,step = setp),
        sliderInput(inputId = "Q_school", label = "school", min = min, max = max, value = 1,step = setp)
        
      ),
      
      mainPanel(
        plotOutput(outputId = "distPlot"),
        plotOutput(outputId = "distPlot2")
      )
    )
  )
  
  # see https://github.com/Atlantis-Ecosystem-Model/ReactiveAtlantis/blob/master/R/compare.R
  server <- function(input, output) {
    
    # specify the model and things that change - reactive
    model<-reactive({
      
      # catchability
      df_param<-df_param %>% 
        mutate(r_max = case_when(spCommon == "lanternfish" ~ r_max*input$Q_lanternfish,
                                 spCommon == "school whiting" ~ r_max*input$Q_whiting,
                                        spCommon == "arrow squid" ~ r_max*input$Q_squid,
                                        spCommon == "bigeye ocean perch" ~ r_max*input$Q_perch,
                                        spCommon == "mackerel" ~ r_max*input$Q_mackerel,
                                        spCommon == "redfish" ~ r_max*input$Q_redfish,
                                        spCommon == "deepwater sharks" ~ r_max*input$Q_deepShark,
                                        spCommon == "jackass morwong" ~ r_max*input$Q_morwong,
                                        spCommon == "tiger flathead" ~ r_max*input$Q_flathead,
                                        spCommon == "dories and oreos" ~ r_max*input$Q_dory,
                                        spCommon == "blue warehou" ~ r_max*input$Q_bluew,
                                        spCommon == "orange roughy" ~ r_max*input$Q_or,
                                        spCommon == "blue grenadier" ~ r_max*input$Q_blueg,
                                        spCommon == "silver warehou" ~ r_max*input$Q_silverw,
                                        spCommon == "gemfish" ~ r_max*input$Q_gemfish,
                                        spCommon == "pink ling" ~ r_max*input$Q_ling,
                                        spCommon == "common sawshark" ~ r_max*input$Q_sawshark,
                                        spCommon == "gummy shark" ~ r_max*input$Q_gummy,
                                        spCommon == "school shark" ~ r_max*input$Q_school))
      # browser()
      # Q
      # these settings are those from sim_fitted in SEAmodel_run.Rmd
      params <- MizerParams(df_param, interaction = interaction, kappa = kappa, kappa_ben = kappa_ben, kappa_alg = kappa_alg, w_pp_cutoff = w_pp_cutoff, min_w_bb = min_w_bb, w_bb_cutoff = w_bb_cutoff, fleetDynamics = fleetDynamics, selectivity_params = selectivity_params, catchability = catchability, target = target)
      
      sim <- project(params, effort = effort, dt = dt, fleetDynamics = fleetDynamics, management = management, price = price, cost = cost, diet_steps = diet_steps)
      
      model<-sim
      
    })
    
    # define outputs plots
    output$distPlot <- renderPlot({
      compareTrends(model(),sim_FD, fleetDynamics = fleetDynamics, type = "yield", yieldObs_timeVariant,ssbObs, rescale = rescale, areaEco)$plotYield_sp
    })
    
    # define second output plots
    output$distPlot2<-renderPlot({
      compareTrends(model(),sim_FD, fleetDynamics = fleetDynamics,type = "SSB",yieldObs_timeVariant,ssbObs,rescale = rescale,areaEco)$plotYield_sp
    })
    
  }
  
  shinyApp(ui = ui, server = server)
  
}

## Shiny for FD model ----

fn_shiny_FD<-function(df_param, theta, kappa, kappa_ben, kappa_alg, w_pp_cutoff, min_w_bb, w_bb_cutoff, fleetDynamics, selectivity_params, catchability, target, effort, dt, management, price, cost, diet_steps, ke, initial_n, initial_n_pp, initial_n_bb, initial_effort, scaling_price, Blevel_management, sim, type, yieldObs_timeVariant,ssbObs,rescale, areaEco){

  min_q = 1
  max_q = 5
  step_q = 0.5
  
  ui <- fluidPage(
    
    titlePanel("Fleet dynamics"),
    
    sidebarLayout(
      sidebarPanel(
        
        # cacthability SH
        sliderInput(inputId = "redfish", label = "redfish_SH", min = min_q, max = max_q, value = 1,step = step_q),
        sliderInput(inputId = "morwong", label = "morwong_SH", min = min_q, max = max_q, value = 1,step = step_q),
        sliderInput(inputId = "bluew", label = "bluew_SH", min = min_q, max = max_q, value = 1,step = step_q),
        sliderInput(inputId = "silverw", label = "silverw_SH", min = min_q, max = max_q, value = 1,step = step_q),
        
        # cacthability UP
        sliderInput(inputId = "redfish_UP", label = "redfish_UP", min = min_q, max = max_q, value = 1,step = step_q),
        sliderInput(inputId = "morwong_UP", label = "morwong_UP", min = min_q, max = max_q, value = 1,step = step_q),
        sliderInput(inputId = "bluew_UP", label = "bluew_UP", min = min_q, max = max_q, value = 1,step = step_q),
        sliderInput(inputId = "silverw_UP", label = "silverw_UP", min = min_q, max = max_q, value = 1,step = step_q),
        sliderInput(inputId = "squid_UP", label = "squid_UP", min = min_q, max = max_q, value = 1,step = step_q)
      ),
      
      mainPanel(
        plotOutput(outputId = "distPlot"),
        plotOutput(outputId = "distPlot2", height = "200"),
        plotOutput(outputId = "distPlot3", height = "200")
      )
    )
  )
  
  server <- function(input, output) {
  
    model<-reactive({
      
      # catchability SH
      target[3, "centroberyx affinis"]<- target[3, "centroberyx affinis"]/input$redfish
      target[3, "nemadactylus macropterus"]<- target[3, "nemadactylus macropterus"]/input$morwong
      target[3, "seriolella brama"]<- target[3, "seriolella brama"]/input$bluew
      target[3, "seriolella punctata"]<- target[3, "seriolella punctata"]/input$silverw
      
      # catchability UP
      target[4, "centroberyx affinis"]<- target[4, "centroberyx affinis"]/input$redfish_UP
      target[4, "nemadactylus macropterus"]<- target[4, "nemadactylus macropterus"]/input$morwong_UP
      target[4, "seriolella brama"]<- target[4, "seriolella brama"]/input$bluew_UP
      target[4, "seriolella punctata"]<- target[4, "seriolella punctata"]/input$silverw_UP
      target[4, "nototodarus gouldi"]<- target[4, "nototodarus gouldi"]/input$squid_UP
      
      params <- MizerParams(df_param, interaction = interaction, kappa = kappa, kappa_ben = kappa_ben, kappa_alg = kappa_alg, w_pp_cutoff = w_pp_cutoff, min_w_bb = min_w_bb, w_bb_cutoff = w_bb_cutoff, fleetDynamics = fleetDynamics, selectivity_params = selectivity_params, catchability = catchability, target = target)
      
      sim_FD <- project(params, effort = effort, dt = dt, fleetDynamics = fleetDynamics, management = management, price = price, cost = cost, diet_steps = diet_steps, ke = ke, initial_n = initial_n, initial_n_pp = initial_n_pp, initial_n_bb = initial_n_bb, initial_effort = initial_effort, scaling_price = scaling_price, Blevel_management = Blevel_management)
      
      model <-sim_FD
      
    })
    
    # define outputs plots
    output$distPlot <- renderPlot({
      compareTrends(sim,model(),fleetDynamics = fleetDynamics, type = "yield", yieldObs_timeVariant,ssbObs, rescale = 3, areaEco)$plotYield_sp
    })
    
    # define second output plots
    output$distPlot2<-renderPlot({
      compareTrends(sim, model(), fleetDynamics = fleetDynamics, type = "yield", yieldObs_timeVariant,ssbObs, rescale = rescale, areaEco)$plotYield_fl
    })
    
    output$distPlot3<-renderPlot({
      compareTrends(sim, model(), fleetDynamics = fleetDynamics, type = "yield", yieldObs_timeVariant,ssbObs, rescale = rescale, areaEco)$plotEffort_fl
    })
    
  }
  
  shinyApp(ui = ui, server = server)
  
}

