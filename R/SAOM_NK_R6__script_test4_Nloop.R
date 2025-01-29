rm(list=ls())  ## uncomment to clear environ before run
########################
##
##   SAOM-NK Runscript:
##
##   1. Test 1
##
##
#######################

## Directories
dir_proj <- 'C://Users//sdr8y//OneDrive - University of Missouri//Research//Search_networks//SaoMNK//R'
dir_data <- 'D://Search_networks'

## default settings: Users do not change; TODO: implment within restricted class attributes
DV_NAME <- 'self$bipartite_rsienaDV'


outlist <- list()
Ns <- 2 * 3^c(1:5)  ## 6 12 24 48 96
for (i in 1:length(Ns)) 
{
  
  
  ##****************************************##
  ## I. SIM SETUP 
  ##****************************************##
  ##----------------------
  ## 1. Environment:  determines what types of entities and strategies are permissible
  ##----------------------
  environ_params <- list(
    M = 3,        ## Actors
    N = Ns[ i ],       ## Components
    BI_PROB = 0, ## Environmental Density (DGP hyperparameter)
    component_matrix_start = 'rand', ##**TODO** Implement: 'rand','modular','semi-modular',...
    rand_seed = 123,
    visualize_init = F,
    name = sprintf('_testSAOMNK_4_loopN%s', Ns[ i ] )#,
    ## use starting matrix param that is character ('modular',etc) or random (random seed 12345)
  )
  
  
  ##----------------------
  ## 2. ACTOR STRATEGY:  constant (for now)
  ##----------------------
  
  # 2.a. Actor Strategies List (covariates X of SAOM structural covariate effecs like outActX,inPopX,simX,... )
  strategies <- list(
    egoX   =  c(-1,0, 1), #c(0),
    inPopX =  c(1,0, -1)  #c(0),
  )
  # strategies <- list(
  #   egoX   =  c(0),
  #   inPopX =  c(0)
  # )
  
  ## 2.b. Component Payoffs vector
  component_payoffs <-  runif(environ_params$N, min = 0, max = 1)
  
  
  ## 2. Strategies sets the objective function as a linear combination of network stats across DVs
  #
  actor_strats <- lapply(strategies, function(strat) rep(strat,  environ_params$M/length(strat)) )
  #
  structure_model <- list(
    dv_bipartite = list(
      name = 'self$bipartite_rsienaDV',
      effects = list( ##**STRUCTURAL EFFECTS -- dyadic/network endogeneity sources**
        list(effect='density', parameter= -10, fix=T, dv_name=DV_NAME), ##interaction1 = NULL
        list(effect='inPop',   parameter=  -1,  fix=T, dv_name=DV_NAME), #interaction1 = NUL
        list(effect='outAct',  parameter= -1, fix=T, dv_name=DV_NAME)#, #interaction1 = NULL
        # list(effect='outInAss', parameter=0, fix=F, dv_name=DV_NAME), #interaction1 = NULL
        # list(effect='cycle4', parameter=.5, fix=T, dv_name=DV_NAME)#, #interaction1 = NULL
      ),
      coCovars = list( ##**STRATEGY -- MONADIC CONSTANT COVARIATE EFFECTS **
        list(effect='altX',   parameter= 15, fix=T,dv_name=DV_NAME,interaction1='self$component_1_coCovar', x= component_payoffs ),
        # list(effect='outActX',parameter= .5,fix=T,dv_name=DV_NAME,interaction1='self$component_1_coCovar', x= component_payoffs ), #interaction1 = NULL
        list(effect='egoX',   parameter= 1,fix=T,dv_name=DV_NAME,interaction1='self$strat_1_coCovar', x= actor_strats[[1]] ), #interaction1 = NULL
        list(effect='inPopX', parameter= 1, fix=T,dv_name=DV_NAME,interaction1='self$strat_2_coCovar', x= actor_strats[[2]] )#, #interaction1 = NULL
        # list(effect=c('egoX','inPopX'), parameter= -2, fix=T,dv_name=DV_NAME,interaction1=c('self$strat_1_coCovar','self$strat_2_coCovar'), x= actor_strats[[2]] )#, #interaction1 = NULL
        # list(effect=c('outActX','inPopX'),parameter= -.3, fix=T, dv_name=DV_NAME,interaction1=c('self$strat_1_coCovar','self$strat_2_coCovar'),
        #      x = c(1,0,-1,1,0,-1,1,0,-1) )#, #interaction1 = NULL
      ),
      varCovars = list() ##**MONADIC TIME-VARYING COVARIATE EFFECTS -- DYNAMIC STRATEGY PROGRAMS**
    )
  )
  
  
  
  
  
  
  
  ##****************************************##
  ## II. SIM ANALYSIS 
  ##****************************************##
  ###############  Load R6 Class DEPENDENCIES ############################
  ## Biparite Environment Search Simulation Class
  SaomNkRSienaBiEnv <- source(file.path(dir_proj, 'SAOM_NK_R6_model.R'))$value
  # ## RSiena search Class
  # SaomNkRSienaBiEnv_search_rsiena <- source(file.path(dir_proj, 'SAOM_NK_R6_search_rsiena_model.R'))$value
  ###########
  ## Working director
  setwd(dir_proj)
  
  
  
  
  ## INIT SIM ENVIRONMENT: 
  env1 <- SaomNkRSienaBiEnv$new(environ_params)
  ## 1.1. Search 1. SHORT Run
  env1$search_rsiena_multiwave_run(
    structure_model, 
    waves=1, ##id='_short_run',
    iterations = env1$M * env1$N, 
    rand_seed = 12345
  )
  
  
  
  ## 1.2. Search 2. LONG Run
  env2 <- SaomNkRSienaBiEnv$new(environ_params)
  env2$search_rsiena_multiwave_run(
    structure_model, 
    waves=1,
    iterations = 8000,  ## 50 * env2$M * round(sqrt(env2$N)), 
    rand_seed = 12345
  )
  
  # Process results
  env1$search_rsiena_multiwave_process_results()
  env2$search_rsiena_multiwave_process_results()
  
  # Plot Utility
  env1$search_rsiena_multiwave_plot('utility_strategy_summary', 
                                    plot_file = as.character(as.numeric(Sys.time())) )  ## thin_factor = 1
  env2$search_rsiena_multiwave_plot('utility_strategy_summary',  thin_factor = 5, 
                                    plot_file = as.character(as.numeric(Sys.time())) ) ## thin_factor = 1
  
  # Plot K 4panel (degree evolution)
  env1$search_rsiena_multiwave_plot('K_4panel', 
                                    plot_file = as.character(as.numeric(Sys.time())) )  ## thin_factor = 1
  env2$search_rsiena_multiwave_plot('K_4panel', thin_factor = 5, 
                                    plot_file = as.character(as.numeric(Sys.time())) ) ## thin_factor = 1
  
  
  
  outlist[[ sprintf('N%s',Ns[i]) ]] <- list(env1=env1, env2=env2, N=Ns[i])
}



saveRDS(outlist, 
        file = file.path(dir_data, sprintf('%s_%s.rds',environ_params$name,as.character(as.integer(Sys.time())))) )
## end




