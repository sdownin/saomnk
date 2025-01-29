rm(list=ls())  ## uncomment to clear environ before run
TIME_START <- Sys.time()
########################
##
##   SAOM-NK Runscript:
##
##   1. Test 1
##
##
#######################
# library(furrr)
# library(progressr)
# library(foreach)
# library(doParallel)
library(dplyr)


###############  Load R DEPENDENCIES ############################
## Directories
dir_proj <- 'C://Users//sdr8y//OneDrive - University of Missouri//Research//Search_networks//SaoMNK//R'
dir_proj_base <- 'C://Users//sdr8y//OneDrive - University of Missouri//Research//Search_networks'
dir_data <- 'D://Search_networks'

## INPUT FILE
grid.input.file <- 'SAOMNK_gridsearch_input_basic_M3.xlsx'


## Biparite Environment Search Simulation Class
SaomNkRSienaBiEnv      <- source(file.path(dir_proj, 'SAOM_NK_R6_model.R'))$value
# SaomNkRSienaBiEnv_base <- source(file.path(dir_proj, 'SAOM_NK_R6_base_model.R'))$value


## Working director
setwd(dir_proj)


#################### FUNCTIONS ####################################
getGridRunNameFromScenario <- function(scenario, environ_seed_params) {
  varnames <- names(scenario)
  vals <- unname(scenario)
  outvec <- c()
  for (i in 1:length(varnames)) {
    if (varnames[i] %in% environ_seed_params) {
      outvec <- c(outvec, paste(varnames[i], vals[i], sep='='))
    }
  }
  outstr <- paste(outvec, collapse='_')
  return(sprintf('_%s_', outstr))
}



##=========== GET LIST OF GRIDSEARCH INPUT PARAMETERS ==================================
# grid.input.file <- 'SAOMNK_gridsearch_params_base_v1.xlsx'
# grid.input.file <- file.choose()
inputs <- readxl::read_excel(file.path(dir_proj_base, grid.input.file))

##
environ_seed_params <- unique(inputs$name[inputs$type %in% c('environment','seed')])
eff_cov_params      <- unique(inputs$name[inputs$type %in% c('effects','coCovars','varCovars','dyCovars')])
eff_struct_params   <- unique(inputs$name[inputs$type %in% c('effects')])
eff_component_params<- unique(inputs$name[inputs$type %in% c('coCovars','varCovars','dyCovars') & grepl('self\\$compon.+',inputs$interaction1) ])
eff_strategy_params <- unique(inputs$name[inputs$type %in% c('coCovars','varCovars','dyCovars') & grepl('self\\$strat.+',inputs$interaction1) ])



##-------------- Get ALL Parameters --------------
parameters <- lapply(strsplit(inputs$parameter, '[|]', perl=T), as.numeric)
names(parameters) <- inputs$name

##-------------- Get Strategies --------------
strats <- inputs %>% filter(type=='coCovars', grepl('strat(egy|egies)?_\\d{1,2}',  interaction1) ) %>% as.list()
if (length(strats$x)!=0 & !any(is.na(strats$x))) {
  strats$x <- lapply(strsplit(strats$x, '[|]'), function(item){
    lapply(strsplit(item, '_'), as.numeric)
  })
}
strats <- lapply(strats, unique)
#
names(strats$x) <- strats$name
#
for (i in 1:length(strats$name)) {
  stratname <- strats$name[ i ]
  strat_interactation1 <- strats$interaction1[ i ]
  parameters[[stratname]] <- list()
  for (j in 1:length(strats$x[[stratname]])) {
    for (k in 1:length(strats$parameter)) {
      .id <- length(parameters[[stratname]]) + 1
      parameters[[stratname]][[ .id ]] <- list( 
        parameter    = as.numeric(strats$parameter[ k ] ) , 
        x            = strats$x[[stratname]][[ j ]],
        interaction1 = strat_interactation1
      )
    }
  }
}

##-------------- Get Components --------------
components <- inputs %>% filter(type=='coCovars', grepl('component(s)?_\\d{1,2}',  interaction1) ) %>% as.list()
#
for (i in 1:length(components$name)) {
  componame <- components$name[ i ]
  compo_interactation1 <- components$interaction1[ i ]
  paramvals <- parameters[[componame]]
  parameters[[componame]] <- list()
  for (j in 1:length(paramvals)) {
    .id <- length( parameters[[componame]]  ) + 1
    parameters[[componame]][[.id]] <- list(
      parameter    = paramvals[ j ],
      x            = components$x,
      interaction1 = compo_interactation1
    )
  }
}


##====================================================================



#########################################################
## FILE MANAGEMENT
#########################################################
grid.input.filebase <- gsub('\\.(xls|xlsx|csv)','',grid.input.file )
grid.input.ts <- gsub('\\.','', as.character(as.numeric(Sys.time())) )
dir_job <- file.path(dir_data, sprintf('__%s__%s',grid.input.filebase, grid.input.ts))
if (!dir.exists(dir_job)) {
  dir.create(dir_job)
}

## RDS files saved in grid_results directory within the gridsearch job directory
dir_job_results <- file.path(dir_job, 'grid_results')
if (!dir.exists(dir_job_results)) {
  dir.create(dir_job_results)
}



#########################################################
## GRIDSEARCH INDEX
#########################################################
## GET PARAMETER GRID dataframe FROM PARAMETER LIST CREATED ABOVE
paramgrid <- expand.grid(parameters)
paramgrid$time_id <- NA ## This field will hold TIME_START IDs for each model run in parallel
paramgrid$time_duration <- NA ## holds the difference of end time - start time (duration of 1 run)
paramgrid$ID <- 1:nrow(paramgrid)

## WRITE INDEX TO JOB DIR
write.csv(paramgrid[1,c('ID','time_id','time_duration')], file=file.path(dir_job, 'paramgrid_index.csv'), row.names = F)

## COPY INPUTS dataframe to JOB DIR
write.csv(inputs, file = file.path(dir_job, grid.input.file))




#########################################################
## BEGIN OUTPUT TXT FILE
#########################################################
sink(file = file.path(dir_job, sprintf('output_%s_%s.txt',grid.input.filebase, grid.input.ts) ) , type='output' )

cat(sprintf('\nJob start time: %s\n', as.character(TIME_START)))
cat(sprintf('\nMain loop start time: %s\n', as.character(Sys.time())))


cat(sprintf('\ndir_job: %s\n', dir_job))


## check parameter inputs
cat('\nParameters:\n')
print(parameters)

cat('\nparamgrid dimensions:\n')
print(dim(paramgrid))



#########################################################
##
## MAIN LOOP
##
#########################################################
# for (i in 1:nrow(paramgrid))
for (i in 1:10) ##**TEST**

  ## ADD TIME_START FOR SCENARIO ID
  paramgrid$time_id[i] <- gsub('\\.','', as.character(as.numeric(Sys.time())) )

  ## Simulation scenario i (inputs list for next simulation) 
  scenario <- paramgrid[i, ]
  
  ## get output file name base of this scenario for log and data output, etc.
  scenario_filebase <-  sprintf('%s_%s_%s',grid.input.filebase, grid.input.ts, scenario$time_id)
  
  ##****************************************##
  ## I. SIM SETUP 
  ##****************************************##
  
  ##----------------------
  ## 1. Environment:  determines what types of entities and strategies are permissible
  ##----------------------
  environ_params <- list(
    M = scenario$M,             ## Actors
    N = scenario$N,             ## Components
    BI_PROB = scenario$BI_PROB, ## Environmental Density (DGP hyperparameter)
    name = getGridRunNameFromScenario(scenario, environ_seed_params),
    rand_seed = scenario$environ_seed,
    component_matrix_start = 'rand',
    visualize_init = F
  )
  
  ##----------------------
  ## 2. ACTOR STRATEGY:  constant (for now)
  ##----------------------
  structure_model <- list(
    dv_bipartite = list(
      name = 'self$bipartite_rsienaDV',
      effects = list( ),  ## structural effects
      coCovars = list( ), ## strategy and component covariates (constant within 1 sim chain)
      varCovars = list() ## ignore for now
    )#**TODO** generalize for multiple DVs b
    ##         currently only works for 1 DV network (dv_bipartite)
  )
  ## 1. ADD EFFECTS
  for (j in 1:length(eff_struct_params)) {
    effname <- eff_struct_params[ j ]
    effval <- scenario[[ effname ]]
    structure_model$dv_bipartite$effects[[ j ]] <- list(
      effect= effname, 
      parameter= scenario[[ effname ]], 
      fix=TRUE, ## restricted from user changes currently
      dv_name='self$bipartite_rsienaDV' ## restricted from user changes currently
    )
  }
  ## 2. ADD COMPONENT
  for (j in 1:length(eff_component_params)) {
    effname <- eff_component_params[j]
    effvals <- scenario[[ effname ]]
    for (k in 1:length(effvals)) {
      ## replace function name with corresponding random vector from that function
      x_input_func <- components$x ## keep input function name
      if (components$x == 'runif') {
        set.seed(scenario$environ_seed) ## environ seed or run seed here ?
        cov_x <- runif(scenario$N, min=0, max=1)
      } else {
        stop('only runif() curently implemented for component_payoffs.')
      }
      structure_model$dv_bipartite$coCovars[[j]] <- list(
        effect= effname, 
        parameter= effvals[[ k ]]$parameter, 
        interaction1 = effvals[[ k ]]$interaction1,
        x = cov_x,
        fix=TRUE, ## restricted from user changes currently
        dv_name='self$bipartite_rsienaDV' ## restricted from user changes currently
      )
    }
  }
  ## 3. ADD STRATEGY
  for (j in 1:length(eff_strategy_params)) {
    effname <- eff_strategy_params[j]
    effvals <- scenario[[ effname ]]
    for (k in 1:length(effvals)) {
      .id <- length( structure_model$dv_bipartite$coCovars ) + 1
      structure_model$dv_bipartite$coCovars[[ .id ]] <- list(
        effect= effname, 
        parameter= effvals[[ k ]]$parameter, 
        interaction1 = effvals[[ k ]]$interaction1,
        x =  rep( effvals[[ k ]]$x,  times= scenario$M / length(effvals[[ k ]]$x) ),  ## repeats [1,2,3,...,s,  1,2,3,...,s,  ...] for M/#strats
        fix=TRUE, ## restricted from user changes currently
        dv_name='self$bipartite_rsienaDV' ## restricted from user changes currently
      )
    }
  }
  

  
  ##****************************************##
  ## II. SIM ANALYSIS 
  ##****************************************##
  ## INIT Environment 
  env1 <- SaomNkRSienaBiEnv$new(environ_params)
  ## Search LONG run
  env1$search_rsiena_multiwave_run(
    structure_model, 
    waves=1,
    iterations = 10000,  
    rand_seed = scenario$run_seed,
    dir_output = dir_job_results,
    file_output = scenario_filebase
  )
  # Process results
  env1$search_rsiena_multiwave_process_results()
  

  ##****************************************##
  ## III. I/O & FILE MANAGEMENT 
  ##****************************************##
  cat("\nSaving simulation scenario model results to RDS file...")
  saveRDS(env1, file = file.path(dir_job_results, sprintf('%s.rds',scenario_filebase) ) )
  cat("done.\n")
  ##
  
  ## UDPATE JOB INDEX
  paramgrid$time_duration[i] <-   Sys.time() - paramgrid$time_id[i]
  write.csv(paramgrid[i,c('ID','time_id','time_duration')], file=file.path(dir_job, 'paramgrid_index.csv'), row.names = F )
  # write.csv(paramgrid[i,c('ID','time_id','time_duration')], file=file.path(dir_job, sprintf('paramgrid_index_%s.csv',i)), row.names = F )
  ## end
} ##**end main loop**




## Ending time
TIME_END <- Sys.time()
cat(sprintf('\nEnd time: %s\n', TIME_END))
cat(sprintf('\nElapsed time: %s\n', TIME_END - TIME_START ))
print( TIME_END - TIME_START )


cat('\n\nSIMULATION COMPLETED SUCCESSFULLY.\n\n')


## return output to console
sink()


##**END**




