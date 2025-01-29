rm(list=ls())  ## uncomment to clear environ before run
########################
##
##   SAOM-NK Runscript:
##
##   1. Estimation Recovery 
##
##
#######################
library(plyr)
library(dplyr)
library(uuid)
library(ggpubr)
library(data.table)


##****************************************************##
## ITERATIONS
MAX_ITERATIONS <- 3000
## Replications
nreps <- 20
##****************************************************##
#
.ts <- gsub('\\.','',as.character(as.numeric(Sys.time())))
#
# .uuid <- UUIDgenerate()
##  DIRECTORY
grid.dir <- sprintf('_SAOMNK_sensitivity_start0_rand20__%s', .ts )


###############  Load R DEPENDENCIES ############################
## Directories
dir_ext <- 'D:\\Search_networks'
dir_r <- 'C:\\Users\\sdr8y\\OneDrive - University of Missouri\\Research\\Search_networks\\SaoMNK\\R'
# dir_proj <- file.path( dir_ext, grid.dir )
# dir_proj_results <- file.path(dir_proj, 'sensitivity_results')

setwd(dir_ext)

# ## script run directory
# if (!dir.exists(dir_proj)) {
#   dir.create(dir_proj)
# }

## default settings: Users do not change; TODO: implment within restricted class attributes
DV_NAME <- 'self$bipartite_rsienaDV'


###############  Load R6 Class DEPENDENCIES ############################
## Biparite Environment Search Simulation Class
SaomNkRSienaBiEnv <- source(file.path(dir_r, 'SAOM_NK_R6_model.R'))$value
# ## RSiena search Class
# SaomNkRSienaBiEnv_search_rsiena <- source(file.path(dir_proj, 'SAOM_NK_R6_search_rsiena_model.R'))$value
###########




#########################################################
## FILE MANAGEMENT
#########################################################
# ## RDS files saved in grid_results directory within the gridsearch job directory
# if (!dir.exists(dir_proj_results)) {
#   dir.create(dir_proj_results)
# } else {
#   # ## if img dir exists, remove simufiles from this run list by UUIDs in existing plots
#   # plotted <- dir(dir_proj_results, pattern = '\\.png')
#   # skipuuids <- unique(sapply(strsplit(plotted, '_'), function(x) x[1] ))
#   # dropids <- unname(sapply(skipuuids, function(x) grep(x, simfiles)))
#   # if(length(dropids))
#   #   simfiles <- simfiles[ -dropids ] ## drop files with UUIDs already plotted
# }







# ## Working director
# setwd(dir_proj)




# #################### FUNCTIONS ####################################
# getGridRunNameFromScenario <- function(scenario, environ_seed_params) {
#   varnames <- names(scenario)
#   vals <- unname(scenario)
#   outvec <- c()
#   for (i in 1:length(varnames)) {
#     if (varnames[i] %in% environ_seed_params) {
#       outvec <- c(outvec, paste(varnames[i], vals[i], sep='='))
#     }
#   }
#   outstr <- paste(outvec, collapse='_')
#   return(sprintf('_%s_', outstr))
# }





strategies <- list(
  egoX   =   c(0), #c(-1,0, 1),
  inPopX =   c(0) #c(1,0, -1)
)



##_-----

environ_params <- list(
  M = 2,        ## Actors
  N = 3,       ## Components
  BI_PROB = 0, ## Environmental Density (DGP hyperparameter)
  component_matrix_start = 'rand', ##**TODO** Implement: 'rand','modular','semi-modular',...
  rand_seed = 1234,
  visualize_init = F,
  name = '_test_sensitivity_'
)

## 2.b. Component Payoffs vector
component_payoffs <-  runif(environ_params$N, min = 0, max = 1)
## 2. Strategies sets the objective function as a linear combination of network stats across DVs
#
actor_strats_list <- lapply(strategies, function(strat) rep(strat,  environ_params$M/length(strat)) )

#
structure_model <- list(
  dv_bipartite = list(
    name = 'self$bipartite_rsienaDV',
    effects = list( ##**STRUCTURAL EFFECTS -- dyadic/network endogeneity sources**
      list(effect='density', parameter= 0, fix=T, dv_name=DV_NAME), ##interaction1 = NULL
      list(effect='inPop',   parameter= 0,  fix=T, dv_name=DV_NAME), #interaction1 = NUL
      list(effect='outAct',  parameter= 0, fix=T, dv_name=DV_NAME)
    ),
    coCovars = list( ##**STRATEGY -- MONADIC CONSTANT COVARIATE EFFECTS **
      list(effect='altX',   parameter= 0.0000000000000001, fix=T,dv_name=DV_NAME,interaction1='self$component_1_coCovar', x= component_payoffs ),
      list(effect='egoX',   parameter= 0,fix=T,dv_name=DV_NAME, interaction1='self$strat_1_coCovar', x= actor_strats_list[[1]] ), #interaction1 = NULL
      list(effect='inPopX', parameter= 0, fix=T,dv_name=DV_NAME,interaction1='self$strat_2_coCovar', x= actor_strats_list[[2]] )
    ),
    varCovars = list() ##**MONADIC TIME-VARYING COVARIATE EFFECTS -- DYNAMIC STRATEGY PROGRAMS**
  )
)

# ##****************************************##
# ## II. SIM ANALYSIS 
# ##****************************************##
# # ## INIT SIM ENVIRONMENT: 
# # env1 <- SaomNkRSienaBiEnv$new(environ_params)
# 
# ## INIT SIM ENVIRONMENT: 
env1 <- SaomNkRSienaBiEnv$new(environ_params)
# 
# # pd_steps <- env1$M * env1$N
# pd_steps <- 2 * env1$M 
# 
# ITERATIONS <- min( 40 * pd_steps, MAX_ITERATIONS)
# 
# ## 1.1. Search 1. SHORT Run
env1$search_rsiena_multiwave_run(
  structure_model,
  waves=1, ##id='_short_run',
  iterations = 500, #4 * env1$M * env1$N,
  rand_seed = 12345
)
# # Process results
env1$search_rsiena_multiwave_process_results()



##############################################################
self <- env1
# nwaves <- 5
# step_ids <- round(seq(2*self$M*self$N, self$rsiena_model$n3, by=1*self$M))[1:nwaves] ##, length=nwaves))
# step_ids <-  self$M + seq(1, self$rsiena_model$n3, by=3)[1:nwaves] ##, length=nwaves))
rsiena_phase2_nsub <- 1
rsiena_n2start_scale <- 1
iterations <- 1000
rand_seed <- 54321
returnDeps <- T
returnChains <- T


# apply(self$bi_env_arr[,,1:5],3, sum)

#############################
#
##--1. INIT RSiena Model--------
set.seed(rand_seed)
# cat('\n\nDEBUG: called init_rsiena_model_from_structure_model_bipartite_matrix()\n\n')
#
ACTORS     <- sienaNodeSet(self$M, nodeSetName="ACTORS")
COMPONENTS <- sienaNodeSet(self$N, nodeSetName="COMPONENTS")
#
# self$config_structure_model <- structure_model
structure_model <- self$config_structure_model
structure_model_dvs <- names(structure_model)
#
# hasVarCovars <- 'varCovars' %in% names(structure_model$dv_bipartite)
hasCoCovars <- 'coCovars' %in% names(structure_model$dv_bipartite)
#
print('DEBUG:  hasCoCovars: ')
print(hasCoCovars)

## init networks to duplicate for the init arrays (two network waves)
# array_bi_net <- array(c(bipartite_matrix1, bipartite_matrix2), dim=c(self$M, self$N, 2) )
#
# array_bi_net <-  self$bi_env_arr[,,step_ids] ##**MAIN STEP**











array_bi_net <-  self$bi_env_arr[,,c(1:15)] ##**MAIN STEP**

# self$bi_env_arr[,,1] + self$bi_env_arr[,,2] + self$bi_env_arr[,,3] 
# 
# self$bi_env_arr[,,4] + self$bi_env_arr[,,5] + self$bi_env_arr[,,6] 


## Drop information above binary ties for RSiena DVs
# array_bi_net[ array_bi_net > 1 ] <- 1
self$bipartite_rsienaDV <- sienaDependent(array_bi_net, type='bipartite',
                                          nodeSet =c('ACTORS', 'COMPONENTS'), allowOnly = F)








##---------------------------------------------
self$rsiena_data <- self$get_rsiena_data_from_structure_model(structure_model)
print('self$rsiena_data : ')
print(self$rsiena_data)

##  2. Init effects
self$rsiena_effects <- getEffects(self$rsiena_data)
## UPDATE PARAMETERS TO FIXED=FALSE so we can estimate them
##  3. Add effects from model objective function list
self$add_rsiena_effects(structure_model)
self$rsiena_effects$fix <- rep(FALSE, length(self$rsiena_effects$fix))

params_norates <- which(self$rsiena_effects$type[self$rsiena_effects$include] != 'rate')
names_theta_in <- self$rsiena_effects$effectName[self$rsiena_effects$include][params_norates]
theta_in <- self$rsiena_effects$parm[self$rsiena_effects$include][params_norates] ## skip the rate parameter _effect[1]

# nrows <- self$rsiena_model$n3
nrows <- 4 * self$M 
ncols <- length(theta_in)
shock.splits <- 8
theta_matrix <- matrix(NA,  nrow = nrows, ncol=ncols  )
for (j in 1:ncols) {
  nportion   <- floor( nrow(theta_matrix) / shock.splits ) 
  beta_seq_1 <- rep(  1 * theta_in[j], nportion * 8 )  ## How many pd_chunks at what theta value?
  beta_seq_2 <- c() # rep(  1 * theta_in[j], nportion * 1 ) ## number of pd_chunks at what theta value?
  beta_seq_3 <- c() # rep(  1 * theta_in[j], nportion * 1  ) ## number of pd_chunks at what theta value?
  vec <- c( beta_seq_1, beta_seq_2, beta_seq_3 )
  rowdiff <- nrow(theta_matrix) - length(vec)
  if (rowdiff != 0) {
    vec <- c( vec, rep(vec[length(vec)], round(abs(rowdiff))) ) ## add the last elements again if too short
  }
  #
  theta_matrix[, j] <- vec
}

# self$rsiena_model$theta
  
  
##  4. RSiena Algorithm
self$rsiena_algorithm <- sienaAlgorithmCreate(projname=sprintf('%s_%s',self$SIM_NAME,self$TIMESTAMP),
                                              simOnly = T,
                                              # nsub = rsiena_phase2_nsub * 1,
                                              nsub = 0,
                                              # n2start = rsiena_n2start_scale * 2.52 * (7+sum(self$rsiena_effects$include)),
                                              n3 = nrow(theta_matrix),
                                              seed = rand_seed)
## 5. Run RSiena simulation
self$rsiena_model <- siena07(self$rsiena_algorithm,
                             data = self$rsiena_data, 
                             effects = self$rsiena_effects,
                             batch = TRUE,
                             returnDeps = TRUE, 
                             returnChains = TRUE,
                             returnThetas = TRUE,
                             returnDataFrame = TRUE, ##**TODO** CHECK
                             returnLoglik = TRUE ,  ##**TODO** CHECK
                             thetaValues = theta_matrix 
)   # returnChains = returnChains

# # Summarize and plot results
# mod_summary <- summary(self$rsiena_model)
# if(!is.null(mod_summary))
#   print(mod_summary)
# # plot(self$rsiena_model)
# digits <- 3
# print(screenreg(list(self$rsiena_model), single.row = T, digits = digits))


# print(self$rsiena_model$tconv)
# print( abs(self$rsiena_model$tconv) <  0.1 )
# print(self$rsiena_model$tconv.max)
# print(self$rsiena_model$tconv.max[1] < 0.25)
# all(
#   abs(self$rsiena_model$tconv) <  0.1 &
#     self$rsiena_model$tconv.max[1] < 0.25
# )


##
self$search_rsiena_multiwave_process_results()



# xints <- (1:shock.splits) * floor(max(self$actor_util_df$chain_step_id)/shock.splits)
xints <- c( length(beta_seq_1), length(beta_seq_1) + length(beta_seq_2))
thin_factor <- 1
actthin <- self$actor_util_df %>%  filter(chain_step_id %% thin_factor == 0)
print(dim(actthin))
tmpdf <-  actthin %>% mutate(actor_id=actor_id)
plt.act <- tmpdf %>%       #ungroup() %>%
  ggplot(aes(x=chain_step_id, y=utility, color=strategy)) + #
  geom_point(pch=16, alpha=.8, size=3) +
  # geom_point(pch=5, alpha=.5) +
  geom_smooth(aes(fill=strategy), method='loess', span=.4, alpha=.15) +
  geom_smooth(aes(x=chain_step_id, y=mean), method='loess', color='black', span=.4, alpha=.05, linewidth=1.1,
              data=actthin %>% group_by(chain_step_id) %>% summarize(mean=mean(utility, na.rm=T)) %>% 
                mutate(strategy=NA)) +
  # geom_smooth(aes(x=chain_step_id, y=utility_mean), method = 'loess', span=loess_span, 
  #             data=actutilstrat%>%group_by(chain_step_id)%>%summarize(utility_mean=mean(utility)),
  #             color='black', linetype=1, alpha=.1) +
  geom_vline(xintercept =xints, linetype=1, 
             color='black') +
  geom_hline(yintercept = 0, linetype=4, color='black') +
  theme_bw()
print(plt.act)





## Individual actors
# xints <- (1:shock.splits) * floor(max(self$actor_util_df$chain_step_id)/shock.splits)
xints <- c( length(beta_seq_1), length(beta_seq_1) + length(beta_seq_2))
thin_factor <- 1
actthin <- self$actor_util_df %>%  filter(chain_step_id %% thin_factor == 0)
print(dim(actthin))
tmpdf <-  actthin %>% mutate(actor_id=actor_id)
plt.act <- tmpdf %>%       #ungroup() %>%
  ggplot(aes(x=chain_step_id, y=utility, color=strategy)) + #
  geom_point(pch=16, alpha=.8, size=3) +
  # geom_point(pch=5, alpha=.5) +
  geom_smooth(aes(linetype=actor_id), method='loess', span=.4, alpha=.15) +
  geom_smooth(aes(x=chain_step_id, y=mean), method='loess', color='black', span=.4, alpha=.05, linewidth=1.1,
               data=actthin %>% group_by(chain_step_id) %>% summarize(mean=mean(utility, na.rm=T)) %>% 
                mutate(strategy=NA)) +
  # geom_smooth(aes(x=chain_step_id, y=utility_mean), method = 'loess', span=loess_span, 
  #             data=actutilstrat%>%group_by(chain_step_id)%>%summarize(utility_mean=mean(utility)),
  #             color='black', linetype=1, alpha=.1) +
  geom_vline(xintercept =xints, linetype=1, 
             color='black') +
  geom_hline(yintercept = 0, linetype=4, color='black') +
  theme_bw()
print(plt.act)








self$rsiena_model$x




# print( self$search_rsiena_multiwave_plot_actor_utility_strategy_summary(thin_factor = 3) )
# 
# 
# self$
# 
# ##





# gof.od <- RSiena::sienaGOF(self$rsiena_model, OutdegreeDistribution, levls=0:self$N, varName = 'self$bipartite_rsienaDV')
# # gof.od <- RSiena::sienaGOF(self$rsiena_model, OutdegreeDistribution, varName = 'self$bipartite_rsienaDV')
# plot(gof.od)
# 
# gof.id <- RSiena::sienaGOF(self$rsiena_model, IndegreeDistribution, levls=0:self$M, varName = 'self$bipartite_rsienaDV')
# plot(gof.id)














