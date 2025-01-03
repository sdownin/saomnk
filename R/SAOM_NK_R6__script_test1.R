########################
##
##   SAOM-NK Runscript:
##
##   1. Test 1
##
##
#######################

rm(list=ls())

## Directories
dir_proj <- 'C://Users//sdr8y//OneDrive - University of Missouri//Research//Search_networks//SaoMNK//R'
dir_data <- 'D://Search_networks'

setwd(dir_proj)

###############   DEPENDENCIES ###############################
SaomNkRSienaBiEnv <- source(file.path(dir_proj, 'SAOM_NK_R6_model.R'))$value


## default settings: Users do not change; TODO: implment within restricted class attributes
dv_name='self$bipartite_rsienaDV'






##****************************************##
## I. SIM SETUP 
##****************************************##


##----------------------
## 1. Environment:  determines what types of entities and strategies are permissible
##----------------------
environ_params <- list(
  M = 12,        ## Actors
  N = 24,       ## Components
  BI_PROB = 1, ## Environmental Density (DGP hyperparameter)
  component_matrix_start = 'rand', ##**TODO** Implement: 'rand','modular','semi-modular',...
  rand_seed = 123,
  visualize_init = F,
  name = '_TESTrsienaPayoffs3_'#,
  ## use starting matrix param that is character ('modular',etc) or random (random seed 12345)
)


##----------------------
## 2. ACTOR STRATEGY:  constant (for now)
##----------------------

## 2.a. Actor Strategies List (covariates X of SAOM structural covariate effecs like outActX,inPopX,simX,... )
actor_strats <- list(
  egoX   = rep( c(-1,0, 1),  4),
  inPopX = rep( c(1,0, -1),   4)#,
)

## 2.b. Component Payoffs vector
component_payoffs <- 1 * runif(environ_params$N, min = 0, max = 1)

## 2. Strategies sets the objective function as a linear combination of network stats across DVs
structure_model <- list(
  dv_bipartite = list(
    name = 'self$bipartite_rsienaDV',
    effects = list( ##**STRUCTURAL EFFECTS -- dyadic/network endogeneity sources**
      list(effect='density', parameter= -1, fix=T, dv_name=dv_name), ##interaction1 = NULL
      list(effect='inPop',   parameter= .1,  fix=T, dv_name=dv_name), #interaction1 = NUL
      list(effect='outAct',  parameter= .1, fix=T, dv_name=dv_name)#, #interaction1 = NULL
      # list(effect='outInAss', parameter=0, fix=F, dv_name=dv_name), #interaction1 = NULL
      # list(effect='cycle4', parameter=.5, fix=T, dv_name=dv_name)#, #interaction1 = NULL
    ),
    coCovars = list( ##**STRATEGY -- MONADIC CONSTANT COVARIATE EFFECTS **
      list(effect='altX',   parameter= 1, fix=T,dv_name=dv_name,interaction1='self$component_1_coCovar', x= component_payoffs ),
      list(effect='outActX',parameter= .5,fix=T,dv_name=dv_name,interaction1='self$component_1_coCovar', x= component_payoffs ), #interaction1 = NULL
      list(effect='egoX',   parameter= 1,fix=T,dv_name=dv_name,interaction1='self$strat_1_coCovar', x= actor_strats[[1]] ), #interaction1 = NULL
      list(effect='inPopX', parameter= .3, fix=T,dv_name=dv_name,interaction1='self$strat_2_coCovar', x= actor_strats[[2]] )#, #interaction1 = NULL
      # list(effect=c('egoX','inPopX'), parameter= -2, fix=T,dv_name=dv_name,interaction1=c('self$strat_1_coCovar','self$strat_2_coCovar'), x= actor_strats[[2]] )#, #interaction1 = NULL
      # list(effect=c('outActX','inPopX'),parameter= -.3, fix=T, dv_name=dv_name,interaction1=c('self$strat_1_coCovar','self$strat_2_coCovar'),
      #      x = c(1,0,-1,1,0,-1,1,0,-1) )#, #interaction1 = NULL
    ),
    varCovars = list() ##**MONADIC TIME-VARYING COVARIATE EFFECTS -- DYNAMIC STRATEGY PROGRAMS**
  )
)
# strat  <- c(rep(-1, round(environ_params$M / 2)), rep(1, round(environ_params$M / 2)))
# strat <- c(-1, 0, 1, -1, 0, 1, -1, 0)




##****************************************##
## II. SIM ANALYSIS 
##****************************************##

## 1. Short Run: INIT SIM object: 
m1 <- SaomNkRSienaBiEnv$new(environ_params)

# 1.1. SHORT RUN
m1$search_rsiena_multiwave_run(structure_model, waves=1, 
                               iterations =500, #1 * m1$M*m1$N, 
                               rand_seed = 12345)
# 1.2. process the chain
m1$search_rsiena_multiwave_process_results()
# 1.3.1. plot the actor utility timeseries
m1$search_rsiena_multiwave_plot('utility_strategy_summary', thin_factor = 1)
# 1.3.2. plot the degree summary (K4) panel 
m1$search_rsiena_multiwave_plot('K_4panel', thin_factor = 1)
 

## 2. Long Run
m1$search_rsiena_multiwave_run(structure_model, waves=1, iterations =2500, rand_seed = 12345)
m1$search_rsiena_multiwave_process_results()
m1$search_rsiena_multiwave_plot('utility_strategy_summary', thin_factor = 1)
m1$search_rsiena_multiwave_plot('K_4panel', thin_factor = 1)


##


################## SIM ANALYSIS ########################################
## 1. INIT SIM object
m1 <- SaomNkRSienaBiEnv$new(environ_params)
## 2. RUN SIM
m1$search_rsiena_run(structure_model, iterations=1000, get_eff_doc = F, run_seed=123)
## 3. Plot
m1$search_rsiena_plot()
#
m1$search_rsiena_plot('utility_by_strategy', smooth_method = 'loess', show_utility_points = F)
m1$search_rsiena_plot('utility_density_by_strategy')


#
m1$search_rsiena_plot('utility_by_strategy', smooth_method = 'loess', show_utility_points = T)
m1$search_rsiena_plot('utility_by_strategy', smooth_method = 'loess', show_utility_points = FALSE)
m1$search_rsiena_plot('utility_density_by_strategy')
m1$search_rsiena_plot('utility_histogram_by_strategy')


##
m1$search_rsiena_plot('utility', smooth_method = 'loess', show_utility_points = FALSE)


m1$search_rsiena_plot('utility', smooth_method = 'loess', show_utility_points = FALSE)
m1$search_rsiena_plot('utility', smooth_method = 'gam', show_utility_points = FALSE)
m1$search_rsiena_plot('utility', smooth_method = 'lm', show_utility_points = FALSE)
#
m1$search_rsiena_plot('utility', smooth_method = 'loess', show_utility_points = TRUE)
m1$search_rsiena_plot('utility', smooth_method = 'gam', show_utility_points = TRUE)
m1$search_rsiena_plot('utility', smooth_method = 'lm', show_utility_points = TRUE)
#
m1$search_rsiena_plot('utility_density')
#

##
m1$search_rsiena_plot('utility')
m1$search_rsiena_plot_actor_utility()
# ## 4. Plot ministep chain progress 
# m1$search_rsiena_plot('ministep_count') ## plot sequence of choices by actor
# m1$search_rsiena_plot('utility', rolling_window = 1)
# m1$search_rsiena_plot('utility', rolling_window = round(.1*iterations) )
# m1$search_rsiena_plot('utility', rolling_window = round(.02*iterations) )
## 5. Plot simulation system stability
m1$search_rsiena_plot('stability')
#################### SIM END #############################################
































































































# 
# self = m1
# 
# forebear_rate <- ( nrow(self$chain_stats) - length(self$rsiena_model$chain) ) / length(self$rsiena_model$chain)
# 
# 
# # getStructModNetStatsFromBipartiteMat <- function(bi_env_mat) {
# #   # eff <- m1$rsiena_effects[m1$rsiena_effects$include, ]
# #   efflist <- self$config_structure_model$dv_bipartite$effects
# #   # set the bipartite environment network matrix
# #   # bi_env_mat <- self$bipartite_matrix
# #   ### empty matrix to hold actor network statistics 
# #   mat <- matrix(rep(0, m1$M * length(efflist) ), nrow=self$M, ncol=length(efflist) )
# #   effnames <- sapply(efflist, function(x) x$effect, simplify = T)
# #   effparams <- sapply(efflist, function(x) x$parameter, simplify = T)
# #   colnames(mat) <- effnames
# #   rownames(mat) <- 1:self$M
# #   #
# #   for (i in 1:length(efflist)) 
# #   {
# #     eff_name <- efflist[[ i ]]$effect
# #     #
# #     xActorDegree  <- rowSums(bi_env_mat, na.rm=T)
# #     xComponentDegree  <- colSums(bi_env_mat, na.rm=T)
# #     ## network statistics dataframe
# #     #
# #     if (eff_name == 'density' ) 
# #     {
# #       mat[ , i] <- c( xActorDegree )
# #     } 
# #     else if (eff_name == 'outAct' )
# #     {
# #       mat[ , i]  <- c( xActorDegree^2  )
# #     } 
# #     else if (eff_name == 'inPop' ) 
# #     {
# #       mat[ , i]  <- c( bi_env_mat %*% (xComponentDegree + 1) )  ## MxN %*% N
# #     } 
# #     else  
# #     {
# #       print(sprintf('Effect not yet implemented: `%s`', eff_name))
# #     }
# #   }
# #   return(mat)
# # }
# 
# # get_struct_mod_net_stats_list_from_bi_mat <- function(bi_env_mat, type='all') {
# #   # eff <- m1$rsiena_effects[m1$rsiena_effects$include, ]
# #   efflist <- self$config_structure_model$dv_bipartite$effects
# #   # set the bipartite environment network matrix
# #   # bi_env_mat <- self$bipartite_matrix
# #   ### empty matrix to hold actor network statistics
# #   df <- data.frame() #nrow=self$M, ncol=length(efflist)
# #   effnames <- sapply(efflist, function(x) x$effect, simplify = T)
# #   effparams <- sapply(efflist, function(x) x$parameter, simplify = T)
# #   mat <- matrix(rep(0, m1$M * length(efflist) ), nrow=self$M, ncol=length(efflist) )
# #   colnames(mat) <- effnames
# #   rownames(mat) <- 1:self$M
# #   #
# #   for (i in 1:length(efflist))
# #   {
# #     eff_name <- efflist[[ i ]]$effect
# #     #
# #     xActorDegree  <- rowSums(bi_env_mat, na.rm=T)
# #     xComponentDegree  <- colSums(bi_env_mat, na.rm=T)
# #     ## network statistics dataframe
# #     #
# #     if (eff_name == 'density' )
# #     {
# #       stat <- c( xActorDegree )
# #     }
# #     else if (eff_name == 'outAct' )
# #     {
# #       stat <- c( xActorDegree^2 )
# #     }
# #     else if (eff_name == 'inPop' )
# #     {
# #       stat <- c( bi_env_mat %*% (xComponentDegree + 1) )
# #     }
# #     else
# #     {
# #       cat(sprintf('\n\nEffect not yet implemented: `%s\n\n`', eff_name))
# #     }
# #     #
# #     mat[ , i] <- stat
# #     #
# #     df <- rbind(df, data.frame(
# #       statistic = stat, 
# #       actor_id = factor(1:self$M), 
# #       effect_id = factor(i), 
# #       effect_name = effnames[i] 
# #     ))
# #   }
# #   #
# #   if (type %in% c('df','data.frame'))   return(df)
# #   if (type %in% c('mat','matrix'))      return(mat)
# #   if (type %in% c('all','both','list')) return(list(df=df, mat=mat))
# #   cat(sprintf('specified return type %s not found', type))
# # }
# 
# 
# ##**TODO** Check if updating preallocated matrix is faster then rbind to data.frame rows ?
# # get_chain_stats_list <- function() { ## 'all','statdf','bi_env_arr', 'util', 'util_diff')
# #   if (is.null(self$chain_stats))
# #     stop('chain_stats missing; run simulation with returnChains=TRUE before computing actor utility')
# #   # ## actor Utility vector
# #   # au <- c( mat %*% self$rsiena_model$theta )
# #   # hist( util )
# #   effnames <- sapply(self$config_structure_model$dv_bipartite$effects, function(x) x$effect)
# #   ## remove chain entries where no change was made (i.e., keep if !stability )
# #   tiechdf <- self$chain_stats[ !self$chain_stats$stability, ]
# #   
# #   ## get matrix timeseries and network statistics timeseries
# #   nchains <- length(self$rsiena_model$chain)
# #   theta   <- self$rsiena_model$theta
# #   ntheta <- length(theta)
# #   #
# #   bi_env_arr <- array(NA, dim=c(self$M, self$N, nchains))
# #   # bi_env_long <- data.frame()
# #   statdf <- data.frame()
# #   utildf <- data.frame()
# #   util_diff <- data.frame()
# #   #
# #   # utility = util,
# #   # chain_step_id = i, 
# #   # actor_id = factor(1:self$M )
# #   #
# #   nrows_1step_stat <- ntheta * self$M
# #   nrows_1step_util <- self$M
# #   #
# #   statmat_long  <- matrix(NA, nrow= (nchains * nrows_1step_stat), ncol=4) ## chain_step_id, actor_id, stat, value
# #   utilmat_long  <- matrix(NA, nrow= (nchains * nrows_1step_util), ncol=3) ## chain_step_id, actor_id, utility
# #   util_diff_long <- matrix(NA, nrow=(nchains * nrows_1step_util), ncol=3)
# #   #
# #   for (i in 1:nrow(tiechdf)) {
# #     mstep <- tiechdf[i,]
# #     cat(sprintf('\n %.2f%s i=%s, j=%s', 100*i/nrow(tiechdf),'%', mstep$id_from,  mstep$id_to))
# #     ## update bipartite environment matrix for one step (toggle one dyad)
# #     actor_row <- mstep$id_from
# #     component_col <- (mstep$id_to - self$M)
# #     bi_env_mat <- toggleBiMat(bi_env_mat, actor_row,  component_col)
# #     #
# #     # bi_env_mat
# #     #
# #     # Get New Statistics
# #     statsobj <- get_struct_mod_net_stats_list_from_bi_mat( bi_env_mat )
# #     #
# #     statsdf  <- rbind(statsdf,  statsobj$df )
# #     #
# #     # stat_step_fill <- matrix(
# #     #   c(
# #     #     rep(i, nrows_1step_stat), 
# #     #     rep(1:self$M, each=ntheta),
# #     #     rep(1:ntheta, self$M), 
# #     #     rep(NA, nrows_1step_stat)
# #     #   ), 
# #     #   ncol = 4 ##nrow = nrows_1step_stat
# #     # )
# #     stat_step_mat <-  expand.grid(chain_step_id=i, actor_id=1:self$M, stat_eff_id=1:ntheta, value=NA)
# #     util_step_fill <- c()
# #     
# #     #
# #     c(statsobj$mat)
# #     #
# #     statrow_ids <- ( 1 + (i-1)*nrows_1step_stat ):( i*nrows_1step_stat )
# #     statmat_long[ statrow_ids, ] <- matrix(c(), nrow=nrows_1step_stat, ncol=4)
# #     #
# #     statsmat <- statsobj$mat  #get_n_by_m_mat_from_long_df
# #     #
# #     ## Add ministep updated bipartite environment to array
# #     bi_env_arr[ , , i]  <- bi_env_mat
# #     #
# #     # stat_arr[ , , i] <- statsmat
# #     ## Add utilities to array
# #     util <- c( statsmat %*% theta )
# #     utildf <- rbind(utildf, data.frame(
# #       utility = util,
# #       chain_step_id = i, 
# #       actor_id = factor(1:self$M )
# #     ))
# #     util_diff <- rbind(util_diff, data.frame(
# #       utility = if(i == 1){ NA } else {util - util_lag }, ## diff 
# #       chain_step_id = i, 
# #       actor_id = factor(1:self$M )
# #     ))
# #     ######
# #     ## update utility lag for next period difference
# #     util_lag <- util
# #     ######
# #   }
# #   #####
# #   return(list(
# #     utildf = utildf,
# #     util_diff = util_diff,
# #     statsdf = statsdf,
# #     bi_env_arr = bi_env_arr
# #   ))
# # }
# 
# 
# 
# 
# # get_chain_stats_list <- function() { ## 'all','statdf','bi_env_arr', 'util', 'util_diff')
# #   if (is.null(self$chain_stats))
# #     stop('chain_stats missing; run simulation with returnChains=TRUE in order to compute actor utility')
# #   # ## actor Utility vector
# #   # au <- c( mat %*% self$rsiena_model$theta )
# #   # hist( util )
# #   effnames <- sapply(self$config_structure_model$dv_bipartite$effects, function(x) x$effect)
# #   ## remove chain entries where no change was made (keep if !stability )
# #   tiechdf <- self$chain_stats[ !self$chain_stats$stability, ]
# #   
# #   ## get matrix timeseries and network statistics timeseries
# #   nchains <- length(self$rsiena_model$chain)
# #   theta   <- self$rsiena_model$theta
# #   #
# #   bi_env_arr <- array(NA, dim=c(self$M, self$N, nchains))
# #   #
# #   stats_li <- list()
# #   util_li  <- list()
# #   util_diff_li <- list()
# #   # bi_env_long <- data.frame()
# #   # statdf <- data.frame()
# #   # utildf <- data.frame()
# #   # util_diff <- data.frame()
# #   for (i in 1:nrow(tiechdf)) {
# #     mstep <- tiechdf[i,]
# #     cat(sprintf('\n %.2f%s i=%s, j=%s', 100*i/nrow(tiechdf),'%', mstep$id_from,  mstep$id_to))
# #     ## update bipartite environment matrix for one step (toggle one dyad)
# #     bi_env_mat <- toggleBiMat(bi_env_mat, mstep$id_from,  (mstep$id_to - self$M)  )
# #     #
# #     # bi_env_mat
# #     #
# #     # Get New Statistics
# #     # statsobj <- get_struct_mod_net_stats_list_from_bi_mat( bi_env_mat )
# #     statmat <- get_struct_mod_stats_mat_from_bi_mat( bi_env_mat )
# #     #
# #     step_statgrid <- expand.grid(chain_step_id=i, actor_id=1:self$M, effect_id=1:ntheta)
# #     step_statgrid$value <- c( statmat )
# #     #
# #     stats_li[[i]] <- step_statgrid   ##**CHECK list**
# #     #
# #     # statsdf  <- rbind(statsdf,  statsobj$df )
# #     #
# #     # statsmat <- statsobj$mat  #get_n_by_m_mat_from_long_df
# #     #
# #     
# #     #
# #     # tmpstatdf <- as.data.frame( statsmat )
# #     # tmpstatdf$chain_step_id <- i
# #     # tmpstatdf$actor_id <- factor(1:self$M)
# #     # statdf <- rbind(statdf,  tmpstatdf )
# #     ## network statistics matrix added to array
# #     # stat_long <- rbind(stat_long, statdf)
# #     ## Add ministep updated bipartite environment to array
# #     bi_env_arr[ , , i]  <- bi_env_mat
# #     #
# #     # stat_arr[ , , i] <- statsmat
# #     ## Add utilities to array
# #     util <- c( statsmat %*% theta )
# #     # utildf <- rbind(utildf, data.frame(
# #     #   utility = util,
# #     #   chain_step_id = i, 
# #     #   actor_id = factor(1:self$M )
# #     # ))
# #     ##
# #     step_utilgrid <- expand.grid(chain_step_id=i, actor_id=1:self$M)
# #     step_utilgrid$utility <- util
# #     util_li[[i]] <- step_utilgrid  ##**CHECK list**
# #     ##
# #     ##
# #     # util_diff <- rbind(util_diff, data.frame(
# #     #   utility = if(i == 1){ NA } else {util - util_lag }, ## diff 
# #     #   chain_step_id = i, 
# #     #   actor_id = factor(1:self$M )
# #     # ))
# #     ##
# #     step_util_diffgrid <- expand.grid(chain_step_id=i, actor_id=1:self$M)
# #     step_util_diffgrid$utility <- if(i == 1){ NA } else { util - util_lag }
# #     util_diff_li[[i]] <- step_util_diffgrid  ##**CHECK list**
# #     ##
# #     ######
# #     ## update utility lag for next period difference
# #     util_lag <- util
# #     ######
# #   }
# #   #####
# #   return(list(
# #     stats_li = stats_li,
# #     util_li = util_li,
# #     util_diff_li = util_diff_li,
# #     bi_env_arr = bi_env_arr
# #   ))
# # }
# # 
# # chain_stats_li <- self$get_chain_stats_list()
# # 
# # ## Actor Network Statistics long dataframe
# # stats_long <- data.table::rbindlist( chain_stats_li$stats_li )
# # stats_long$actor_id <- as.factor(stats_long$actor_id)
# # stats_long$effect_name <- as.factor(effnames[ stats_long$effect_id ])
# # ## Actor Utility  long dataframe
# # util_long <- data.table::rbindlist( chain_stats_li$util_li ) 
# # util_long$actor_id <- as.factor(util_long$actor_id)
# # ## Actor Utility Difference long dataframe
# # util_diff_long <- data.table::rbindlist( chain_stats_li$util_diff_li ) 
# # util_diff_long$actor_id <- as.factor(util_diff_long$actor_id)
# # 
# # 
# # csl <- get_chain_stats_list()
# 
# 
# # utildf %>% group_by(actor_id) %>% mutate(rollmean= zoo::rollmean(utility, k=10)) %>%
# 
# # avgutil <- utildf %>% group_by(chain_step_id ) %>% summarize(mean=mean(utility)) %>% pull(mean)
# # actutil1 <- utildf %>% filter(actor_id==1) %>% pull(utility)
# # 
# # 
# # ff <- fft(avgutil, inverse = F)
# # freq <- seq(0, 1, length.out = (length(avgutil)/2) + 1)
# # ampl <- Mod(ff[1:((length(avgutil)/2)+1)])
# # plot(freq, ampl, type='l', log='y')
# # # mvfft() ## multivariate fft 
# # 
# # ## Compare 2 actors utilty
# # util_df <- self$actor_stats
# # ggplot(util_df%>%filter(actor_id %in% c(1,2)), aes(x=chain_step_id, y=utility, color=actor_id, shape=actor_id)) + 
# #   geom_point(alpha=.2, size=3) + # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
# #   theme_bw()
# # 
# # 
# # ## Compare 2 actors utilty
# # ggplot(utildf%>%filter(actor_id %in% c(1,2)), aes(x=chain_step_id, y=utility, color=actor_id, shape=actor_id)) + 
# #   geom_point(alpha=.2, size=3) + # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
# #   theme_bw()
# 
# ## plot actors utility points and smoothed loess lines
# thin_factor <- 1
# ggplot(utildf%>%filter(chain_step_id%%thin_factor==0), aes(x=chain_step_id, y=utility, color=actor_id)) + 
#   geom_point(alpha=.2, shape=1, size=2) +
#   geom_smooth(aes(linetype=actor_id), method = 'loess', linewidth=.8, alpha=.15) + 
#   theme_bw()
# ## thin factor = 10 (takes every 10th chain step)
# thin_factor <- 10
# ggplot(utildf%>%filter(chain_step_id%%thin_factor==0), aes(x=chain_step_id, y=utility, color=actor_id)) + 
#   geom_point(alpha=.5, shape=1, size=3) +
#   geom_smooth(aes(linetype=actor_id), method = 'loess', linewidth=.8, alpha=.15) + 
#   theme_bw()
# 
# 
# ## plot actors utility smoothed loess lines (no points)
# thin_factor <- 1
# ggplot(utildf%>%filter(chain_step_id%%thin_factor==0), aes(x=chain_step_id, y=utility, color=actor_id)) + 
#   geom_smooth(aes(linetype=actor_id), method = 'loess', linewidth=1, alpha=.1) + 
#   geom_hline(yintercept =  mean(utildf%>%filter(chain_step_id%%thin_factor==0)%>%pull(utility))) +
#   theme_bw()
# #
# thin_factor <- 10
# ggplot(utildf%>%filter(chain_step_id%%thin_factor==0), aes(x=chain_step_id, y=utility, color=actor_id)) + 
#   geom_smooth(aes(linetype=actor_id), method = 'loess', linewidth=1, alpha=.05) + 
#   geom_hline(yintercept =  mean(utildf%>%filter(chain_step_id%%thin_factor==0)%>%pull(utility))) +
#   theme_bw()
# 
# # ggplot(utildf, aes(x=chain_step_id, y=utility, color=actor_id)) + 
# #   geom_line(alpha=.3, shape=1) + 
# #   theme_bw()
# 
# ## Densities of actor utilities
# ggplot(utildf, aes(x=utility, color=actor_id, fill=actor_id, linetype=actor_id)) + 
#   geom_density(alpha=.1, linewidth=1)  + 
#   # facet_wrap(~chain_step_id)+ 
#   theme_bw()
# 
# ## Actor density fact plots comparing H1 to H2 utility distribution
# utildf$chain_half <- factor(1 + 1*(utildf$chain_step_id >= median(utildf$chain_step_id)), levels = c(2,1)) ## reverse order for linetype_1 used for Half_2
# ggplot(utildf, aes(x=utility, color=actor_id, fill=actor_id, linetype=chain_half)) + 
#   geom_density(alpha=.1, linewidth=1)  + 
#   facet_wrap(~actor_id)+
#   theme_bw()
# 
# ## Actor utility difference timeseries
# ggplot(util_diff, aes(x=chain_step_id, y=utility, color=actor_id)) + 
#   geom_line(alpha=.9, shape=1) +  facet_grid(actor_id ~ .) +
#   geom_hline(yintercept = 0, linetype=2) +
#   theme_bw()
# 
# ## Actor utility difference timeseries
# ggplot(util_diff, aes(x=utility, color=actor_id, fill=actor_id)) + 
#   geom_histogram(alpha=.3, bins=35) + 
#   facet_grid(.~actor_id ) +
#   scale_y_log10() +
#   # geom_hline(yintercept = 0, linetype=2) +
#   theme_bw()
# 
# 
# # ggplot(util_diff, aes(x=chain_step_id, y=utility, color=actor_id, fill=actor_id, linetype=actor_id)) + 
# #   geom_line()  + 
# #   # facet_wrap(~chain_step_id)+ 
# #   theme_bw()
# 
# util_diff_wide <- tidyr::pivot_wider(util_diff, 
#     values_from='utility', 
#     names_from = 'chain_step_id',  
#     id_cols = 'actor_id' 
#   ) %>% 
#   select( !c(actor_id) )
# 
# util_diff_t    <- t(util_diff_wide )
# util_diff_m1_t <- t(util_diff_wide %>% select( !c(`1`) ))  ## drop NA's in first columns representing first chain step
# 
# util_diff_cor <- cor(util_diff_m1_t)
# diag(util_diff_cor) <- 0
# #
# heatmap(util_diff_cor)
# #
# rownames(util_diff_cor) <- paste0('A',1:self$M)
# colnames(util_diff_cor) <- paste0('A',1:self$M)
# #
# cor_long <- as.data.frame(util_diff_cor) %>%
#   tibble::rownames_to_column(var = "Actor_1") %>%
#   pivot_longer(cols = -Actor_1, names_to = "Actor_2", values_to = "correlation")
# 
# ggplot(cor_long, aes(x=Actor_1, y=Actor_2, fill=correlation)) +
#   geom_tile(color = "white") +
#   scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0) +
#   geom_text(aes(label = round(correlation, 2)), size = 4, color='white') +  # Text labels
#   labs(title = "Actor Utility Correlation Heatmap", x = "Actor  i", y = "Actor  j") +
#   theme_bw()
# 
# ## Get dataframe of nchains (rows) by n_summary_stats (columns) for summarying actor utility changes
# util_diff_t_summary <- data.frame( 
#   ministep_count_id = 1:nrow(util_diff_t),
#   sum = rowSums( util_diff_t ),
#   mean = rowMeans( util_diff_t ),
#   sd = apply( util_diff_t, 1, sd ),
#   # skew = apply( util_diff_t, 1, skewness ),
#   # kurt = apply( util_diff_t, 1, kurtosis ),
#   min = apply( util_diff_t, 1, min ),
#   max = apply( util_diff_t, 1, max ),
#   range_magnitude = apply( util_diff_t, 1, function(x) abs(diff(range(x))) ),
#   cnt_change = apply( util_diff_t, 1, function(x) sum(x != 0) ),  ##sum of logical counts the TRUEs
#   cnt_stable = apply( util_diff_t, 1, function(x) sum(x == 0) ),  ##sum of logical counts the TRUEs
#   cnt_pos    = apply( util_diff_t, 1, function(x) sum(x > 0) ),   ##sum of logical counts the TRUEs
#   cnt_neg  = apply( util_diff_t, 1, function(x) sum(x < 0) )   ##sum of logical counts the TRUEs
# )
#  ## convert to long for plotting  
# util_diff_t_summary_long <- util_diff_t_summary %>% 
#   pivot_longer(names_to = 'utility_stat', 
#                values_to = 'utility_stat_value', 
#                cols = c('sum', 'mean', 'sd', 'min', 'max', 'range_magnitude', 
#                         'cnt_change', 'cnt_stable', 'cnt_pos', 'cnt_neg'))
# 
# # ggplot(util_diff_t_summary_long, aes(x=utility_stat_value, color=utility_stat, fill=utility_stat)) + 
# #   geom_density(alpha=.2) + facet_wrap(~utility_stat, scales = 'free') + 
# #   theme_bw()
# ### Histograph facet_wrap of statistics of actor_utilty change from one ministep of any other actor
# ggplot(util_diff_t_summary_long, aes(x=utility_stat_value, color=utility_stat, fill=utility_stat)) + 
#   geom_histogram(alpha=.2) + facet_wrap(~utility_stat, scales = 'free_x') + 
#   geom_vline(xintercept = 0) +
#   theme_bw()
# 
# 
# # statmat_arr_long <- melt(statmat_arr, varnames = effnames, )
# 
# statmat_arr_long <- statmat_arr #%>% rename(actor = Var1, statistic = Var2, period = Var3, value = Freq)
# 
# # Preview the result
# head(long_df_tidyr)
# 
# # melt(wavemat, varnames = c("Simulation", "Wave", "NetworkEffect"), value.name = "Mean")
# 
# 
# # statdf <- ldply(statlist) 
# 
# 
# # # eff <- m1$rsiena_effects[m1$rsiena_effects$include, ]
# # efflist <- m1$config_structure_model$dv_bipartite$effects
# # # mat <- matrix(NA, nrow=m1$M, ncol=length(efflist$efflist))
# # statlist <- list()
# # for (i in 1:length(efflist)) {
# #   #
# #   eff <- efflist[[ i ]]
# #   bi_env_mat <- m1$bipartite_matrix
# #   
# #   xActorDegree  <- rowSums(bi_env_mat, na.rm=T)
# #   xComponentDegree  <- colSums(bi_env_mat, na.rm=T)
# #   ## network statistics dataframe
# #   # xmat <- matrix(NA, )
# #   #
# #   if (eff$effect == 'density') {
# #     statlist[['density']] <- xActorDegree 
# #   } else if (eff$effect == 'outAct'){
# #     statlist[['outAct']]  <- xActorDegree^2 
# #   } else if (eff$effect == 'inPop') {
# #     statlist[['inPop']]  <- c( bi_env_mat %*% (xComponentDegree + 1) )  ## MxN %*% N
# #   } else  {
# #     print(sprintf('Effect not yet implemented: `%s`', eff$effect))
# #   }
# # }
# # statdf <- ldply(statlist) 
# 
# 
# 
# # ## 2. Strategies sets the objective function as a linear combination of network stats across DVs
# # structure_model <- list(
# #   dv_social = list(
# #     dv_name = 'self$social_rsienaDV',
# #     effects = list(
# #       list(effect='density', parameter= 0, fix=T, dv_name='self$social_rsienaDV')#, #interaction1 = NULL
# #       # list(effect='transTriads', parameter=1, fix=F, dv_name='self$social_rsienaDV')#, #interaction1 = NULL
# #       # list(effect='density', parameter=0, fix=T),
# #     )
# #   ),
# #   dv_search = list(
# #     dv_name = 'self$search_rsienaDV',
# #     effects = list(
# #       list(effect='density', parameter= 0, fix=T, dv_name='self$search_rsienaDV')#, #interaction1 = NULL
# #       # list(effect='transTriads', parameter=1, fix=F, dv_name='self$search_rsienaDV')#, #interaction1 = NULL
# #       # list(effect='density', parameter=0, fix=T),
# #     )
# #   ),
# #   dv_bipartite = list(
# #     name = 'self$bipartite_rsienaDV',
# #     effects = list(
# #       list(effect='density', parameter= -2, fix=T, dv_name='self$bipartite_rsienaDV'), ##interaction1 = NULL
# #       list(effect='inPop', parameter= 1, fix=F, dv_name='self$bipartite_rsienaDV'), #interaction1 = NUL
# #       list(effect='outAct', parameter= -1, fix=F, dv_name='self$bipartite_rsienaDV')#, #interaction1 = NULL
# #       # list(effect='outInAss', parameter=0, fix=F, dv_name='self$bipartite_rsienaDV'), #interaction1 = NULL
# #       # list(effect='cycle4', parameter=.5, fix=T, dv_name='self$bipartite_rsienaDV')#, #interaction1 = NULL
# #       # list(effect='density', parameter=0, fix=T),
# #     )
# #   )
# # )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # getChangeContributions <- function(algorithm, data, effects)
# # {
# #   ## Gets the simulated statistics.
# #   ## The following initializations data, effects, and model
# #   ## for calling "getTargets" in "siena07.setup.h"
# #   ## is more or less copied from "getTargets" in "getTargets.r".
# #   ## However, some modifications have been necessary to get it to work.
# #   f <- unpackData(data,algorithm)
# #   
# #   effects <- effects[effects$include,]
# #   if (!is.null(algorithm$settings))
# #   {
# #     stop('not implemented: RI together with settings')
# #     # effects <- addSettingsEffects(effects, algorithm)
# #   }
# #   else
# #   {
# #     effects$setting <- rep("", nrow(effects))
# #   }
# #   pData <- .Call(C_setupData, PACKAGE=pkgname,
# #                  list(as.integer(f$observations)),
# #                  list(f$nodeSets))
# #   ## register a finalizer
# #   ans <- reg.finalizer(pData, clearData, onexit = FALSE)
# #   ans<- .Call(C_OneMode, PACKAGE=pkgname,
# #               pData, list(f$nets))
# #   ans <- .Call(C_Bipartite, PACKAGE=pkgname, # added 1.1-299
# #                pData, list(f$bipartites))
# #   ans<- .Call(C_Behavior, PACKAGE=pkgname, pData,
# #               list(f$behavs))
# #   ans<-.Call(C_ConstantCovariates, PACKAGE=pkgname,
# #              pData, list(f$cCovars))
# #   ans<-.Call(C_ChangingCovariates,PACKAGE=pkgname,
# #              pData,list(f$vCovars))
# #   ans<-.Call(C_DyadicCovariates,PACKAGE=pkgname,
# #              pData,list(f$dycCovars))
# #   ans<-.Call(C_ChangingDyadicCovariates,PACKAGE=pkgname,
# #              pData, list(f$dyvCovars))
# #   
# #   storage.mode(effects$parm) <- 'integer'
# #   storage.mode(effects$group) <- 'integer'
# #   storage.mode(effects$period) <- 'integer'
# #   
# #   effects$effectPtr <- rep(NA, nrow(effects))
# #   depvarnames <- names(data$depvars)
# #   tmpeffects <- split(effects, effects$name)
# #   myeffectsOrder <- match(depvarnames, names(tmpeffects))
# #   ans <- .Call(C_effects, PACKAGE=pkgname, pData, tmpeffects)
# #   pModel <- ans[[1]][[1]]
# #   for (i in 1:length(ans[[2]]))
# #   {
# #     effectPtr <- ans[[2]][[i]]
# #     tmpeffects[[i]]$effectPtr <- effectPtr
# #   }
# #   myeffects <- tmpeffects
# #   for(i in 1:length(myeffectsOrder)){
# #     myeffects[[i]]<-tmpeffects[[myeffectsOrder[i]]]
# #   }
# #   ans <- .Call(C_getTargets, PACKAGE=pkgname, pData, pModel, myeffects,
# #                parallelrun=TRUE, returnActorStatistics=FALSE,
# #                returnStaticChangeContributions=TRUE)
# #   # See getTargets in siena07setup.cpp; also see rTargets in StatisticsSimulation.cpp
# #   ans
# # }
# # 
# # unpackData <- function(data, x)
# # {
# #   f <- NULL
# #   observations<- data$observations
# #   types <- sapply(data$depvars, function(x) attr(x, "type"))
# #   f$nDepvars <- length(data$depvars)
# #   oneModes <- data$depvars[types == "oneMode"]
# #   Behaviors <- data$depvars[types == "behavior"]
# #   continuousBehaviors <- data$depvars[types == "continuous"]
# #   bipartites <- data$depvars[types == "bipartite"]
# #   ## add the settings
# #   # oneModes <- lapply(oneModes, function(depvar) {
# #   #                    name <- attr(depvar, "name")
# #   #                    if (name %in% names(x$settings)) {
# #   #                      # attr(depvar, "settings") <- c("universal", "primary", x$settings[[name]])
# #   #                      attr(depvar, "settings") <- c(x$settings[[name]])
# #   #                    }
# #   #                    depvar
# #   #                    })
# #   f$nets <- lapply(oneModes, function(x, n, comp)
# #     unpackOneMode(x, n, comp),
# #     n = observations, comp=data$compositionChange)
# #   names(f$nets) <- names(oneModes)
# #   f$bipartites <- lapply(bipartites, function(x, n, comp)
# #     unpackBipartite(x, n, comp),
# #     n = observations, comp=data$compositionChange)
# #   names(f$bipartites) <- names(bipartites)
# #   f$behavs <-  lapply(Behaviors, function(x, n) unpackBehavior(x, n),
# #                       n = observations)
# #   names(f$behavs) <- names(Behaviors)
# #   f$contbehavs <- lapply(continuousBehaviors, function(x, n)
# #     unpackBehavior(x, n), n = observations)
# #   names(f$contbehavs) <- names(continuousBehaviors)
# #   f$observations <- observations
# #   f$seed<- vector("list", observations - 1)
# #   f$depvars <- data$depvars
# #   f$nodeSets <- data$nodeSets
# #   f$oneModes <- oneModes
# #   f$Behaviors <- Behaviors
# #   f$continuousBehaviors <- continuousBehaviors
# #   f$oneModeUpOnly <- sapply(oneModes, function(x) attr(x, "uponly"))
# #   f$oneModeDownOnly <- sapply(oneModes, function(x) attr(x, "downonly"))
# #   f$behaviorUpOnly <- sapply(Behaviors, function(x) attr(x, "uponly"))
# #   f$behaviorDownOnly <- sapply(Behaviors, function(x) attr(x,
# #                                                            "downonly"))
# #   f$distances <- sapply(data$depvars, function(x) attr(x, "distance"))
# #   f$cCovars <- data$cCovars
# #   f$vCovars <- data$vCovars
# #   ## dyadic covars need to be edgelists
# #   f$dycCovars <- lapply(data$dycCovars, function(x) unpackCDyad(x))
# #   f$dyvCovars <- lapply(data$dyvCovars, function(x,n) unpackVDyad(x,n),
# #                         n=observations)
# #   ## create the composition change event lists
# #   f$exog <- lapply(data$compositionChange, function(x)
# #     unpackCompositionChange(x))
# #   f
# # }
# # 
# # unpackBipartite <- function(depvar, observations, compositionChange)
# # {
# #   edgeLists <- vector("list", observations)
# #   networks <- vector("list", observations)
# #   actorSet <- attr(depvar, "nodeSet")
# #   compActorSets <- sapply(compositionChange, function(x)attr(x, "nodeSet"))
# #   thisComp <- match(actorSet, compActorSets)
# #   compChange <- any(!is.na(thisComp))
# #   if (compChange)
# #   {
# #     #  stop("Composition change is not yet implemented for bipartite",
# #     #       "networks")
# #     action <- attr(compositionChange[[thisComp]], "action")
# #     ccOption <- attr(compositionChange[[thisComp]], "ccOption")
# #   }
# #   else
# #   {
# #     ccOption <- 0
# #     action <- matrix(0, nrow=attr(depvar, "netdims")[1], ncol=observations)
# #   }
# #   sparse <- attr(depvar, "sparse")
# #   allowOnly <- attr(depvar, "allowOnly")
# #   if (sparse)
# #   {
# #     ## require(Matrix)
# #     ## have a list of sparse matrices in triplet format
# #     ## with missings and structurals embedded and 0 based indices!
# #     netmiss <- vector("list", observations)
# #     for (i in 1:observations)
# #     {
# #       ## extract this matrix
# #       networks[[i]] <- depvar[[i]]
# #       nActors <- nrow(depvar[[i]])
# #       nReceivers <- ncol(depvar[[i]])
# #       ## stop if any duplicates
# #       netmat <- cbind(networks[[i]]@i+1, networks[[i]]@j+1,
# #                       networks[[i]]@x)
# #       if (any(duplicated(netmat[, 1:2])))
# #       {
# #         stop("duplicate entries in sparse matrix")
# #       }
# #       ## extract missing entries
# #       netmiss[[i]] <- netmat[is.na(netmat[, 3]), , drop = FALSE]
# #       ## carry forward missing values if any
# #       if (i == 1) # set missings to zero
# #       {
# #         netmat <- netmat[!is.na(netmat[,3]), ]
# #         networks[[i]] <- spMatrix(nActors, nReceivers, netmat[, 1],
# #                                   netmat[, 2], netmat[,3])
# #       }
# #       else
# #       {
# #         netmiss1 <- netmiss[[i]][, 1:2]
# #         storage.mode(netmiss1) <- "integer"
# #         networks[[i]][netmiss1[, 1:2]] <-
# #           networks[[i-1]][netmiss1[, 1:2]]
# #       }
# #     }
# #     for (i in 1:observations)
# #     {
# #       mat1 <- networks[[i]]
# #       mat1 <- cbind(mat1@i + 1, mat1@j + 1, mat1@x)
# #       ##missing edgelist
# #       mat2 <- netmiss[[i]]
# #       mat2[, 3] <- 1
# #       ## rows of mat1 with structural values
# #       struct <- mat1[, 3] %in% c(10, 11)
# #       ## reset real data
# #       mat1[struct, 3] <- mat1[struct, 3] - 10
# #       ## copy reset data to structural edgelist
# #       mat3 <- mat1[struct, , drop = FALSE]
# #       ## now remove the zeros from reset data
# #       mat1 <- mat1[!mat1[, 3] == 0, ]
# #       ## do comp change
# #       if (compChange)
# #       {
# #         ## revert to sparse matrices temporarily
# #         mat1 <- spMatrix(nrow=nActors, ncol=nReceivers, i = mat1[, 1],
# #                          j=mat1[, 2], x=mat1[, 3])
# #         mat2 <- spMatrix(nrow=nActors, ncol=nReceivers, i = mat2[, 1],
# #                          j=mat2[, 2], x=mat2[, 3])
# #         mat3 <- spMatrix(nrow=nActors, ncol=nReceivers, i = mat3[, 1],
# #                          j=mat3[, 2], x=mat3[, 3])
# #         ones <- which(action[, i] == 1)
# #         twos <- which(action[, i] == 2)
# #         threes <- which(action[, i] == 3)
# #         for (j in ones) ## False data is not preceded by anything real
# #         {
# #           if (ccOption %in% c(1, 2))
# #           {
# #             ## find missing values for this actor
# #             use <- mat2[j, ] > 0
# #             ## remove from real data (i.e. zero)
# #             mat1[j, use] <- 0
# #             ## remove from missing data
# #             mat2[j, use] <- 0
# #             ## remove from raw data for distances later
# #             depvar[[i]][j, use] <- 0 ## zero
# #           }
# #           else if (ccOption == 3)
# #           {
# #             ## add the row  to the missing data
# #             mat2[j, ] <- 1
# #             ## set to missing in raw data for distances later
# #             depvar[[i]][j, ] <- NA
# #           }
# #         }
# #         for (j in threes) ## False data is preceded and followed by real
# #         {
# #           if (ccOption %in% c(1, 2))
# #           {
# #             ## find missing values for this actor
# #             use <- mat2[j, ] > 0
# #             ## remove these from mat2, the missing data
# #             mat2[j, use] <- 0
# #             ## carry forward
# #             if (i == 1)
# #             {
# #               ## 0 any matches from mat1, the real data
# #               mat1[j, use] <- 0
# #             }
# #             else
# #             {
# #               mat1[j, use] <- networks[[i-1]][j, use]
# #             }
# #             depvar[[i]][j, use] <- 0 ##  not missing
# #           }
# #           else if (ccOption == 3)
# #           {
# #             ## add the row to the missing data
# #             mat2[j, ] <- 1
# #             depvar[[i]][j, ] <- NA
# #           }
# #         }
# #         for (j in twos) ## False data is not followed by anything real
# #         {
# #           if (ccOption == 1)
# #           {
# #             ## find missing values for this actor
# #             use <- mat2[j, ] > 0
# #             ## remove these from mat2, the missing data
# #             mat2[j, use] <- 0
# #             depvar[[i]][j, use] <- 0 ##  not missing
# #             ## carry forward
# #             if (i == 1)
# #             {
# #               ## 0 any matches from mat1, the real data
# #               mat1[j, use] <- 0
# #             }
# #             else
# #             {
# #               mat1[j, use] <- networks[[i-1]][j , use]
# #             }
# #           }
# #           else if (ccOption %in% c(2, 3))
# #           {
# #             ## add the row  to the missing data
# #             mat2[j, ] <- 1
# #             depvar[[i]][j, ] <- NA
# #           }
# #         }
# #         
# #         ## now revert to triplet matrices, after updating networks
# #         networks[[i]] <- mat1
# #         mat1 <- cbind(mat1@i + 1, mat1@j + 1, mat1@x)
# #         mat2 <- cbind(mat2@i + 1, mat2@j + 1, mat2@x)
# #         mat3 <- cbind(mat3@i + 1, mat3@j + 1, mat3@x)
# #         if (any (mat1[, 3] == 0) || any (mat2[, 3] == 0) ||
# #             any (mat3[, 3] == 0))
# #         {
# #           stop("zero values in sparse matrices")
# #         }
# #         if (any (duplicated(mat1[, -3])) ||
# #             any (duplicated(mat2[, -3])) ||
# #             any (duplicated(mat3[, -3])))
# #         {
# #           stop("duplicate values in sparse matrices")
# #         }
# #       }
# #       ##fix up storage mode to be integer
# #       storage.mode(mat1) <- "integer"
# #       storage.mode(mat2) <- "integer"
# #       storage.mode(mat3) <- "integer"
# #       ## add attribute of size
# #       attr(mat1,"nActors") <- c(nActors, nReceivers)
# #       attr(mat2,"nActors") <- c(nActors, nReceivers)
# #       attr(mat3,"nActors") <- c(nActors, nReceivers)
# #       if (i < observations)
# #       {
# #         ## recreate the distance etc
# #         mymat1 <- depvar[[i]]
# #         mymat2 <- depvar[[i + 1]]
# #         ##remove structural values
# #         x1 <- mymat1@x
# #         x2 <- mymat2@x
# #         x1[x1 %in% c(10, 11)] <- NA
# #         x2[x2 %in% c(10, 11)] <- NA
# #         mymat1@x <- x1
# #         mymat2@x <- x2
# #         mydiff <- mymat2 - mymat1
# #         attr(depvar, "distance")[i] <- sum(mydiff != 0,
# #                                            na.rm = TRUE)
# #         if (allowOnly)
# #         {
# #           if (all(mydiff@x >= 0, na.rm=TRUE))
# #           {
# #             attr(depvar, "uponly")[i] <- TRUE
# #           }
# #           if (all(mydiff@x <= 0, na.rm=TRUE))
# #           {
# #             attr(depvar, "downonly")[i] <- TRUE
# #           }
# #         }
# #       }
# #       edgeLists[[i]] <- list(mat1 = t(mat1), mat2 = t(mat2),
# #                              mat3 = t(mat3))
# #     }
# #   }
# #   else
# #   {
# #     for (i in 1:observations) ## carry missings forward  if exist
# #     {
# #       networks[[i]] <- depvar[, , i]
# #       if (i == 1)
# #         networks[[i]][is.na(depvar[, , i])] <-0
# #       else ##carry missing forward!
# #         networks[[i]][is.na(depvar[, , i])] <-
# #           networks[[i-1]][is.na(depvar[, , i])]
# #     }
# #     for (i in 1:observations)
# #     {
# #       ones <- which(action[, i] == 1)
# #       twos <- which(action[, i] == 2)
# #       threes <- which(action[, i] == 3)
# #       for (j in ones) ## False data is not preceded by anything real
# #       {
# #         if (ccOption %in% c(1, 2))
# #         {
# #           use <- is.na(depvar[j, , i])
# #           depvar[j, use, i] <- 0 ## not missing
# #           networks[[i]][j, use] <- 0 ## zero
# #         }
# #         else if (ccOption == 3)
# #         {
# #           depvar[j, , i] <- NA ## missing
# #         }
# #       }
# #       for (j in threes) ## False data is preceded and followed by real
# #       {
# #         
# #         if (ccOption %in% c(1, 2))
# #         {
# #           use <- is.na(depvar[j, , i])
# #           depvar[j, use, i] <- 0 ##  not missing
# #           ## carry forward already done
# #           if (i == 1)
# #           {
# #             networks[[i]][j, use] <- 0
# #           }
# #           else
# #           {
# #             networks[[i]][j, use] <- networks[[i-1]][j, use]
# #           }
# #         }
# #         else if (ccOption == 3)
# #         {
# #           depvar[j, , i] <- NA ## missing
# #         }
# #       }
# #       for (j in twos) ## False data is not followed by anything real
# #       {
# #         if (ccOption == 1)
# #         {
# #           use <- is.na(depvar[j, , i])
# #           depvar[j, use, i] <- 0 ##  not missing
# #           ## carry forward already done
# #           if (i == 1)
# #           {
# #             networks[[i]][j, use] <- 0
# #           }
# #           else
# #           {
# #             networks[[i]][j, use] <- networks[[i-1]][j, use]
# #             
# #           }
# #         }
# #         else if (ccOption %in% c(2, 3))
# #         {
# #           depvar[j, , i] <- NA ## missing
# #         }
# #       }
# #     }
# #     for (i in 1:observations)
# #     {
# #       if (i < observations)
# #       {
# #         ## recreate distances, as we have none in c++. (no longer true)
# #         mymat1 <- depvar[,,i, drop=FALSE]
# #         mymat2 <- depvar[,,i + 1,drop=FALSE]
# #         ##remove structural values
# #         mymat1[mymat1 %in% c(10,11)] <- NA
# #         mymat2[mymat2 %in% c(10,11)] <- NA
# #         mydiff <- mymat2 - mymat1
# #         attr(depvar, "distance")[i] <- sum(mydiff != 0,
# #                                            na.rm = TRUE)
# #         if (allowOnly)
# #         {
# #           if (all(mydiff >= 0, na.rm=TRUE))
# #           {
# #             attr(depvar, "uponly")[i] <- TRUE
# #           }
# #           if (all(mydiff <= 0, na.rm=TRUE))
# #           {
# #             attr(depvar, "downonly")[i] <- TRUE
# #           }
# #         }
# #       }
# #       
# #       edgeLists[[i]] <- createEdgeLists(networks[[i]], depvar[, , i], TRUE)
# #     }
# #   }
# #   ## add attribute of nodeset
# #   attr(edgeLists, "nodeSet") <- attr(depvar, "nodeSet")
# #   ## add attribute of name
# #   attr(edgeLists, "name") <- attr(depvar, "name")
# #   ## add attribute of distance
# #   attr(edgeLists, "distance") <- attr(depvar, "distance")
# #   ## attr uponly and downonly
# #   attr(edgeLists, "uponly") <- attr(depvar, "uponly")
# #   attr(edgeLists, "downonly") <- attr(depvar, "downonly")
# #   ## attr symmetric
# #   attr(edgeLists, "symmetric") <- attr(depvar, "symmetric")
# #   ## attr balmean
# #   attr(edgeLists, "balmean") <- attr(depvar, "balmean")
# #   ## attr structmean
# #   attr(edgeLists, "structmean") <- attr(depvar, "structmean")
# #   attr(edgeLists, "averageOutDegree") <- attr(depvar, "averageOutDegree")
# #   return(edgeLists = edgeLists)
# # }
# # 
# # createEdgeLists<- function(mat, matorig, bipartite)
# # {
# #   ## mat1 is basic values, with missings and structurals replaced
# #   tmp <- lapply(1 : nrow(mat), function(x, y)
# #   {
# #     mymat <- matrix(0, nrow = sum(y[x, ] > 0), ncol = 3)
# #     mymat[, 1] <- x
# #     mymat[, 2] <- which(y[x, ] != 0)
# #     mymat[, 3] <- y[x, mymat[, 2]]
# #     mymat
# #   }, y = mat)
# #   mat1 <- do.call(rbind, tmp)
# #   ## mat2 reverts to matorig to get the missing values
# #   tmp <- lapply(1 : nrow(matorig), function(x, y)
# #   {
# #     mymat <- matrix(0, nrow = sum(is.na(y[x, ])), ncol = 3)
# #     mymat[, 1] <- x
# #     mymat[, 2] <- which(is.na(y[x, ]))
# #     mymat[, 3] <- 1
# #     mymat
# #   }, y = matorig)
# #   mat2 <- do.call(rbind, tmp)
# #   ## remove the diagonal if not bipartite
# #   if (!bipartite)
# #   {
# #     mat2 <- mat2[mat2[, 1] != mat2[, 2], , drop=FALSE]
# #   }
# #   ## mat3 structurals
# #   struct <- mat1[,3] %in% c(10, 11)
# #   mat1[struct, 3] <- mat1[struct,3] - 10
# #   mat3 <- mat1[struct, , drop=FALSE]
# #   mat3[, 3] <- 1
# #   mat1 <- mat1[!mat1[,3] == 0, , drop=FALSE] ##remove any zeros just created
# #   ##fix up storage mode to be integer
# #   storage.mode(mat1) <- "integer"
# #   storage.mode(mat2) <- "integer"
# #   storage.mode(mat3) <- "integer"
# #   ## add attribute of size
# #   if (bipartite)
# #   {
# #     attr(mat1, "nActors") <- c(nrow(mat), ncol(mat))
# #     attr(mat2, "nActors") <- c(nrow(mat), ncol(mat))
# #     attr(mat3, "nActors") <- c(nrow(mat), ncol(mat))
# #   }
# #   else
# #   {
# #     attr(mat1, "nActors") <- nrow(mat)
# #     attr(mat2, "nActors") <- nrow(mat)
# #     attr(mat3, "nActors") <- nrow(mat)
# #   }
# #   
# #   list(mat1 = t(mat1), mat2 = t(mat2), mat3 = t(mat3))
# # }
# # 
# # ##@createCovarEdgeLists siena07 Reformat data for C++
# # createCovarEdgeList<- function(mat, matorig)
# # {
# #   tmp <- lapply(1 : nrow(mat), function(x, y)
# #   {
# #     mymat <- matrix(0, nrow = sum(y[x, ] != 0), ncol = 3)
# #     mymat[, 1] <- x
# #     mymat[, 2] <- which(y[x, ] != 0)
# #     mymat[, 3] <- y[x, mymat[, 2]]
# #     mymat
# #   }, y = mat)
# #   mat1 <- do.call(rbind, tmp)
# #   ##mat2 reverts to matorig to get the missing values
# #   tmp <- lapply(1 : nrow(matorig), function(x, y)
# #   {
# #     mymat <- matrix(0, nrow = sum(is.na(y[x, ])), ncol = 3)
# #     mymat[, 1] <- x
# #     mymat[, 2] <- which(is.na(y[x, ]))
# #     mymat[, 3] <- 1
# #     mymat
# #   }, y = matorig)
# #   mat2 <- do.call(rbind, tmp)
# #   ## add attribute of size
# #   attr(mat1, "nActors1") <- nrow(mat)
# #   attr(mat1, "nActors2") <- ncol(mat)
# #   list(mat1=t(mat1), mat2=t(mat2))
# # }
# 
# 
# 
# #
# # getChangeContributions(m1$rsiena_algorithm, m1$rsiena_data, m1$rsiena_effects)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # ##**TODO** Action count data frames (explore, exploit) 
# # ##          rows (actors) by columns (iterations), counts: 0,1,2,3,...
# # 
# # rolling_window = 20
# # 
# # self <- m1
# # 
# # m1$chain_stats %>% group_by(iteration_id) %>% count() %>% plot(type='l')
# # 
# # m1$chain_stats %>% group_by(X__dv_name, X__id_from) %>% count() %>% print(n=60)
# # 
# # util_after <- m1$chain_stats$X__utility_after
# # plot(util_after, type='l')
# # plot(cummean(util_after), type='l')
# # 
# # 
# # ## Get Actor-period count of ministep decisions
# # step_summary_df <- m1$chain_stats %>%
# #   filter( X__dv_name=='self$bipartite_rsienaDV') %>% 
# #   group_by(X__id_from)
# # iters_with_change <- unique(step_summary_df$iteration_id)
# # cnt_mat <- matrix(0, nrow=self$M, ncol=length(self$rsiena_model$chain) )
# # for (iter in iters_with_change) {
# #   cat(sprintf(' %s ', iter))
# #   iter_df <- step_summary_df %>% filter(iteration_id == iter)
# #   actor_cnts <- iter_df %>% group_by(X__id_from) %>% count()
# #   id_not_0s <- ( ! actor_cnts$X__id_from %in% c(0,'0') )
# #   if ( any( ! id_not_0s ) ) {
# #     print('zero in id_from ?:')
# #     print(id_not_0s)
# #     actor_cnts <- actor_cnts %>% filter( ! X__id_from %in% c(0,'0') )
# #   }
# #   cnt_mat[ actor_cnts$X__id_from, iter ] <- actor_cnts$n
# # }
# # # View(.)
# # 
# # cnt_df <- as.data.frame(cnt_mat)
# # colnames(cnt_df) <- paste0("Period_", 1:length(self$rsiena_model$chain) )
# # cnt_df$Actor <- 1:self$M
# # rownames(cnt_df) <- paste0("Actor_", cnt_df$Actor)
# # 
# # long_data <- cnt_df %>% 
# #     pivot_longer(cols = starts_with("Period_"),
# #                  names_to = "period",
# #                  values_to = "count") %>% 
# #   mutate(
# #     period = as.numeric(gsub("Period_", "", period)),
# #     actor = as.numeric(gsub("Actor_", "", Actor))
# #   ) %>%
# #   select(actor, period, count)  # Reorder columns
# # 
# # ## plot actor event sequences
# # ggplot(long_data, aes(y = factor(actor), x=period, fill = count)) + 
# #   geom_tile() + 
# #   scale_fill_gradient(low = "white", high = "blue", name = "Ministeps Count") +
# #   theme_bw() + labs('y' = 'Actor')
# # 
# # 
# # # ggplot(long_data, aes(x=period, y=count, color=factor(actor))) + geom_point() + geom_smooth()
# # 
# # 
# # 
# # 
# # ## Get timeseries of actor utility function value
# # util_step_df <- m1$chain_stats %>% filter( X__dv_name=='self$bipartite_rsienaDV')
# # n_util_steps <- nrow(util_step_df)
# # util_mat <- matrix(NA, nrow=self$M, ncol=length(self$rsiena_model$chain) )
# # for (i in 1:n_util_steps) {
# #   cat(sprintf(' %s ', i))
# #   x <- util_step_df[i, ]
# #   util_mat[ x$X__id_from , x$iteration_id ] <- x$X__utility_after
# #   
# # }
# # # Function to fill forward values along rows
# # fill_forward_func <- function(row) {
# #   for (i in 2:length(row)) {
# #     if (is.na(row[i])) {
# #       row[i] <- row[i - 1]
# #     }
# #   }
# #   return(row)
# # }
# # # Apply the fill-forward function to each row
# # util_mat_filled <- t(apply(util_mat, 1, fill_forward_func))
# # 
# # util_df_fill <- as.data.frame(util_mat_filled)
# # colnames(util_df_fill) <- paste0("Period_", 1:length(self$rsiena_model$chain) )
# # util_df_fill$Actor <- 1:self$M
# # rownames(util_df_fill) <- paste0("Actor_", util_df_fill$Actor)
# # 
# # long_util <- util_df_fill %>% 
# #   pivot_longer(cols = starts_with("Period_"),
# #                names_to = "period",
# #                values_to = "count") %>% 
# #   mutate(
# #     period = as.numeric(gsub("Period_", "", period)),
# #     actor = as.numeric(gsub("Actor_", "", Actor))
# #   ) %>%
# #   select(actor, period, count)  # Reorder columns
# # 
# # ggplot(long_util, aes(x=period, y=count, color=factor(actor))) + 
# #   geom_point() + geom_smooth() + theme_bw()
# # 
# # ggplot(long_util, aes(x=period, y=count, color=factor(actor))) +  # geom_point(shape=1, alpha=.5) + 
# #   geom_line(linewidth=.8) + # geom_smooth() + 
# #   theme_bw() + labs(y='Utility', legend='Actor')
# # 
# # # compute utility moving average
# # long_util_ma <- long_util %>%
# #   group_by(actor) %>%
# #   mutate(rolling_avg = zoo::rollmean(count, k = rolling_window, fill = NA, align = "right"))
# # #
# # util_avg <- long_util_ma %>% group_by(period) %>% summarize(mean=mean(rolling_avg, na.rm=T))
# # # plot utility moving average
# # ggplot(long_util_ma, aes(x=period, y=rolling_avg, color=factor(actor))) +  
# #   geom_point(shape=1, alpha=.1) + 
# #   geom_line(linewidth=.9, alpha=.85) + # geom_smooth() + 
# #   geom_line(data = util_avg, linewidth=2, color='black', aes(x=period, y=mean, color='Mean')) + 
# #   theme_bw() + labs(y='Utility', color='Actor')
# # 
# # 
# # ##
# # plot(colMeans(util_mat_filled, na.rm = T), type='l')
# # 
# # 
# # # %>%  mutate(Period = as.numeric(gsub("Period_", "", Period)))  
# # 
# # ## The chain has the structure chain[[run]][[depvar]][[period]][[ministep]].
# # ##ch[[5]] = '[[depvar]][[period]][[ministep]]'
# # 
# # ch <- m1$rsiena_model$chain
# # id_ministep <- 5
# # run <- ch[[id_ministep]]
# # depvar <- run[[1]]
# # period <- depvar[[1]]
# # ministep <- period[[1]]
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # # 1. INIT SIM
# # m2 <- SaomNkRSienaBiEnv$new(environ_params, rand_seed=54321)
# # # 2. RUN SIM
# # m2$search_rsiena_run(structure_model, payoff_formulas, iterations=1000, 
# #                            returnDeps = F,  n_snapshots =1, 
# #                            rsiena_phase2_nsub=1, rsiena_n2start_scale = 1, 
# #                            get_eff_doc = FALSE)
# # # 3. PLOT SIM (EXPLORE RESULTS)
# # m2$search_rsiena_plot_stability()
# # 
# # 
# # 
# # binet <- network(m1$bipartite_matrix, directed = F)
# # ergstats <- ergmMPLE(binet ~ edges + cycle(4), output = 'array')
# # bi_edges <- ergpreds[,, 'edges']
# # bi_cyc4  <- ergpreds[,, 'cycle4']
# # 
# # 
# # 
# # 
# # 
# #  
# # id_wave <- 1
# # ### ## m1$rsiena_model$sf2[ <sims> , id_wave, <effects> ]
# # 
# # m1$rsiena_model$sf2
# # 
# # mean <-  m1$rsiena_model$sf2 
# # rownames(wavemat) <- as.character( 1:nrow(wavemat) )
# # colnames(wavemat) <- gsub('\\s', '_' , m1$rsiena_model$effects$effectName )
# # # dat[,,3] <- dat[1,1,][ m1$rsiena_model$effects$effectName ]
# # 
# # df_long <- melt(wavemat, varnames = c("Simulation", "Wave", "NetworkEffect"), value.name = "Mean")
# # 
# # # Plot using ggplot2
# # ggplot(df_long, aes(x = Mean, color=NetworkEffect, fill=NetworkEffect)) +
# #   # geom_boxplot(aes(group = Wave), outlier.shape = NA) +
# #   geom_density(alpha=.1, bins=11) +
# #   facet_wrap(~ NetworkEffect, scales = "free_x") +
# #   theme_minimal() +
# #   labs(title = "Distribution of Simulations by Effects",
# #        x = "X",
# #        y = "Y")
# # 
# # 
# # 
# # ##$sf2  -> [sims,  waves,  effects]
# # ##$sf   -> [sims, effects]
# # 
# # apply(m1$rsiena_model$sf2[,1,], 2, mean)
# # 
# # colMeans(m1$rsiena_model$sf)
# # 
# # 
# # 
# # ##------------------------------------
# # 
# # # m1 <-
# # 
# # # We can also format the predictor matrix into an array:
# # mplearray <- ergmMPLE(formula, output="array")
# # 
# # #
# # dvnet <- 
# # formula <- as.formula(sprintf('%s ~ edges',dvnet))
# # 
# # # The resulting matrices are big, so only print the first 8 actors:
# # mplearray$response[1:8,1:8]
# # mplearray$predictor[1:8,1:8, ]
# # mplearray$weights[1:8,1:8]
# # 
# # ##------------------------------------
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ###--------------------------------------------------------------------------
# # 
# # m1$search_rsiena_run(objective_list, iterations=2000, returnDeps = T, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
# # m1$search_rsiena_plot_stability()
# # 
# # 
# # # ##
# # # m1$search_rsiena_run(iterations=4000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# # # m1$search_rsiena_plot_stability()
# # # m1$search_rsiena_run(iterations=400, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# # # m1$search_rsiena_plot_stability()
# # 
# # 
# # 
# # 
# # m1 <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .11, sim_name = '_TESTrsiena_')
# # m1$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# # m1$search_rsiena_plot_stability()
# # #
# # m2 <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .11, sim_name = '_TESTrsiena_')
# # m2$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
# # m2$search_rsiena_plot_stability()
# # 
# # #
# # ms1 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .11, sim_name = '_TESTrsiena_')
# # ms1$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# # ms1$search_rsiena_plot_stability()
# # #
# # ms2 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .11, sim_name = '_TESTrsiena_')
# # ms2$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
# # ms2$search_rsiena_plot_stability()
# # 
# # 
# # 
# # md1 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .2, sim_name = '_TESTrsiena_')
# # md1$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# # md1$search_rsiena_plot_stability()
# # #
# # md2 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .2, sim_name = '_TESTrsiena_')
# # md2$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
# # md2$search_rsiena_plot_stability()
# # 
# # # #
# # # mc1 <- SAOM_NK_RSiena$new(M = 12, N = 30, BI_PROB = .05, sim_name = '_TESTrsiena_')
# # # mc1$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# # # mc1$search_rsiena_plot_stability()
# # # #
# # # mc2 <- SAOM_NK_RSiena$new(M = 12, N = 30, BI_PROB = .05, sim_name = '_TESTrsiena_')
# # # mc2$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
# # # mc2$search_rsiena_plot_stability()
# # 
# # # #
# # # ms2 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .15, sim_name = '_TESTrsiena_')
# # # ms2$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = .1)
# # # ms2$search_rsiena_plot_stability()
# # 
# # 
# # #
# # # m3 <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .15, sim_name = '_TESTrsiena_')
# # # m3$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=3, rsiena_n2start_scale = 1)
# # # m3$search_rsiena_plot_stability()
# # #
# # # m4 <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .15, sim_name = '_TESTrsiena_')
# # # m4$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=4, rsiena_n2start_scale = 1)
# # # m4$search_rsiena_plot_stability()
# # 
# # 
# # saomnkrsiena <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .15, sim_name = '_TESTrsiena_')
# # ##
# # # saomnkrsiena$search_rsiena_run(iterations=10, rsiena_phase2_nsub=1) 
# # # saomnkrsiena$search_rsiena_run(iterations=100, rsiena_phase2_nsub=2)
# # # saomnkrsiena$search_rsiena_run(iterations=5000, rsiena_phase2_nsub=3)
# # saomnkrsiena$search_rsiena_run(iterations=1000, n_snapshots =2, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# # saomnkrsiena$search_rsiena_plot_stability()
# # # saomnkrsiena$search_rsiena_run(iterations=1000, rsiena_phase2_nsub=5)
# # 
# # # saom_nk_enhanced$search_rsiena(1)
# # saomnkrsiena$search_rsiena_run(iterations=100) 
# # 
# # saomnkrsiena$search_rsiena_run(iterations=10000)
# # saomnkrsiena$search_rsiena_run(iterations=100000)
# # # saomnkrsiena$search_rsiena_run(iterations=20, overwrite = F) 
# # 
# # saomnkrsiena2 <- SAOM_NK_RSiena$new(M = 15, N = 24, BI_PROB = .15, sim_name = '_TESTrsiena_')
# # # saom_nk_enhanced$search_rsiena(1)
# # saomnkrsiena2$search_rsiena_run(iterations=1000) 
# # 
# # 
# # 
# # 
# # ########################################
# # self <- saomnkrsiena$clone(deep=T)
# # sims = self$rsiena_model$sims
# # n <- length(sims)
# # outlist <- list()
# # difflist <- list()
# # jaccardlist <- list()
# # K_env_list <- list()
# # K_soc_list <- list()
# # for(i in 1:length(sims)) {
# #   cat(sprintf(' %s ', i))
# #   ##
# #   el_bi_env <- sims[[ i ]][[1]][[3]]$`1`
# #   ## update numbering of second mode (the comonent integer names shift upward by the number of actors)
# #   el_bi_env[,2] <- el_bi_env[,2] + self$M
# #   ##
# #   MplusN <- self$M + self$N
# #   #############
# #   ## Bipartite matrix space (N+M by N+M)
# #   ## Undirected --> Upper right rectangle of full bipartite matrix
# #   ##   M N
# #   ## M[0,X], for X in [0,1]
# #   ## N[0,0]
# #   bi_env_mat_sp <- sparseMatrix(i = el_bi_env[,1],
# #                                 j = el_bi_env[,2],
# #                                 x = el_bi_env[,3],
# #                                 dims = c(MplusN, MplusN))
# #   ## get the top-right triangle for the bipartite environment
# #   bi_env_mat_new  <- as.matrix(bi_env_mat_sp)[ 1:self$M, (self$M+1):(MplusN) ]
# #   ## Add ACTOR NAMES on rows
# #   rownames(bi_env_mat_new) <- as.character( 1:self$M )
# #   ## Add COMPONENT NAMES on columns (N+1 ... N+M)
# #   colnames(bi_env_mat_new) <- as.character( (1:self$N) + self$M  )
# #   ##
# #   # self$plot_bipartite_system_from_mat(bi_env_mat_new, i, plot_save = TRUE)
# #   ##
# #   outlist[[ sprintf('sim%d',i) ]] <- bi_env_mat_new
# #   
# #   if (i == 1) {
# #     difflist[[ sprintf('diff%d-%d', i-1, i) ]] <- 0
# #     jaccardlist[[sprintf('jac%d-%d', i-1, i)]] <- 0
# #   } else {
# #     difflist[[ sprintf('diff%d-%d', i-1, i) ]] <-  bi_env_mat_new - outlist[[ (i-1) ]]
# #     jaccardlist[[sprintf('jac%d-%d', i-1, i)]] <- get_jaccard_index(m0 = outlist[[ (i-1) ]], m1 = bi_env_mat_new )
# #   }
# #   
# #   new_bi_g <- igraph::graph_from_biadjacency_matrix(bi_env_mat_new, 
# #                                                       directed = F, mode = 'all', 
# #                                                       multiple = T, weighted = T, 
# #                                                       add.names = T)
# #   projections <- igraph::bipartite_projection(new_bi_g, multiplicity = T, which = 'both')
# #   new_g_soc <- projections$proj1
# #   new_g_env <- projections$proj2
# # 
# #   K_soc_list[[i]] <- igraph::degree(new_g_soc)
# #   K_env_list[[i]] <- igraph::degree(new_g_env)
# #   
# #   
# #   # print(sim_i[[1]][[3]]`1`)
# # }
# # 
# # saomnkrsiena$plot_degree_progress_from_sims(K_soc_list = K_soc_list, K_env_list = K_env_list, plot_save = T)
# # 
# # 
# # # par(mfrow=c(1,3))
# # # stability_vec <- cumsum(plyr::ldply(jaccardlist)[,2])/(1:n)
# # # plot(stability_vec , type='l' , main='Change in Stability  [t-1, t]')
# # # 
# # # 
# # degree_summary <- do.call(rbind, lapply(1:length(K_soc_list), function(iter) {
# #   data.frame(
# #     Iteration = iter,
# #     Mean_K_S = mean(K_soc_list[[iter]]),
# #     Q25_K_S = quantile(K_soc_list[[iter]], 0.25),
# #     Q75_K_S = quantile(K_soc_list[[iter]], 0.75),
# #     Mean_K_E = mean(K_env_list[[iter]]),
# #     Q25_K_E = quantile(K_env_list[[iter]], 0.25),
# #     Q75_K_E = quantile(K_env_list[[iter]], 0.75)
# #   )
# # }))
# # (plt <- ggplot(degree_summary, aes(x = Iteration)) +
# #     geom_line(aes(y = Mean_K_S, color = "Mean K_S")) +
# #     geom_ribbon(aes(ymin = Q25_K_S, ymax = Q75_K_S, fill = "K_S"), alpha = 0.1) +
# #     geom_line(aes(y = Mean_K_E, color = "Mean K_E")) +
# #     geom_ribbon(aes(ymin = Q25_K_E, ymax = Q75_K_E, fill = "K_E"), alpha = 0.1) +
# #     scale_color_manual(values = c("Mean K_S" = "blue", "Mean K_E" = "red")) +
# #     scale_fill_manual(values = c("K_S" = "blue", "K_E" = "red")) +
# #     labs(title = "Degree Progress Over Iterations",
# #          x = "Iteration",
# #          y = "Degree",
# #          color = "Mean Degree",
# #          fill = "IQR (Mid-50%)") +
# #     theme_minimal()
# # )
# # 
# # 
# # 
# # ########################################
# # 
# # # sims = self$rsiena_model$sims
# # # n <- length(sims)
# # # outlist <- list()
# # # difflist <- list()
# # # jaccardlist <- list()
# # # for(i in 1:length(sims)) {
# # #   cat(sprintf(' %s ', i))
# # #   ##
# # #   el_bi_env <- sims[[ i ]][[1]][[3]]$`1`
# # #   ## update numbering of second mode (the comonent integer names shift upward by the number of actors)
# # #   el_bi_env[,2] <- el_bi_env[,2] + self$M
# # #   ##
# # #   MplusN <- self$M + self$N
# # #   #############
# # #   ## Bipartite matrix space (N+M by N+M)
# # #   ## Undirected --> Upper right rectangle of full bipartite matrix
# # #   ##   M N
# # #   ## M[0,X], for X in [0,1]
# # #   ## N[0,0]
# # #   bi_env_mat_sp <- sparseMatrix(i = el_bi_env[,1],
# # #                                 j = el_bi_env[,2],
# # #                                 x = el_bi_env[,3],
# # #                                 dims = c(MplusN, MplusN))
# # #   ##
# # #   bi_env_mat_new  <- as.matrix(bi_env_mat_sp)[ 1:self$M, (self$M+1):(MplusN) ]
# # #   ##
# # #   # self$plot_bipartite_system_from_mat(bi_env_mat_new, i, plot_save = TRUE)
# # #   ##
# # #   outlist[[ sprintf('sim%d',i) ]] <- bi_env_mat_new
# # # 
# # #   if (i == 1) {
# # #     difflist[[ sprintf('diff%d-%d', i-1, i) ]] <- 0
# # #     jaccardlist[[sprintf('jac%d-%d', i-1, i)]] <- 0
# # #   } else {
# # #     difflist[[ sprintf('diff%d-%d', i-1, i) ]] <-  bi_env_mat_new - outlist[[ (i-1) ]]
# # #     jaccardlist[[sprintf('jac%d-%d', i-1, i)]] <- get_jaccard_index(m0 = outlist[[ (i-1) ]], m1 = bi_env_mat_new )
# # #   }
# # # 
# # #   # print(sim_i[[1]][[3]]`1`)
# # # }
# # # 
# # # par(mfrow=c(1,3))
# # # stability_vec <- cumsum(plyr::ldply(jaccardlist)[,2])/(1:n)
# # # plot(stability_vec , type='l' , main='Change in Stability  [t-1, t]')
# # 
# # # stability_delta <- stability_vec / (1:n)
# # # plot(stability_delta, type='o', log='y',
# # #      xlab='t', 
# # #      # ylim=c( stability_delta[n-1],  1 ),
# # #      ylab='Ln Stability Change [t-1, t]', main='Stabilization Rate (Ln Change in Stability)' 
# # # ); abline(h = tol, col='pink', lty=2)
# # # 
# # # burn_prop <- 0.2
# # # iter_postburn <- round(c( 1-burn_prop, burn_prop) * n )
# # # hist(stability_vec[ iter_postburn[1]:iter_postburn[2] ], main='Stability (post-burn)')
# # 
# # 
# # # ########################################################################
# # # sims = saomnkrsiena$rsiena_model$sims
# # # outlist <- list()
# # # difflist <- list()
# # # jaccardlist <- list()
# # # for(i in 1:length(sims)) {
# # #   cat(sprintf(' %s ', i))
# # #   ##
# # #   el_bi_env <- sims[[ i ]][[1]][[3]]$`1`
# # #   ## update numbering of second mode (the comonent integer names shift upward by the number of actors)
# # #   el_bi_env[,2] <- el_bi_env[,2] + saomnkrsiena$M
# # #   ##
# # #   MplusN <- saomnkrsiena$M + saomnkrsiena$N
# # #   #############
# # #   ## Bipartite matrix space (N+M by N+M)
# # #   ## Undirected --> Upper right rectangle of full bipartite matrix
# # #   ##   M N
# # #   ## M[0,X], for X in [0,1]
# # #   ## N[0,0]
# # #   bi_env_mat_sp <- sparseMatrix(i = el_bi_env[,1],
# # #                                 j = el_bi_env[,2],
# # #                                 x = el_bi_env[,3],
# # #                                 dims = c(MplusN, MplusN))
# # #   ##
# # #   bi_env_mat_new  <- as.matrix(bi_env_mat_sp)[ 1:saomnkrsiena$M, (saomnkrsiena$M+1):(MplusN) ]
# # #   ##
# # #   # saomnkrsiena$plot_bipartite_system_from_mat(bi_env_mat_new, i, plot_save = TRUE) 
# # #   ##
# # #   outlist[[ sprintf('sim%d',i) ]] <- bi_env_mat_new
# # #   
# # #   if (i > 1) {
# # #     difflist[[ sprintf('diff%d-%d', i-1, i) ]] <-  bi_env_mat_new - outlist[[ (i-1) ]] 
# # #     
# # #     jaccardlist[[sprintf('jac%d-%d', i-1, i)]] <- get_jaccard_index(bi_env_mat_new, outlist[[ (i-1) ]] )
# # #   }
# # #   
# # #   # print(sim_i[[1]][[3]]`1`)
# # # }
# # # 
# # # par(mfrow=c(1,3))
# # # stability_vec <- cumsum(plyr::ldply(jaccardlist)[,2])/(1:length(jaccardlist))
# # # plot(stability_vec , type='l' )
# # # 
# # # stability_delta <- stability_vec / 1:length(stability_vec)
# # # plot(stability_delta, type='o', log='y',
# # #      xlab='t', ylab='Ln Stability Change [t-1, t]', main='Change in Stability' )
# # # 
# # # hist(stability_delta[100:length(stability_delta)])
# # # 
# # # ########################################################################
# # 
# # 
# # saom_nk_enhanced$local_search_saom_batchrun(batches = 2, iterations = 3, plot_save = T)
# # saom_nk_enhanced$plot_fitness_progress()
# # saom_nk_enhanced$plot_degree_progress()
# # # saom_nk_enhanced$visualize_networks()
# # 
# # 
# # saom_nk_enhanced$local_search_saom_batchrun(batches = 4, iterations = 10, 
# #                                             plot_save=TRUE, filename='__TEST_PLOT2__') 
# # # saom_nk_enhanced$visualize_networks(plot_save = T)
# # 
# # 
# # 
# # #### SANDBOX ########
# # 
# # saomnkrsiena$rsiena_model$sims[[1]][[1]][[2]]
# # 
# # M <- saomnkrsiena$M
# # N <- saomnkrsiena$N
# # actors <- 1:saomnkrsiena$M 
# # components <- 1:saomnkrsiena$N
# # el_bi_env <- saomnkrsiena$rsiena_model$sims[[1]][[1]][[3]]$`1`
# # el_bi_env[,2] <- el_bi_env[,2] + M
# # # bi_env_el <- el_bi_env[,1:2]
# # ##
# # el_proj1  <- saomnkrsiena$rsiena_model$sims[[1]][[1]][[1]]$`1`
# # el_proj2  <- saomnkrsiena$rsiena_model$sims[[1]][[1]][[2]]$`1`
# # el_proj2 <- el_proj2 + M
# # 
# # ## Undirected --> Upper right rectangle of full bipartite matrix
# # ##   M N
# # ## M[0,X] for X=[0,1]
# # ## N[0,0]
# # bi_env_mat_sp <- sparseMatrix(i = el_bi_env[,1],
# #                               j = el_bi_env[,2],
# #                               x = el_bi_env[,3],
# #                               dims = c(M+N, M+N))
# # bi_env_mat <- as.matrix(bi_env_mat_sp)[1:M,(M+1):(M+N)]
# # 
# # # self$set_system_from_bipartite_igraph( self$random_bipartite_igraph() )
# # 
# # 
# # #####################
# # 
# # 
# # 
# # #########################################
# # 
# # ## Simple environment 
# # env_simple <- SAOM_NK_Enhanced$new(M = 12, N = 8, BI_PROB = .15, sim_name = '_ENV_SIMPLE_')
# # env_simple$local_search_saom_batchrun(batches = 10, iterations = 5, plot_save=TRUE) 
# # 
# # ## Mid environment 
# # env_mid <- SAOM_NK_Enhanced$new(M = 12, N = 12, BI_PROB = .15, sim_name = '_ENV_MID_')
# # env_mid$local_search_saom_batchrun(batches = 10, iterations = 5, plot_save=TRUE)  
# # 
# # ## Complex environment 
# # env_complex <- SAOM_NK_Enhanced$new(M = 12, N = 18, BI_PROB = .15, sim_name = '_ENV_COMPLEX_')
# # env_complex$local_search_saom_batchrun(batches = 10, iterations = 5, plot_save=TRUE) 
# # 
# # ## Very Complex environment 
# # env_verycomplex <- SAOM_NK_Enhanced$new(M = 12, N = 30, BI_PROB = .15, sim_name = '_ENV_VERYCOMPLEX_')
# # env_verycomplex$local_search_saom_batchrun(batches = 15, iterations = 5, plot_save=TRUE) 
# # 
# # 
# # #########################################
# 
# 
# 
# 
# # library(RSiena) # or RSienaTest
# # 
# # ###############################################################################
# # ###                                                                         ###
# # ###         First main function: SimulateNetworks                           ###
# # ###                                                                         ###
# # ###############################################################################
# # 
# # 
# # SimulateNetworks <- function(n, M, rate, dens, rec, tt, c3,
# #                              Vaego, Vaalt, Vasim, Vbego, Vbalt, Vbsim){
# #   # Simulates M consecutive network waves, with n actors,
# #   # according to a stochastic actor-oriented model
# #   # with parameter values rate for rate,
# #   # dens for outdegree, rec for reciprocity,
# #   # tt for transitive triplets, c3 for 3-cycles,
# #   # an actor covariate Va with values alternating between 0 and 1,
# #   # with parameter values Vaego, Vaalt, Vasim
# #   # for egoX, altX, and simX with respect to Va,
# #   # and an actor covariate Vb with a standard normal distribution,
# #   # with parameter values Vbego, Vbalt, Vbsim
# #   # for egoX, altX, and simX with respect to Vb.
# #   ##
# #   # Create actor covariates
# #   V0 <- rep(0, n)
# #   V0[2*(1:(n %/% 2))] <- 1 # equal binary
# #   V1 <- rnorm(n, 0, 1)
# #   # Create initial 2-wave data to get a suitable data structure.
# #   # arbitrarily, this initial network has an expected average degree of 3
# #   X0 <- matrix(rbinom(n*n,1,3/(n-1)),n,n)
# #   diag(X0) <- 0
# #   X1 <- X0
# #   # but X0 and X1 should not be identical for use in sienaDependent
# #   X0[1,2] <- 0
# #   X0[2,1] <- 1
# #   X1[1,2] <- 1
# #   X1[2,1] <- 0
# #   XX <- array(NA,c(n,n,2))
# #   XX[,,1] <- X0
# #   XX[,,2] <- X1
# #   # With this data structure, we now can create the data.
# #   Va <- coCovar(V0)
# #   Vb <- coCovar(V1)
# #   X   <- sienaDependent(XX, allowOnly = FALSE)
# #   InitData <- sienaDataCreate(X, Va, Vb)
# #   InitEff0 <- getEffects(InitData)
# #   # sink to avoid printing to the screen
# #   sink("eff.txt")
# #   # Specify the parameters.
# #   # The rate parameter is first multiplied by 10,
# #   # which will be used only to get from the totally random network XX[,,1] = X0
# #   # to the network that will be the simulated first wave.
# #   InitEff0 <- setEffect(InitEff0, Rate, type="rate", initialValue = 10*rate)
# #   InitEff0 <- setEffect(InitEff0, density, initialValue = dens)
# #   InitEff0 <- setEffect(InitEff0, recip, initialValue = rec)
# #   InitEff0 <- setEffect(InitEff0, transTrip, initialValue = tt)
# #   InitEff0 <- setEffect(InitEff0, cycle3, initialValue = c3)
# #   InitEff0 <- setEffect(InitEff0, egoX, interaction1="Va", initialValue = Vaego)
# #   InitEff0 <- setEffect(InitEff0, altX, interaction1="Va", initialValue = Vaalt)
# #   InitEff0 <- setEffect(InitEff0, simX, interaction1="Va", initialValue = Vasim)
# #   InitEff0 <- setEffect(InitEff0, egoX, interaction1="Vb", initialValue = Vbego)
# #   InitEff0 <- setEffect(InitEff0, altX, interaction1="Vb", initialValue = Vbalt)
# #   InitEff0 <- setEffect(InitEff0, simX, interaction1="Vb", initialValue = Vbsim)
# #   # The parameter given for n3 should be larger than sum(InitEff0$include)
# #   nthree <- sum(InitEff0$include)	+ 5
# #   InitAlg <- sienaAlgorithmCreate(projname="Init", useStdInits=FALSE,
# #                                   cond=FALSE, nsub=0, n3=nthree, simOnly=TRUE)
# #   # Simulate the first wave.
# #   InitSim   <- siena07(InitAlg, data=InitData, eff=InitEff0,
# #                        returnDeps=TRUE, batch=TRUE, silent=TRUE)
# #   # Now prepare for simulating waves 2 to M.
# #   # Create empty result network.
# #   Xs <- array(0, dim=c(n,n,M))
# #   # The rate parameter value from the function call is reinstated in InitEff.
# #   InitEff <- InitEff0
# #   InitEff <- setEffect(InitEff, Rate, type="rate", initialValue = rate)
# #   sink()
# #   for (m in 1:M){
# #     # Note that we start this loop with a previously simulated network.
# #     # Transform the previously simulated network
# #     # from edge list into adjacency matrix
# #     XXsim <- matrix(0,n,n)
# #     nsim  <- InitAlg$n3
# #     XXsim[InitSim$sims[[nsim]][[1]]$X[[1]][,1:2]]  <- InitSim$sims[[nsim]][[1]]$X[[1]][,3]
# #     # Put simulated network into the result matrix.
# #     Xs[,,m] <- XXsim
# #     # Put simulated network in desired places for the next simulation
# #     XX[,,2] <- XX[,,1] # used only to get the data structure
# #     XX[,,1] <- XXsim
# #     if (m < M){
# #       # The following is only to prevent the error that would occur
# #       # in the very unlikely event XX[,,1] == XX[,,2].
# #       if (identical(XX[,,1], XX[,,2])){XX[1,2,1] <- 1 - XX[1,2,2]}
# #       # Specify the two-wave network data set starting with XX[,,1].
# #       X <- sienaDependent(XX, allowOnly = FALSE)
# #       # Simulate wave m+1 starting at XX[,,1] which is the previous XXsim
# #       InitData  <- sienaDataCreate(X, Va, Vb)
# #       InitSim <- siena07(InitAlg, data=InitData, eff=InitEff,
# #                          returnDeps=TRUE, batch=TRUE, silent=TRUE)
# #     }
# #   }
# #   # Present the average degrees to facilitate tuning the outdegree parameter
# #   # to achieve a desired average value for the average degrees.
# #   cat("Average degrees ", round(colSums(Xs,dims=2)/n, digits=2), "\n")
# #   # Result: simulated data set; covara and covarb are vectors of length n;
# #   # networks is an array of dimension nxnxM
# #   list(covara = V0, covarb = V1, networks = Xs)
# # }
# # 
# # 
# # ###############################################################################
# # ###                                                                         ###
# # ###         Examples                                                        ###
# # ###                                                                         ###
# # ###############################################################################
# # 
# # 
# # # Trial values:
# # n <- 20
# # M <- 4
# # rate <- 2
# # dens <- -1.9
# # rec <- 2
# # tt <- 0.3
# # c3 <- -0.3
# # Vaego <- 0
# # Vaalt <- 0
# # Vasim <- 0.6
# # Vbego <- 0.5
# # Vbalt <- -0.5
# # Vbsim <- 0.5
# # 
# # # Example call:
# # SN <- SimulateNetworks(n, M, rate, dens, rec, tt, c3, Vaego, Vaalt, Vasim,
# #                        Vbego, Vbalt, Vbsim)
# # # You can repeat this call a few times, and then see the varying values
# # # reported for the average degrees.
# # # You can also experiment this with other values for dens,
# # # keeping everything else the same,
# # # varying dens by values of +/- 0.05 to +/- 0.2, for example.
# # 
# # # For larger n, slightly lower values of dens are required
# # # to achieve roughly the same average degrees. For example:
# # n <- 30
# # dens <- -2.0
# # SN <- SimulateNetworks(n, M, rate, dens, rec, tt, c3, Vaego, Vaalt, Vasim,
# #                        Vbego, Vbalt, Vbsim)
# # 
# # # Results:
# # SN[[1]]
# # # the same as
# # SN$covara
# # 
# # SN[[2]]
# # # the same as
# # SN$covarb
# # 
# # SN[[3]]
# # # the same as
# # SN$networks
# 
# 
# 
# 
