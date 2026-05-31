rm(list=ls())  ## uncomment to clear environ before run
########################
##
##   SAOM-NK Runscript:
##
##   1. Test 1
##
##
#######################
library(plyr)
library(dplyr)

## Directories
dir_proj <- 'C://Users//sdr8y//OneDrive - University of Missouri//Research//Search_networks//SaoMNK//R'
dir_proj_base <- 'C://Users//sdr8y//OneDrive - University of Missouri//Research//Search_networks'
dir_data <- 'D://Search_networks'

## default settings: Users do not change; TODO: implment within restricted class attributes
DV_NAME <- 'self$bipartite_rsienaDV'


###############  Load R6 Class DEPENDENCIES ############################
## Biparite Environment Search Simulation Class
SaomNkRSienaBiEnv <- source(file.path(dir_proj, 'SAOM_NK_R6_model.R'))$value
# ## RSiena search Class
# SaomNkRSienaBiEnv_search_rsiena <- source(file.path(dir_proj, 'SAOM_NK_R6_search_rsiena_model.R'))$value
###########
## Working director
setwd(dir_proj)

# readxl::read_excel(file.path(dir_proj_base, 'SAOMNK_gridsearch_params_base_v1.xlsx'))
# 
# 
# read.csv(file.path(dir_proj, 'SAOM_NK_gridsearch_params_base_v1.csv'), stringsAsFactors = F)



# self <- readRDS('D:\\Search_networks\\__SAOMNK_gridsearch_input_basic_M3__173666069463852\\grid_results\\2_SAOMNK_gridsearch_input_basic_M3_173666069463852_1601c9ba-f770-4bd9-a8b2-62bb926a4f7f.rds')




environ_params <- list(
  M = 27,        ## Actors
  N = 16,       ## Components
  BI_PROB = .01, ## Environmental Density (DGP hyperparameter)
  component_matrix_start = 'rand', ##**TODO** Implement: 'rand','modular','semi-modular',...
  rand_seed = 1234,
  visualize_init = F,
  name = '_test_fitness_'
)
strategies <- list(
  egoX   =  c(-1,0, 1), #c(0),
  inPopX =  c(1,0, -1)  #c(0),
)
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
      list(effect='density', parameter= -1, fix=T, dv_name=DV_NAME), ##interaction1 = NULL
      list(effect='inPop',   parameter= 2,  fix=T, dv_name=DV_NAME), #interaction1 = NUL
      list(effect='outAct',  parameter= .1, fix=T, dv_name=DV_NAME)
    ),
    coCovars = list( ##**STRATEGY -- MONADIC CONSTANT COVARIATE EFFECTS **
      list(effect='altX',   parameter= 30, fix=T,dv_name=DV_NAME,interaction1='self$component_1_coCovar', x= component_payoffs ),
      list(effect='egoX',   parameter= .1,fix=T,dv_name=DV_NAME,interaction1='self$strat_1_coCovar', x= actor_strats[[1]] ), #interaction1 = NULL
      list(effect='inPopX', parameter= 1, fix=T,dv_name=DV_NAME,interaction1='self$strat_2_coCovar', x= actor_strats[[2]] )
    ),
    varCovars = list() ##**MONADIC TIME-VARYING COVARIATE EFFECTS -- DYNAMIC STRATEGY PROGRAMS**
  )
)

##****************************************##
## II. SIM ANALYSIS 
##****************************************##
# ## INIT SIM ENVIRONMENT: 
# env1 <- SaomNkRSienaBiEnv$new(environ_params)

set.seed(123)
seeds <- sample(1:999999, 5, replace = F)
plist <- list()
for (i in 1:length(seeds)) {
  
  ## INIT SIM ENVIRONMENT: 
  env1 <- SaomNkRSienaBiEnv$new(environ_params)
  
  ## 1.1. Search 1. SHORT Run
  env1$search_rsiena_multiwave_run(
    structure_model, 
    waves=1, ##id='_short_run',
    iterations = 4 * env1$M * env1$N, 
    rand_seed = seeds[i]
  )
  # Process results
  env1$search_rsiena_multiwave_process_results()
  #
  # env1$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T) #+  geom_vline(xintercept=seq(0,440,8*9))
  
  plist[[as.character(seeds[i])]] <- list(
    plot = env1$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T, thin_factor=3, loess_span=.25), #+  geom_vline(xintercept=seq(0,440,8*9))
    sim = env1,
    seed=seeds[i]
  )
  #
  # #
  # env1$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T, thin_factor=10, loess_span=.3)
}

sapply(plist,print)
# env1$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T, show_utility_points = T, thin_factor = 8)
# # env1$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T, show_utility_points = F, thin_factor = 1)
# 
# env1$search_rsiena_multiwave_plot_K_4panel(return_plot = T, thin_factor = 1)
# 
# env1$search_rsiena_multiwave_plot_utility_ridge_density_by_strategy(return_plot = T, show_strategy_means = F, plot_periods=6)
# 
# 
# 
# # env0 = env1
# 
# #
# env1$
# env1$bi_env_arr[,,1:5]


##############################################################
self <- env1
nwaves <- 2
# step_ids <- round(seq(2*self$M*self$N, self$rsiena_model$n3, by=1*self$M))[1:nwaves] ##, length=nwaves))
step_ids <-  1*self$M*self$N + seq(1, self$rsiena_model$n3, by=round(self$M * log2(self$N)/2))[1:nwaves] ##, length=nwaves))
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
array_bi_net <-  self$bi_env_arr[,,step_ids] ##**MAIN STEP**
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
##  4. RSiena Algorithm
self$rsiena_algorithm <- sienaAlgorithmCreate(projname=sprintf('%s_%s',self$SIM_NAME,self$TIMESTAMP),
                                              # simOnly = F,
                                              nsub = rsiena_phase2_nsub * 1,
                                              # n2start = rsiena_n2start_scale * 2.52 * (7+sum(self$rsiena_effects$include)),
                                              n3 = iterations*15,
                                              seed = rand_seed)
## 5. Run RSiena simulation
self$rsiena_model <- siena07(self$rsiena_algorithm,
                             data = self$rsiena_data, 
                             effects = self$rsiena_effects,
                             batch = TRUE,
                             returnDeps = returnDeps, 
                             # returnChains = returnChains,
                             returnDataFrame = TRUE, ##**TODO** CHECK
                             returnLoglik = TRUE #,  ##**TODO** CHECK
)   # returnChains = returnChains

# Summarize and plot results
mod_summary <- summary(self$rsiena_model)
if(!is.null(mod_summary))
  print(mod_summary)
# plot(self$rsiena_model)
digits <- 3
print(screenreg(list(self$rsiena_model), single.row = T, digits = digits))


print(self$rsiena_model$tconv)
print( abs(self$rsiena_model$tconv) <  0.1 )
print(self$rsiena_model$tconv.max)
print(self$rsiena_model$tconv.max[1] < 0.25)
all(
  abs(self$rsiena_model$tconv) <  0.1 &
  self$rsiena_model$tconv.max[1] < 0.25
)

gof.od <- RSiena::sienaGOF(self$rsiena_model, OutdegreeDistribution, levls=0:self$N, varName = 'self$bipartite_rsienaDV')
# gof.od <- RSiena::sienaGOF(self$rsiena_model, OutdegreeDistribution, varName = 'self$bipartite_rsienaDV')
plot(gof.od)

gof.id <- RSiena::sienaGOF(self$rsiena_model, IndegreeDistribution, levls=0:self$M, varName = 'self$bipartite_rsienaDV')
plot(gof.id)












############################

actor_strats <- self$get_actor_strategies()
madf <- self$actor_util_df %>% mutate(
    strategy=actor_strats[actor_id],
    step_zoo = zoo::as.zoo(chain_step_id)
  ) %>%
  group_by(step_zoo, strategy) %>%
  mutate(
    utility_ma5 = zoo::rollmean(utility, k=5, fill = NA, align = 'center'),
    utility_cumsum = cumsum(utility),
    utility_ma  = cumsum(utility)/unique(chain_step_id)
  ) %>%
  ungroup() %>%
  as.data.table()

# library(data.table)
# # madf <- as.data.table(madf)
# madf <- madf[ , cusum_utility := cumsum(utility), by=.(step_zoo, actor_id) ]



loess_span <- .4
alpha_val <- .25

madf %>% ggplot(aes(x=chain_step_id, y=utility)) +
  geom_smooth(method='loess', span=loess_span) + geom_point(alpha=alpha_val) + 
  theme_bw() +  ggtitle(sprintf('Loess span = %.2f',loess_span))
 



madf %>% ggplot(aes(x=chain_step_id, y=utility, color=strategy, linetype=actor_id)) +
  geom_smooth(method='loess', span=loess_span) + 
  geom_smooth(aes(x=chain_step_id, y=utility_mean), method = 'loess', span=loess_span, 
              data=madf %>% group_by(chain_step_id)%>%summarize(utility_mean=mean(utility)),
              color='black', linetype=1, alpha=.1) +
  geom_point(pch=1, alpha=alpha_val) + 
  # scale_y_log10() +
  theme_bw() +  ggtitle(sprintf('Loess span = %.2f',loess_span)) # + xlim(0,1000)




madf %>% filter(actor_id %in% c(1:3)) %>%
  ggplot(aes(x=chain_step_id, y=utility, color=strategy, linetype=actor_id)) +
  geom_smooth(method='loess', span=loess_span) + 
  geom_smooth(aes(x=chain_step_id, y=utility_mean), method = 'loess', span=loess_span, 
              data=madf %>% group_by(chain_step_id)%>%summarize(utility_mean=mean(utility)),
              color='black', linetype=1, alpha=.1) +
  geom_point(pch=1, alpha=alpha_val) + 
  # scale_y_log10() +
  theme_bw() + # xlim(0,1000)
  ggtitle(sprintf('Loess span = %.2f',loess_span)) #+ xlim(0,1000)


## Utility Transitions Panel by Strategy with population average black line
madf %>% ggplot(aes(x=chain_step_id, y=utility, color=strategy, linetype=actor_id)) +
  geom_smooth(method='loess', span=loess_span) + 
  geom_smooth(aes(x=chain_step_id, y=utility_mean), method = 'loess', span=loess_span, 
              data=madf %>% group_by(chain_step_id)%>%summarize(utility_mean=mean(utility)),
              color='black', linetype=1, alpha=.1) +
  geom_point(pch=1, alpha=alpha_val) + 
  facet_grid(.~strategy) +
  # scale_y_log10() +
  theme_bw() + # xlim(0,1000)
  ggtitle(sprintf('Loess span = %.2f',loess_span)) #+ xlim(0,1000)


# madf %>% #filter(actor_id %in% c(1,3,4,6,7,9)) %>% 
#   ggplot(aes(x=chain_step_id, y=utility_cumsum, color=strategy)) +
#   geom_smooth(method='loess') + geom_point() + theme_bw()




## 6. Update System
## update simulation object environment from RSiena simulation model
new_bi_env_igraph <- self$get_bipartite_igraph_from_rsiena_model()
#
self$set_system_from_bipartite_igraph( new_bi_env_igraph )
############################


mf <- Mod(fft( madf %>% filter(actor_id %in% c(1)) %>% pull(utility)  ))
plot(mf)
##

# tiechdf <- env1$chain_stats %>% filter( !stability )

# > act
# # A tibble: 9 × 8
#      id_from n_decision n_nochange forbearance n_comps_explored exploitation_rate  entropy tie_seq                                 
#        <dbl>      <int>      <int>       <dbl>            <int>             <dbl>   <dbl> <chr>                                   
#   1       1        101          2      0.0198               33              3.06    17.6 7|5|31|20|14|11|4|4|2|6|16|18|1|29|27|2…
#   2       2         72          5      0.0694               31              2.32    18.0 32|26|4|26|14|25|17|4|8|31|3|33|31|13|1…
#   3       3         97          5      0.0515               31              3.13    16.8 17|11|8|33|33|6|14|29|23|7|20|7|2|14|18…
#   4       4         99          5      0.0505               32              3.09    17.3 12|8|12|28|25|31|27|29|15|32|17|28|27|2…
#   5       5        109          3      0.0275               31              3.52    16.6 13|25|5|24|30|3|23|33|21|30|9|24|5|1|12…
#   6       6         94          2      0.0213               33              2.85    18.1 32|32|33|27|22|16|27|19|15|2|14|9|5|31|…
#   7       7        109          0      0                    32              3.41    16.8 1|22|18|13|27|2|20|4|13|18|31|28|22|7|2…
#   8       8        114          2      0.0175               30              3.8     15.6 26|17|9|22|14|9|28|31|4|17|8|18|6|19|31…
#   9       9         96          3      0.0312               31              3.10    16.9 2|8|23|33|33|27|1|15|15|32|23|28|2|9|20…
#
# > round(cor(act[,-c(1,ncol(act))]), 3)
#   n_decision n_nochange forbearance n_comps_explored exploitation_rate entropy
#   n_decision             1.000     -0.590      -0.767           -0.111             0.970  -0.749
#   n_nochange            -0.590      1.000       0.962           -0.285            -0.487   0.221
#   forbearance           -0.767      0.962       1.000           -0.262            -0.659   0.352
#   n_comps_explored      -0.111     -0.285      -0.262            1.000            -0.349   0.726
#   exploitation_rate      0.970     -0.487      -0.659           -0.349             1.000  -0.883
#   entropy               -0.749      0.221       0.352            0.726            -0.883   1.000


#####
#
self <- env1
#
#
act <- self$chain_stats %>% 
  mutate(
    chain_step_id = 1:nrow(self$chain_stats),
    actor_id = as.factor(id_from),
    strategy = as.factor(actor_strats$egoX[id_from])
  ) %>%
  mutate(chain_h1 = ifelse(chain_step_id < median(chain_step_id, na.rm=T), 'First Half','Second Half')) %>%
  group_by(actor_id, chain_h1, strategy)  %>% 
  dplyr::summarize(
   n_nochange=sum(stability), ## sum logical = count TRUE
   n_decision=n(), 
   forbearance = sum(stability)/n(),
   n_comps_explored =length(unique(id_to)),
   strategy = as.factor( paste(unique(self$strat_1_coCovar[ actor_id ]), collapse="|")),
   # explore = length(unique(id_to)) / n(), ## more components per decision = explore
   exploitation_rate = n()/length(unique(id_to)),  ## more decisions among fewer number of components = exploit
   entropy= -sum( (plyr::count(id_to)/n()) * log2(plyr::count(id_to)/n()) ), 
   tie_seq=paste(id_to, collapse="|"),
   step_seq = paste(chain_step_id, collapse = '|'),
   step_precedence_sum = sum(1 - chain_step_id/nrow(self$chain_stats), na.rm=T),
   step_precedence_mean=mean(1 - chain_step_id/nrow(self$chain_stats), na.rm=T),
   step_precedence_var=var(1 - chain_step_id/nrow(self$chain_stats), na.rm=T),
  ) %>% left_join(
    self$actor_util_df %>% 
      mutate(actor_id=as.factor(actor_id), 
             chain_h1 = ifelse(chain_step_id < median(chain_step_id, na.rm=T), 'First Half','Second Half')) %>% 
      mutate(strategy = as.factor(actor_strats$egoX[actor_id])) %>%
      group_by(actor_id, chain_h1, strategy ) %>%
      dplyr::summarize(utility_mean=mean(utility, na.rm=T)), 
    by=c('actor_id', 'chain_h1','strategy')
  )

act

## levels order
coveffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
covparams <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$parameter)
covfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)ifelse(x$fix,'','(var)'))
#
structeffs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$effect)
structparams <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$parameter)
structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)ifelse(x$fix,'','(var)'))
# structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)substr(as.character(x$fix),1,1))
covDvTypes <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$interaction1)
componentDV_ids <- grep('self\\$component_\\d{1,2}_coCovar', covDvTypes) ## ex: "self$component_1_coCovar"
stratDV_ids     <- grep('self\\$strat_\\d{1,2}_coCovar', covDvTypes ) ## ex: "self$strat_1_coCovar" 
#
compoeffs   <- coveffs[ componentDV_ids ]
compoparams <- covparams[ componentDV_ids ]
compofixs   <- covfixs[ componentDV_ids ]
#
strateffs   <- coveffs[ stratDV_ids ]
stratparams <- covparams[ stratDV_ids ]
stratfixs   <- covfixs[ stratDV_ids ]
efflvls <- c('utility', strateffs, compoeffs, structeffs)
facet_strip_colors <- c('black',
                        rep('darkgray', length(strateffs)), 
                        rep('lightgray', length(compoeffs)),
                        rep('white', length(structeffs)))
facet_strip_text_colors <- c('white',rep('black', length(efflvls)-1))

env1$actor_stats_df$strategy <- as.factor( env1$strat_1_coCovar[ env1$actor_stats_df$actor_id ] )
## Add utility as extra 'effect' 
act_effs <- env1$actor_stats_df %>% bind_rows( 
    env1$actor_util_df %>% mutate(
      effect_id=NA, 
      value=utility, 
      effect_name='utility', 
      strategy=factor(env1$strat_1_coCovar[actor_id]), 
      utility=NULL
    )
  ) 
act_effs$effect_name <- factor(act_effs$effect_name, levels=efflvls)

# library(gridExtra)  
# library(ggh4x)
# gg_col_strip_colors <- ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = facet_strip_colors),
#                                            text_y=ggh4x::elem_list_text(color = facet_strip_text_colors)
#                                           )

##**SIGNALS PANEL**
act_effs %>%  ggplot(aes(x=chain_step_id, y=value, color=strategy,fill=strategy, linetype=strategy)) + 
  geom_point(alpha=.04, shape=1)  + 
  geom_smooth(method='loess', alpha=.1) +
  # ggh4x::facet_wrap2(~effect_name , strip = gg_col_strip_colors, scales='free_y') +
  facet_grid(effect_name ~ ., scales='free_y') +
  theme_bw() + theme(legend.position = 'left') +
  ggtitle('Actor Time Series By Strategy: Utility and Component Effects')
# +  theme(  
#     strip.background = element_rect(fill = facet_strip_colors, color = facet_strip_colors),  # Background color  
#     strip.text = element_text(color = facet_strip_text_colors, size = 12)  # Text color and size  
  # ) 

# env1$actor_stats_df %>% 
#   ggplot(aes(x=chain_step_id, y=value, color=strategy,fill=strategy, linetype=strategy)) + 
#   geom_point(alpha=.1, shape=1)  + 
#   geom_smooth(method='lm', alpha=.1) +
#   facet_grid(effect_name ~ ., scales='free_y') +
#   theme_cowplot()

 

act %>% ggplot(aes(x=entropy, fill=strategy, color=strategy)) + 
  geom_density(alpha=.2, linewidth=1) + 
  facet_wrap(~chain_h1) +
  geom_vline(data=act%>%group_by(chain_h1,strategy)%>%dplyr::summarize(mean=mean(entropy)), 
             aes(xintercept = mean, color=strategy, linetype=strategy),  linewidth=1.1) +
  theme_bw()

act %>% ggplot(aes(x=utility_mean, fill=strategy, color=strategy)) +
  geom_density(alpha=.2, linewidth=1) + 
   geom_vline(data=act%>%group_by(chain_h1,strategy)%>%dplyr::summarize(mean=mean(utility_mean)), 
             aes(xintercept = mean, color=strategy, linetype=strategy),  linewidth=1.1) +
  facet_wrap(~chain_h1) +
  theme_bw()

act %>% ggplot(aes(x=forbearance, fill=strategy, color=strategy)) +
  geom_density(alpha=.2, linewidth=1) + 
  geom_vline(data=act%>%group_by(chain_h1,strategy)%>%dplyr::summarize(mean=mean(forbearance)), 
             aes(xintercept = mean, color=strategy, linetype=strategy), linewidth=1.1) +
  facet_wrap(~chain_h1) +
  theme_bw()

act %>% ggplot(aes(x=exploitation_rate, fill=strategy, color=strategy)) +
  geom_density(alpha=.2, linewidth=1) + 
  geom_vline(data=act%>%group_by(chain_h1,strategy)%>%dplyr::summarize(mean=mean(exploitation_rate)), 
             aes(xintercept = mean, color=strategy, linetype=strategy),  linewidth=1.1) +
  facet_wrap(~chain_h1) +
  theme_bw()

act %>% ggplot(aes(x=n_comps_explored, fill=strategy, color=strategy)) +
  geom_density(alpha=.2, linewidth=1) +
  # geom_histogram(alpha=.1, position = 'jitter') +
  geom_vline(data=act%>%group_by(chain_h1,strategy)%>%dplyr::summarize(mean=mean(n_comps_explored)), 
             aes(xintercept = mean, color=strategy, linetype=strategy),  linewidth=1.1) +
  facet_wrap(~chain_h1) +
  theme_bw()




.##----------------------------------------
# act %>% ggplot(aes(x=n_decision, y=utility_mean, fill=strategy, color=strategy, shape=strategy)) +
#   # geom_density(alpha=.2, linewidth=1) +
#   geom_point(size=3) +
#   # geom_histogram(alpha=.1, position = 'jitter') +
#   geom_hline(data=act%>%group_by(chain_h1,strategy)%>%dplyr::summarize(mean=mean(utility_mean)), 
#              aes(yintercept = mean, color=strategy, linetype=strategy),  linewidth=1.1) +
#   facet_wrap(~chain_h1) +
#   theme_bw()


actaggr <- act %>% mutate(aggr_f = ifelse(n_decision-n_nochange > median(n_decision-n_nochange, na.rm=T),'High Aggr.','Low Aggr.')) 
actaggr %>%
  ggplot(aes(x=entropy, y=utility_mean, fill=strategy, color=strategy, shape=strategy)) +
  # geom_density(alpha=.2, linewidth=1) +
  geom_point(size=2) +
  # geom_density_2d(aes(colour  = strategy), size = 0.5, linetype=2) +  # Density contours
  stat_density2d(data = actaggr, geom='polygon', contour=TRUE,
                 aes(fill=strategy, color=strategy),bins=9, linewidth=.5, alpha=.05) +
  geom_smooth(method = 'lm', alpha=.01, linewidth=.9, linetype=1) +
  # geom_histogram(alpha=.1, position = 'jitter') +
  # geom_hline(data=actaggr%>%group_by(chain_h1,strategy,aggr_f)%>%dplyr::summarize(mean=mean(utility_mean)), 
  #            aes(yintercept = mean, color=strategy, linetype=strategy),  linewidth=1.1) +
  facet_grid(aggr_f ~ chain_h1) +
  theme_pubclean()

# actaggr <- act %>% mutate(aggr_f = ifelse(n_decision-n_nochange > median(n_decision-n_nochange, na.rm=T),'High Aggr.','Low Aggr.')) 
actaggr %>%
  ggplot(aes(x=forbearance, y=utility_mean, fill=strategy, color=strategy, shape=strategy)) +
  # geom_density(alpha=.2, linewidth=1) +
  geom_point(size=2) +
  # geom_density_2d(aes(colour  = strategy), size = 0.5, linetype=2) +  # Density contours
  stat_density2d(data = actaggr, geom='polygon', contour=TRUE,
                 aes(fill=strategy, color=strategy),bins=6, linewidth=.7, alpha=.05) +
  geom_smooth(method = 'lm', alpha=.01, linewidth=.9, linetype=1) +
  # geom_histogram(alpha=.1, position = 'jitter') +
  # geom_hline(data=actaggr%>%group_by(chain_h1,strategy,aggr_f)%>%dplyr::summarize(mean=mean(utility_mean)), 
  #            aes(yintercept = mean, color=strategy, linetype=strategy),  linewidth=1.1) +
  facet_grid(aggr_f~chain_h1) +
  theme_pubclean()


actaggr %>%
  ggplot(aes(x=n_comps_explored, y=utility_mean, fill=strategy, color=strategy, shape=strategy)) +
  # geom_density(alpha=.2, linewidth=1) +
  geom_point(size=2) +
  # geom_density_2d(aes(colour  = strategy), size = 0.5, linetype=2) +  # Density contours
  stat_density2d(data = actaggr, geom='polygon', contour=TRUE,
                 aes(fill=strategy, color=strategy),bins=6, linewidth=.7, alpha=.05) +
  geom_smooth(method = 'lm', alpha=.01, linewidth=.9, linetype=1) +
  # geom_histogram(alpha=.1, position = 'jitter') +
  # geom_hline(data=actaggr%>%group_by(chain_h1,strategy,aggr_f)%>%dplyr::summarize(mean=mean(utility_mean)), 
  #            aes(yintercept = mean, color=strategy, linetype=strategy),  linewidth=1.1) +
  facet_grid(aggr_f~chain_h1) +
  theme_pubclean()

actaggr %>%
  ggplot(aes(x=exploitation_rate, y=utility_mean, fill=strategy, color=strategy, shape=strategy)) +
  facet_grid(aggr_f~chain_h1) +
  # geom_density(alpha=.2, linewidth=1) +
  geom_point(size=2) +
  # geom_density_2d(aes(colour  = strategy), size = 0.5, linetype=2) +  # Density contours
  stat_density2d(data = actaggr, geom='polygon', contour=TRUE,
                 aes(fill=strategy, color=strategy),bins=7, linewidth=.7, alpha=.1) +
  geom_smooth(method = 'lm', alpha=.01, linewidth=.9, linetype=1) +
  # geom_histogram(alpha=.1, position = 'jitter') +
  # geom_hline(data=actaggr%>%group_by(chain_h1,strategy,aggr_f)%>%dplyr::summarize(mean=mean(utility_mean)), 
  #            aes(yintercept = mean, color=strategy, linetype=, strategy),  linewidth) +
  theme_pubclean()

# plot(y=act$utility_mean, x=act$n_decision)




######################
##**TODO**

self <- env1


#
config_param_vals <- c(
  unlist(sapply(self$config_structure_model$dv_bipartite$effects, function(x) x$parameter)),
  unlist(sapply(self$config_structure_model$dv_bipartite$coCovars, function(x) x$parameter)),
  unlist(sapply(self$config_structure_model$dv_bipartite$varCovars, function(x) x$parameter))
)
fixed_params <- c(
  unlist(sapply(self$config_structure_model$dv_bipartite$effects, function(x) x$fix)),
  unlist(sapply(self$config_structure_model$dv_bipartite$coCovars, function(x) x$fix)),
  unlist(sapply(self$config_structure_model$dv_bipartite$varCovars, function(x) x$fix))
)
#
#UTILITY
theta   <- self$rsiena_model$theta
## replace fixed theta param with given param values (instead of zero default when effect is fixed)
if (any(fixed_params)) theta[ fixed_params ] <- config_param_vals[ fixed_params ]



## PLOT ACTOR UTILITY LANDSCAPE
utilist <- list()
statsl <- list()

# step_ids <- round(c(1,8,64,512, self$rsiena_model$n3/2, self$rsiena_model$n3))
acpds <- seq(1, self$rsiena_model$n3-1, by= self$M*self$N )
step_ids <- c(acpds[1:6])
for (step_id in step_ids) {
  # step_id <- 1
  cat(sprintf('\nstep %s: actors ', step_id))
  
  bi_env_mat_step <- self$bi_env_arr[, , step_id]
  
  act_counterfacts <- list()
  
  pltlist <- list()
  for (i in 1:9) {
    cat(sprintf(' %s ', i))
    
    
    ##**ACTOR i DECISION PERSPECTIVE** 0000000000000000000000000000000000
    
    ##**TODO**
    ## ALL actor-component counterfactual configuations for Actor i  (2^N rows by N cols)
    iland = expand.grid(lapply(1:self$N, function(x) 0:1 ))
    iland_config_step_row_id <- which(apply(iland, 1, function(x) all(x == bi_env_mat_step[i,]) ))
    
    tmpmat <- bi_env_mat_step
    ## ifit dimensions [ M, 2^N ]
    ifit <- apply(iland, 1, function(x){
      tmpmat[i,] <- x  ## set counterfactual actor-component configuration
      return( self$get_struct_mod_stats_mat_from_bi_mat( tmpmat ) %*% theta ) ## compute utility vector
    }) 
    
    # apply
    
    ids.max <- which(ifit == max(ifit), arr.ind = TRUE) 
    fits.max <- ifit[ ids.max ]
    nmax <- length(fits.max)
    
    # ifit
    
    # ids.max[1,]
    
    act_counterfacts[[i]] <- list(
      ifit = ifit, #Given all other ties, Actor i's configurations applied to utility func for all actors
      ids.max = ids.max,
      fits.max = fits.max,
      nmax = nmax
    )
    
    ## distances of each counterfactual configuration 
    dist_counterfac <- iland - bi_env_mat_step[i, ]
    z <- apply( dist_counterfac, 1, function(x) {
      c(dist=sum(x!=0),`drop`=sum(x==-1), `nochange`=sum(x==0), `add`=sum(x==1), 
        counterfac=paste(x, collapse = '|'), start=paste(bi_env_mat_step[i, ], collapse = '|') )
    })
    # z <- z[,order(z['d',], decreasing = F)]
    ##**TODO**
    ##**INFORMATION LOSS**
    ##  This uses only ego's own counterfactual fits (landscape)
    ##  but the corresponding other actor's affected fits (landscapes are not currently used )
    z <- rbind(z, utility_ego=ifit[ i , ] )
    z <- rbind(z, utility_alter_mean= colMeans(ifit[ -i , ], na.rm = T) )
    z <- rbind(z, utility_alter_sd  = apply(ifit[ -i , ], 2, function(x) sd(x, na.rm = T)) )
    
    wdf <- as.data.frame( t(z) )
    
    actfit_long <- wdf %>% 
      mutate(
        dist = as.numeric(dist),
        drop = as.numeric(drop),
        nochange = as.numeric(nochange),
        add = as.numeric(add) ,
        utility_ego=as.numeric(utility_ego),
        utility_alter_mean = as.numeric(utility_alter_mean),
        utility_alter_sd = as.numeric(utility_alter_sd)
      ) %>%
      pivot_longer(cols = c( drop:add, utility_ego:utility_alter_sd )) 
     
    
    plt <- actfit_long %>% ggplot(aes(x=value))+ 
      # geom_density() + 
      geom_histogram() +
      facet_grid(dist ~ name, scales='free_x') + theme_bw() + 
      ggtitle(sprintf('Utility Transition Paths: Actor %s, Step %s', i, step_id))
    
    pltlist[[i]] <- plt
    
    # z <- apply(iland, 1, function(x) {
    #   x <- plyr::count( x - bi_env_mat_step[i, ]) 
    #   cnts <- merge(x$freq, data.frame(x=c())
    #   names(cnts) <- x$x
    #   return(cnts)
    # }, simplify = T)
    # 
    # stringdist::stringdist( , , method='hamming')
    
    step_actor_key <- sprintf('%s|%s',step_id,i)
    
    statsl[[step_actor_key]] <- list(
      actfit_long = actfit_long,
      act_counterfacts = act_counterfacts,
      plt = plt,
      actor_id=i, chain_step_id = step_id,
      strategy=as.factor(self$strat_1_coCovar[i])
    )
    
    utilist[[step_actor_key]] <- actfit_long %>% 
      filter(name %in% c('utility_ego','utility_alter_mean')) %>%
      mutate(actor_id=i, chain_step_id = step_id, strategy=as.factor(self$strat_1_coCovar[i]))
    
    # ## STEP 1
    # # actor_step_land_arr <- array(NA, dim=c(self$N, self$N))
    # for (j1 in 1:self$N) {
    #   
    #   bi_env_mat_step_togj <- self$toggleBiMat( bi_env_mat_step,  i ,  j )
    #   
    #   stats_step_togj <- self$get_struct_mod_stats_mat_from_bi_mat( bi_env_mat_step_togj )
    #    
    #   
    #   stringdist::seq_dist()
    #   
    # 
    #   dfij <- as.data.frame( stats_step_togj )
    #   ## UTILITY
    #   dfij$utility  <- c( stats_step_togj %*% theta ) ## c() converts 1-col matrix to numeric vector
    #   ##
    #   dfij$ij_toggled_on <- ( bi_env_mat_step_togj[i,j] == 1 )
    #   dfij$ij_added <- bi_env_mat_step_togj[i,j] > bi_env_mat_step[i,j]
    #   dfij$ij_dropped <- bi_env_mat_step_togj[i,j] < bi_env_mat_step[i,j]
    #   dfij$id_from <- i
    #   dfij$id_to   <- j
    #   dfij$id_actor_counterfac <- 1:nrow(dfij)
    #   dfij$chain_step_id <- step_id
    #   dfij$decider_strategy   <- as.factor(actor_strats[[1]])[i]
    #   dfij$others_strategy <- as.factor(actor_strats[[1]])
    #   dfij$chain_half <- ifelse(step_id < (self$rsiena_model$n3 / 2), 'H1', 'H2')
    #   
    #   #
    #   statsl[[sprintf('%s|%s-%s',step_id,i,j)]] <- dfij
    # }
    
  }
  
  # Arrange plots across multiple pages
  multi_page_plot <- marrangeGrob(pltlist, nrow = 1, ncol = 1)
  
  # Save to a PDF
  ggsave("multi_page_plots.pdf", multi_page_plot, width = 8, height = 8)
  
}

#
utildf <- data.table::rbindlist(utilist, use.names = T, idcol = 'step_actor_key')
#
# statsdf <- data.table::rbindlist(statsl, use.names = T, idcol = 'toggle_id')


#
varname <- 'utility_ego'
utildfvar <- utildf%>%filter(name==varname)
ggplot(utildfvar, aes(x=value, color=strategy, fill=strategy, linetype=strategy)) + 
  # geom_density(alpha=.05, size=.9) +
  geom_histogram(alpha=.15, position = 'identity')+
  geom_vline(aes(xintercept = mean, color=strategy, linetype=strategy),
             data = utildfvar %>%
               group_by(dist,strategy,chain_step_id)%>%
               dplyr::summarize(mean=mean(value))
             ) +
  facet_grid(dist ~ chain_step_id)+ theme_cowplot() + 
  # ggtitle('Fitness Landscapes by Strategy over Time (Actor Counterfactual Utility)')
  ggtitle(sprintf('How far are the peaks over time?\nEvolution of counterfactual distributions by rewiring distance: %s', varname))

## ACTOR POINT SERIES BY PERIOD (lines by actor; colors by strategy)
actutilstrat <- self$actor_util_df%>%
  mutate(strategy=as.factor(self$strat_1_coCovar[actor_id]),
         period = ceiling(chain_step_id / (self$M*self$N)) )
ggplot(actutilstrat, 
       aes(x=chain_step_id, y=utility, color=strategy, linetype = actor_id)) + 
  # geom_density(alpha=.05, size=.9) +
  geom_point(alpha=.4, shape=1)+
  # facet_grid(actor_id ~ ., scales=)+
  geom_vline(xintercept = step_ids) +
  geom_line() +
  # geom_smooth(method='loess') +
  xlim(c(0,self$M * self$N * 6.05)) +
  theme_bw() + 
  ggtitle('ACTOR UTILITY POINT-LINE SERIES colored by STRATEGY vline by PERIOD')

glist <- list()
for (pd_i in 1:6) {
  steps_in_pd <- self$M * self$N
  glist[[pd_i]] <-  ggplot(self$actor_util_df%>%
           mutate(strategy=as.factor(self$strat_1_coCovar[actor_id]),
                  period=ceiling(chain_step_id / steps_in_pd)) %>%
           filter(period <= 6), 
         aes(x=chain_step_id, y=utility, color=strategy, fill=strategy)
    ) + 
    # geom_density(alpha=.05, size=.9) +
    geom_point(alpha=.4, shape=1)+
    geom_vline(xintercept = step_ids) +
    geom_smooth(method='loess') +
    # facet_grid(period~.) +
    xlab('') +
    xlim(c(0,self$M * self$N * pd_i + .05)) +
    theme_bw() #+ ggtitle('SMOOTHED UTILITY by STRATEGY vline by PERIOD')
}
# grid.arrange(glist)

ml <- marrangeGrob(glist, ncol=1, nrow=length(glist), 
                   top='SMOOTHED UTILITY colored by STRATEGY vline by PERIOD',
                   bottom='Chain Step ID') ## nrow=2, 




NtimesM <- self$M * self$N
nstrats <- length(unique(self$strat_1_coCovar))
library(ggrain)


plt.str.rain <- ggplot(actutilstrat %>% filter(period <= 6) %>%
                      mutate(period=factor(period, levels=c('1','2','3','4','5','6'))), #%>%   mutate(strategy_pd=unique(paste(self$strat_1_coCovar[actor_id],period,collapse = '_'))),
                      aes(x=period, y=utility, fill = strategy, color=strategy)) + #shape=actor_id, 
  geom_rain(alpha = .15, rain.side = 'r',  # id.long.var = 'actor_pd',
            boxplot.args.pos = list(width = .1,
                                    position = ggpp::position_dodgenudge(width = .1,
                                                                         x = c(.13, .13, # t1 old, t1 young
                                                                               .13, .13,
                                                                               .13, .13,
                                                                               .13, .13,
                                                                               .13, .13,
                                                                               .13, .13)
                                                                         )
                                    ),
            violin.args.pos = list(width = 1.1,
                                   position = position_nudge(
                                     x = c(rep(.2, 256*2), rep(.2, 256*2),# t1
                                           rep(.2, 256*2), rep(.2, 256*2), # t2
                                           rep(.2, 256*2), rep(.2, 256*2),
                                           rep(.2, 256*2), rep(.2, 256*2),
                                           rep(.2, 256*2), rep(.2, 256*2),
                                           rep(.2, 256*2), rep(.2, 256*2)
                                           )
                                                             )
                                   )
            ) +
  stat_summary(fun = mean, geom = "line",
               aes(group = strategy, color = strategy)) +
  stat_summary(fun = mean, geom = "point",
               aes(group = strategy, color = strategy)) +
  geom_hline(yintercept = 0, linetype=2) +
  guides(fill = 'none', shape='none') + 
  theme_classic() + 
  ggtitle(sprintf('Strategy Transition Paths: Utility Distributions by Period (%s steps)',self$M*self$N))


# varname <- 'utility_alter_mean'
# ggplot(utildf%>%filter(name==varname), 
#        aes(x=value, color=strategy, fill=strategy, linetype=strategy)) + 
#   # geom_density(alpha=.05, size=.9) +
#   geom_histogram(alpha=.2, position = 'identity')+
#   # geom_vline(aes(xintercept = mean, color=alter_strategy, linetype=alter_strategy), 
#   #            data=statsdf%>%group_by(alter_strategy,chain_step_id)%>%dplyr::summarize(mean=mean(utility))) +
#   facet_grid(dist ~ chain_step_id)+ theme_cowplot() + 
#   # ggtitle('Fitness Landscapes by Strategy over Time (Actor Counterfactual Utility)')
#   ggtitle(sprintf('%s', varname))

##
ggplot(statsdf, aes(x=utility,color=fct_rev(factor(chain_step_id)),linetype=fct_rev(factor(chain_step_id)))) + 
  geom_density(alpha=.1) +
  geom_vline(aes(xintercept = mean, color=fct_rev(factor(chain_step_id)), linetype=fct_rev(factor(chain_step_id))), 
             data=statsdf%>%group_by(id_from,chain_step_id)%>%dplyr::summarize(mean=mean(utility))) +
  facet_wrap(~id_from)+ theme_bw() +
  ggtitle('Actor Utility Landscape Distributions over Time (Actor Counterfactual Utility)')

##########

ggplot(statsdf, aes(x=inPopX,
                    color=alter_strategy,
                    fill=alter_strategy,
                    linetype=alter_strategy)) + 
  geom_density(alpha=.1) +
  geom_vline(aes(xintercept = mean, 
                 color=alter_strategy, 
                 linetype=alter_strategy),
             data=statsdf%>%
               group_by(alter_strategy,chain_half)%>%
               dplyr::summarize(mean=mean(inPopX))
             ) +
  facet_grid(chain_half ~.)+
  theme_bw() +
  ggtitle('Stat Distribution')
# ##
# ggplot(statsdf, aes(x=utility,color=fct_rev(factor(chain_step_id)),linetype=fct_rev(factor(chain_step_id)))) + 
#   geom_density(alpha=.1) +
#   facet_wrap(~id_from)+ theme_bw()


# ggplot(statsdf, aes(x=utility, color=chain_step_id)) + 
#   geom_density( alpha=.1) +
#   scale_color_discrete() +
#   facet_wrap(~id_from)+ theme_bw()



# library(parallel)
# library(foreach)


statsdf %>% mutate(actor=id_from) %>% group_by(actor) %>% mutate(is_ego=id_from==actor)
for (i in 1:length(self$M)) {
  actor_id <- i
  # statsi <- statsdf %>% filter(id_from==i)
  # statsi$
  adfi <- self$actor_util_df %>% filter(actor_id == actor_id )
  
  adfi$is_ego <-  (adfi$actor_id == actor_id)
}



##########################################


env1$bi_env_arr



# Create a function to plot the network at each time step  
plot_network <- function(data, time_step) {  
  # Create a network object  
  net <- network(data[data$time == time_step, ], directed = TRUE)  
  
  # Plot the network  
  plot(net, displaylabels = TRUE, main = paste("Time Step:", time_step))  
}  

# Create a list of plots for each time step  
plots <- lapply(1:time_steps, function(t) {  
  plot_network(network_data, t)  
})  

library(gganimate)
library(ggraph)
library(network)

# Combine plots into a single animation  
animation <- ggplot() +  
  geom_edge_link(data = env1$chain_stats, aes(x = id_from, y = id_to)) +  
  geom_node_point(data = env1$chain_stats, aes(x = id_from, y = id_to), size = 5) +  
  labs(title = 'Network Animation: Time Step {chain_step_id}', x = 'From', y = 'To') +  
  transition_time(time) +  
  ease_aes('linear')  

# Render the animation  
animate(animation, nframes = 100, fps = 10)  

#






# 
# #################

actaggr %>%
  ggplot(aes(x=n_decision, y=utility_mean, fill=strategy, color=strategy, shape=strategy)) +
  geom_point(alpha = 0.6, size=3) +  # Scatter points
  # stat_ellipse(aes(color = strategy), alpha = 0.2, level = 0.95) +  # Correlation ellipses
  geom_density_2d(aes(color = strategy), size = .5) +  # Density contours
  geom_hline(data=actaggr%>%group_by(chain_h1,strategy,aggr_f)%>%dplyr::summarize(mean=mean(utility_mean)),
             aes(yintercept = mean, color=strategy, linetype=strategy),  linewidth=1.1) +
  facet_grid(aggr_f~chain_h1) +
  theme_bw()


actaggr %>%
  ggplot(aes(x=entropy, y=utility_mean, fill=strategy, color=strategy, shape=strategy)) +
  geom_point(alpha = 0.6, size=3) +  # Scatter points
  # stat_ellipse(aes(color = strategy), alpha = 0.2, level = 0.95) +  # Correlation ellipses
  geom_density_2d(aes(color = strategy), size = 0.5) +  # Density contours
  geom_hline(data=actaggr%>%group_by(chain_h1,strategy,aggr_f)%>%dplyr::summarize(mean=mean(utility_mean)),
             aes(yintercept = mean, color=strategy, linetype=strategy),  linewidth=1.1) +
  facet_grid(aggr_f~chain_h1) +
  theme_bw()

actaggr %>%
  ggplot(aes(x=forbearance, y=utility_mean, fill=strategy, color=strategy, shape=strategy)) +
  geom_point(alpha = 0.6, size=3) +  # Scatter points
  # stat_ellipse(aes(color = strategy), alpha = 0.2, level = 0.95) +  # Correlation ellipses
  geom_density_2d(aes(color = strategy), size = 0.5) +  # Density contours
  geom_hline(data=actaggr%>%group_by(chain_h1,strategy,aggr_f)%>%dplyr::summarize(mean=mean(utility_mean)),
             aes(yintercept = mean, color=strategy, linetype=strategy),  linewidth=1.1) +
  facet_grid(aggr_f~chain_h1) +
  theme_bw()


# plot(x=act$step_precedence_var, y=act$utility_mean)

# act %>% ggplot(aes(x=step_precedence_var, y=utility_mean, fill=strategy, color=strategy)) +
#   geom_density(alpha=.2, linewidth=1) +
#   # geom_histogram(alpha=.1, position = 'jitter') +
#   geom_vline(data=act%>%group_by(chain_h1,strategy)%>%dplyr::summarize(mean=mean(step_precedence_var)), 
#              aes(xintercept = mean, color=strategy, linetype=strategy),  linewidth=1.1) +
#   facet_wrap(~chain_h1) +
#   theme_bw()

  # %>% left_join(
#   env1$actor_util_df %>% group_by(actor_id) %>% summarize(utility_mean=mean(utility, na.rm=T), utility_sd=sd(utility,na.rm=T), utility_max=max(utility, na.rm=T)),
#   x=actor_id, y=actor_id
# )
# act

env1$actor

names(act)

names(env1$actor_stats_df)

actj <- merge(x=act, y=env1$actor_stats_df, x='actor_id', y='actor_id', all=T)


reg.cols <- c('n_decision',
              'forbearance',
              'n_comps_explored',
              'exploitation_rate',
              'entropy',
              'step_precedence_sum',
              'step_precedence_mean',
              'step_precedence_var')
reg.cols.strat <- c('strategy', reg.cols)
reg.cols.actor <- c('id_from', reg.cols.strat)

actcor <- act[,reg.cols]
actreg <- act[,reg.cols.strat]
act2 <- act[,reg.cols.actor]

round( cor(actcor), 3)

folm <- lm(forbearance ~ . , data=actreg)
summary(folm)

library(sandwich)
library(lmtest)
# vcovHC(folm)
coeftest(folm, vcov. = vcovHC)

# explorelm <- lm(n_comps_explored ~ . , data=actreg)
# summary(explorelm)
# 
# exploitlm <- lm(exploitation_rate ~ . , data=actreg)
# summary(exploitlm)

plot(y=actreg$forbearance, x= actreg$entropy  ); abline(lm(forbearance ~ entropy, data=actreg))

# actreg %>% tidyr::pivot_longer(cols = everything(), names_to = 'variable', values_to='value')
###############################
#######


















tiechdf <- env1$chain_stats %>% filter( !stability )

tiechdf[1:5,]

bi_env_mat <- env1$bi_env_arr[,,1]


for (i in 2:length()) {
  # i <- 2
  # bi_env_mat_i <- env1$bi_env_arr[,,1:5]
  
  mstep <- tiechdf[i,]
  
  ( bi_env_mat <- env1$bi_env_arr[,,i] )
  
  ## update bipartite environment matrix for one step (self$toggle one dyad)
  bi_env_mat <- self$toggleBiMat(bi_env_mat, mstep$id_from,  (mstep$id_to - self$M)  )
  
  
}

# i <- 2
# # bi_env_mat_i <- env1$bi_env_arr[,,1:5]
# 
# mstep <- tiechdf[i,]
# 
# ( bi_env_mat <- env1$bi_env_arr[,,i] )
# 
# ## update bipartite environment matrix for one step (self$toggle one dyad)
# bi_env_mat <- self$toggleBiMat(bi_env_mat, mstep$id_from,  (mstep$id_to - self$M)  )
# 

###################################
























library(ggridges) 
library(dplyr)
library(forcats)

# ## https://wilkelab.org/ggridges/articles/gallery.html
# ##  EXAMPLE OF density ridges like a movie of evolving density
# # text labels
# ggplot(lincoln_weather, aes(x = `Mean Wind Speed[MPH]`, y = `Month`, fill = ..x..)) +
#   geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
#   scale_y_discrete(expand = c(0, 0, 0.2, 0)) + # add top margin
#   scale_x_continuous(expand = c(0, 0)) +
#   theme(legend.position = "none") +
#   scale_fill_viridis(option = "C", direction = -1) +
#   labs(title = "Daily Mean Wind Speeds in Lincoln, NE in 2016") +
#   xlab("Mean Wind Speed (MPH)")
# 
# ##  EXAMPLE OF multi-group (multicolor) density ridges like a movie of evolving group density



inputs <- read_excel(file.path(dir_proj_base, grid.input.file))

##
environ_seed_params <- unique(inputs$name[inputs$type %in% c('environment','seed')])
eff_cov_params      <- unique(inputs$name[inputs$type %in% c('effects','coCovars','varCovars','dyCovars')])
eff_struct_params   <- unique(inputs$name[inputs$type %in% c('effects')])
eff_component_params<- unique(inputs$name[inputs$type %in% c('coCovars','varCovars','dyCovars') & grepl('self\\$compon.+',inputs$interaction1) ])
eff_strategy_params <- unique(inputs$name[inputs$type %in% c('coCovars','varCovars','dyCovars') & grepl('self\\$strat.+',inputs$interaction1) ])




# get_actor_strategies = function() {
#   coveffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
#   covDvTypes <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$interaction1)
#   stratDV_ids     <- grep('self\\$strat_\\d{1,2}_coCovar', covDvTypes ) ## ex: "self$strat_1_coCovar" 
#   strat_effs <- coveffs[ stratDV_ids ]
#   nstrat <- length(strat_effs)
#   if (nstrat == 1) return(self$strat_1_coCovar)
#   if (nstrat == 2) return(paste(self$strat_1_coCovar, 
#                                 self$strat_2_coCovar, sep='_'))
#   if (nstrat == 3) return(paste(self$strat_1_coCovar, 
#                                 self$strat_2_coCovar, 
#                                 self$strat_3_coCovar, sep='_'))
#   if (nstrat == 4) return(paste(self$strat_1_coCovar, 
#                                 self$strat_2_coCovar, 
#                                 self$strat_3_coCovar, 
#                                 self$strat_4_coCovar, sep='_'))
#   if (nstrat == 5) return(paste(self$strat_1_coCovar, 
#                                 self$strat_2_coCovar, 
#                                 self$strat_3_coCovar, 
#                                 self$strat_4_coCovar, 
#                                 self$strat_5_coCovar, sep='_'))
#   # if (nstrat == 1) return(paste(strat_effs, self$strat_1_coCovar, sep='= '))
#   # if (nstrat == 2) return(paste(strat_effs, paste(self$strat_1_coCovar, 
#   #                                                 self$strat_2_coCovar, sep='_'), sep='= '))
#   # if (nstrat == 3) return(paste(strat_effs, paste(self$strat_1_coCovar, 
#   #                                                 self$strat_2_coCovar, 
#   #                                                 self$strat_3_coCovar, sep='_'), sep='= '))
#   # if (nstrat == 4) return(paste(strat_effs, paste(self$strat_1_coCovar, 
#   #                                                 self$strat_2_coCovar, 
#   #                                                 self$strat_3_coCovar, 
#   #                                                 self$strat_4_coCovar, sep='_'), sep='= '))
#   # if (nstrat == 5) return(paste(strat_effs, paste(self$strat_1_coCovar, 
#   #                                                 self$strat_2_coCovar, 
#   #                                                 self$strat_3_coCovar, 
#   #                                                 self$strat_4_coCovar, 
#   #                                                 self$strat_5_coCovar, sep='_'), sep='= '))
# }




actor_ids=c() 
wave_ids=c()
thin_factor=1
thin_wave_factor=1
smooth_method='loess'  ##"lm", "glm", "gam", "loess","auto"
show_utility_points=T
scale_utility=TRUE
return_plot=FALSE
plot_file=NA
plot_dir=NA


actor_strat <- as.factor( self$get_actor_strategies() )
#
nstep <- sum(!self$chain_stats$stability)
#
coveffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
covparams <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$parameter)
covfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)ifelse(x$fix,'','(var)'))
# stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)substr(as.character(x$fix),1,1))
#
structeffs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$effect)
structparams <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$parameter)
structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)ifelse(x$fix,'','(var)'))
# structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)substr(as.character(x$fix),1,1))
covDvTypes <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$interaction1)
componentDV_ids <- grep('self\\$component_\\d{1,2}_coCovar', covDvTypes) ## ex: "self$component_1_coCovar"
stratDV_ids     <- grep('self\\$strat_\\d{1,2}_coCovar', covDvTypes ) ## ex: "self$strat_1_coCovar" 
#
compoeffs   <- coveffs[ componentDV_ids ]
compoparams <- covparams[ componentDV_ids ]
compofixs   <- covfixs[ componentDV_ids ]
#
strateffs   <- coveffs[ stratDV_ids ]
stratparams <- covparams[ stratDV_ids ]
stratfixs   <- covfixs[ stratDV_ids ]
#
actor_component_period <- self$M * self$N
#
transition_pds <- 2
density_ridges_rel_min_height = 1e-07  ## prevents density ridges colored lines from covering full x-axis (clarifies group separation)
## Compare 2 actors utilty
dat <- self$actor_wave_util %>% 
  filter(chain_step_id %% thin_factor == 0) %>% 
  filter(wave_id %% thin_wave_factor == 0 ) %>% 
  mutate(
    strategy = actor_strat[ actor_id ],
    chain_below_med =  chain_step_id < median(chain_step_id),
    actor_component_period = 1 + floor( chain_step_id / actor_component_period )
  ) %>% 
  mutate(
    stabilization_summary_period = ifelse(actor_component_period <= transition_pds, 
                                          actor_component_period, 
                                          sprintf('%s+\n(%s-%s)',transition_pds+1,transition_pds+1,max(actor_component_period)))
  )
# dat$chain_half <- factor(ifelse(dat$chain_below_med, '1st Half', '2nd Half'))

#
util_lab <- 'Actor Utility'
if(scale_utility) {
  util_sc <- scale(dat$utility)
  util_lab <- sprintf('Actor Utility\n(Standardized Center = %.2f; Scale = %.2f)',
                      attr(util_sc, 'scaled:center'), 
                      attr(util_sc, 'scaled:scale'))
  if (!all(dat$utility == 0))
    dat <- dat %>% mutate(utility = c(scale(utility)))
}
#
density_rng <- range(dat$utility, na.rm=T)
density_absdiff_scale <- abs(diff(density_rng)) * 0.15
util_lim <- c(density_rng[1] - density_absdiff_scale,  density_rng[2] + density_absdiff_scale)
#
point_size <- 10 / log( nstep )
point_alpha <- min( 1,  1/log10( nstep ) )
#
if(length(actor_ids))
  dat <- dat %>% filter(actor_id %in% actor_ids)
if(length(wave_ids))
  dat <- dat %>% filter(wave_id %in% wave_ids)
dat_acp_stabil_means <- dat %>% group_by(stabilization_summary_period, strategy) %>% 
  dplyr::summarize(mean=mean(utility, na.rm=T)) %>%
  mutate(PeriodFct = fct_rev(as.factor(stabilization_summary_period)))
#
# plt <- ggplot(dat, aes(x=chain_step_id, y=utility)) + 
#   geom_hline(data=dat_acp_means, aes(yintercept=mean), linetype=3, col='black' ) +
#   facet_grid(wave_id ~ .) 


##==============================================
strat_legend_title <- sprintf("Strategy (%s) :  ", paste(strateffs, collapse = '_'))
strat_break <- levels(actor_strat) 
strat_labs <- sapply(1:length(levels(actor_strat)), function(i) {
  a <- levels(actor_strat)[i]
  names(a) <- a ## # names(a) <- sprintf('%s: %s', i, a)
  return(a)
}) 
dat_dens_rigde <- dat %>%
  mutate(PeriodFct = fct_rev(as.factor(stabilization_summary_period))) 
  # mutate(PeriodFct = as.factor(stabilization_summary_period))

group_dens_means <- dat_dens_rigde %>% ungroup() %>% group_by(strategy) %>% 
  dplyr::summarize(n=n(),mean=mean(utility,na.rm=T))

# dat_ridge_lines_density <-  dat_dens_rigde %>% ungroup() %>%
#   group_by(stabilization_summary_period, strategy) %>%
#   summarise(density = density(utility), 
#             strategy=paste(unique(strategy),collapse = "|"),
#             stabilization_summary_period=paste(unique(stabilization_summary_period),collapse = "|") ) #%>%
#   # unnest_wider(density) %>%
#   # rename(density_value = y, group = y) %>%
#   # unnest(cols = c(x, y)) %>%
#   # mutate(y = as.numeric(factor(group))) # Assign numeric y-values for alignment


plt.dr <- ggplot(dat_dens_rigde, aes(y = PeriodFct, x = utility, color=strategy, fill=strategy)) +
  # geom_density_ridges(stat = "binline", alpha=.15,
  #                     bins = 90, draw_baseline = FALSE) +
  ######
  stat_density_ridges(aes(point_color = strategy, point_fill = strategy, point_shape = strategy),
                      quantile_lines = TRUE, alpha = .3, rel_min_height = density_ridges_rel_min_height,
                      point_size=.4,
                      jittered_points = T, 
                      position = position_raincloud(adjust_vlines = FALSE, ygap = -.1, height = .15),# "raincloud",
                      # position = position_points_sina(rel_min = 0.1, rel_max = 0.9, seed = NULL),
                      quantiles = c(0.5), linewidth=.75 ) +
  # stat_density_ridges(aes(point_color = strategy, point_fill = strategy, point_shape = strategy),
  #                     quantile_lines = TRUE, alpha = 0.25,
  #                     point_size=.2,
  #                     quantiles = c(0.5), linewidth=.7 ) +
  ######
  ## geom_density_ridges(
  ##   aes(x = utility , fill = strategy),
  ##   # jittered_points = TRUE, scale = .95, rel_min_height = 0.2,
  ##   # point_shape = "|", point_size=1,
  ##   # position =  position_points_jitter(height = 0),
  ##   alpha = .4, #color = "darkgray",
  ##   from = min(dat$utility)-0.5, to = max(dat$utility)+0.5
  ## ) +
  # geom_line(data = dat_ridge_lines_density, aes(x = x, y = y + ..dat_ridge_lines_density..), 
  #           inherit.aes = FALSE, size = 0.8) +
  ###
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_cyclical(
    breaks = strat_break,
    labels = strat_labs,
    values = hue_pal()(length(levels(actor_strat))), # c("#ff0000", "#0000ff", "#ff8080", "#8080ff"),
    # name = sprintf("strategy\n(%s)", paste(strateffs, collapse = '_')), 
    guide = "legend"
  ) +
  labs(
    x = util_lab,
    y = sprintf(" Time Period \n(Actor-Component-Period = %s decision steps)", actor_component_period ),
    title = "Actor Utility Stabilization Paths",
    subtitle = sprintf("(Decision Chain Iterations: %s )", nstep), 
    # caption = "(plot design adapted from source: Marc Belzunces (@marcbeldata))",
    color = strat_legend_title,
    fill = strat_legend_title,
    point_color = strat_legend_title,
    point_fill = strat_legend_title,
    point_shape = strat_legend_title
  ) +
  geom_vline(xintercept = 0, linetype=1) +
  # geom_vline(data = group_dens_means, 
  #            aes(xintercept = mean, color=strategy),
  #            linetype=3, linewidth=1.3) +
  # geom_text(data = group_dens_means, 
  #           aes(x = mean, y = 1+length(unique(dat_dens_rigde$stabilization_summary_period)), label = round(mean, 2), color=strategy), 
  #           inherit.aes = FALSE, size = 5, nudge_x=-.1, nudge_y=.35  ) +
  coord_cartesian(clip = "off") +
  theme_ridges(grid = T, center=T) + theme(legend.position = 'bottom')
print(plt.dr) 

# #
# plt.dr.dat <- ggplot_build(plt.dr)$data[[1]]
# plt.dr.b$data
##==============================================

self$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T)

# plyr::count(plt.dr.dat$group)
# nstrats <- length(levels(actor_strat))
# npds <- length(levels(as.factor(dat$stabilization_summary_period)))
# 
# getGrpId <- function(strat, pd, nstrats, npds) (strat-1)*npds + pd #(npds-pd+1)
# getGrpId(strat=2, pd=1, nstrats=nstrats, npds=npds)
# 
# getPdFromPlotDf <- function(strat, pd, nstrats, npds)  
# getStratFromPlotDf <- function(strat, pd, nstrats, npds)  
# 
# plt.dr.dat %>% 
#   mutate(group_fct = as.factor(group)) %>%
#   filter(!is.na(density)) %>%   
#   group_by(group = group_fct) %>% 
#   dplyr::summarize(
#     n=n(), x_min= min(x, na.rm=T), x_max=max(x, na.rm=T), y_pos=first(y)
#   )
# plt.dr.dat$strategy <- NA
# plt.dr.dat$pd <- NA
# for (i in 1:nrow(plt.dr.dat)) {
#   plt.dr.dat$strategy[i] <- plt.dr.dat$group[i]
#   plt.dr.dat$pd[i]       <- 
# }


head(plt.dr.dat)


# 
# line_data <- dat_dens_rigde %>%
#   group_by(group, as.factor(y) ) %>%
#   summarise(
#     x_min = min(x[!is.na(density)]),  # Minimum x with density > 0
#     x_max = max(x[!is.na(density)]),  # Maximum x with density > 0
#     y_pos = first(y)                 # y-position for the group
#   )
# # Add connecting lines to the plot
# (plt.dr.lined <- plt.dr +
#   geom_segment(data = line_data,
#                aes(x = x_min, xend = x_max, y = y_pos, yend = y_pos),
#                inherit.aes = FALSE, size = 0.8)
# )




###






















##---------------------------------------------------------------------------------------------
##
search_rsiena_multiwave_plot_actor_utility_strategy_summary = function(actor_ids=c(), 
                                                                       wave_ids=c(),
                                                                       thin_factor=1, 
                                                                       thin_wave_factor=1,
                                                                       smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                                       show_utility_points=T,
                                                                       scale_utility=TRUE,
                                                                       return_plot=FALSE,
                                                                       plot_file=NA, plot_dir=NA
) {
  ## actor strategy
  if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
    stop("Actor Strategy self$strat_1_coCovar not set.")
  actor_strat <- as.factor( self$get_actor_strategies() )
  #
  nstep <- sum(!self$chain_stats$stability)
  #
  strateffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
  stratparams <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$parameter)
  stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)ifelse(x$fix,'','(var)'))
  # stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)substr(as.character(x$fix),1,1))
  #
  structeffs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$effect)
  structparams <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$parameter)
  structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)ifelse(x$fix,'','(var)'))
  # structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)substr(as.character(x$fix),1,1))
  covDvTypes <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$interaction1)
  componentDV_ids <- grep('self\\$component_\\d{1,2}_coCovar', covDvTypes) ## ex: "self$component_1_coCovar"
  stratDV_ids     <- grep('self\\$strat_\\d{1,2}_coCovar', covDvTypes ) ## ex: "self$strat_1_coCovar" 
  #
  actor_component_period <- self$M * self$N
  ## Compare 2 actors utilty
  dat <- self$actor_wave_util %>% 
    filter(chain_step_id %% thin_factor == 0) %>% 
    filter(wave_id %% thin_wave_factor == 0 ) %>% 
    mutate(
      strategy = actor_strat[ actor_id ],
      chain_below_med =  chain_step_id < median(chain_step_id),
      actor_component_period = 1 + floor( chain_step_id / actor_component_period )
    )
  dat$chain_half <- factor(ifelse(dat$chain_below_med, '1st Half', '2nd Half'))

  #
  util_lab <- 'Actor Utility'
  if(scale_utility) {
    util_sc <- scale(dat$utility)
    util_lab <- sprintf('Actor Utility\n(Standardized Center = %.2f; Scale = %.2f)',
                        attr(util_sc, 'scaled:center'), 
                        attr(util_sc, 'scaled:scale'))
    if (!all(dat$utility == 0))
      dat <- dat %>% mutate(utility = c(scale(utility)))
  }
  #
  density_rng <- range(dat$utility, na.rm=T)
  density_absdiff_scale <- abs(diff(density_rng)) * 0.15
  util_lim <- c(density_rng[1] - density_absdiff_scale,  density_rng[2] + density_absdiff_scale)
  #
  point_size <- 10 / log( nstep )
  point_alpha <- min( 1,  1/log10( nstep ) )
  #
  if(length(actor_ids))
    dat <- dat %>% filter(actor_id %in% actor_ids)
  if(length(wave_ids))
    dat <- dat %>% filter(wave_id %in% wave_ids)
  dat_wave_means <- dat %>% group_by(wave_id) %>% 
    dplyr::summarize(mean=mean(utility, na.rm=T))
  #
  plt <- ggplot(dat, aes(x=chain_step_id, y=utility)) + 
    geom_hline(data=dat_wave_means, aes(yintercept=mean), linetype=3, col='black' ) +
    facet_grid(wave_id ~ .) 
  if(show_utility_points)
    plt <- plt + geom_point(aes(color=strategy), alpha=point_alpha, shape=1, size=point_size)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
  if(self$exists(smooth_method))
    plt <- plt + geom_smooth(aes(linetype=actor_id, color=strategy, fill=strategy), method = smooth_method, linewidth=1, alpha=.09)
  #
  plt <- plt + geom_vline(xintercept = which(1:max(dat$chain_step_id) %% actor_component_period==0) , linetype = 3, color='gray')
  # 
  plt <- plt + theme_bw() + 
    ylim(util_lim) + 
    ylab(util_lab) +
    xlab('Actor Decision Chain Ministep') +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.key.width = unit(0.8, "cm"),  # Adjust legend key width
      legend.spacing.x = unit(0.3, "cm")   # Adjust spacing between keys
    )
  
  ACTORS     <- sienaNodeSet(self$M, nodeSetName="ACTORS")
  COMPONENTS <- sienaNodeSet(self$N, nodeSetName="COMPONENTS")
  # strat_mat_varCovar <- varCovar(strat_eff$x, nodeSet = 'ACTORS')
  if (! 'coCovars' %in% names(self$config_structure_model$dv_bipartite) ) {
    rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV), nodeSets = list(ACTORS, COMPONENTS))
    return(rsiena_data)
  }
  
  plt <- plt + ggtitle(sprintf('Environment: Actors (M) = %s, Components (N) = %s, Init.Prob. = %.2f\nActor Strategy:  %s\nComponent Payoff:  %s\nStructure:  %s', 
                               self$M, self$N, self$BI_PROB,
                               paste( paste(paste(strateffs[stratDV_ids], stratparams[stratDV_ids], sep='= '), stratfixs[stratDV_ids], sep='' ), collapse = ';  '),
                               paste( paste(paste(strateffs[componentDV_ids], stratparams[componentDV_ids], sep='= '), stratfixs[componentDV_ids], sep=''), collapse = ';  '),
                               paste( paste(paste(structeffs, structparams, sep='= '), structfixs, sep=''), collapse = ';  ')
  ))
  plt <- plt +  guides(color = guide_legend(nrow = 1))
  
  #### Density
  stratmeans <- dat %>% group_by(strategy, wave_id) %>% 
    dplyr::summarize(n=n(), mean=mean(utility, na.rm=T), sd=sd(utility, na.rm=T))
  ## Actor density fact plots comparing H1 to H2 utility distribution
  plt2 <- ggplot(dat, aes(x=utility, color=strategy, fill=strategy)) + ##linetype=chain_half
    geom_density(alpha=.1, linewidth=1)  +
    # geom_histogram(alpha=.1, position = 'dodge') +
    geom_vline(data = stratmeans, aes(xintercept = mean, color=strategy), linetype=2, linewidth=.9) +
    geom_vline(data = dat_wave_means,  aes(xintercept=mean), linetype=3, col='black' ) +
    labs(y='', x='') +
    # xlim(c(ggplot_build(plt)$layout$panel_params[[1]]$y.range)) + 
    xlim(util_lim) +
    coord_flip() +
    facet_grid(wave_id ~ .) +
    ylab('Actor Utility Density') +
    # geom_vline(xintercept = 0, linetype=2, linewidth=.9, color='gray')+
    theme_bw() + theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = 'none', 
      plot.margin=unit(c(5.5, 5.5, 5.5, -23), 'pt'),
      axis.text.y = element_blank(),
      axis.ticks.y=element_blank()#,
      # axis.text.x = element_blank(),
      # axis.ticks.x=element_blank()#,
    ) + ggtitle('\n\n\n')
  
  # combined_plot <- plot_grid(plt ,#+ theme(panel.spacing = unit(0, "lines")), 
  #                            plt2 ,# + theme(panel.spacing = unit(0, "lines")), 
  #                            ncol = 2, rel_widths = c(4, 1)) 
  # print(combined_plot)
  
  combined_plot <- ggarrange(
    plt, plt2, 
    ncol = 2, 
    widths = c(4.1,0.9), # Adjust column widths
    common.legend = TRUE, # Share a common legend if needed
    legend = "bottom"#,     # Place legend at the bottom
    # align = 'h'
  ) 
  # print(combined_plot)
  #
  # print(plt)
  #
  if(!is.na(plot_file))
    ggsave(file = file.path(ifelse(is.na(plot_dir),getwd(),plot_dir), sprintf("%s_%s.png", self$config_environ_params$name, plot_file)), 
           combined_plot, 
           width = 10, height = 8, units = 'in', dpi = 600)
  
  if(return_plot)
    return(combined_plot)
}
##---------------------------------------------------------------------------------------------











































# outlist <- list()
# Ns <- 2 * 3^c(1:5)  ## 6 12 24 48 96
# for (i in 1:length(Ns)) 
# {
  
  
  ##****************************************##
  ## I. SIM SETUP 
  ##****************************************##
  ##----------------------
  ## 1. Environment:  determines what types of entities and strategies are permissible
  ##----------------------
  environ_params <- list(
    M = 3,        ## Actors
    N = 12,       ## Components
    BI_PROB = 0, ## Environmental Density (DGP hyperparameter)
    component_matrix_start = 'rand', ##**TODO** Implement: 'rand','modular','semi-modular',...
    rand_seed = 123,
    visualize_init = F,
    name = sprintf('_testSAOMNK_4_loopN%s', 12 )#,
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
  
#   
#   
#   outlist[[ sprintf('N%s',Ns[i]) ]] <- list(env1=env1, env2=env2, N=Ns[i])
# }
# 
# 
# 
# saveRDS(outlist, 
#         file = file.path(dir_data, sprintf('%s_%s.rds',environ_params$name,as.character(as.integer(Sys.time())))) )
# ## end




