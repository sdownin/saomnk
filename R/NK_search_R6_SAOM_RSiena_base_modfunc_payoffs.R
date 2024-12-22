##
##
# #
# get_rsiena_model_snapshots = function(n=4, plot_save = TRUE, return_list=TRUE) {
#   #
#   snapshot_sim_ids <- seq(0,1,length.out= n) * length(self$rsiena_model$sims)
#   snapshot_sim_ids[1] <- 1
#   #
#   bi_env_dv_id <- which( names(self$structure_model) == 'dv_bipartite' )
#   wave_id <- 1
#   ##
#   MplusN <- self$M + self$N
#   ##
#   snapshots <- list()
#   for (sim_id in snapshot_sim_ids) {
#     ##
#     el_bi_env <- self$rsiena_model$sims[[ sim_id ]][[ wave_id ]][[ bi_env_dv_id ]]$`1`
#     ## update numbering of second mode (the comonent integer names shift upward by the number of actors)
#     el_bi_env[,2] <- el_bi_env[,2] + self$M
#     #############
#     ## Bipartite matrix space (N+M by N+M)
#     ## Undirected --> Upper right rectangle of full bipartite matrix
#     ##   M N
#     ## M[0,X], for X in [0,1]
#     ## N[0,0]
#     bi_env_mat_sp <- sparseMatrix(i = el_bi_env[,1],
#                                   j = el_bi_env[,2],
#                                   x = el_bi_env[,3],
#                                   dims = c(MplusN, MplusN))
#     
#       
#     bi_env_mat_new  <- as.matrix(bi_env_mat_sp)[ 1:self$M, (self$M+1):(MplusN) ]
#     
#     self$plot_bipartite_system_from_mat(bi_env_mat_new, sim_id, plot_save = TRUE) 
#     
#     snapshots[[ sprintf('sim%d',sim_id) ]] <- bi_env_mat_new
#   }
#  
#   ##
#   if(return_list) return( snapshots )
#   
# },

#
# chain_stats = data.frame(),  ## 
# sims_stats = list(
#   'K_A'=list(), ##**TODO** Actor social space [degree = central position in social network]
#   'K_B'=list(), #  Bipartite space  [ degree = resource/affiliation density ]
#   'K_C'=list()  #  Component space  [ degree = interdependence/structuration : epistatic interactions ]
# ),
# #
# self$chain_stats <-

# ##
# search_rsiena_process_ministep_chain_PREVIOUS = function() {
#   if (is.null(self$rsiena_model$chain)) {
#     stop("Chain not available. Ensure returnChains=TRUE was set in the siena07 call.")
#   }
# 
#   depvar <- 1
#   period <- 1
# 
#   .getTieChange <- function(dvName, ministep) {
#     if(dvName == 'self$bipartite_rsienaDV') {
#       return(ifelse(as.numeric(ministep[[5]])==self$N), FALSE, TRUE)
#     } else if (dvName %in% c('self$social_rsienaDV','self$search_rsienaDV')) {
#       return(ifelse(ministep[[4]]==ministep[[5]], FALSE, TRUE))
#     } else {
#       stop(sprintf('dvName %s not implemented in .getTieChange()', dvName))
#     }
#   }
# 
#   ## Bipatite network chain --> value Alter=m means that no change has occurred.
#   ##** "chain[[run]][[depvar]][[period]][[ministep]]"**
#   self$chain_stats <- ldply(seq_along(self$rsiena_model$chain), function(run_id) {
#     run <- self$rsiena_model$chain[[ run_id ]]
#     period_ministeps <- run[[depvar]][[period]]
#     n_ministeps <- length(period_ministeps)
#     dfsteps <- data.frame()
#     if (  n_ministeps ) {
#       for (i in 1:n_ministeps) {
#         ministep <- period_ministeps[[ i ]]
#         dfsteps <- rbind(dfsteps,  data.frame(
#           `iteration_id` = run_id,
#           `ministep_id` = i,
#           `_period_mu` = attr(period_ministeps, 'mu'),
#           `_period_sigma2` = attr(period_ministeps, 'sigma2'),
#           `_period_finalReciprocalRate` = attr(period_ministeps, 'finalReciprocalRate'),
#           `_period_ministep_cnt` = n_ministeps,
#           `_tie_change` = ifelse(),
#           `__aspect` = ministep[[1]],  ## aspect = network or behavior function
#           `__Var`    = ministep[[2]],  ## Var  (same as aspect, but binary [0=Network,1=Behavior])
#           `__VarName` = ministep[[3]],      ## VarName (same as aspect|Var, but name character)
#           `__id_from` = as.numeric( ministep[[4]] ) + 1, ## Ego (+1 updates C++ 0-index to R 1-index)
#           `__id_to` = as.numeric( ministep[[5]] ) + (1 + self$M), ## Alter  (+1 updates C++ 0-index to R 1-index)
#           `__difference` = ministep[[6]],  ## Difference (0 == no change)
#           `__reciprocal_rate` = ministep[[7]], ## reciprocal rate
#           `__LogOptionSetProb` = ministep[[8]],  ##**TODO** CHANGE NAME:  LogOptionSetProb
#           `__LogChoiceProb` = ministep[[9]],   ##**TODO** CHANGE NAME:  LogChoiceProb
#           `__diagnoal` = ifelse(is.null(ministep[[10]]), NA, ministep[[10]]),  ## Diagonal
#           `__stability` = ifelse(is.null(ministep[[11]]), NA, ministep[[11]]), ## Stability: TRUE=no change; FALSE=change
#           `__bool1` = ministep[[12]], ##??
#           `__bool2` = ministep[[13]]  ##??
#         ))
#       }
#     }
# 
#     return(dfsteps)
# 
#   })
# 
# 
#   print(rbind(head(self$chain_stats),tail(self$chain_stats)))
#   print(dim(self$chain_stats))
# },

##

# Load necessary libraries
library(R6)
library(Matrix)
library(network)
library(igraph)
library(visNetwork)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(ggraph)
library(ggpubr)
library(RSiena)
library(texreg)
library(grid) # Load the necessary library for grid.text
library(plyr)
library(dplyr)
library(tidyr)


##
##
##
get_jaccard_index <- function(m0, m1) {
  diffvec <-  c(m1 - m0)   ## new m1 - old m0
  cnt_maintain <- sum( m1 * m0 )
  cnt_change <- sum( diffvec != 0 )  ## sum = count true cases (added + dropped)
  return( cnt_maintain / (cnt_maintain + cnt_change) )
}

##
exists <- function(x){
  return(!is.null(x) && !is.na(x) && !is.nan(x))
}

# Define the one-mode network matrix toggle function:
toggle <- function(m,i,j){
  if (i != j) {
    # print(sprintf('test i %s != j %s',i,j))
    m[i,j] <-  ( 1 - m[i,j] )
  }
  return(m)
}

# Define the bipartite (2-mode) network matrix toggle function:
toggleBiMat <- function(m,i,j){
  m[i,j] <-  ( 1 - m[i,j] )
  return(m)
}



##**TODO**
## 1. [list] environ_params     Environmental parameters list
## 2. [list] structure_model    SAOM for environmental structural evolution (network ties: social, search, bipartite)
## 4. [list] < other search specific attributes, like environ shock from regulatory change ? >
##**/**


##
##
##
SaomNkRSienaBiEnv_base <- R6Class( 
  ##
  "SaomNkRSienaBiEnv_base",
  ##
  
  public = list(
    #
    ITERATION = 0,
    TIMESTAMP = NULL,
    SIM_NAME = NULL,
    #
    M = NULL,                # Number of actors
    N = NULL,                # Number of components
    BI_PROB = NULL,          # probability of tie in bipartite network dyads (determines bipartite density)
    #
    config_environ_params = list(),   ## environment of simulation (DGP) **TODO** add shocks or landscape changes **
    config_structure_model = list(),  ## list of SAOM model DVs and effects
    ## config_payoff_formulas = list(),  ## list of formulae for each DV's network statistics predictors to be used in 
    ## K = NULL,             # Avg. degree of component-[actor]-component interactions
    bipartite_igraph = NULL, # Bipartite network object 
    social_igraph = NULL,   # Social space projection (network object)
    search_igraph = NULL, # Search space projection (network object)
    #
    bipartite_matrix = NULL,
    social_matrix = NULL,
    search_matrix = NULL,
    #
    bipartite_matrix_waves = list(), ## list of simulated networks for extended multi-wave simulation; implements dynamic strategy choice/changes
    #
    bipartite_rsienaDV = NULL,
    social_rsienaDV = NULL,
    search_rsienaDV = NULL,
    #
    strat_1_coCovar = NULL,
    strat_2_coCovar = NULL,
    strat_3_coCovar = NULL,
    strat_4_coCovar = NULL,
    strat_5_coCovar = NULL,
    # strat_6_coCovar = NULL,
    strat_mat_varCovar = NULL,
    #
    component_1_coCovar = NULL,
    component_2_coCovar = NULL,
    # component_3_coCovar = NULL,
    #
    plots = list(),
    multiwave_plots = list(),
    #
    fitness_landscape = NULL, # Fitness landscape matrix
    progress_scores = c(),  # Track scores over iterations
    P_change = NULL,          # Probability of changing the component interaction matrix
    #
    chain_stats = NULL,  ##
    #
    actor_stats_df = NULL,
    actor_util_df = NULL,
    actor_util_diff_df = NULL,
    #
    component_stats_df = NULL,
    actor_proj_stats_df = NULL,
    component_proj_stats_df = NULL,
    #
    actor_wave_stats = NULL,
    actor_wave_util = NULL,
    actor_wave_util_diff = NULL,
    #
    component_wave_stats =  NULL,
    #
    actor_proj_wave_stats = list(),
    component_proj_wave_stats = list(),

    ## Coupled Degree Statistics
    K_A_df = NULL, ## degree distribution statistics of Actor projected social network (common components)
    K_B1_df = NULL, ## degree distribution statistics of Bipartite social network - [1]Actors
    K_B2_df = NULL, ## degree distribution statistics of Bipartite social network - [2]Components
    K_C_df = NULL, ## degree distribution statistics of Component projected network (common actors)
    #
    K_wave_A = NULL,
    K_wave_B1 = NULL,
    K_wave_B2 = NULL,
    K_wave_C = NULL,
    ##
    sims_stats = NULL,
    #   'K_A'=list(), ##**TODO** Actor social space [degree = central position in social network]
    #   'K_B'=list(), #  Bipartite space  [ degree = resource/affiliation density ]
    #   'K_C'=list()  #  Component space  [ degree = interdependence/structuration : epistatic interactions ]
    #   ),
    # #
    rsiena_model = NULL,
    rsiena_data = NULL,
    rsiena_effects = NULL,
    rsiena_algorithm = NULL,
    #
    rsiena_model_waves = list(), ## list of models from multiwave simulation; implements dynamic strategy choice/changes
    #
    rsiena_run_seed = NULL,
    rsiena_env_seed = NULL,
    
    # Constructor to initialize the SAOM-NK model
    initialize = function(config_environ_params) {
      cat('\nCALLED _BASE_ INIT\n')
      ## ----- prevent clashes with sna package-----------
      sna_err_check <-tryCatch(expr = { detach('package:sna') }, error=function(e)e )
      ## -------------------------------------------------
      ##**TODO:  LOAD DEPENDENCY FUNCTIONS ETC**
      ##
      self$config_environ_params = config_environ_params
      self$M <- config_environ_params[['M']]
      self$N <- config_environ_params[['N']]
      self$BI_PROB <- config_environ_params[['BI_PROB']]
      # self$P_change <- config_environ_params[['P_change']]
      #
      # self$bipartite_igraph <- self$generate_bipartite_igraph()
      # self$social_network <- self$project_social_space()
      # self$search_landscape <- self$project_search_space()
      default_seed <- 123
      self$rsiena_env_seed <- ifelse(is.null(config_environ_params[['rand_seed']]), 
                                      default_seed,
                                      config_environ_params[['rand_seed']])
      start_bipartite_igraph <- self$random_bipartite_igraph(self$rsiena_env_seed)
      self$set_system_from_bipartite_igraph( start_bipartite_igraph )
      #
      self$TIMESTAMP <- round( as.numeric(Sys.time())*100 )
      self$SIM_NAME <- config_environ_params[['name']]
    },
    
    #
    set_system_from_bipartite_igraph = function(bipartite_igraph) {
      if( ! 'igraph' %in% class(bipartite_igraph)) 
        stop(sprintf('\nbipartite_igraph is not an igraph object; class %s\n', class(bipartite_igraph)))
      if(!length(E(bipartite_igraph)$weight)) 
        E(bipartite_igraph)$weight <- 1
      #
      self$bipartite_igraph <- bipartite_igraph
      self$social_igraph <- self$project_social_space(bipartite_igraph)
      self$search_igraph <- self$project_search_space(bipartite_igraph)
      #
      self$bipartite_matrix  <- igraph::as_biadjacency_matrix(bipartite_igraph,attr = 'weight', sparse = F)
      self$social_matrix <- igraph::as_adjacency_matrix(self$social_igraph, attr = 'weight', sparse = F)
      self$search_matrix <- igraph::as_adjacency_matrix(self$search_igraph, attr = 'weight', sparse = F)
    },
    
    # Generate a random bipartite network object
    random_bipartite_igraph = function(rand_seed = 123) {
      set.seed(rand_seed)  # For reproducibility
      probs <- c( 1 - self$BI_PROB, self$BI_PROB )
      bipartite_matrix <- matrix(sample(0:1, self$M * self$N, replace = TRUE, prob = probs ),
                                 nrow = self$M, ncol = self$N)
      # bipartite_net <- network(bipartite_matrix, bipartite = TRUE, directed = TRUE)
      bipartite_igraph <- igraph::graph_from_biadjacency_matrix(bipartite_matrix, directed = TRUE, mode = 'all')
      return(bipartite_igraph)
    },
    
    ##
    get_bipartite_projections = function(ig_bipartite) {
      projs <- igraph::bipartite_projection(ig_bipartite, multiplicity = T, which = 'both')
      return(projs)    
    },

    # Project the social space (actor network) from the bipartite network
    project_social_space = function(ig_bipartite) {
      projs <- self$get_bipartite_projections(ig_bipartite)
      social_igraph <- projs$proj1
      if(!length(E(social_igraph)$weight)) 
        E(social_igraph)$weight <- 1
      return(social_igraph)
    },

    # Project the search landscape space (component interaction network)
    project_search_space = function(ig_bipartite) {
      projs <- self$get_bipartite_projections(ig_bipartite)
      search_igraph <- projs$proj2
      if(!length(E(search_igraph)$weight)) 
        E(search_igraph)$weight <- 1
      return(search_igraph)
    },
    
    # Update Iteration progress
    increment_sim_iter = function(val = 1) {
      self$ITERATION <- self$ITERATION + val
    }, 
    
    
    ##
    include_rsiena_effect_from_eff_list = function(eff) {
      
      if (eff$effect == 'density') 
      {
        self$rsiena_effects <- includeEffects(self$rsiena_effects,  density, ## get network statistic function from effect name (character)
                                             name = eff$dv_name,  # interaction1 = eff$interaction1,
                                             fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  density, 
                                         name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      }
      else if (eff$effect == 'inPop') 
      {
        self$rsiena_effects <- includeEffects(self$rsiena_effects, inPop, ## get network statistic function from effect name (character)
                                             name = eff$dv_name, # interaction1 = eff$interaction1,
                                             fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  inPop, 
                                          name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      }
      else if (eff$effect == 'outAct') 
      {
        self$rsiena_effects <- includeEffects(self$rsiena_effects, outAct, ## get network statistic function from effect name (character)
                                             name = eff$dv_name, # interaction1 = eff$interaction1,
                                             fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  outAct, 
                                          name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      }
      else if (eff$effect == 'cycle4') 
      {
        self$rsiena_effects <- includeEffects(self$rsiena_effects, cycle4, ## get network statistic function from effect name (character)
                                             name = eff$dv_name, # interaction1 = eff$interaction1,
                                             fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  cycle4, 
                                          name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      }
      else if (eff$effect == 'transTriads') 
      {
        self$rsiena_effects <- includeEffects(self$rsiena_effects,  transTriads, ## get network statistic function from effect name (character)
                                              name = eff$dv_name, # interaction1 = eff$interaction1,
                                              fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  transTriads, 
                                         name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      }
      # else if (eff$effect == 'cycle4ND') 
      # {
      #   self$rsiena_effects <- includeEffects(self$rsiena_effects, cycle4ND, ## get network statistic function from effect name (character)
      #                                         name = eff$dv_name, # interaction1 = eff$interaction1,
      #                                         fix = eff$fix)
      #   self$rsiena_effects <- setEffect(self$rsiena_effects,  cycle4ND, 
      #                                    name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      # }
      else if (eff$effect == 'egoX')
      {
        
        # activity_covar <- varCovar(activity_data)
        # effects <- includeEffects(effects, egoX, interaction1 = "activity_covar")
        
        self$rsiena_effects <- includeEffects(self$rsiena_effects,  egoX, ## get network statistic function from effect name (character)
                                              name = eff$dv_name, 
                                              interaction1 = eff$interaction1,
                                              fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  egoX, 
                                         interaction1 = eff$interaction1,
                                         name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      }
      else if (eff$effect == 'altX')
      {
        self$rsiena_effects <- includeEffects(self$rsiena_effects,  altX, ## get network statistic function from effect name (character)
                                              name = eff$dv_name, 
                                              interaction1 = eff$interaction1,
                                              fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  altX, 
                                         interaction1 = eff$interaction1,
                                         name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      }
      else if (eff$effect == 'outActX')
      {
        self$rsiena_effects <- includeEffects(self$rsiena_effects,  outActX, ## get network statistic function from effect name (character)
                                              name = eff$dv_name, 
                                              interaction1 = eff$interaction1,
                                              fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  outActX, 
                                         interaction1 = eff$interaction1,
                                         name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      }
      else if (eff$effect == 'inPopX')
      {
        self$rsiena_effects <- includeEffects(self$rsiena_effects,  inPopX, ## get network statistic function from effect name (character)
                                              name = eff$dv_name, 
                                              interaction1 = eff$interaction1,
                                              fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  inPopX, 
                                         interaction1 = eff$interaction1,
                                         name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      }
      else 
      {
        ##
        print(eff)
        stop('effect printed above was not included in rsiena_effects.')
      }
      
    },
    
    # ##
    # get_struct_mod_net_stats_list_from_bi_mat = function(bi_env_mat, type='all') {
    #   # eff <- m1$rsiena_effects[m1$rsiena_effects$include, ]
    #   efflist <- self$config_structure_model$dv_bipartite$effects
    #   # set the bipartite environment network matrix
    #   # bi_env_mat <- self$bipartite_matrix
    #   ### empty matrix to hold actor network statistics
    #   df <- data.frame() #nrow=self$M, ncol=length(efflist)
    #   effnames <- sapply(efflist, function(x) x$effect, simplify = T)
    #   effparams <- sapply(efflist, function(x) x$parameter, simplify = T)
    #   mat <- matrix(rep(0, m1$M * length(efflist) ), nrow=self$M, ncol=length(efflist) )
    #   colnames(mat) <- effnames
    #   rownames(mat) <- 1:self$M
    #   #
    #   for (i in 1:length(efflist))
    #   {
    #     eff_name <- efflist[[ i ]]$effect
    #     #
    #     xActorDegree  <- rowSums(bi_env_mat, na.rm=T)
    #     xComponentDegree  <- colSums(bi_env_mat, na.rm=T)
    #     ## network statistics dataframe
    #     #
    #     if (eff_name == 'density' )
    #     {
    #       stat <- c( xActorDegree )
    #     }
    #     else if (eff_name == 'outAct' )
    #     {
    #       stat <- c( xActorDegree^2 )
    #     }
    #     else if (eff_name == 'inPop' )
    #     {
    #       stat <- c( bi_env_mat %*% (xComponentDegree + 1) )
    #     }
    #     # else if (eff_name == 'transTriads' )
    #     # {
    #     #   stat <- 
    #     # }
    #     # else if (eff_name == 'cycle4' )
    #     # {
    #     #   stat <- 
    #     # }
    #     else
    #     {
    #       cat(sprintf('\n\nEffect not yet implemented: `%s\n\n`', eff_name))
    #     }
    #     #
    #     mat[ , i] <- stat
    #     #
    #     df <- rbind(df, data.frame(
    #       statistic = stat, 
    #       actor_id = factor(1:self$M), 
    #       effect_id = factor(i), 
    #       effect_name = effnames[i] 
    #     ))
    #   }
    #   #
    #   if (type %in% c('df','data.frame'))   return(df)
    #   if (type %in% c('mat','matrix'))      return(mat)
    #   if (type %in% c('all','both','list',NA)) return(list(df=df, mat=mat))
    #   cat(sprintf('specified return type %s not found', type))
    # }, 
    
    
    ##
    get_struct_mod_stats_mat_from_bi_mat = function(bi_env_mat, type='all') {
      #
      xActorDegree      <- rowSums(bi_env_mat, na.rm=T)
      xComponentDegree  <- colSums(bi_env_mat, na.rm=T)
      #
      # eff <- m1$rsiena_effects[m1$rsiena_effects$include, ]
      efflist <- self$config_structure_model$dv_bipartite$effects
      #
      efflist <- c(efflist, self$config_structure_model$dv_bipartite$coCovars )
      #
      efflist <- c(efflist, self$config_structure_model$dv_bipartite$varCovars )
      #
      neffs <- length(efflist)
      ### empty matrix to hold actor network statistics
      effnames <- sapply(efflist, function(x) x$effect, simplify = T)
      effparams <- sapply(efflist, function(x) x$parameter, simplify = T)
      #
      mat <- matrix(rep(0, m1$M * neffs ), nrow=self$M, ncol=neffs )
      colnames(mat) <- effnames
      rownames(mat) <- 1:self$M
      #  
      for (i in 1:neffs)
      {
        item <- efflist[[ i ]]
        # print('DEBUG  get_struct_mod_stats_mat_from_bi_mat() ')
        # print(item)
        
        ## network statistics dataframe
        #
        if (item$effect == 'density' )
        {
          stat <- c( xActorDegree )
        }
        else if (item$effect == 'outAct' )
        {
          stat <- c( xActorDegree^2 )
        }
        else if (item$effect == 'inPop' )
        {
          stat <- c( bi_env_mat %*% (xComponentDegree + 1) )
        }
        # else if (item$effect == 'transTriads' )
        # {
        #   stat <- 
        # }
        # else if (item$effect == 'cycle4' )
        # {
        #   stat <- 
        # }
        else if (item$effect == 'egoX')
        {
          covar <- item$x
          checkConform <-  all(
            (  ## A or B
              class(covar) %in% c('array','matrix') & nrow(covar)==self$M
            ) | ( 
              length(covar)==self$M  ## array, matrix; vector
            )
          )
          if( ! checkConform )  
            stop('egoX covar not conformable for multiplication given number of actors')
          stat <- c( covar * xActorDegree ) ##**vector element-wise multiplication by rows of covar matrix, or elements of covar array
            
        }
        else if (item$effect == 'altX')
        {
          covar <- item$x
          checkConform <-  all(
            (  ## A or B
              class(covar) %in% c('array','matrix') & nrow(covar)==self$M
            ) | ( 
              length(covar)==self$N  ## COMPONENT array, for ACTOR statistic 
            )
          )
          if( ! checkConform )
            stop('altX covar not conformable for multiplication given number of components or actors')
          covarComponentMat <- matrix(rep(covar, self$M), nrow=self$M, ncol=self$N, byrow = TRUE)
          stat <- rowSums( covarComponentMat * bi_env_mat, na.rm=T ) ##**vector element-wise multiplication by rows of covar matrix, or elements of covar array
          
        }
        ###--------------------------------
        else if (item$effect == 'outActX') ## interaction1 component_coCovar
        {
          ## N-vector of component covariate
          covar <- item$x 
          # MxN matrix of row-stacked component covariate (repeated for each actor)
          covarComponentMat <- matrix(rep(covar, self$M), nrow=self$M, ncol=self$N, byrow = TRUE)
          ## M-vector of actor's squared sum of component-covariate-weighted component connections (weighted version of the squared degree)
          stat <- rowSums( covarComponentMat * bi_env_mat, na.rm = T)^2
        }
        else if (item$effect == 'inPopX') ## interaction1 strat_coCovar
        {
          covar <- item$x ## M-vector of actor strategy covariate
          ## MxN matrix holding actor strategy covariate as columns stacked for each component
          covarActorMat <- matrix(rep(covar, self$N), nrow=self$M, ncol=self$N,  byrow = FALSE)
          ## N-vector of component weights = sum of actor covariate for the component's connected actors
          component_weights_from_actor_stats <- colSums(covarActorMat * bi_env_mat, na.rm=T)
          ## M-vector of actor sum of it's connected component weights (which are computed as the sum of the connected actor covariates)
          stat <- rowSums( bi_env_mat * component_weights_from_actor_stats, na.rm=T ) ##**vector element-wise multiplication by rows of covar matrix, or elements of covar array
        }
        else
        {
          cat(sprintf('\n\nEffect not yet implemented: `%s\n\n`', eff_name))
        }
        
        #
        mat[ , i] <- stat
        
      }
      
      return(mat)
    }
    
    
    
  )
    
    
)




###############################################################################
###############################################################################

# Define the SAOM_NK_Enhanced class with local search, complex search strategies, network metrics, and visualization
SaomNkRSienaBiEnv <- R6Class( 
  ##
  "SaomNkRSienaBiEnv",
  ##
  inherit = SaomNkRSienaBiEnv_base,
  ##
  
  public = list(
    #

    # Constructor to initialize the SAOM-NK model
    initialize = function(config_environ_params) {
      ##
      cat('\nTEST FROM CALLED CLASS INIT *BEFORE* BASE INIT\n')
      ##
      super$initialize(config_environ_params)
      ##
      cat('\nTEST FROM CALLED CLASS INIT *AFTER* BASE INIT\n')
      ##
      #
      visualize_init <- ifelse(is.null(config_environ_params[['visualize_init']]), 
                               TRUE, 
                               config_environ_params[['visualize_init']])
      if (visualize_init)
        self$visualize_system_bi_env_rsiena(plot_save = TRUE)
    },
    
    
    
    # Record fitness for each actor at the current iteration
    record_fitness = function() {
      fitness_snapshot <- sapply(1:self$M, function(actor_id) {
        self$objective_func(actor_id, 1) # Assuming a generic component ID for uniformity 
      })  
      self$fitness_history[[self$ITERATION]] <- fitness_snapshot
    },
    
    record_Ks = function() {
      ## Proj1: K_S Social space degree 
      self$degree_history_K_S[[self$ITERATION]] <- degree( self$get_social_igraph() )
      ## Proj2: K_E Component Environment space interaction (epistasis) degree
      self$degree_history_K_E[[self$ITERATION]] <- degree( self$get_component_igraph() )
    },
    
    ##
    get_rsiena_data_from_structure_model = function(structure_model) {
      ACTORS     <- sienaNodeSet(self$M, nodeSetName="ACTORS")
      COMPONENTS <- sienaNodeSet(self$N, nodeSetName="COMPONENTS")
      # strat_mat_varCovar <- varCovar(strat_eff$x, nodeSet = 'ACTORS')
      if (! 'coCovars' %in% names(structure_model$dv_bipartite) ) {
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      dv_name <- structure_model$dv_bipartite$name
        
      covDvTypes <- sapply(structure_model$dv_bipartite$coCovars, function(x)x$interaction1)
      componentDV_ids <- grep('self\\$component_\\d{1,2}_coCovar', covDvTypes) ## ex: "self$component_1_coCovar"
      stratDV_ids     <- grep('self\\$strat_\\d{1,2}_coCovar', covDvTypes ) ## ex: "self$strat_1_coCovar" 
      
      nstrat <- length(stratDV_ids)
      ncomponent <- length(componentDV_ids)
    
      ## 1 Component ------------------------------
      component_1_eff <- structure_model$dv_bipartite$coCovars[[ componentDV_ids[1] ]]
      self$component_1_coCovar <- coCovar(component_1_eff$x, nodeSet = 'COMPONENTS')
      ##-------------------------------------------
      #
      strat_1_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[1] ]]
      self$strat_1_coCovar <- coCovar(strat_1_eff$x, nodeSet = 'ACTORS')
      if (ncomponent==1 & nstrat==1 ) 
      {
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV, self$component_1_coCovar, self$strat_1_coCovar), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      if (ncomponent==1 & nstrat==2 ) 
      {
        strat_2_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[2] ]]
        self$strat_2_coCovar <- coCovar(strat_2_eff$x, nodeSet = 'ACTORS')
        #
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV, self$component_1_coCovar, self$strat_1_coCovar,self$strat_2_coCovar), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      if (ncomponent==1 & nstrat==3 ) 
      {
        strat_2_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[2] ]]
        self$strat_2_coCovar <- coCovar(strat_2_eff$x, nodeSet = 'ACTORS')
        strat_3_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[3] ]]
        self$strat_3_coCovar <- coCovar(strat_3_eff$x, nodeSet = 'ACTORS')
        #
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV, self$component_1_coCovar, self$strat_1_coCovar,self$strat_2_coCovar,self$strat_3_coCovar), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      if (ncomponent==1 & nstrat==4 ) 
      {
        strat_2_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[2] ]]
        self$strat_2_coCovar <- coCovar(strat_2_eff$x, nodeSet = 'ACTORS')
        strat_3_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[3] ]]
        self$strat_3_coCovar <- coCovar(strat_3_eff$x, nodeSet = 'ACTORS')
        strat_4_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[4] ]]
        self$strat_4_coCovar <- coCovar(strat_4_eff$x, nodeSet = 'ACTORS')
        #
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV, self$component_1_coCovar, self$strat_1_coCovar,self$strat_2_coCovar,self$strat_3_coCovar,self$strat_4_coCovar), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      if (ncomponent==1 & nstrat==5 ) 
      {
        strat_2_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[2] ]]
        self$strat_2_coCovar <- coCovar(strat_2_eff$x, nodeSet = 'ACTORS')
        strat_3_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[3] ]]
        self$strat_3_coCovar <- coCovar(strat_3_eff$x, nodeSet = 'ACTORS')
        strat_4_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[4] ]]
        self$strat_4_coCovar <- coCovar(strat_4_eff$x, nodeSet = 'ACTORS')
        strat_5_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[5] ]]
        self$strat_5_coCovar <- coCovar(strat_5_eff$x, nodeSet = 'ACTORS')
        #
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV, self$component_1_coCovar, self$strat_1_coCovar,self$strat_2_coCovar,self$strat_3_coCovar,self$strat_4_coCovar,self$strat_5_coCovar), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      #
      ## 2 Component --------------------------------------
      component_2_eff <- structure_model$dv_bipartite$coCovars[[ componentDV_ids[2] ]]
      self$component_2_coCovar <- coCovar(component_2_eff$x, nodeSet = 'COMPONENTS')
      ##----------------------------------------------------
      #
      if (ncomponent>=2 & nstrat==1 ) 
      {
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV, self$component_1_coCovar,self$component_2_coCovar, self$strat_1_coCovar), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      if (ncomponent>=2 & nstrat==2 ) 
      {
        strat_2_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[2] ]]
        self$strat_2_coCovar <- coCovar(strat_2_eff$x, nodeSet = 'ACTORS')
        #
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV, self$component_1_coCovar,self$component_2_coCovar,  self$strat_1_coCovar,self$strat_2_coCovar), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      if (ncomponent>=2 & nstrat==3 ) 
      {
        strat_2_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[2] ]]
        self$strat_2_coCovar <- coCovar(strat_2_eff$x, nodeSet = 'ACTORS')
        strat_3_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[3] ]]
        self$strat_3_coCovar <- coCovar(strat_3_eff$x, nodeSet = 'ACTORS')
        #
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV, self$component_1_coCovar,self$component_2_coCovar,  self$strat_1_coCovar,self$strat_2_coCovar,self$strat_3_coCovar), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      if (ncomponent>=2 & nstrat==4 ) 
      {
        strat_2_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[2] ]]
        self$strat_2_coCovar <- coCovar(strat_2_eff$x, nodeSet = 'ACTORS')
        strat_3_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[3] ]]
        self$strat_3_coCovar <- coCovar(strat_3_eff$x, nodeSet = 'ACTORS')
        strat_4_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[4] ]]
        self$strat_4_coCovar <- coCovar(strat_4_eff$x, nodeSet = 'ACTORS')
        #
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV, self$component_1_coCovar,self$component_2_coCovar,  self$strat_1_coCovar,self$strat_2_coCovar,self$strat_3_coCovar,self$strat_4_coCovar), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      if (ncomponent>=2 & nstrat==5 ) 
      {
        strat_2_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[2] ]]
        self$strat_2_coCovar <- coCovar(strat_2_eff$x, nodeSet = 'ACTORS')
        strat_3_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[3] ]]
        self$strat_3_coCovar <- coCovar(strat_3_eff$x, nodeSet = 'ACTORS')
        strat_4_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[4] ]]
        self$strat_4_coCovar <- coCovar(strat_4_eff$x, nodeSet = 'ACTORS')
        strat_5_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[5] ]]
        self$strat_5_coCovar <- coCovar(strat_5_eff$x, nodeSet = 'ACTORS')
        #
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV, self$component_1_coCovar,self$component_2_coCovar,  self$strat_1_coCovar,self$strat_2_coCovar,self$strat_3_coCovar,self$strat_4_coCovar,self$strat_5_coCovar), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      #
    },
    
    #
    init_rsiena_model_from_structure_model_bipartite_matrix = function(structure_model, bipartite_matrix, 
                                                                       rand_seed=123) {
      set.seed(rand_seed)
      # cat('\n\nDEBUG: called init_rsiena_model_from_structure_model_bipartite_matrix()\n\n')
      #
      ACTORS     <- sienaNodeSet(self$M, nodeSetName="ACTORS")
      COMPONENTS <- sienaNodeSet(self$N, nodeSetName="COMPONENTS")
      #
      self$config_structure_model <- structure_model
      structure_model_dvs <- names(structure_model)

      ## Simulation baseline nets should not be same; make one small change (toggle one dyad)
      ## @see https://www.stats.ox.ac.uk/~snijders/siena/NetworkSimulation.R
      bipartite_matrix1 <- bipartite_matrix
      bipartite_matrix2 <- bipartite_matrix
      .i <- sample(1:self$M, 1)
      .j <- sample(1:self$N, 1)
      bipartite_matrix2 <-  toggleBiMat(bipartite_matrix2, .i, .j )
      ##
      social_matrix1 <- bipartite_matrix1 %*% t(bipartite_matrix1)
      search_matrix1 <- t(bipartite_matrix1) %*% bipartite_matrix1
      ##
      social_matrix2 <- bipartite_matrix2 %*% t(bipartite_matrix2)
      search_matrix2 <- t(bipartite_matrix2) %*% bipartite_matrix2
      
      ## init networks to duplicate for the init arrays (two network waves)
      array_bi_net <- array(c(bipartite_matrix1, bipartite_matrix2), dim=c(self$M, self$N, 2) )
      array_social <- array(c(social_matrix1, social_matrix2), dim=c(self$M, self$M, 2) )
      array_search <- array(c(search_matrix1, search_matrix2), dim=c(self$N, self$N, 2) )
      ## Drop information above binary ties for RSiena DVs
      ##**TODO** Simulate potential influence / bias from this information loss
      array_bi_net[ array_bi_net > 1 ] <- 1
      array_social[ array_social > 1 ] <- 1
      array_search[ array_search > 1 ] <- 1
    
      
      #
      # sienaDepVars <- list()
      #
      if ('dv_social' %in% structure_model_dvs ) {
        self$social_rsienaDV <- sienaDependent(array_social, type='oneMode', nodeSet = 'ACTORS', allowOnly = F)
      }
      if ('dv_search' %in% structure_model_dvs) {
        self$search_rsienaDV <- sienaDependent(array_search, type='oneMode', nodeSet = 'COMPONENTS', allowOnly = F)
      }
      if ('dv_bipartite' %in% structure_model_dvs) {
        self$bipartite_rsienaDV <- sienaDependent(array_bi_net, type='bipartite', nodeSet =c('ACTORS', 'COMPONENTS'), allowOnly = F)
      }
      ##---------------------------------------------
      
      ##---------------------------------------------
      ##**TODO** debug/redo for arbitrary DVs inputs
      # if (all( c('dv_social','dv_search','dv_bipartite') %in% structure_model_dvs ))
      # {
      #   self$rsiena_data <- sienaDataCreate(list(self$social_rsienaDV, self$search_rsienaDV, self$bipartite_rsienaDV), 
      #                                       nodeSets = list(ACTORS, COMPONENTS))
      # }
      # else if ( all( c('dv_social','dv_bipartite') %in% structure_model_dvs ))
      # {
      #   self$rsiena_data <- sienaDataCreate(list(self$social_rsienaDV, self$bipartite_rsienaDV), 
      #                                       nodeSets = list(ACTORS, COMPONENTS))
      # }
      # else if ( all( c('dv_search','dv_bipartite') %in% structure_model_dvs ))
      # {
      #   self$rsiena_data <- sienaDataCreate(list(self$search_rsienaDV, self$bipartite_rsienaDV), 
      #                                       nodeSets = list(ACTORS, COMPONENTS))
      # } 
      # else if ( all( c('dv_social','dv_search') %in% structure_model_dvs ))
      # {
      #   self$rsiena_data <- sienaDataCreate(list(self$social_rsienaDV, self$search_rsienaDV), 
      #                                       nodeSets = list(ACTORS, COMPONENTS))
      # }
      # else if ('dv_social' %in% structure_model_dvs )
      # {
      #   self$rsiena_data <- sienaDataCreate(list(self$social_rsienaDV), 
      #                                       nodeSets = list(ACTORS, COMPONENTS))
      # }
      # else if ('dv_search' %in% structure_model_dvs )
      # {
      #   self$rsiena_data <- sienaDataCreate(list(self$search_rsienaDV), 
      #                                       nodeSets = list(ACTORS, COMPONENTS))
      # }
      if ('dv_bipartite' %in% structure_model_dvs)
      {
        # strat_mat_varCovar <- varCovar(strat_eff$x, nodeSet = 'ACTORS')
        self$rsiena_data <- self$get_rsiena_data_from_structure_model(structure_model)
      }
      else 
      {
        Stop('structural model has no dependent variables.')
      }
      ##---------------------------------------------
      ##
      # print('sienaDepVars')
      # print(sienaDepVars)
      ###
      # self$rsiena_data <- sienaDataCreate(list(self$social_rsienaDV, 
      #                                          self$search_rsienaDV, 
      #                                          self$bipartite_rsienaDV), 
      #                                     nodeSets = list(ACTORS, COMPONENTS))
      # self$rsiena_data <- sienaDataCreate(sienaDepVars, nodeSets = list(ACTORS, COMPONENTS))
      ###
      # self$rsiena_data <- rsienaDataCreate(list(self$bipartite_rsienaDV), 
      #                                    nodeSets = list(ACTORS, COMPONENTS))
      
      # print(self$rsiena_data )
      # print('END init_rsiena_model_from_bipartite_matrix ')
      # print(self)
      
    },
    
    #
    init_multiwave_rsiena_model_from_structure_model_bipartite_matrix = function(structure_model, 
                                                                                 bipartite_matrix1, 
                                                                                 bipartite_matrix2,
                                                                                 rand_seed=123) {
      set.seed(rand_seed)
      # cat('\n\nDEBUG: called init_rsiena_model_from_structure_model_bipartite_matrix()\n\n')
      #
      ACTORS     <- sienaNodeSet(self$M, nodeSetName="ACTORS")
      COMPONENTS <- sienaNodeSet(self$N, nodeSetName="COMPONENTS")
      #
      self$config_structure_model <- structure_model
      structure_model_dvs <- names(structure_model)
      #
      # hasVarCovars <- 'varCovars' %in% names(structure_model$dv_bipartite)
      hasCoCovars <- 'coCovars' %in% names(structure_model$dv_bipartite)
      #
      print('DEBUG:  hasCoCovars: ')
      print(hasCoCovars)
      ## Simulation baseline nets should not be same; make one small change (toggle one dyad)
      ## @see https://www.stats.ox.ac.uk/~snijders/siena/NetworkSimulation.R
      #
      if(identical(bipartite_matrix1, bipartite_matrix2)) {
        .i <- sample(1:self$M, 1)
        .j <- sample(1:self$N, 1)
        bipartite_matrix2 <-  toggleBiMat(bipartite_matrix2, .i, .j )
      }
      ##
      social_matrix1 <- bipartite_matrix1 %*% t(bipartite_matrix1)
      search_matrix1 <- t(bipartite_matrix1) %*% bipartite_matrix1
      ##
      social_matrix2 <- bipartite_matrix2 %*% t(bipartite_matrix2)
      search_matrix2 <- t(bipartite_matrix2) %*% bipartite_matrix2
      ## init networks to duplicate for the init arrays (two network waves)
      array_bi_net <- array(c(bipartite_matrix1, bipartite_matrix2), dim=c(self$M, self$N, 2) )
      array_social <- array(c(social_matrix1, social_matrix2), dim=c(self$M, self$M, 2) )
      array_search <- array(c(search_matrix1, search_matrix2), dim=c(self$N, self$N, 2) )
      ## Drop information above binary ties for RSiena DVs
      array_bi_net[ array_bi_net > 1 ] <- 1
      array_social[ array_social > 1 ] <- 1
      array_search[ array_search > 1 ] <- 1
      if ('dv_social' %in% structure_model_dvs ) {
        self$social_rsienaDV <- sienaDependent(array_social, type='oneMode', nodeSet = 'ACTORS', allowOnly = F)
      }
      if ('dv_search' %in% structure_model_dvs) {
        self$search_rsienaDV <- sienaDependent(array_search, type='oneMode', nodeSet = 'COMPONENTS', allowOnly = F)
      }
      if ('dv_bipartite' %in% structure_model_dvs) {
        self$bipartite_rsienaDV <- sienaDependent(array_bi_net, type='bipartite', nodeSet =c('ACTORS', 'COMPONENTS'), allowOnly = F)
      }
      ##---------------------------------------------
      self$rsiena_data <- self$get_rsiena_data_from_structure_model(structure_model)
    },
    
    
    ##
    ##**TODO**
    ##**Create custom RSiena interaction functions for only bipartite DV, **
    ##**but objective function includes statistics of the projections (social net, search landscape)**
    ##
    add_rsiena_effects = function(structure_model) {
      
      if (is.null(self$rsiena_effects))
        stop('initiate self$rsiena_effects before adding effects.')

      
      for (i in 1:length(structure_model)) {
        
        dv <- structure_model[[ i ]]
        
        for (j in 1:length(dv$effects)) {
          cat(sprintf('\n structural effects i=%s, j=%s\n', i, j))
          
          eff <- dv$effects[[ j ]]
          
          print(eff)
          # print(eff)
          # stop('DEBUG')
          
          self$include_rsiena_effect_from_eff_list(eff)
          # self$rsiena_effects <- setEffect(self$rsiena_effects, density, 
          #                                 name = 'self$bipartite_rsienaDV', parameter = 0, fix=T)  ## effect parameter 
        }
        
        for (j in 1:length(dv$coCovars)) {
          cat(sprintf('\n coCovars i=%s, j=%s\n', i, j))
          
          eff <- dv$coCovars[[ j ]]
          
          print(eff)
          # print(eff)
          # stop('DEBUG')
          
          self$include_rsiena_effect_from_eff_list(eff)
        }
        
        for (j in 1:length(dv$varCovars)) {
          ##**TODO** Implement variable covariates 
          ## (requires custom script to iteratively simulate next network and update covariate one wave at a time)
        }
        
      }
      
    },
  
    
    # RSIENA
    search_rsiena_init = function(structure_model, get_eff_doc = FALSE) {
      ##--1. RSiena Model
      ##  1.1. INIT: bipartite matrix --> RSiena model
      self$init_rsiena_model_from_structure_model_bipartite_matrix(self$bipartite_matrix, structure_model, self$rsiena_env_seed)
      
      print('self$rsiena_data : ')
      print(self$rsiena_data)
      
      ##  1.2. INIT effects list in simulation model (in RSiena model)
      self$rsiena_effects <- getEffects(self$rsiena_data)
      
      # Effects Documentation
      if(get_eff_doc)
        effectsDocumentation(self$rsiena_effects)
      ##-----------------------------
      
      # stop('DEBUG XXXXXXXXXXXXX')
      ##--2. NETWORK: STRUCTURE EVOLUTION (structure_Model)-----
      ##  2.1. Add effects from model objective function list
      self$add_rsiena_effects(structure_model)
      
      # ##--3. SEARCH (FITNESS): PAYOFFS (payoff_formulas) ----------
      # ##  3.1. set payoff rules (formulas)
      # self$config_payoff_formulas <- payoff_formulas
    },
    
    # #
    # search_rsiena_multiwave_init = function(structure_model, bipartite_matrix_0, bipartite_matrix_1, get_eff_doc = FALSE) {
    #   ##--1. RSiena Model
    #   ##  1.1. INIT: bipartite matrix --> RSiena model
    #   self$init_multiwave_rsiena_model_from_structure_model_bipartite_matrix(structure_model, 
    #                                                                            bipartite_matrix_0, 
    #                                                                            bipartite_matrix_1, 
    #                                                                            self$rsiena_env_seed)
    #   
    #   print('self$rsiena_data : ')
    #   print(self$rsiena_data)
    #   
    #   ##  1.2. INIT effects list in simulation model (in RSiena model)
    #   self$rsiena_effects <- getEffects(self$rsiena_data)
    #   
    #   # Effects Documentation
    #   if(get_eff_doc)
    #     effectsDocumentation(self$rsiena_effects)
    #   ##-----------------------------
    #   
    #   # stop('DEBUG XXXXXXXXXXXXX')
    #   ##--2. NETWORK: STRUCTURE EVOLUTION (structure_Model)-----
    #   ##  2.1. Add effects from model objective function list
    #   self$add_rsiena_effects(structure_model)
    #   
    #   # ##--3. SEARCH (FITNESS): PAYOFFS (payoff_formulas) ----------
    #   # ##  3.1. set payoff rules (formulas)
    #   # self$config_payoff_formulas <- payoff_formulas
    # },
    
    
    #
    search_rsiena_execute_sim = function(iterations, 
                                         returnDeps=T, 
                                         returnChains=T,
                                         rsiena_phase2_nsub=1, rsiena_n2start_scale=1, 
                                         digits=3,
                                         seed=123) {
      if (is.null(self$rsiena_effects))
        stop('Set rsiena_effects before running simulation.')
      # if (is.null(self$config_payoff_formulas))
      #   stop('Set search payoff formulas before running search simulation.')
      
      
      # cat(sprintf('\nDEBUG CHECK algorithm created, about to run siena07...\n'))
      
      
      ##-----------------------------
      ## RSiena Algorithm
      self$rsiena_algorithm <- sienaAlgorithmCreate(projname=sprintf('%s_%s',self$SIM_NAME,self$TIMESTAMP),
                                                    simOnly = T,
                                                    nsub = rsiena_phase2_nsub,
                                                    n2start = rsiena_n2start_scale * 2.52 * (7+sum(self$rsiena_effects$include)),
                                                    n3 = iterations,
                                                    seed = seed)
      
       
      # Run RSiena simulation
      self$rsiena_model <- siena07(self$rsiena_algorithm,
                                   data = self$rsiena_data, 
                                   effects = self$rsiena_effects,
                                   batch = TRUE,
                                   returnDeps = returnDeps, 
                                   returnChains = returnChains,
                                   returnDataFrame = TRUE, ##**TODO** CHECK
                                   returnLoglik = TRUE #,  ##**TODO** CHECK
                                   )   # returnChains = returnChains

      # Summarize and plot results
      mod_summary <- summary(self$rsiena_model)
      if(!is.null(mod_summary))
        print(mod_summary)
      # plot(self$rsiena_model)


      print(screenreg(list(self$rsiena_model), single.row = T, digits = digits))

      ## update simulation object environment from RSiena simulation model
      new_bi_env_igraph <- self$get_bipartite_igraph_from_rsiena_model()
      #
      self$set_system_from_bipartite_igraph( new_bi_env_igraph )

      #
      # self$plot_bipartite_system_from_mat(self$bipartite_matrix, iterations, plot_save = T)
      self$visualize_system_bi_env_rsiena(plot_save = TRUE)

      # stop('OK_DEBG')
    },  
    
    
    
    # # Extend existing simulation environment
    # search_rsiena_extend = function(iterations, returnDeps = TRUE) {
    #   self$rsiena_algorithm$n3 <- iterations
    #   self$rsiena_model <- siena07(self$rsiena_algorithm, data = self$rsiena_data, effects = self$rsiena_effects,
    #                                prevAns = self$rsiena_model,
    #                                batch = TRUE, returnDeps = returnDeps)
    #   ##
    #   mod_summary <- summary(self$rsiena_model)
    #   if(!is.null(mod_summary)) 
    #     print(mod_summary)
    #   ##
    #   print(screenreg(list(self$rsiena_model), single.row = T, digits = 3))
    # },
    
    
    ## Single Simulation Run
    search_rsiena_run = function(structure_model, 
                                 iterations=1000, 
                                 returnDeps=TRUE, returnChains=TRUE, ## TRUE=simulation only
                                 n_snapshots=1, 
                                 plot_save = TRUE, overwrite=TRUE, 
                                 rsiena_phase2_nsub=1, rsiena_n2start_scale=1,
                                 get_eff_doc = FALSE, digits=3,
                                 run_seed=123
                                 ) {
      set.seed(run_seed)
      self$rsiena_run_seed <- run_seed
      if( overwrite | is.null(self$rsiena_model) ) {
        ## 1. Init simulation
        # cat(sprintf('\nDEBUG RUN: 1. search_rsiena_init\n'))
        self$search_rsiena_init(structure_model, get_eff_doc)
        ## 2. Execute simulation
        # cat(sprintf('\nDEBUG RUN: 2. search_rsiena_execute_sim\n'))
        self$search_rsiena_execute_sim(iterations,
                                       returnDeps=returnDeps,
                                       returnChains=returnChains,
                                       rsiena_phase2_nsub=rsiena_phase2_nsub,
                                       rsiena_n2start_scale=rsiena_n2start_scale, 
                                       digits=digits,
                                       rand_seed=run_seed)
        ## 3. Process chain of simulation ministeps
        # cat(sprintf('\nDEBUG RUN: 3. search_rsiena_process_ministep_chain\n'))
        self$search_rsiena_process_ministep_chain()
        ## 4. Process actor statistics (e.g., utility)
        # cat(sprintf('\nDEBUG RUN: 4. search_rsiena_process_actor_stats\n'))
        self$search_rsiena_process_stats()
      } else {
        # self$search_rsiena_extend(objective_list, iterations)
        stop('_extend() method not yet implemented.')
      }
      # # self$visualize_system_bi_env_rsiena(plot_save = plot_save)
      # self$get_rsiena_model_snapshots(n=n_snapshots)
    },
    
    # ## Single Simulation Run
    # search_rsiena_multiwave_run = function(structure_model, 
    #                                        waves=4,
    #                                        iterations=1000, 
    #                                        returnDeps=TRUE, returnChains=TRUE, ## TRUE=simulation only
    #                                        n_snapshots=1, 
    #                                        plot_save = TRUE, overwrite=TRUE, 
    #                                        rsiena_phase2_nsub=1, rsiena_n2start_scale=1,
    #                                        get_eff_doc = FALSE, digits=3,
    #                                        run_seed=123
    # ) {
    #   
    #   set.seed(run_seed)
    #   self$rsiena_run_seed <- run_seed
    #   ## 1. Init simulation
    #   # cat(sprintf('\nDEBUG RUN: 1. search_rsiena_init\n'))
    #   self$search_rsiena_init(structure_model, get_eff_doc=FALSE)
    #   ## 2. Execute simulation
    #   # cat(sprintf('\nDEBUG RUN: 2. search_rsiena_execute_sim\n'))
    #   self$search_rsiena_execute_multiwave_sim(iterations,
    #                                            returnDeps=returnDeps,
    #                                            returnChains=returnChains,
    #                                            rsiena_phase2_nsub=rsiena_phase2_nsub,
    #                                            rsiena_n2start_scale=rsiena_n2start_scale, 
    #                                            digits=digits,
    #                                            rand_seed=run_seed)
    #   ## 3. Process chain of simulation ministeps
    #   # cat(sprintf('\nDEBUG RUN: 3. search_rsiena_process_ministep_chain\n'))
    #   self$search_rsiena_process_ministep_chain()
    #   ## 4. Process actor statistics (e.g., utility)
    #   # cat(sprintf('\nDEBUG RUN: 4. search_rsiena_process_stats\n'))
    #   self$search_rsiena_process_stats()
    #   
    #   self$rsiena_waves <- list()
    #   
    #   for (w in 1:waves) {
    #     ## DYNAMIC DECISIONS (STRATEGY COVARIATE CHANGES)
    #     # structure_model$coCovars[[1]]$x <- c(0,0,0,0,1,1,1,1) ## new strategy or function for time-based/decision-rule strategy choice
    #     #
    #     ## 1. Init simulation
    #     # cat(sprintf('\nDEBUG RUN: 1. search_rsiena_init\n'))
    #     self$search_rsiena_init(structure_model, get_eff_doc=FALSE)
    #     ## 2. Execute simulation
    #     # cat(sprintf('\nDEBUG RUN: 2. search_rsiena_execute_sim\n'))
    #     self$search_rsiena_execute_multiwave_sim(iterations,
    #                                              returnDeps=returnDeps,
    #                                              returnChains=returnChains,
    #                                              rsiena_phase2_nsub=rsiena_phase2_nsub,
    #                                              rsiena_n2start_scale=rsiena_n2start_scale, 
    #                                              digits=digits,
    #                                              rand_seed=run_seed)
    #     ## 3. Process chain of simulation ministeps
    #     # cat(sprintf('\nDEBUG RUN: 3. search_rsiena_process_ministep_chain\n'))
    #     self$search_rsiena_process_ministep_chain()
    #     ## 4. Process actor statistics (e.g., utility)
    #     # cat(sprintf('\nDEBUG RUN: 4. search_rsiena_process_stats\n'))
    #     self$search_rsiena_process_stats()
    #     
    #     self$rsiena_waves[[w]] <- list(
    #       rsiena_model = self$rsiena_model,
    #       actor_stats = self$actor_stats
    #     )
    #     
    #   }
    #   
    #   cat(sprintf('\nSimulated %s network waves.\n', waves))
    #   
    # },
    
    
    #**TODO**
    search_rsiena_multiwave_run = function(structure_model,
                                           waves=2, 
                                           iterations=1000, 
                                           returnDeps=T, 
                                           returnChains=T,
                                           rsiena_phase2_nsub=1, rsiena_n2start_scale=1, 
                                           digits=3,
                                           rand_seed=123) {
      bipartite_matrix_0 <- self$bipartite_matrix
      ##--1. INIT RSiena Model: set $rsiena_data --------
      self$init_rsiena_model_from_structure_model_bipartite_matrix(structure_model, bipartite_matrix_0, self$rsiena_env_seed)
      #
      print('self$rsiena_data : ')
      print(self$rsiena_data)
      ##  2. Init effects
      self$rsiena_effects <- getEffects(self$rsiena_data)
      ##  3. Add effects from structure_model list
      self$add_rsiena_effects(structure_model)
      ##  4. RSiena Algorithm 
      self$rsiena_algorithm <- sienaAlgorithmCreate(projname=sprintf('%s_%s',self$SIM_NAME,self$TIMESTAMP),
                                                    simOnly = T,
                                                    nsub = rsiena_phase2_nsub,
                                                    # n2start = rsiena_n2start_scale * 2.52 * (7+sum(self$rsiena_effects$include)),
                                                    n3 = iterations,
                                                    seed = rand_seed)
      ## 5. Run RSiena simulation
      self$rsiena_model <- siena07(self$rsiena_algorithm,
                                   data = self$rsiena_data, 
                                   effects = self$rsiena_effects,
                                   batch = TRUE,
                                   returnDeps = returnDeps, 
                                   returnChains = returnChains,
                                   returnDataFrame = TRUE, ##**TODO** CHECK
                                   returnLoglik = TRUE #,  ##**TODO** CHECK
      )   # returnChains = returnChains
      
      # Summarize and plot results
      mod_summary <- summary(self$rsiena_model)
      if(!is.null(mod_summary))
        print(mod_summary)
      # plot(self$rsiena_model)
      print(screenreg(list(self$rsiena_model), single.row = T, digits = digits))
      
      # 6. Update System
      ## update simulation object environment from RSiena simulation model
      new_bi_env_igraph <- self$get_bipartite_igraph_from_rsiena_model()
      self$set_system_from_bipartite_igraph( new_bi_env_igraph )
      
      # sink() ## write output text to file
      for (w in 1:waves){
        #
        bipartite_matrix_previous <- if(w == 1){ bipartite_matrix_0 }else{ self$bipartite_matrix_waves[[w-1]] }
        ##--1. INIT RSiena Model--------
        self$init_multiwave_rsiena_model_from_structure_model_bipartite_matrix(structure_model,
                                                                               bipartite_matrix_previous, ## matrix1 (previous)
                                                                               self$bipartite_matrix,     ## matrix2 (latest)
                                                                               self$rsiena_env_seed)
        print('self$rsiena_data : ')
        print(self$rsiena_data)
        ##  2. Init effects
        self$rsiena_effects <- getEffects(self$rsiena_data)
        ##  3. Add effects from model objective function list
        self$add_rsiena_effects(structure_model)
        ##  4. RSiena Algorithm
        self$rsiena_algorithm <- sienaAlgorithmCreate(projname=sprintf('%s_%s',self$SIM_NAME,self$TIMESTAMP),
                                                      simOnly = T,
                                                      nsub = rsiena_phase2_nsub,
                                                      # n2start = rsiena_n2start_scale * 2.52 * (7+sum(self$rsiena_effects$include)),
                                                      n3 = iterations,
                                                      seed = rand_seed)
        ## 5. Run RSiena simulation
        self$rsiena_model <- siena07(self$rsiena_algorithm,
                                     data = self$rsiena_data, 
                                     effects = self$rsiena_effects,
                                     batch = TRUE,
                                     returnDeps = returnDeps, 
                                     returnChains = returnChains,
                                     returnDataFrame = TRUE, ##**TODO** CHECK
                                     returnLoglik = TRUE #,  ##**TODO** CHECK
        )   # returnChains = returnChains
        
        # Summarize and plot results
        mod_summary <- summary(self$rsiena_model)
        if(!is.null(mod_summary))
          print(mod_summary)
        # plot(self$rsiena_model)
        print(screenreg(list(self$rsiena_model), single.row = T, digits = digits))
        
        ## 6. Update System
        ## update simulation object environment from RSiena simulation model
        new_bi_env_igraph <- self$get_bipartite_igraph_from_rsiena_model()
        #
        self$set_system_from_bipartite_igraph( new_bi_env_igraph )
        #
        # # self$plot_bipartite_system_from_mat(self$bipartite_matrix, iterations, plot_save = T)
        # self$visualize_system_bi_env_rsiena(plot_save = TRUE)
        
        ##---------- 2. Waves 2,3,4,... in Multiwave Simulation -------------------
        
        self$rsiena_model_waves[[w]] <- self$rsiena_model
        self$bipartite_matrix_waves[[w]] <- self$bipartite_matrix
        
      }
      
    },
    
    #**TODO**
    search_rsiena_multiwave_extend = function(waves=1, 
                                               iterations=1000, 
                                               returnDeps=T, 
                                               returnChains=T,
                                               rsiena_phase2_nsub=1, rsiena_n2start_scale=1, 
                                               digits=3,
                                               rand_seed=123) {
      ##  4. RSiena Algorithm 
      self$rsiena_algorithm <- sienaAlgorithmCreate(projname=sprintf('%s_%s',self$SIM_NAME,self$TIMESTAMP),
                                                    simOnly = T,
                                                    nsub = rsiena_phase2_nsub,
                                                    # n2start = rsiena_n2start_scale * 2.52 * (7+sum(self$rsiena_effects$include)),
                                                    n3 = iterations,
                                                    seed = rand_seed)
      ## 5. Run RSiena simulation
      self$rsiena_model <- siena07(self$rsiena_algorithm,
                                   data = self$rsiena_data, 
                                   effects = self$rsiena_effects,
                                   batch = TRUE,
                                   returnDeps = returnDeps, 
                                   returnChains = returnChains,
                                   returnDataFrame = TRUE, ##**TODO** CHECK
                                   returnLoglik = TRUE ,  ##**TODO** CHECK
                                   prevAns= self$rsiena_model
      )   # returnChains = returnChains
      
      # Summarize and plot results
      mod_summary <- summary(self$rsiena_model)
      if(!is.null(mod_summary))
        print(mod_summary)
      # plot(self$rsiena_model)
      print(screenreg(list(self$rsiena_model), single.row = T, digits = digits))
      
      # 6. Update System
      ## update simulation object environment from RSiena simulation model
      new_bi_env_igraph <- self$get_bipartite_igraph_from_rsiena_model()
      self$set_system_from_bipartite_igraph( new_bi_env_igraph )
      
      # # sink() ## write output text to file
      # for (w in 1:waves){
      #   #
      #   bipartite_matrix_previous <- if(w == 1){ bipartite_matrix_0 }else{ self$bipartite_matrix_waves[[w-1]] }
      #   ##--1. INIT RSiena Model--------
      #   self$init_multiwave_rsiena_model_from_structure_model_bipartite_matrix(structure_model,
      #                                                                          bipartite_matrix_previous, ## matrix1 (previous)
      #                                                                          self$bipartite_matrix,     ## matrix2 (latest)
      #                                                                          self$rsiena_env_seed)
      #   print('self$rsiena_data : ')
      #   print(self$rsiena_data)
      #   ##  2. Init effects
      #   self$rsiena_effects <- getEffects(self$rsiena_data)
      #   ##  3. Add effects from model objective function list
      #   self$add_rsiena_effects(structure_model)
      #   ##  4. RSiena Algorithm
      #   self$rsiena_algorithm <- sienaAlgorithmCreate(projname=sprintf('%s_%s',self$SIM_NAME,self$TIMESTAMP),
      #                                                 simOnly = T,
      #                                                 nsub = rsiena_phase2_nsub,
      #                                                 # n2start = rsiena_n2start_scale * 2.52 * (7+sum(self$rsiena_effects$include)),
      #                                                 n3 = iterations,
      #                                                 seed = rand_seed)
      #   ## 5. Run RSiena simulation
      #   self$rsiena_model <- siena07(self$rsiena_algorithm,
      #                                data = self$rsiena_data, 
      #                                effects = self$rsiena_effects,
      #                                batch = TRUE,
      #                                returnDeps = returnDeps, 
      #                                returnChains = returnChains,
      #                                returnDataFrame = TRUE, ##**TODO** CHECK
      #                                returnLoglik = TRUE #,  ##**TODO** CHECK
      #   )   # returnChains = returnChains
      #   
      #   # Summarize and plot results
      #   mod_summary <- summary(self$rsiena_model)
      #   if(!is.null(mod_summary))
      #     print(mod_summary)
      #   # plot(self$rsiena_model)
      #   print(screenreg(list(self$rsiena_model), single.row = T, digits = digits))
      #   
      #   ## 6. Update System
      #   ## update simulation object environment from RSiena simulation model
      #   new_bi_env_igraph <- self$get_bipartite_igraph_from_rsiena_model()
      #   #
      #   self$set_system_from_bipartite_igraph( new_bi_env_igraph )
      #   #
      #   # # self$plot_bipartite_system_from_mat(self$bipartite_matrix, iterations, plot_save = T)
      #   # self$visualize_system_bi_env_rsiena(plot_save = TRUE)
      #   
      #   ##---------- 2. Waves 2,3,4,... in Multiwave Simulation -------------------
      #   
      #   self$rsiena_model_waves[[w]] <- self$rsiena_model
      #   self$bipartite_matrix_waves[[w]] <- self$bipartite_matrix
      #   
      # }
      
    },
    
    
    #
    search_rsiena_multiwave_process_results = function() {
      actor_wave_stats <- list()
      actor_wave_util <- list()
      actor_wave_util_diff <- list()
      #
      K_wave_A <- list()
      K_wave_B1 <- list()
      K_wave_B2 <- list()
      K_wave_C <- list()
      #
      for (w in 1:length(self$rsiena_model_waves)) {
        #
        rsiena_model_w <- self$rsiena_model_waves[[ w ]]
        #
        bipartite_igraph_w <- self$get_bipartite_igraph_from_rsiena_model(rsiena_model_w)
        #
        self$set_system_from_bipartite_igraph( bipartite_igraph_w )
        ## 3. Process chain of simulation ministeps
        # cat(sprintf('\nDEBUG RUN: 3. search_rsiena_process_ministep_chain\n'))
        self$search_rsiena_process_ministep_chain()
        ## 4. Process actor statistics (e.g., utility)
        # cat(sprintf('\nDEBUG RUN: 4. search_rsiena_process_stats\n'))
        self$search_rsiena_process_stats()
        ##**TODO**
        actor_wave_stats[[w]]     <- self$actor_stats_df %>% mutate(wave_id=w)
        actor_wave_util[[w]]      <- self$actor_util_df %>% mutate(wave_id=w)
        actor_wave_util_diff[[w]] <- self$actor_util_diff_df %>% mutate(wave_id=w)
        # #
        # component_wave_stats[[w]] <- self$component_stats %>% mutate(wave_id=w)
        # # self$actor_proj_wave_stats[[w]]     <- self$actor_proj_stats %>% mutate(wave_id=w)
        # # self$component_proj_wave_stats[[w]] <- self$component_proj_stats %>% mutate(wave_id=w)
        
        K_wave_A[[w]]  <- self$K_A_df  %>% mutate(wave_id=w)
        K_wave_B1[[w]] <- self$K_B1_df %>% mutate(wave_id=w)
        K_wave_B2[[w]] <- self$K_B2_df %>% mutate(wave_id=w)
        K_wave_C[[w]]  <- self$K_C_df  %>% mutate(wave_id=w)
      }
      
      ##**TODO**
      self$actor_wave_stats     <- data.table::rbindlist( actor_wave_stats )
      self$actor_wave_util      <- data.table::rbindlist( actor_wave_util )
      self$actor_wave_util_diff <- data.table::rbindlist( actor_wave_util_diff )
      
      # ##
      # self$component_wave_stats <- data.table::rbindlist( component_wave_stats )
      
      self$K_wave_A  <- data.table::rbindlist( K_wave_A ) 
      self$K_wave_B1 <- data.table::rbindlist( K_wave_B1 ) 
      self$K_wave_B2 <- data.table::rbindlist( K_wave_B2 ) 
      self$K_wave_C  <- data.table::rbindlist( K_wave_C ) 
    },
    
    # #
    # chain_stats = NULL,  ## 
    # actor_stats = NULL,
    # component_stats = NULL,
    # actor_proj_stats = NULL,
    # component_proj_stats = NULL,
    # ##
    # wave_stats_list = list(),
    
      
    #
    get_bipartite_igraph_from_rsiena_model = function(rsiena_model = NULL, sim_iteration = NULL) {
      rsiena_model <- if(is.null(rsiena_model)) { self$rsiena_model } else { rsiena_model }
      new_bi_env_mat <- self$get_bipartite_matrix_from_rsiena_model(rsiena_model, sim_iteration)
      new_bi_env_igraph <- igraph::graph_from_biadjacency_matrix(new_bi_env_mat, directed = T, weighted = T, mode = 'out') ##**TODO: CHECK** all vs. out
      return(new_bi_env_igraph)
    },
    
    #
    get_bipartite_matrix_from_rsiena_model = function(rsiena_model = NULL, sim_iteration = NULL, wave_id=1) {
      rsiena_model <- if(is.null(rsiena_model)) { self$rsiena_model } else { rsiena_model }
      ######### UPDATE SYSTEM ENVIRONMENT SNAPSHOT FROM EVOLVED Bipartite Network DV #####################
      if (! length(rsiena_model$sims) )
        stop('Run RSiena simulations to set rsiena_model$sims before get_bipartite_matrix_from_rsiena_model.')
      # sim_id_last <- length( rsiena_model$sims ) ##
      sim_id <- ifelse(is.null(sim_iteration), 
                       length(rsiena_model$sims), ## default current state is the last simulation in sims list
                       sim_iteration)
      bi_env_dv_id <- which( names(self$config_structure_model) == 'dv_bipartite' )
      ##
      # actors <- 1:self$M 
      # components <- 1:self$N
      MplusN <- self$M + self$N
      ## Get DV (bi-partite network) from simulation iteration=sim_id
      el_bi_env <- rsiena_model$sims[[ sim_id ]][[ wave_id ]][[ bi_env_dv_id ]]$`1`
      ## update numbering of second mode (the N component integer names shift upward by the number of actors M to match bipartite naming)
      el_bi_env[,2] <- el_bi_env[,2] + self$M
      ### get networks from other prjected space ties 
      # el_proj1  <- rsiena_model$sims[[1]][[1]][[1]]$`1` ## Social network
      # el_proj2  <- rsiena_model$sims[[1]][[1]][[2]]$`1` ## 
      # el_proj2 <- el_proj2 + self$M
      #############
      ## Undirected --> Upper right rectangle of full bipartite matrix
      ##   M N
      ## M[0,X], for X in [0,1]
      ## N[0,0]
      bi_env_mat_sp <- sparseMatrix(i = el_bi_env[,1],
                                    j = el_bi_env[,2],
                                    x = el_bi_env[,3],
                                    dims = c(MplusN, MplusN))
      new_bi_env_mat <- as.matrix(bi_env_mat_sp)[ 1:self$M, (self$M+1):(MplusN) ]
      return( new_bi_env_mat )
    },
    
    ##
    get_chain_stats_list = function() { ## 'all','statdf','bi_env_arr', 'util', 'util_diff')
      if (is.null(self$chain_stats))
        stop('chain_stats missing; run simulation with returnChains=TRUE in order to compute actor utility')
      # ## actor Utility vector
      # au <- c( mat %*% self$rsiena_model$theta )
      # hist( util )
      effnames <- unlist(sapply(self$config_structure_model$dv_bipartite$effects, function(x) x$effect))
      #
      coCovar_effnames  <- unlist(sapply(self$config_structure_model$dv_bipartite$coCovars, function(x) x$effect))
      varCovar_effnames <- unlist(sapply(self$config_structure_model$dv_bipartite$varCovars, function(x) x$effect))
      ## ALL effect names
      effnames <- c(effnames, coCovar_effnames, varCovar_effnames)
      
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
      
      ## Use system bipartite matrix as start
      bi_env_mat <- self$bipartite_matrix
      
      ## remove chain entries where no change was made (keep if !stability )
      tiechdf <- self$chain_stats[ !self$chain_stats$stability, ]
      
      ## get matrix timeseries and network statistics timeseries
      nchains <- length(self$rsiena_model$chain)
      ## paramters
      theta   <- self$rsiena_model$theta
      ## replace fixed theta param with given param values (instead of zero default when effect is fixed)
      if (any(fixed_params)) 
        theta[ fixed_params ] <- config_param_vals[ fixed_params ]
      #
      ntheta <- length(theta)
      #
      bi_env_arr <- array(NA, dim=c(self$M, self$N, nchains))
      #
      stats_li <- list()
      util_li  <- list()
      util_diff_li <- list()
      K_A_li <- list()
      K_B1_li <- list()
      K_B2_li <- list()
      K_C_li <- list()
      # bi_env_long <- data.frame()
      # statdf <- data.frame()
      # utildf <- data.frame()
      # util_diff <- data.frame()
      for (i in 1:nrow(tiechdf)) {
        mstep <- tiechdf[i,]
        # cat(sprintf('\n %.2f%s i=%s, j=%s', 100*i/nrow(tiechdf),'%', mstep$id_from,  mstep$id_to))
        if(i %% 100 == 0) cat(sprintf('\n %.2f%s', 100*i/nrow(tiechdf),'%'))
        ## update bipartite environment matrix for one step (toggle one dyad)
        bi_env_mat <- toggleBiMat(bi_env_mat, mstep$id_from,  (mstep$id_to - self$M)  )
        #
        # bi_env_mat
        #
        # Get New Statistics
        # statsobj <- self$get_struct_mod_net_stats_list_from_bi_mat( bi_env_mat )
        statmat <- self$get_struct_mod_stats_mat_from_bi_mat( bi_env_mat )
        #
        step_statgrid <- expand.grid(chain_step_id=i, actor_id=1:self$M, effect_id=1:ntheta)
        step_statgrid$value <- c( statmat )
        #
        stats_li[[i]] <- step_statgrid   ##**CHECK list**
        #
        # statsdf  <- rbind(statsdf,  statsobj$df )
        #
        # statsmat <- statsobj$mat  #get_n_by_m_mat_from_long_df
        #
        
        #
        # tmpstatdf <- as.data.frame( statsmat )
        # tmpstatdf$chain_step_id <- i
        # tmpstatdf$actor_id <- factor(1:self$M)
        # statdf <- rbind(statdf,  tmpstatdf )
        ## network statistics matrix added to array
        # stat_long <- rbind(stat_long, statdf)
        ## Add ministep updated bipartite environment to array
        bi_env_arr[ , , i]  <- bi_env_mat
        #
        # stat_arr[ , , i] <- statmat
        ## Add utilities to array
        util <- c( statmat %*% theta )
        # utildf <- rbind(utildf, data.frame(
        #   utility = util,
        #   chain_step_id = i, 
        #   actor_id = factor(1:self$M )
        # ))
        ##
        step_utilgrid <- expand.grid(chain_step_id=i, actor_id=1:self$M)
        step_utilgrid$utility <- util
        util_li[[i]] <- step_utilgrid  ##**CHECK list**
        ##
        ##
        # util_diff <- rbind(util_diff, data.frame(
        #   utility = if(i == 1){ NA } else {util - util_lag }, ## diff 
        #   chain_step_id = i, 
        #   actor_id = factor(1:self$M )
        # ))
        ##
        step_util_diffgrid <- expand.grid(chain_step_id=i, actor_id=1:self$M)
        step_util_diffgrid$utility <- if(i == 1){ NA } else { util - util_lag }
        util_diff_li[[i]] <- step_util_diffgrid  ##**CHECK list**
        ##
        ######
        ## update utility lag for next period difference
        util_lag <- util
        ######
        
        #
        K_A_grid <- expand.grid(chain_step_id=i, actor_id=1:self$M)
        # K_A_grid$value <- rowSums( bi_env_mat %*% t(bi_env_mat) , na.rm=T) ## WEIGHTED
        K_A_grid$value <- apply( bi_env_mat %*% t(bi_env_mat), 1, function(x)sum(x>0, na.rm=T) )
        K_A_li[[i]] <- K_A_grid
        #
        K_B1_grid <- expand.grid(chain_step_id=i, actor_id=1:self$M)
        # K_B1_grid$value <- rowSums( bi_env_mat, na.rm=T)  ## WEIGHTED
        K_B1_grid$value <- apply( bi_env_mat, 1, function(x)sum(x>0, na.rm=T))
        K_B1_li[[i]] <- K_B1_grid
        #
        K_B2_grid <- expand.grid(chain_step_id=i, component_id=1:self$N)
        # K_B2_grid$value <- colSums( bi_env_mat, na.rm=T)  ## WEIGHTED
        K_B2_grid$value <- apply( bi_env_mat, 2, function(x)sum(x>0, na.rm=T))
        K_B2_li[[i]] <- K_B2_grid
        #
        K_C_grid <- expand.grid(chain_step_id=i, component_id=1:self$N)
        # K_C_grid$value <- colSums( t(bi_env_mat) %*% bi_env_mat  , na.rm=T)  ## WEIGHTED
        K_C_grid$value <- apply( t(bi_env_mat) %*% bi_env_mat, 2, function(x)sum(x>0, na.rm=T))  ## WEIGHTED
        K_C_li[[i]] <- K_C_grid
      }
      ##-----------------
      ## Actor Network Statistics long dataframe
      stats_df <- data.table::rbindlist( stats_li )
      stats_df$actor_id <- as.factor(stats_df$actor_id)
      stats_df$effect_name <- as.factor(effnames[ stats_df$effect_id ])
      ## Actor Utility  long dataframe
      util_df <- data.table::rbindlist( util_li ) 
      util_df$actor_id <- as.factor(util_df$actor_id)
      ## Actor Utility Difference long dataframe
      util_diff_df <- data.table::rbindlist( util_diff_li ) 
      util_diff_df$actor_id <- as.factor(util_diff_df$actor_id)
      ##---
      
      K_A_df <- data.table::rbindlist( K_A_li ) 
      K_A_df$actor_id <- as.factor(K_A_df$actor_id)
      #
      K_B1_df <- data.table::rbindlist( K_B1_li ) 
      K_B1_df$actor_id <- as.factor(K_B1_df$actor_id)
      #
      K_B2_df <- data.table::rbindlist( K_B2_li ) 
      K_B2_df$component_id <- as.factor(K_B2_df$component_id)
      #
      K_C_df <- data.table::rbindlist( K_C_li ) 
      K_C_df$component_id <- as.factor(K_C_df$component_id)
      
      #####---------------
      return(list(
        stats_df = stats_df,
        util_df = util_df,
        util_diff_df = util_diff_df,
        K_A_df = K_A_df,
        K_B1_df = K_B1_df,
        K_B2_df = K_B2_df,
        K_C_df = K_C_df,
        bi_env_arr = bi_env_arr
      ))
    },
    
    
    ##
    ##
    ## ## @see https://www.stats.ox.ac.uk/~snijders/siena/WorkOnChains.r
    # For the meaning of the 13 fields of a ministep:
    # From siena07utilities.cpp:
    #       SET_STRING_ELT(colnames, 0, mkChar("Aspect")); network or behaviour;
    #       SET_STRING_ELT(colnames, 1, mkChar("Var"));
    #       SET_STRING_ELT(colnames, 2, mkChar("VarName"));
    #       SET_STRING_ELT(colnames, 3, mkChar("Ego"));
    #       SET_STRING_ELT(colnames, 4, mkChar("Alter"));
    #       SET_STRING_ELT(colnames, 5, mkChar("Diff"));
    #               difference in dependent behavior variable: 0 if nothing changes
    #       SET_STRING_ELT(colnames, 6, mkChar("ReciRate"));
    #       SET_STRING_ELT(colnames, 7, mkChar("LogOptionSetProb"));
    #       SET_STRING_ELT(colnames, 8, mkChar("LogChoiceProb"));
    #       SET_STRING_ELT(colnames, 9, mkChar("Diagonal"));
    # This is from C++ code; note that in C++, numbering starts at 0.
    # Therefore, for a one-mode network the values of Ego and Alter run
    # from 0 to n-1, where n is the number of actors.
    # If the mini-step is a network mini-step for a one-mode network,
    # the change made can be inferred by comparing ego and alter.
    # If ego and alter are the same node, then no change has occurred;
    # if ego and alter are different nodes, then the value for the tie
    # from ego to alter is toggled (1 -> 0; 0 -> 1).
    # For a two-mode network the values of Alter run from 0 to m,
    # where m is the number of nodes in the second mode;
    # here the value Alter=m means that no change has occurred.
    ##
    ### Interpretation of the columns on the outputted df - 
    #
    # 1 - network or behavior function 
    # 2 - same as 1, denoted as 0/1 
    # 3 - same as 1/2, denoted using varname 
    # 4 - ego ID (starting @0) 
    # 5 - if network function - alter ID (starting @0) for changes 
    # ego ID for no change 
    # 0 if behavior function 
    # 6 - 0 if network 
    # -1, 0, 1 for change in behavior level behavior 
    # Our aims don't make use of columns 7-10
    # 11 - Designates stability (i.e., “True” means no change, “False” means change) 
    ##
    search_rsiena_process_ministep_chain = function() {
      
      if (is.null(self$rsiena_model$chain)) {
        stop("Chain not available. Ensure returnChains=TRUE was set in the siena07 call.")
      }
      simChain <- self$rsiena_model$chain
      ##
      depvar <- 1
      period <- 1

      ###--------
      ## Bipatite network chain --> value Alter=m means that no change has occurred.
      ##** "chain[[run]][[depvar]][[period]][[ministep]]"**
      chainDatZeroIndex <- ldply(seq_along(simChain), function(iter){
        ncolsIter <- length(simChain[[iter]][[depvar]][[period]])
        t(matrix(unlist(simChain[[iter]][[depvar]][[period]]), nc=ncolsIter))
      })
      ### one-index chain as new object
      chainDat <- chainDatZeroIndex
      # chainDat[,1] <- chainDat[,1]           ## aspect
      chainDat[,2] <- as.integer(chainDat[,2]) ## 0=Network; 1=Behavior
      chainDat[,4] <- as.numeric(chainDat[,4]) + 1 ## Ego (from) 1-indexing from C++ 0-index
      chainDat[,5] <- as.numeric(chainDat[,5]) + 1 ## Alter (to) 1-indexing from C++ 0-index
      chainDat[,6] <- as.numeric(chainDat[,6]) ## Behavior difference
      chainDat[,7] <- as.numeric(chainDat[,7]) ##
      chainDat[,8] <- as.numeric(chainDat[,8]) ##
      chainDat[,9] <- as.numeric(chainDat[,9]) ##
      # chainDat[,10] <-chainDat[,10] ## Diagonal
      chainDat[,11] <-as.logical(chainDat[,11]) ## Stability
      # Set dataframe Names
      
      # print(chainDat)
      # print(head(chainDat))
      
      ###--------
      .getTieChangeAfterOneindexing <- function(x, N) {
        ## already re-indexed id_from and id_to -- changing from C++ 0-index to R 1-index
        dv_name <- x[3]
        id_from <- x[4]
        id_to   <- x[5]
        if(dv_name == 'self$bipartite_rsienaDV') {
          ## bipartite network after oneIndexing the node ids:  id_to==(N+1) means no tie 
          return(ifelse(id_to == (N+1), FALSE, TRUE))
        } else if (dv_name %in% c('self$social_rsienaDV','self$search_rsienaDV')) {
          ## bipartite network after oneIndexing the node ids:  id_to==id_from means no tie
          return(ifelse(id_from == id_to, FALSE, TRUE))
        } else {
          stop(sprintf('dv_name %s not implemented in .getTieChange()', dv_name))
        }
      }
      
      names(chainDat) <- c('dv_type','dv_type_bin','dv_varname','id_from','id_to','beh_difference',
                           'reciprocal_rate','LogOptionSetProb', 'LogChoiceProb', 'diagonal','stability')
      # Set tie change by rules
      chainDat$tie_change <- apply(chainDat, 1, function(x) .getTieChangeAfterOneindexing(x, N=self$N) )
      
      ## set to simulation self
      self$chain_stats <- chainDat
        
      
      print(head(self$chain_stats))
      print(summary(self$chain_stats))
      print(dim(self$chain_stats))
    },
    
    #
    search_rsiena_process_stats = function() {
      ##
      # actor_stats_df = NULL,
      # actor_util_df = NULL,
      # actor_util_diff_df = NULL,
      # ##
      chain_stats_list <- self$get_chain_stats_list()
      #
      self$actor_stats_df      <- chain_stats_list$stats_df
      self$actor_util_df       <- chain_stats_list$util_df
      self$actor_util_diff_df  <- chain_stats_list$util_diff_df
      ###
      # self$component_stats      <- list() ##self$get_component_stats_list()
      # self$actor_proj_stats     <- list() ##self$get_component_stats_list()
      # self$component_proj_stats <- list() ##self$get_component_stats_list()
      ###
      self$K_A_df  <- chain_stats_list$K_A_df
      self$K_B1_df <- chain_stats_list$K_B1_df
      self$K_B2_df <- chain_stats_list$K_B2_df
      self$K_C_df  <- chain_stats_list$K_C_df
    },
    
    
    #
    process_sims = function() {
     
      #   'K_A'=list(), ##**TODO** Actor social space [degree = central position in social network]
      #   'K_B'=list(), #  Bipartite space  [ degree = resource/affiliation density ]
      #   'K_C'=list()  #  Component space  [ degree = interdependence/structuration : epistatic interactions ]
      #   ), 
      
      self$sims_stats <- NULL
    },
    
    
    ##===================== PLOTTING ============================

    
    #
    search_rsiena_plot_stability = function(tol=1e-5, step_size=1, wave_id=1) {
      
      sims = self$rsiena_model$sims
      n <- length(sims)
      bi_env_dv_id <- which( names(self$config_structure_model) == 'dv_bipartite' )
      if (!length(bi_env_dv_id)) stop('dv_bipartite is missing from self$config_structure_model.')
      #
      outlist <- list()
      difflist <- list()
      jaccardlist <- list()
      K_soc_list <- list()
      K_env_list <- list()
      ##
      sim_ids_plot <- seq(1, n, by=step_size) ##**TODO CHECK**
      # ##**TODO: CHECK** if not plotting every simulation, then drop first sim ( (nsims/step_size) - 1)
      # if(step_size>1) sim_ids_plot <- sim_ids_plot[-1]
      ###
      for(i in 1:length(sim_ids_plot)) {
        sim_id <- sim_ids_plot[ i ]
        cat(sprintf(' %s ', i))
        ##
        el_bi_env <- sims[[ sim_id ]][[ wave_id ]][[ bi_env_dv_id ]]$`1`
        ## update numbering of second mode (the comonent integer names shift upward by the number of actors)
        el_bi_env[,2] <- el_bi_env[,2] + self$M
        ##
        MplusN <- self$M + self$N
        #############
        ## Bipartite matrix space (N+M by N+M)
        ## Undirected --> Upper right rectangle of full bipartite matrix
        ##   M N
        ## M[0,X], for X in [0,1]
        ## N[0,0]
        bi_env_mat_sp <- sparseMatrix(i = el_bi_env[,1],
                                      j = el_bi_env[,2],
                                      x = el_bi_env[,3],
                                      dims = c(MplusN, MplusN))
        ##
        bi_env_mat_new  <- as.matrix(bi_env_mat_sp)[ 1:self$M, (self$M+1):(MplusN) ]
        ##
        ## Add ACTOR NAMES on rows
        rownames(bi_env_mat_new) <- as.character( 1:self$M )
        ## Add COMPONENT NAMES on columns (N+1 ... N+M)
        colnames(bi_env_mat_new) <- as.character( (1:self$N) + self$M  )
        ##
        # self$plot_bipartite_system_from_mat(bi_env_mat_new, i, plot_save = TRUE) 
        ##
        outlist[[ sprintf('sim%d',i) ]] <- bi_env_mat_new
        
        if (i == 1) {
          difflist[[ sprintf('diff%d-%d', i-1, i) ]] <- 0
          jaccardlist[[sprintf('jac%d-%d', i-1, i)]] <- 0
        } else {
          difflist[[ sprintf('diff%d-%d', i-1, i) ]] <-  bi_env_mat_new - outlist[[ (i-1) ]] 
          jaccardlist[[sprintf('jac%d-%d', i-1, i)]] <- get_jaccard_index(m0 = outlist[[ (i-1) ]], m1 = bi_env_mat_new )
        }
        
        new_bi_g <- igraph::graph_from_biadjacency_matrix(bi_env_mat_new, 
                                                          directed = T, mode = 'all', 
                                                          multiple = T, weighted = T, 
                                                          add.names = T)
        projections <- igraph::bipartite_projection(new_bi_g, multiplicity = T, which = 'both')
        K_soc_list[[i]] <- igraph::degree(projections$proj1)
        K_env_list[[i]] <- igraph::degree(projections$proj2)
        
        # print(sim_i[[1]][[3]]`1`)
      }
      
      ##_-----------------------------------------
      ## K interdependencies (K_E, K_S)
      
      degree_summary <- do.call(rbind, lapply(1:length(K_soc_list), function(iter) {
         data.frame(
           Iteration = iter,
           Mean_K_S = mean(K_soc_list[[iter]]),
           Q25_K_S = quantile(K_soc_list[[iter]], 0.25),
           Q75_K_S = quantile(K_soc_list[[iter]], 0.75),
           Mean_K_E = mean(K_env_list[[iter]]),
           Q25_K_E = quantile(K_env_list[[iter]], 0.25),
           Q75_K_E = quantile(K_env_list[[iter]], 0.75)
         )
      }))
      K_plt <- ggplot(degree_summary, aes(x = Iteration)) +
         geom_line(aes(y = Mean_K_S, color = "Mean K_S")) +
         geom_ribbon(aes(ymin = Q25_K_S, ymax = Q75_K_S, fill = "K_S"), alpha = 0.05) +
         geom_line(aes(y = Mean_K_E, color = "Mean K_E")) +
         geom_ribbon(aes(ymin = Q25_K_E, ymax = Q75_K_E, fill = "K_E"), alpha = 0.05) +
         scale_color_manual(values = c("Mean K_S" = "blue", "Mean K_E" = "red")) +
         scale_fill_manual(values = c("K_S" = "blue", "K_E" = "red")) +
         labs(title = "Degree Progress Over Iterations",
                   x = "Iteration",
                   y = "Degree",
                   color = "Mean Degree",
                   fill = "IQR (Mid-50%)") +
         theme_minimal()
      
      print(K_plt)
      
      
      #-------------------------------------------
      par(mfrow=c(1,3))
      ##------------------------------------------
      jaccard_vec <- plyr::ldply(jaccardlist)[sim_ids_plot[-1], 2] ## skip first period (no change yet)
      n_changes <- length(jaccard_vec)
      ##
      # sim_ids_plot_steps <- if(step_size > 1) { sim_ids_plot[-c(1:2)] } else { sim_ids_plot[-1] }
      sim_ids_plot_steps <- sim_ids_plot[-1] 
      ##
      stability_vec <- cumsum(jaccard_vec) / (1:n_changes)
      #
      # stability_delta <- cumsum(stability_vec) / (1:n_changes)
      stability_delta <- c(0, abs(diff(stability_vec))) / abs(stability_vec)
      #
      plot(x=sim_ids_plot_steps, y=jaccard_vec, 
           type='l', ylab='Jaccard Index [t-1, t]', xlab='Simulation Iteration',
           main='Inter-Sim Distance\n(Jaccard Index between same-wave sims)')
      #
      plot(x=sim_ids_plot_steps, y=stability_vec , 
           type='l' , ylab='Jaccard Index [t-1, t]', xlab='Simulation Iteration',
           ylim=c( .9, 1),
           main='Stabiliation by Iteration\n(Inter-Sim Distance Moving Average)'
           ); abline(h=1, col='gray', lty=2)
      #
      plot(x = sim_ids_plot_steps, y=stability_delta, 
           type='l', log='y', xlab='Simulation Iteration', 
           # ylim=c( stability_delta[n-1],  1 ),
           ylab='Ln Stability Change [t-1, t]', 
           main='Sufficient Iterations?\n(Inter-Sim Distance Moving Average Change)' 
           ); abline(h = tol, col='pink', lty=2)
      
      # burn_prop <- 0.2
      # iter_postburn <- round(c( 1-burn_prop, 1) * n_changes )
      # idx <- iter_postburn[1]:iter_postburn[2]
      # if (any(idx))
      #   hist(stability_vec[ idx ], main='Stability (post-burn)', breaks=31)
      
    },
    
    # Convenience function for plotting all relevant plots 
    # or one plot by specifying plot 
    search_rsiena_multiwave_plot = function(type=NA, 
                                            rolling_window = 10, 
                                            actor_ids=c(),
                                            component_ids=c(),
                                            wave_ids=c(),
                                            thin_factor=1, 
                                            thin_wave_factor=1,
                                            smooth_method='loess',
                                            show_utility_points=TRUE,
                                            append_plot=FALSE,
                                            histogram_position='identity',
                                            scale_utility=TRUE,
                                            return_plot=TRUE
    ) {
      # if (all(is.na(type)) |  'ministep_count' %in% type)
      #   return(self$search_rsiena_plot_actor_ministep_count(rolling_window))
      # if (all(is.na(type)) | type == 'utility')
      #   return(self$search_rsiena_plot_actor_utility(rolling_window))
      plist <- list()
      # if (all(is.na(type)) |  'utility' %in% type)
      #   plist[['utility']] <- self$search_rsiena_plot_actor_utility(actor_ids, thin_factor, smooth_method, show_utility_points, return_plot=T)
      # 
      if (all(is.na(type)) |  'K_4panel' %in% type)
        plist[['K_4panel']] <- self$search_rsiena_multiwave_plot_K_4panel(actor_ids, component_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, return_plot=T)
      
      if (all(is.na(type)) |  'K_A_strategy_summary' %in% type)
        plist[['K_A_strategy_summary']] <- self$search_rsiena_multiwave_plot_K_A_strategy_summary(actor_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, return_plot=T)
      
      if (all(is.na(type)) |  'K_B1_strategy_summary' %in% type)
        plist[['K_B1_strategy_summary']] <- self$search_rsiena_multiwave_plot_K_B1_strategy_summary(actor_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, return_plot=T)
      
      if (all(is.na(type)) |  'K_B2_strategy_summary' %in% type)
        plist[['K_B2_strategy_summary']] <- self$search_rsiena_multiwave_plot_K_B2_strategy_summary(component_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, return_plot=T)
      
      if (all(is.na(type)) |  'K_C_strategy_summary' %in% type)
        plist[['K_C_strategy_summary']] <- self$search_rsiena_multiwave_plot_K_C_strategy_summary(component_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, return_plot=T)
      
      
      if (all(is.na(type)) |  'utility_strategy_summary' %in% type)
        plist[['utility_strategy_summary']] <- self$search_rsiena_multiwave_plot_actor_utility_strategy_summary(actor_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, scale_utility, return_plot=T)
      
      if (all(is.na(type)) |  'utility_by_strategy' %in% type)
        plist[['utility_by_strategy']] <- self$search_rsiena_multiwave_plot_actor_utility_by_strategy(actor_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, return_plot=T)
      # 
      # if (all(is.na(type)) |  'utility_density' %in% type)
      #   plist[['utility_density']] <- self$search_rsiena_plot_actor_utility_density(return_plot=T)
      # 
      if (all(is.na(type)) |  'utility_density_by_strategy' %in% type)
        plist[['utility_density_by_strategy']] <- self$search_rsiena_multiwave_plot_actor_utility_density_by_strategy(thin_wave_factor, return_plot=T)
      # 
      # if (all(is.na(type)) |  'utility_histogram_by_strategy' %in% type)
      #   plist[['utility_histogram_by_strategy']] <- self$search_rsiena_plot_actor_utility_histogram_by_strategy(histogram_position, return_plot=T)
      # 
      # if (all(is.na(type)) |  'stability' %in% type)
      #   plist[['stability']] <- self$search_rsiena_plot_stability() ##**TODO** Fix return plot
      # 
      #  SET plots 
      self$multiwave_plots <- if(append_plot) { append(self$multiwave_plots, plist) } else { plist }
      
      # stop(sprintf('Plot type not implemented: %s', type))
      if(return_plot)
        return(plist)
    },
    
    #
    search_rsiena_multiwave_plot_K_4panel = function(actor_ids=c(), 
                                                     component_ids=c(), 
                                                     wave_ids=c(), 
                                                     thin_factor=1, 
                                                     thin_wave_factor=1, 
                                                     smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                     show_utility_points=T, 
                                                     return_plot=FALSE
                                                     ) {
      K_A  <- self$search_rsiena_multiwave_plot_K_A_strategy_summary(actor_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, show_legend=T, show_title=F, return_plot=T)
      K_B1 <- self$search_rsiena_multiwave_plot_K_B1_strategy_summary(actor_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, show_legend=F, show_title=T, return_plot=T)
      K_B2 <- self$search_rsiena_multiwave_plot_K_B2_strategy_summary(component_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, show_legend=F, show_title=T, return_plot=T)
      K_C  <- self$search_rsiena_multiwave_plot_K_C_strategy_summary(component_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, show_legend=T, show_title=F, return_plot=T)
      # K_B1_leg <- self$search_rsiena_multiwave_plot_K_B1_strategy_summary(actor_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, show_legend=T, show_title=F, return_plot=T)
      # K_B2_leg <- self$search_rsiena_multiwave_plot_K_B2_strategy_summary(component_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, show_legend=T, show_title=F, return_plot=T)
      combined_plot <- ggarrange(
        K_B1, K_B2, 
        K_A, K_C,
        nrow = 2, ncol = 2, 
        widths = c(7, 7), # Adjust column widths
        heights = c(7, 7),
        common.legend = TRUE, # Share a common legend if needed
        legend = "bottom"#,     # Place legend at the bottom
        # align = 'v'
      ) 
      # # Extract legends
      # legend_B1 <- get_legend(K_B1_leg + theme(legend.position = "bottom"))
      # legend_B2 <- get_legend(K_B2_leg + theme(legend.position = "bottom"))
      # 
      # # Remove legends from individual plots
      # legend_B1 <- legend_B1 + theme(legend.position = "none")
      # legend_B2 <- legend_B2 + theme(legend.position = "none")
      # 
      # # Arrange plots with shared legend
      # plot_grid(
      #   ggarrange(K_A, K_C,
      #             K_B1+theme(legend.position = "none"), K_B2+theme(legend.position = "none"), 
      #             ncol = 2, nrow=2),
      #   ggarrange(legend_B1,legend_B2, ncol=2),
      #   nrow = 2,
      #   rel_heights = c(.95, .05) # Adjust legend height
      # )
      # self$multiwave_plots <- list(K_A=K_A, K_B1=K_B1, K_B2=K_B2, K_C=K_C)
      self$multiwave_plots <- list(combined_plot=combined_plot)
      #
      if(return_plot)
        return(combined_plot)
    },
    
    ##
    search_rsiena_multiwave_plot_K_C_strategy_summary = function(component_ids=c(), 
                                                                  wave_ids=c(),
                                                                  thin_factor=1, 
                                                                  thin_wave_factor=1,
                                                                  smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                                  show_utility_points=T,
                                                                  show_legend=TRUE,
                                                                  show_title=TRUE,
                                                                  return_plot=FALSE
                                                                 ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not are set.")
      if ( attr(self$component_1_coCovar, 'nodeSet') != 'COMPONENTS' )
        stop("Component payoff values in self$component_1_coCovar are not set.")
      # actor_strat <- as.factor( self$strat_coCovar )
      # component_types <- as.factor( self$component_coCovar )
      # component_types <- as.factor( sample(0:1, self$N, replace = T) )
      component_types <- as.factor( ifelse(self$component_1_coCovar > 0.5, 'High', 'Low') )
      #
      nstep <- sum(!self$chain_stats$stability)
      #
      strateffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
      stratparams <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$parameter)
      stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)ifelse(x$fix,'(Fix)',''))
      # stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)substr(as.character(x$fix),1,1))
      #
      structeffs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$effect)
      structparams <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$parameter)
      structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)ifelse(x$fix,'(Fix)',''))
      # structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)substr(as.character(x$fix),1,1))
      ## Compare 2 actors utilty
      dat <- self$K_wave_C %>% 
        filter(chain_step_id %% thin_factor == 0) %>% 
        filter(wave_id %% thin_wave_factor == 0 ) %>% 
        mutate(
          # strategy = actor_strat[ actor_id ],
          component_type = component_types[ component_id ],
          chain_below_med =  chain_step_id < median(chain_step_id)
        )
      dat$chain_half <- factor(ifelse(dat$chain_below_med, '1st Half', '2nd Half'))
      #
      y_lab <- 'K_C: Component Epistasis Degree'
      # if(scale_utility) {
      #   util_sc <- scale(dat$value)
      #   y_lab <- sprintf('K_A: Social Degree\n(Standardized Center = %.2f; Scale = %.2f)',
      #                       attr(util_sc, 'scaled:center'), 
      #                       attr(util_sc, 'scaled:scale'))
      #   dat <- dat %>% mutate(value = c(scale(value)))
      # }
      #
      density_rng <- range(dat$value, na.rm=T)
      density_absdiff_scale <- abs(diff(density_rng)) * 0.1
      y_lim <- c(density_rng[1] - density_absdiff_scale,  density_rng[2] + density_absdiff_scale)
      #
      point_size <- 10 / log( nstep )
      point_alpha <- min( 1,  1/log10( nstep ) )
      #
      if(length(component_ids))
        dat <- dat %>% filter(component_id %in% component_ids)
      if(length(wave_ids))
        dat <- dat %>% filter(wave_id %in% wave_ids)
      dat_wave_means <- dat %>% group_by(wave_id) %>% 
        summarize(mean=mean(value, na.rm=T))
      plt <- ggplot(dat, aes(x=chain_step_id, y=value)) + 
        geom_hline(data=dat_wave_means, aes(yintercept=mean), linetype=3, col='black' ) +
        facet_grid(wave_id ~ .) 
      if(show_utility_points)
        plt <- plt + geom_point(aes(color=component_type), alpha=point_alpha, shape=1, size=point_size)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(exists(smooth_method))
        plt <- plt + geom_smooth(aes(linetype=component_type, color=component_type, fill=component_type), method = smooth_method, linewidth=1, alpha=.09)
      #
      plt <- plt + theme_bw() + 
        # scale_linetype_manual(values = rep(1:8, length.out = length(unique(dat$component_id)))) +
        ylim(y_lim) + 
        ylab(y_lab) +
        xlab('Actor Decision Chain Ministep') +
        theme(
          panel.grid.minor = element_blank(),
          legend.position = "bottom"#,
          # legend.box = "horizontal",
          # legend.box.just = "center",
          # legend.key.width = unit(0.8, "cm"),  # Adjust legend key width
          # legend.spacing.x = unit(0.3, "cm")   # Adjust spacing between keys
        )
      if (show_title)
        plt <- plt + ggtitle(sprintf('Strategy:  %s\nStructure:  %s', 
                                     paste( paste(paste(strateffs, stratparams, sep='= '), stratfixs, sep='' ), collapse = ';  '),
                                     paste( paste(paste(structeffs, structparams, sep='= '), structfixs, sep=''), collapse = ';  ')
        ))
      #
      plt <- plt +  guides(color = guide_legend(nrow = 1))
      
      #### Density
      stratmeans <- dat %>% group_by(component_type, wave_id) %>% 
        summarize(n=n(), mean=mean(value, na.rm=T), sd=sd(value, na.rm=T))
      ## Actor density fact plots comparing H1 to H2 utility distribution
      plt2 <- ggplot(dat, aes(x=value, color=component_type, fill=component_type)) + ##linetype=chain_half
        geom_density(alpha=.1, linewidth=1)  +
        # geom_histogram(alpha=.1, position = 'dodge') +
        geom_vline(data = stratmeans, aes(xintercept = mean, color=component_type), linetype=2, linewidth=.9) +
        geom_vline(data = dat_wave_means,  aes(xintercept=mean), linetype=3, col='black' ) +
        labs(y='', x='') +
        # xlim(c(ggplot_build(plt)$layout$panel_params[[1]]$y.range)) + 
        xlim(y_lim) +
        coord_flip() +
        facet_grid(wave_id ~ .) +
        ylab('K_C Density') +
        # labs(color='component_type', fill='component_type') +
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
        ) 
      if (show_title)
        plt2 <- plt2 + ggtitle('\n')
      
      # combined_plot <- plot_grid(plt ,#+ theme(panel.spacing = unit(0, "lines")), 
      #                            plt2 ,# + theme(panel.spacing = unit(0, "lines")), 
      #                            ncol = 2, rel_widths = c(4, 1)) 
      # print(combined_plot)
      
      combined_plot <- ggarrange(
        plt, plt2, 
        ncol = 2, 
        widths = c(4.1,0.9), # Adjust column widths
        common.legend = TRUE, # Share a common legend if needed
        legend = ifelse(show_legend, "bottom", "none")#,     # Place legend at the bottom
        # align = 'h'
      ) 
      # print(combined_plot)
      #
      # print(plt)
      #
      if(return_plot)
        return(combined_plot)
    },
    
    ##
    search_rsiena_multiwave_plot_K_B2_strategy_summary = function(component_ids=c(), 
                                                                  wave_ids=c(),
                                                                  thin_factor=1, 
                                                                  thin_wave_factor=1,
                                                                  smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                                  show_utility_points=T,
                                                                  show_legend=TRUE,
                                                                  show_title=TRUE,
                                                                  return_plot=FALSE
    ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar are not set.")
      if ( attr(self$component_1_coCovar, 'nodeSet') != 'COMPONENTS' )
        stop("Component payoff values in self$component_1_coCovar are not set.")
      # actor_strat <- as.factor( self$strat_coCovar )
      # component_types <- as.factor( self$component_type_coCovar )
      # component_types <- as.factor( rep(1, self$N ) )
      component_types <- as.factor( ifelse(self$component_1_coCovar > 0.5, 'High', 'Low') )
      #
      nstep <- sum(!self$chain_stats$stability)
      #
      strateffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
      stratparams <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$parameter)
      stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)ifelse(x$fix,'(Fix)',''))
      # stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)substr(as.character(x$fix),1,1))
      #
      structeffs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$effect)
      structparams <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$parameter)
      structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)ifelse(x$fix,'(Fix)',''))
      # structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)substr(as.character(x$fix),1,1))
      ## Compare 2 actors utilty
      dat <- self$K_wave_B2 %>% 
        filter(chain_step_id %% thin_factor == 0) %>% 
        filter(wave_id %% thin_wave_factor == 0 ) %>% 
        mutate(
          component_type = component_types[ component_id ],
          chain_below_med =  chain_step_id < median(chain_step_id)
        )
      dat$chain_half <- factor(ifelse(dat$chain_below_med, '1st Half', '2nd Half'))
      #
      y_lab <- 'K_B2: Component-Actor Degree'
      # if(scale_utility) {
      #   util_sc <- scale(dat$value)
      #   y_lab <- sprintf('K_A: Social Degree\n(Standardized Center = %.2f; Scale = %.2f)',
      #                       attr(util_sc, 'scaled:center'), 
      #                       attr(util_sc, 'scaled:scale'))
      #   dat <- dat %>% mutate(value = c(scale(value)))
      # }
      #
      density_rng <- range(dat$value, na.rm=T)
      density_absdiff_scale <- abs(diff(density_rng)) * 0.1
      y_lim <- c(density_rng[1] - density_absdiff_scale,  density_rng[2] + density_absdiff_scale)
      #
      point_size <- 10 / log( nstep )
      point_alpha <- min( 1,  1/log10( nstep ) )
      #
      if(length(component_ids))
        dat <- dat %>% filter(component_id %in% component_ids)
      if(length(wave_ids))
        dat <- dat %>% filter(wave_id %in% wave_ids)
      dat_wave_means <- dat %>% group_by(wave_id) %>% 
        summarize(mean=mean(value, na.rm=T))
      plt <- ggplot(dat, aes(x=chain_step_id, y=value)) + 
        geom_hline(data=dat_wave_means, aes(yintercept=mean), linetype=3, col='black' ) +
        facet_grid(wave_id ~ .) 
      if(show_utility_points)
        plt <- plt + geom_point(aes(color=component_type), alpha=point_alpha, shape=1, size=point_size)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(exists(smooth_method))
        plt <- plt + geom_smooth(aes(linetype=component_type, color=component_type, fill=component_type), method = smooth_method, linewidth=1, alpha=.09)
      #
      plt <- plt + theme_bw() + 
        # scale_linetype_manual(values = rep(1:8, length.out = length(unique(dat$component_id)))) +
        ylim(y_lim) + 
        ylab(y_lab) +
        xlab('Actor Decision Chain Ministep') +
        theme(
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.box.just = "center",
          legend.key.width = unit(0.8, "cm"),  # Adjust legend key width
          legend.spacing.x = unit(0.3, "cm")   # Adjust spacing between keys
        )
      if(show_title)
        plt <- plt + ggtitle(sprintf('Strategy:  %s\nStructure:  %s', 
                                     paste( paste(paste(strateffs, stratparams, sep='= '), stratfixs, sep='' ), collapse = ';  '),
                                     paste( paste(paste(structeffs, structparams, sep='= '), structfixs, sep=''), collapse = ';  ')
        ))
      #
      plt <- plt +  guides(color = guide_legend(nrow = 1))
      
      #### Density
      stratmeans <- dat %>% group_by(component_type, wave_id) %>% 
        summarize(n=n(), mean=mean(value, na.rm=T), sd=sd(value, na.rm=T))
      ## Actor density fact plots comparing H1 to H2 utility distribution
      plt2 <- ggplot(dat, aes(x=value, color=component_type, fill=component_type)) + ##linetype=chain_half
        geom_density(alpha=.1, linewidth=1)  +
        # geom_histogram(alpha=.1, position = 'dodge') +
        geom_vline(data = stratmeans, aes(xintercept = mean, color=component_type), linetype=2, linewidth=.9) +
        geom_vline(data = dat_wave_means,  aes(xintercept=mean), linetype=3, col='black' ) +
        labs(y='', x='') +
        # xlim(c(ggplot_build(plt)$layout$panel_params[[1]]$y.range)) + 
        xlim(y_lim) +
        coord_flip() +
        facet_grid(wave_id ~ .) +
        ylab('K_B2 Density') +
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
        ) 
      if (show_title)
        plt2 <- plt2 + ggtitle('\n')
      
      # combined_plot <- plot_grid(plt ,#+ theme(panel.spacing = unit(0, "lines")), 
      #                            plt2 ,# + theme(panel.spacing = unit(0, "lines")), 
      #                            ncol = 2, rel_widths = c(4, 1)) 
      # print(combined_plot)
      
      combined_plot <- ggarrange(
        plt, plt2, 
        ncol = 2, 
        widths = c(4.1,0.9), # Adjust column widths
        common.legend = TRUE, # Share a common legend if needed
        legend = ifelse(show_legend, "bottom", "none")#,     # Place legend at the bottom
        # align = 'h'
      ) 
      # print(combined_plot)
      #
      # print(plt)
      #
      if(return_plot)
        return(combined_plot)
    },
    
    
    ##
    search_rsiena_multiwave_plot_K_B1_strategy_summary = function(actor_ids=c(), 
                                                                 wave_ids=c(),
                                                                 thin_factor=1, 
                                                                 thin_wave_factor=1,
                                                                 smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                                 show_utility_points=T,
                                                                 show_legend=TRUE,
                                                                 show_title=TRUE,
                                                                 return_plot=FALSE
    ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <- as.factor( self$strat_1_coCovar )
      #
      nstep <- sum(!self$chain_stats$stability)
      #
      strateffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
      stratparams <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$parameter)
      stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)ifelse(x$fix,'(Fix)',''))
      # stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)substr(as.character(x$fix),1,1))
      #
      structeffs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$effect)
      structparams <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$parameter)
      structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)ifelse(x$fix,'(Fix)',''))
      # structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)substr(as.character(x$fix),1,1))
      ## Compare 2 actors utilty
      dat <- self$K_wave_B1 %>% 
        filter(chain_step_id %% thin_factor == 0) %>% 
        filter(wave_id %% thin_wave_factor == 0 ) %>% 
        mutate(
          strategy = actor_strat[ actor_id ],
          chain_below_med =  chain_step_id < median(chain_step_id)
        )
      dat$chain_half <- factor(ifelse(dat$chain_below_med, '1st Half', '2nd Half'))
      #
      y_lab <- 'K_B1: Actor-Component Degree'
      # if(scale_utility) {
      #   util_sc <- scale(dat$value)
      #   y_lab <- sprintf('K_A: Social Degree\n(Standardized Center = %.2f; Scale = %.2f)',
      #                       attr(util_sc, 'scaled:center'), 
      #                       attr(util_sc, 'scaled:scale'))
      #   dat <- dat %>% mutate(value = c(scale(value)))
      # }
      #
      density_rng <- range(dat$value, na.rm=T)
      density_absdiff_scale <- abs(diff(density_rng)) * 0.1
      y_lim <- c(density_rng[1] - density_absdiff_scale,  density_rng[2] + density_absdiff_scale)
      #
      point_size <- 10 / log( nstep )
      point_alpha <- min( 1,  1/log10( nstep ) )
      #
      if(length(actor_ids))
        dat <- dat %>% filter(actor_id %in% actor_ids)
      if(length(wave_ids))
        dat <- dat %>% filter(wave_id %in% wave_ids)
      dat_wave_means <- dat %>% group_by(wave_id) %>% 
        summarize(mean=mean(value, na.rm=T))
      plt <- ggplot(dat, aes(x=chain_step_id, y=value)) + 
        geom_hline(data=dat_wave_means, aes(yintercept=mean), linetype=3, col='black' ) +
        facet_grid(wave_id ~ .) 
      if(show_utility_points)
        plt <- plt + geom_point(aes(color=strategy), alpha=point_alpha, shape=1, size=point_size)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(exists(smooth_method))
        plt <- plt + geom_smooth(aes(linetype=actor_id, color=strategy, fill=strategy), method = smooth_method, linewidth=1, alpha=.09)
      #
      plt <- plt + theme_bw() + 
        ylim(y_lim) + 
        ylab(y_lab) +
        xlab('Actor Decision Chain Ministep') +
        theme(
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.box.just = "center",
          legend.key.width = unit(0.8, "cm"),  # Adjust legend key width
          legend.spacing.x = unit(0.3, "cm")   # Adjust spacing between keys
        )
      if(show_title)
        plt <- plt + ggtitle(sprintf('Strategy:  %s\nStructure:  %s', 
                                     paste( paste(paste(strateffs, stratparams, sep='= '), stratfixs, sep='' ), collapse = ';  '),
                                     paste( paste(paste(structeffs, structparams, sep='= '), structfixs, sep=''), collapse = ';  ')
        ))
      #
      plt <- plt +  guides(color = guide_legend(nrow = 1))
      
      #### Density
      stratmeans <- dat %>% group_by(strategy, wave_id) %>% 
        summarize(n=n(), mean=mean(value, na.rm=T), sd=sd(value, na.rm=T))
      ## Actor density fact plots comparing H1 to H2 utility distribution
      plt2 <- ggplot(dat, aes(x=value, color=strategy, fill=strategy)) + ##linetype=chain_half
        geom_density(alpha=.1, linewidth=1)  +
        # geom_histogram(alpha=.1, position = 'dodge') +
        geom_vline(data = stratmeans, aes(xintercept = mean, color=strategy), linetype=2, linewidth=.9) +
        geom_vline(data = dat_wave_means,  aes(xintercept=mean), linetype=3, col='black' ) +
        labs(y='', x='') +
        # xlim(c(ggplot_build(plt)$layout$panel_params[[1]]$y.range)) + 
        xlim(y_lim) +
        coord_flip() +
        facet_grid(wave_id ~ .) +
        ylab('K_B1 Density') +
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
        ) 
      if (show_title)
        plt2 <- plt2 + ggtitle('\n')
      
      # combined_plot <- plot_grid(plt ,#+ theme(panel.spacing = unit(0, "lines")), 
      #                            plt2 ,# + theme(panel.spacing = unit(0, "lines")), 
      #                            ncol = 2, rel_widths = c(4, 1)) 
      # print(combined_plot)
      
      combined_plot <- ggarrange(
        plt, plt2, 
        ncol = 2, 
        widths = c(4.1,0.9), # Adjust column widths
        common.legend = TRUE, # Share a common legend if needed
        legend = ifelse(show_legend, "bottom", "none")#,     # Place legend at the bottom
        # align = 'h'
      ) 
      # print(combined_plot)
      #
      # print(plt)
      #
      if(return_plot)
        return(combined_plot)
    },
    
    
    ##
    search_rsiena_multiwave_plot_K_A_strategy_summary = function(actor_ids=c(), 
                                                                 wave_ids=c(),
                                                                 thin_factor=1, 
                                                                 thin_wave_factor=1,
                                                                 smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                                 show_utility_points=T,
                                                                 show_legend=TRUE,
                                                                 show_title=TRUE,
                                                                 return_plot=FALSE
    ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <- as.factor( self$strat_1_coCovar )
      #
      nstep <- sum(!self$chain_stats$stability)
      #
      strateffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
      stratparams <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$parameter)
      stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)ifelse(x$fix,'(Fix)',''))
      # stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)substr(as.character(x$fix),1,1))
      #
      structeffs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$effect)
      structparams <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$parameter)
      structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)ifelse(x$fix,'(Fix)',''))
      # structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)substr(as.character(x$fix),1,1))
      ## Compare 2 actors utilty
      dat <- self$K_wave_A %>% 
        filter(chain_step_id %% thin_factor == 0) %>% 
        filter(wave_id %% thin_wave_factor == 0 ) %>% 
        mutate(
          strategy = actor_strat[ actor_id ],
          chain_below_med =  chain_step_id < median(chain_step_id)
        )
      dat$chain_half <- factor(ifelse(dat$chain_below_med, '1st Half', '2nd Half'))
      #
      y_lab <- 'K_A: Social Degree'
      density_rng <- range(dat$value, na.rm=T)
      density_absdiff_scale <- abs(diff(density_rng)) * 0.1
      y_lim <- c(density_rng[1] - density_absdiff_scale,  density_rng[2] + density_absdiff_scale)
      #
      point_size <- 10 / log( nstep )
      point_alpha <- min( 1,  1/log10( nstep ) )
      #
      if(length(actor_ids))
        dat <- dat %>% filter(actor_id %in% actor_ids)
      if(length(wave_ids))
        dat <- dat %>% filter(wave_id %in% wave_ids)
      dat_wave_means <- dat %>% group_by(wave_id) %>% 
        summarize(mean=mean(value, na.rm=T))
      plt <- ggplot(dat, aes(x=chain_step_id, y=value)) + 
        geom_hline(data=dat_wave_means, aes(yintercept=mean), linetype=3, col='black' ) +
        facet_grid(wave_id ~ .) 
      if(show_utility_points)
        plt <- plt + geom_point(aes(color=strategy), alpha=point_alpha, shape=1, size=point_size)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(exists(smooth_method))
        plt <- plt + geom_smooth(aes(linetype=actor_id, color=strategy, fill=strategy), method = smooth_method, linewidth=1, alpha=.09)
      #
      plt <- plt + theme_bw() + 
        ylim(y_lim) + 
        ylab(y_lab) +
        xlab('Actor Decision Chain Ministep') +
        theme(
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.box.just = "center",
          legend.key.width = unit(0.8, "cm"),  # Adjust legend key width
          legend.spacing.x = unit(0.3, "cm")   # Adjust spacing between keys
        )
      if(show_title)
        plt <- plt + ggtitle(sprintf('Strategy:  %s\nStructure:  %s', 
                                     paste( paste(paste(strateffs, stratparams, sep='= '), stratfixs, sep='' ), collapse = ';  '),
                                     paste( paste(paste(structeffs, structparams, sep='= '), structfixs, sep=''), collapse = ';  ')
        ))
      #
      plt <- plt +  guides(color = guide_legend(nrow = 1))
      
      #### Density
      stratmeans <- dat %>% group_by(strategy, wave_id) %>% 
        summarize(n=n(), mean=mean(value, na.rm=T), sd=sd(value, na.rm=T))
      ## Actor density fact plots comparing H1 to H2 utility distribution
      plt2 <- ggplot(dat, aes(x=value, color=strategy, fill=strategy)) + ##linetype=chain_half
        geom_density(alpha=.1, linewidth=1)  +
        # geom_histogram(alpha=.1, position = 'dodge') +
        geom_vline(data = stratmeans, aes(xintercept = mean, color=strategy), linetype=2, linewidth=.9) +
        geom_vline(data = dat_wave_means,  aes(xintercept=mean), linetype=3, col='black' ) +
        labs(y='', x='') +
        # xlim(c(ggplot_build(plt)$layout$panel_params[[1]]$y.range)) + 
        xlim(y_lim) +
        coord_flip() +
        facet_grid(wave_id ~ .) +
        ylab('K_A Density') +
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
        )  
      if (show_title)
        plt2 <- plt2 + ggtitle('\n')
      
      # combined_plot <- plot_grid(plt ,#+ theme(panel.spacing = unit(0, "lines")), 
      #                            plt2 ,# + theme(panel.spacing = unit(0, "lines")), 
      #                            ncol = 2, rel_widths = c(4, 1)) 
      # print(combined_plot)
      
      combined_plot <- ggarrange(
        plt, plt2, 
        ncol = 2, 
        widths = c(4.1,0.9), # Adjust column widths
        common.legend = TRUE, # Share a common legend if needed
        legend = ifelse(show_legend, "bottom", "none")#,     # Place legend at the bottom
        # align = 'h'
      ) 
      # print(combined_plot)
      #
      # print(plt)
      #
      if(return_plot)
        return(combined_plot)
    },
    
    
    ##
    search_rsiena_multiwave_plot_actor_utility_strategy_summary = function(actor_ids=c(), 
                                                                           wave_ids=c(),
                                                                            thin_factor=1, 
                                                                            thin_wave_factor=1,
                                                                            smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                                            show_utility_points=T,
                                                                           scale_utility=TRUE,
                                                                            return_plot=FALSE
    ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <- as.factor( self$strat_1_coCovar )
      #
      nstep <- sum(!self$chain_stats$stability)
      #
      strateffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
      stratparams <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$parameter)
      stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)ifelse(x$fix,'(Fix)',''))
      # stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)substr(as.character(x$fix),1,1))
      #
      structeffs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$effect)
      structparams <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$parameter)
      structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)ifelse(x$fix,'(Fix)',''))
      # structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)substr(as.character(x$fix),1,1))
      ## Compare 2 actors utilty
      dat <- self$actor_wave_util %>% 
        filter(chain_step_id %% thin_factor == 0) %>% 
        filter(wave_id %% thin_wave_factor == 0 ) %>% 
        mutate(
          strategy = actor_strat[ actor_id ],
          chain_below_med =  chain_step_id < median(chain_step_id)
        )
      dat$chain_half <- factor(ifelse(dat$chain_below_med, '1st Half', '2nd Half'))
      #
      util_lab <- 'Actor Utility'
      if(scale_utility) {
        util_sc <- scale(dat$utility)
        util_lab <- sprintf('Actor Utility\n(Standardized Center = %.2f; Scale = %.2f)',
                            attr(util_sc, 'scaled:center'), 
                            attr(util_sc, 'scaled:scale'))
        dat <- dat %>% mutate(utility = c(scale(utility)))
      }
      #
      density_rng <- range(dat$utility, na.rm=T)
      density_absdiff_scale <- abs(diff(density_rng)) * 0.1
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
        summarize(mean=mean(utility, na.rm=T))
      plt <- ggplot(dat, aes(x=chain_step_id, y=utility)) + 
        geom_hline(data=dat_wave_means, aes(yintercept=mean), linetype=3, col='black' ) +
        facet_grid(wave_id ~ .) 
      if(show_utility_points)
        plt <- plt + geom_point(aes(color=strategy), alpha=point_alpha, shape=1, size=point_size)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(exists(smooth_method))
        plt <- plt + geom_smooth(aes(linetype=actor_id, color=strategy, fill=strategy), method = smooth_method, linewidth=1, alpha=.09)
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
      plt <- plt + ggtitle(sprintf('Strategy:  %s\nStructure:  %s', 
                      paste( paste(paste(strateffs, stratparams, sep='= '), stratfixs, sep='' ), collapse = ';  '),
                      paste( paste(paste(structeffs, structparams, sep='= '), structfixs, sep=''), collapse = ';  ')
      ))
      plt <- plt +  guides(color = guide_legend(nrow = 1))
     
      #### Density
      stratmeans <- dat %>% group_by(strategy, wave_id) %>% 
        summarize(n=n(), mean=mean(utility, na.rm=T), sd=sd(utility, na.rm=T))
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
          ) + ggtitle('\n')
      
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
      if(return_plot)
        return(combined_plot)
    },
    
    ##
    search_rsiena_multiwave_plot_actor_utility_by_strategy = function(actor_ids=c(), 
                                                                      thin_factor=1, 
                                                                      thin_wave_factor=1,
                                                                      smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                                      show_utility_points=TRUE,
                                                                      return_plot=FALSE
    ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <- as.factor( self$strat_1_coCovar )
      ## Compare 2 actors utilty
      dat <- self$actor_wave_util %>% 
        filter(chain_step_id %% thin_factor == 0) %>% 
        filter(wave_id %% thin_wave_factor == 0 ) %>% 
        mutate(strategy = actor_strat[ actor_id ] )
      if(length(actor_ids))
        dat <- dat %>% filter(actor_id %in% actor_ids)
      plt <- ggplot(dat, aes(x=chain_step_id, y=utility)) + 
        geom_hline(data=dat%>%group_by(wave_id)%>%summarize(mean=mean(utility, na.rm=T)), aes(yintercept=mean), linetype=2, col='black' ) +
        facet_wrap( ~ wave_id)
      if(show_utility_points)
        plt <- plt + geom_point(aes(color=strategy), alpha=.25, shape=1, size=2)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(exists(smooth_method))
        plt <- plt + geom_smooth(aes(linetype=actor_id, color=strategy), method = smooth_method, linewidth=1, alpha=.15)
      #
      plt <- plt + theme_bw()
      #
      # Add marginal density plots
      plt <- ggMarginal(plt, type = "density", margins = "y")
      
      #
      print(plt)
      #
      if(return_plot)
        return(plt)
    },
    
    
    ##
    search_rsiena_multiwave_plot_actor_utility_density_by_strategy = function(thin_wave_factor=1, return_plot=FALSE) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <- as.factor( self$strat_1_coCovar )
      #
      strateffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
      stratparams <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$parameter)
      stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)ifelse(x$fix,'(Fix)',''))
      # stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)substr(as.character(x$fix),1,1))
      #
      structeffs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$effect)
      structparams <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$parameter)
      structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)ifelse(x$fix,'(Fix)',''))
      # structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)substr(as.character(x$fix),1,1))
      ## Compare 2 actors utilty
      dat <- self$actor_wave_util %>% 
        filter(wave_id %% thin_wave_factor == 0 ) %>% 
        mutate(
          strategy = actor_strat[ actor_id ] ,
          chain_below_med =  chain_step_id < median(chain_step_id)
        )
      dat$chain_half <- factor(ifelse(dat$chain_below_med, '1st Half', '2nd Half'))
      ##filter(chain_step_id %% thin_factor == 0) %>% 
      # dat$chain_half <- factor(1 + 1*(dat$chain_step_id >= median(dat$chain_step_id)), levels = c(2,1)) ## reverse order for linetype_1 used for Half_2
      # dat <- self$actor_util_df ##%>% filter(chain_step_id %% thin_factor == 0)
      ##
      stratmeans <- dat %>% group_by(strategy, chain_half, wave_id) %>% 
        summarize(n=n(), mean=mean(utility, na.rm=T), sd=sd(utility, na.rm=T))
      ## Actor density fact plots comparing H1 to H2 utility distribution
      plt <- ggplot(dat, aes(x=utility, color=strategy, fill=strategy)) + ##linetype=chain_half
        geom_density(alpha=.1, linewidth=1)  +
        # geom_histogram(alpha=.1, position = 'dodge') +
        facet_grid( wave_id ~ chain_half) +
        geom_vline(data = stratmeans, aes(xintercept = mean, color=strategy), linetype=2, linewidth=.9) +
        # geom_vline(xintercept = 0, linetype=2, linewidth=.9, color='gray')+
        theme_bw() + 
        ggtitle(sprintf('Strategy: %s\nStructure: %s', 
                        paste( paste(paste(strateffs, stratparams, sep='='), stratfixs, sep='' ), collapse = '; '),
                        paste( paste(paste(structeffs, structparams, sep='='), structfixs, sep=''), collapse = '; ')
                        ))
      print(plt)
      if(return_plot)
        return(plt)
    },
    
   
    
    # Convenience function for plotting all relevant plots 
    # or one plot by specifying plot 
    search_rsiena_plot = function(type=NA, 
                                  rolling_window = 10, 
                                  actor_ids=c(), 
                                  thin_factor=1, 
                                  smooth_method='loess',
                                  show_utility_points=TRUE,
                                  return_plot=FALSE,
                                  append_plot=FALSE,
                                  histogram_position='identity'
                                  ) {
      # if (all(is.na(type)) |  'ministep_count' %in% type)
      #   return(self$search_rsiena_plot_actor_ministep_count(rolling_window))
      # if (all(is.na(type)) | type == 'utility')
      #   return(self$search_rsiena_plot_actor_utility(rolling_window))
      plist <- list()
      if (all(is.na(type)) |  'utility' %in% type)
        plist[['utility']] <- self$search_rsiena_plot_actor_utility(actor_ids, thin_factor, smooth_method, show_utility_points, return_plot=T)
      
      if (all(is.na(type)) |  'utility_by_strategy' %in% type)
        plist[['utility_by_strategy']] <- self$search_rsiena_plot_actor_utility_by_strategy(actor_ids, thin_factor, smooth_method, show_utility_points, return_plot=T)
      
      if (all(is.na(type)) |  'utility_density' %in% type)
        plist[['utility_density']] <- self$search_rsiena_plot_actor_utility_density(return_plot=T)
      
      if (all(is.na(type)) |  'utility_density_by_strategy' %in% type)
        plist[['utility_density_by_strategy']] <- self$search_rsiena_plot_actor_utility_density_by_strategy(return_plot=T)
      
      if (all(is.na(type)) |  'utility_histogram_by_strategy' %in% type)
        plist[['utility_histogram_by_strategy']] <- self$search_rsiena_plot_actor_utility_histogram_by_strategy(histogram_position, return_plot=T)
      
      if (all(is.na(type)) |  'stability' %in% type)
        plist[['stability']] <- self$search_rsiena_plot_stability() ##**TODO** Fix return plot
      
      #  SET plots 
      self$plots <- if(append_plot) { append(self$plots, plist) } else { plist }
      
      # stop(sprintf('Plot type not implemented: %s', type))
      if(return_plot)
        return(plist)
    },
    # ##**TODO** Action count data frames (explore, exploit) 
    # ##          rows (actors) by columns (iterations), counts: 0,1,2,3,...
    # m1$chain_stats %>% group_by(iteration_id) %>% count() %>% plot(type='l')
    # m1$chain_stats %>% group_by(X__dv_name, X__id_from) %>% count() %>% print(n=60)
    # ##
    # util_after <- m1$chain_stats$X__utility_after
    # plot(util_after, type='l')
    # plot(cummean(util_after), type='l')
    
    ##
    search_rsiena_plot_actor_utility = function(actor_ids=c(), thin_factor=1, 
                                                smooth_method='loess', ##"lm", "glm", "gam", "loess","auto"
                                                show_utility_points=TRUE,
                                                return_plot=FALSE
                                                ) {
      ## Compare 2 actors utilty
      dat <- self$actor_util_df %>% filter(chain_step_id %% thin_factor == 0)
      if(length(actor_ids))
        dat <- dat%>%filter(actor_id %in% actor_ids)
      plt <- ggplot(dat, aes(x=chain_step_id, y=utility, color=actor_id)) + 
        geom_hline(yintercept = mean(dat$utility, na.rm=T), linetype=2, col='black' )
      if(show_utility_points)
        plt <- plt + geom_point(alpha=.25, shape=1, size=2)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(exists(smooth_method))
        plt <- plt + geom_smooth(aes(linetype=actor_id), method = smooth_method, linewidth=1, alpha=.15)
      #
      plt <- plt + theme_bw()
      #
      print(plt)
      #
      if(return_plot)
        return(plt)
    },
    
    ##
    search_rsiena_plot_actor_utility_by_strategy = function(actor_ids=c(), thin_factor=1, 
                                                        smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                        show_utility_points=TRUE,
                                                        return_plot=FALSE
                                                        ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <- as.factor( self$strat_1_coCovar )
      ## Compare 2 actors utilty
      dat <- self$actor_util_df %>% 
        filter(chain_step_id %% thin_factor == 0) %>% 
        mutate(strategy = actor_strat[ actor_id ] )
      if(length(actor_ids))
        dat <- dat %>% filter(actor_id %in% actor_ids)
      plt <- ggplot(dat, aes(x=chain_step_id, y=utility)) + 
        geom_hline(yintercept = mean(dat$utility, na.rm=T), linetype=2, col='black' )
      if(show_utility_points)
        plt <- plt + geom_point(aes(color=strategy), alpha=.25, shape=1, size=2)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(exists(smooth_method))
        plt <- plt + geom_smooth(aes(linetype=actor_id, color=strategy), method = smooth_method, linewidth=1, alpha=.15)
      #
      plt <- plt + theme_bw()
      #
      print(plt)
      #
      if(return_plot)
        return(plt)
    },
    
    ##
    search_rsiena_plot_actor_utility_density = function(return_plot=FALSE) {
      #
      dat <- self$actor_util_df ##%>% filter(chain_step_id %% thin_factor == 0)
      ## Actor density fact plots comparing H1 to H2 utility distribution
      dat$chain_half <- factor(1 + 1*(dat$chain_step_id >= median(dat$chain_step_id)), levels = c(2,1)) ## reverse order for linetype_1 used for Half_2
      plt <- ggplot(dat, aes(x=utility, color=actor_id, fill=actor_id, linetype=chain_half)) + 
         geom_density(alpha=.1, linewidth=1)  + 
         facet_wrap( ~ actor_id)+
         theme_bw()
      print(plt)
      if(return_plot)
        return(plt)
    },
    
    ##
    search_rsiena_plot_actor_utility_density_by_strategy = function(return_plot=FALSE) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <- as.factor( self$strat_1_coCovar )
      ## Compare 2 actors utilty
      dat <- self$actor_util_df %>%  
        mutate(
          strategy = actor_strat[ actor_id ] ,
          chain_below_med =  chain_step_id < median(chain_step_id)
        )
      dat$chain_half <- factor(ifelse(dat$chain_below_med, '1st Half', '2nd Half'))
      ##filter(chain_step_id %% thin_factor == 0) %>% 
      # dat$chain_half <- factor(1 + 1*(dat$chain_step_id >= median(dat$chain_step_id)), levels = c(2,1)) ## reverse order for linetype_1 used for Half_2
      # dat <- self$actor_util_df ##%>% filter(chain_step_id %% thin_factor == 0)
      ##
      stratmeans <- dat %>% group_by(strategy, chain_half) %>% 
        summarize(n=n(), mean=mean(utility, na.rm=T), sd=sd(utility, na.rm=T))
      ## Actor density fact plots comparing H1 to H2 utility distribution
      plt <- ggplot(dat, aes(x=utility, color=strategy, fill=strategy)) + ##linetype=chain_half
        geom_density(alpha=.1, linewidth=1)  +
        # geom_histogram(alpha=.1, position = 'dodge') +
        facet_grid( chain_half ~ .) +
        geom_vline(data = stratmeans, aes(xintercept = mean, color=strategy), linetype=2, linewidth=.9) +
        # geom_vline(xintercept = 0, linetype=2, linewidth=.9, color='gray')+
        theme_bw() 
      print(plt)
      if(return_plot)
        return(plt)
    },
    
    ##
    search_rsiena_plot_actor_utility_histogram_by_strategy = function(histogram_position='identity', 
                                                                      return_plot=FALSE
                                                                      ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <- as.factor( self$strat_1_coCovar )
      ## Compare 2 actors utilty
      dat <- self$actor_util_df %>%  
        mutate(
          strategy = actor_strat[ actor_id ] ,
          chain_below_med =  chain_step_id < median(chain_step_id)
        )
      dat$chain_half <- factor(ifelse(dat$chain_below_med, '1st Half', '2nd Half'))
      ##
      stratmeans <- dat %>% group_by(strategy, chain_half) %>% 
        summarize(n=n(), mean=mean(utility, na.rm=T), sd=sd(utility, na.rm=T))
      ## Actor density fact plots comparing H1 to H2 utility distribution
      plt <- ggplot(dat, aes(x=utility, color=strategy, fill=strategy)) + ##linetype=chain_half
        geom_histogram(alpha=.1, position = histogram_position) +
        facet_grid( chain_half ~ strategy ) +
        geom_vline(data = stratmeans, aes(xintercept = mean, color=strategy), linetype=2, linewidth=.9) +
        # geom_vline(xintercept = 0, linetype=2, linewidth=.9, color='gray')+
        theme_bw()
      print(plt)
      if(return_plot)
        return(plt)
    },
    
    ##
    search_rsiena_plot_actor_ministep_count = function(rolling_window = 20) {
      if (is.null(self$chain_stats))
        stop("Run RSiena to set chain before plotting chain stats.")
      ## Get Actor-period count of ministep decisions
      step_summary_df <- self$chain_stats %>% mutate(chain_step_id = row_number()) %>%
        filter( dv_varname =='self$bipartite_rsienaDV') %>% 
        group_by(id_from)
      #
      iters_with_change <- unique(step_summary_df$chain_step_id)
      #
      cnt_mat <- matrix(0, nrow=self$M, ncol=length(self$rsiena_model$chain) )
      for (iter in iters_with_change) {
        cat(sprintf(' %s ', iter))
        iter_df <- step_summary_df %>% filter(chain_step_id == iter) ##%>% mutate(X__id_from=X__id_from+1)
        actor_cnts <- iter_df %>% group_by(id_from) %>% count()
        cnt_mat[ actor_cnts$id_from, iter ] <- actor_cnts$n
      }
      # View(.)
      
      cnt_df <- as.data.frame(cnt_mat)
      colnames(cnt_df) <- paste0("Period_", 1:length(self$rsiena_model$chain) )
      cnt_df$Actor <- 1:self$M
      rownames(cnt_df) <- paste0("Actor_", cnt_df$Actor)
      
      long_data <- cnt_df %>% 
        pivot_longer(cols = starts_with("Period_"),
                     names_to = "period",
                     values_to = "count") %>% 
        mutate(
          period = as.numeric(gsub("Period_", "", period)),
          actor = as.numeric(gsub("Actor_", "", Actor))
        ) %>%
        select(actor, period, count)  # Reorder columns
      
      ## plot actor event sequences
      plt <- ggplot(long_data, aes(y = factor(actor), x=period, fill = count)) + 
        geom_tile() + 
        scale_fill_gradient(low = "white", high = "blue", name = "Ministeps Count") +
        theme_bw() + labs('y' = 'Actor')
      
      print(plt)
      
      return(plt)
      # ggplot(long_data, aes(x=period, y=count, color=factor(actor))) + geom_point() + geom_smooth()
    },
    

    
    # ##
    # search_rsiena_plot_actor_utility = function(rolling_window = 20) {
    #   
    #   if (is.null(self$chain_stats))
    #     stop("Run RSiena to set chain before plotting chain stats.")
    #   
    #   ## Get timeseries of actor utility function value
    #   util_step_df <- self$chain_stats %>% filter( X__dv_name=='self$bipartite_rsienaDV')
    #   n_util_steps <- nrow(util_step_df)
    #   util_mat <- matrix(NA, nrow=self$M, ncol=length(self$rsiena_model$chain) )
    #   for (i in 1:n_util_steps) {
    #     cat(sprintf(' %s ', i))
    #     x <- util_step_df[i, ]
    #     ## x$X__id_from <- as.numeric(x$X__id_from) + 1
    #     util_mat[ x$X__id_from , x$iteration_id ] <- x$X__utility_after
    #     
    #   }
    #   # Function to fill forward values along rows
    #   fill_forward_func <- function(row) {
    #     for (i in 2:length(row)) {
    #       if (is.na(row[i])) {
    #         row[i] <- row[i - 1]
    #       }
    #     }
    #     return(row)
    #   }
    #   # Apply the fill-forward function to each row
    #   util_mat_filled <- t(apply(util_mat, 1, fill_forward_func))
    #   
    #   util_df_fill <- as.data.frame(util_mat_filled)
    #   colnames(util_df_fill) <- paste0("Period_", 1:length(self$rsiena_model$chain) )
    #   util_df_fill$Actor <- 1:self$M
    #   rownames(util_df_fill) <- paste0("Actor_", util_df_fill$Actor)
    #   
    #   long_util <- util_df_fill %>% 
    #     pivot_longer(cols = starts_with("Period_"),
    #                  names_to = "period",
    #                  values_to = "count") %>% 
    #     mutate(
    #       period = as.numeric(gsub("Period_", "", period)),
    #       actor = as.numeric(gsub("Actor_", "", Actor))
    #     ) %>%
    #     select(actor, period, count)  # Reorder columns
    #   
    #   # ggplot(long_util, aes(x=period, y=count, color=factor(actor))) + 
    #   #   geom_point() + geom_smooth() + theme_bw()
    #   # 
    #   # ggplot(long_util, aes(x=period, y=count, color=factor(actor))) +  # geom_point(shape=1, alpha=.5) + 
    #   #   geom_line(linewidth=.8) + # geom_smooth() + 
    #   #   theme_bw() + labs(y='Utility', legend='Actor')
    #   
    #   # compute utility moving average
    #   long_util_ma <- long_util %>%
    #     group_by(actor) %>%
    #     mutate(rolling_avg = zoo::rollmean(count, k = rolling_window, fill = NA, align = "right"))
    #   #
    #   util_avg <- long_util_ma %>% group_by(period) %>% summarize(mean=mean(rolling_avg, na.rm=T))
    #   # plot utility moving average
    #   plt <- ggplot(long_util_ma, aes(x=period, y=rolling_avg, color=factor(actor))) +  
    #     geom_point(shape=1, alpha=.1) + 
    #     geom_line(linewidth=.9, alpha=.85) + # geom_smooth() + 
    #     geom_line(data = util_avg, linewidth=2, color='black', aes(x=period, y=mean, color='Mean')) + 
    #     theme_bw() + labs(y='Utility', color='Actor')
    #   
    #   print(plt)
    #   return(plt)
    #   
    #   
    # },
    
    

    # ## Batch of Multiple Simulation Extensions to 
    # search_rsiena_batchrun = function(batches=10, iterations=100, plot_save = TRUE, overwrite=FALSE, rsiena_phase2_nsub=1) {
    #   
    #   for (i in 1:batches) {
    #   
    #     overwrite_by_batch <- ifelse(i == 1, TRUE, overwrite)
    #     if(overwrite_by_batch) {
    #       self$search_rsiena_extend(iterations)
    #     } else {
    #       self$search_rsiena_init(iterations, rsiena_phase2_nsub=rsiena_phase2_nsub)
    #     }
    # 
    #    
    #     #rsiena_run(iterations=iterations, plot_save=plot_save, overwrite=overwrite_by_batch)
    #     
    #     ## Take snapshot of evolving system at this iteration of progress
    #     self$visualize_system_bi_env_rsiena(plot_save = plot_save)
    #     
    #   }
    #   
    # },
    
    
    ##
    ##
    ##
    ##**TODO**
    ##**CREATE VISUALIZE_NETWORKS plots for the RSiena approach search_rsiena_batchrun **
    ##
    ##
    ##
    # ===================================================================
    #   Model 1             
    # -------------------------------------------------------------------
    #   basic rate parameter self$social_rsienaDV          0.000 (0.000) ***
    #   self$social_rsienaDV: degree (density)            -2.000 (0.000) *** [suppressed K_S]
    #   basic rate parameter self$search_rsienaDV          0.045 (0.000) ***
    #   self$search_rsienaDV: degree (density)            -0.557 (0.000) *** [suppressed K_E]
    #   basic rate parameter self$bipartite_rsienaDV       0.000 (0.000) ***
    #   self$bipartite_rsienaDV: outdegree (density)      57.398 (0.000) ***
    #   self$bipartite_rsienaDV: indegree - popularity     8.086 (0.000) ***
    #   self$bipartite_rsienaDV: outdegree - activity     -2.389 (0.000) ***
    #   -------------------------------------------------------------------
    #   Iterations                                     2307                
    # ===================================================================
    
    # Visualization of the search progress
    ####
    # plot_search_progress = function() {
    #   if (length(self$progress_scores) == 0) {
    #     stop("No search progress data available. Run a search method first.")
    #   }
    #   plot.new()
    #   plot(self$progress_scores, type = "l", col = "blue", lwd = 2,
    #        main = "Search Progress",
    #        xlab = "Iteration", ylab = "Interaction Score",
    #        ylim = range(self$progress_scores))
    # },
    ####
    # Plot fitness history over iterations for each actor
    plot_fitness_progress = function(plot_save = FALSE) {
      if (length(self$fitness_history) == 0) {
        stop("No fitness history recorded. Run the simulation first.")
      }
      
      fitness_df <- do.call(rbind, lapply(1:length(self$fitness_history), function(iter) {
        data.frame(Iteration = iter, Actor = 1:self$M, Fitness = self$fitness_history[[iter]])
      }))
      
      (plt <- ggplot(fitness_df, aes(x = Iteration, y = Fitness, color = factor(Actor))) +
        geom_line() +
        labs(title = "Fitness Progress Over Iterations",
             x = "Iteration",
             y = "Fitness",
             color = "Actor") +
        theme_minimal()
      )
      
      # Calculate average degree (K) for the social space and component interaction space
      avg_degree_social <- mean(degree(self$get_social_igraph()))
      avg_degree_component <- mean(degree(self$get_component_igraph()))
      
      if (plot_save) {
        keystring <- sprintf("%s_sim%.0f_iter%.0f_N%d_M%d_BI_PROB_%.2f_K_S_%.2f_K_C_%.2f", 
                             self$SIM_NAME, self$TIMESTAMP, self$ITERATION, self$N, self$M, self$BI_PROB, avg_degree_social, avg_degree_component)
        ggsave(filename = sprintf('fitness_SAOM-NK_networks_%s_%s.png', keystring, self$TIMESTAMP), 
               plot = plt, height = 4.5, width = 9, units = 'in', dpi = 300)
      }
    },
  
  
    # Plot time series of mean and 95% quantile of degrees for K_S and K_E
    plot_degree_progress_from_sims = function(K_soc_list, K_env_list, plot_save = FALSE) {
      
      n <- length(K_soc_list)
      
      degree_summary <- do.call(rbind, lapply(1:n, function(iter) {
        data.frame(
          Iteration = iter,
          Mean_K_S = mean(K_soc_list[[iter]]),
          Q25_K_S = quantile(K_soc_list[[iter]], 0.25),
          Q75_K_S = quantile(K_soc_list[[iter]], 0.75),
          Mean_K_E = mean(K_env_list[[iter]]),
          Q25_K_E = quantile(K_env_list[[iter]], 0.25),
          Q75_K_E = quantile(K_env_list[[iter]], 0.75)
        )
      }))
      
      (plt <- ggplot(degree_summary, aes(x = Iteration)) +
          geom_line(aes(y = Mean_K_S, color = "Mean K_S")) +
          geom_ribbon(aes(ymin = Q25_K_S, ymax = Q75_K_S, fill = "K_S"), alpha = 0.1) +
          geom_line(aes(y = Mean_K_E, color = "Mean K_E")) +
          geom_ribbon(aes(ymin = Q25_K_E, ymax = Q75_K_E, fill = "K_E"), alpha = 0.1) +
          scale_color_manual(values = c("Mean K_S" = "blue", "Mean K_E" = "red")) +
          scale_fill_manual(values = c("K_S" = "blue", "K_E" = "red")) +
          labs(title = "Degree Progress Over Iterations",
               x = "Iteration",
               y = "Degree",
               color = "Mean Degree",
               fill = "IQR (Mid-50%)") +
          theme_minimal()
      )
      
      # Calculate average degree (K) for the social space and component interaction space
      avg_degree_social <- mean(degree_summary$Mean_K_S)
      avg_degree_component <- mean(degree_summary$Mean_K_E)
      
      if (plot_save) {
        keystring <- sprintf("K_S_%.2f_K_C_%.2f", 
                              avg_degree_social, avg_degree_component)
        ggsave(filename = sprintf('Ks_degree_SAOM-NK_networks_%s_%s.png',
                                  keystring, round(as.numeric(Sys.time())*100)), 
               plot = plt, height = 4.5, width = 9, units = 'in', dpi = 300)
      }
    },
    
    
    # Plot time series of mean and 95% quantile of degrees for K_S and K_E
    plot_degree_progress = function(plot_save = FALSE) {
      if (length(self$degree_history_K_S) == 0 || length(self$degree_history_K_E) == 0) {
        stop("No degree history recorded. Run the simulation first.")
      }
      
      degree_summary <- do.call(rbind, lapply(1:self$ITERATION, function(iter) {
        data.frame(
          Iteration = iter,
          Mean_K_S = mean(self$degree_history_K_S[[iter]]),
          Q25_K_S = quantile(self$degree_history_K_S[[iter]], 0.25),
          Q75_K_S = quantile(self$degree_history_K_S[[iter]], 0.75),
          Mean_K_E = mean(self$degree_history_K_E[[iter]]),
          Q25_K_E = quantile(self$degree_history_K_E[[iter]], 0.25),
          Q75_K_E = quantile(self$degree_history_K_E[[iter]], 0.75)
        )
      }))
      
      (plt <- ggplot(degree_summary, aes(x = Iteration)) +
        geom_line(aes(y = Mean_K_S, color = "Mean K_S")) +
        geom_ribbon(aes(ymin = Q25_K_S, ymax = Q75_K_S, fill = "K_S"), alpha = 0.1) +
        geom_line(aes(y = Mean_K_E, color = "Mean K_E")) +
        geom_ribbon(aes(ymin = Q25_K_E, ymax = Q75_K_E, fill = "K_E"), alpha = 0.1) +
        scale_color_manual(values = c("Mean K_S" = "blue", "Mean K_E" = "red")) +
        scale_fill_manual(values = c("K_S" = "blue", "K_E" = "red")) +
        labs(title = "Degree Progress Over Iterations",
             x = "Iteration",
             y = "Degree",
             color = "Mean Degree",
             fill = "IQR (Mid-50%)") +
        theme_minimal()
      )
      
      # Calculate average degree (K) for the social space and component interaction space
      avg_degree_social <- mean(degree(self$get_social_igraph()))
      avg_degree_component <- mean(degree(self$get_component_igraph()))
      
      if (plot_save) {
        keystring <- sprintf("%s_sim%.0f_iter%.0f_N%d_M%d_BI_PROB_%.2f_K_S_%.2f_K_C_%.2f", 
                             self$SIM_NAME, self$TIMESTAMP, self$ITERATION, self$N, self$M, self$BI_PROB, avg_degree_social, avg_degree_component)
        ggsave(filename = sprintf('Ks_degree_SAOM-NK_networks_%s_%s.png', keystring, self$TIMESTAMP), 
               plot = plt, height = 4.5, width = 9, units = 'in', dpi = 300)
      }
    },
    
    
    
    # Visualization using ggplot2 and ggraph: bipartite network, social network, and component interaction matrix heatmap
    visualize_system_bi_env_rsiena = function(plot_save = FALSE) {
      #
      RSIENA_ITERATION <- length(self$rsiena_model$sims)
      
      # Generate labels for components in letter-number sequence
      generate_component_labels <- function(n) {
        letters <- LETTERS  # Uppercase alphabet letters
        if (n <= 26)
          return(letters[1:n])
        labels <- c()
        repeat_count <- ceiling(n / length(letters))
        for (i in 1:repeat_count) {
          labels <- c(labels, paste0(letters, i))
        }
        return(labels[1:n])  # Return only as many as needed
      }
      
      component_labels <- generate_component_labels(self$N)
      
      # 1. Bipartite network plot using ggraph with vertex labels
      ig_bipartite <- self$bipartite_igraph
      
      # Set node attributes for shape and label
      # V(ig_bipartite)$type <- V(ig_bipartite)$type
      V(ig_bipartite)$shape <- ifelse(V(ig_bipartite)$type, "square", "circle")
      V(ig_bipartite)$color <- ifelse(V(ig_bipartite)$type, "lightblue", "darkorange")
      V(ig_bipartite)$label <- c( 1:self$M, component_labels)  # Numeric labels for actors
      V(ig_bipartite)$node_text_color <- ifelse(V(ig_bipartite)$type, 'black', 'lightblue') 
      
      bipartite_plot <- ggraph(ig_bipartite, layout = "fr") +
        geom_edge_link(color = "gray") +
        geom_node_point(aes(shape = shape, color = color), size = 5) +
        geom_node_text(aes(label = label, color=node_text_color), vjust = 0.5, hjust = 0.5, size = 3) +  # Center vertex labels over vertices
        scale_shape_manual(values = c("circle" = 16, "square" = 15)) +
        scale_color_identity() +
        labs(title = "[DGP] Bipartite Environment\n(Actors and Components)") +  # Split title into two lines
        theme_minimal() +
        theme(
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()
        )
      
      # 2. Social network plot using ggraph with different statistics for node color and size
      ig_social <- self$social_igraph
      
      # Calculate network statistics
      node_size <- degree(ig_social)  # Degree centrality for node size
      node_color <- eigen_centrality(ig_social)$vector  # Eigenvector centrality for node color
      node_text <- 1:vcount(ig_social)
      
      # V(ig_social)$label <- ifelse(V(ig_social)$type, 1:self$M, 1:self$M)
      
      social_plot <- ggraph(ig_social, layout = "fr") +
        geom_edge_link(color = "gray") +
        geom_node_point(aes(size = node_size, color = node_color)) +
        geom_node_text(aes(label = node_text), vjust = 0.5, hjust = 0.5, size = 3, color='white') +  # Center vertex labels over vertices
        scale_color_gradient(low = "green", high = "red") +
        labs(title = "[Proj1] Actor Social Network\n(Component Overlaps)", 
             color = "Eigenvector\nCentrality", size = "Degree\nCentrality") +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          legend.box = "vertical",  # Stack legends vertically
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank()
        ) +
        guides(
          color = guide_legend(order = 1, nrow = 2),  # Stacks color legend on two rows if needed
          size = guide_legend(order = 2, nrow = 2)    # Stacks size legend on two rows if needed
        )
      
      
      # 3. Component interaction matrix heatmap
      component_matrix <- self$search_matrix
      component_df <- melt(component_matrix)
      colnames(component_df) <- c("Component1", "Component2", "Interaction")
      
      ##**DEBUG**
      # stop('DEBUG print')

      # Replace component numbers with labels in the heatmap data frame
      component_df$Component1 <- factor(component_df$Component1, labels = component_labels)
      component_df$Component2 <- factor(component_df$Component2, labels = component_labels)
      
      heatmap_plot <- ggplot(component_df, aes(x = Component1, y = Component2, fill = Interaction)) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "red") +
        labs(title = "[Proj2] Component Interaction\n(Epistasis) Heatmap", x = "Component 1", y = "Component 2") +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      # Calculate average degree (K) for the social space and component interaction space
      avg_degree_social <- mean(degree(ig_social))
      avg_degree_component <- mean(degree(graph_from_adjacency_matrix(component_matrix)))
      
      # Create a main title using sprintf with simulation parameters
      main_title <- sprintf("Simulation Parameters:\niter = %d, N = %d, M = %d, BI_PROB = %.2f, K_S = %.2f, K_C = %.2f", 
                            RSIENA_ITERATION, self$N, self$M, self$BI_PROB, avg_degree_social, avg_degree_component)
      
      # Arrange plots with a main title
      (plt <- grid.arrange(
        social_plot, bipartite_plot, heatmap_plot, ncol = 3,
        top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold"))
      ))
      if(plot_save){
        keystring <- sprintf("%s_sim%.0f_iter%.0f_N%d_M%d_BI_PROB_%.2f_K_S_%.2f_K_C_%.2f", 
                             self$SIM_NAME, self$TIMESTAMP, RSIENA_ITERATION, self$N, self$M, self$BI_PROB, avg_degree_social, avg_degree_component)
        ggsave(filename = sprintf('3plot_SAOM-NK_networks_%s_%s.png', keystring, self$TIMESTAMP), 
               plot = plt, height = 4.5, width = 9, units = 'in', dpi = 300)
      }
    },
    
  
  
  
  
  # Visualization using ggplot2 and ggraph: bipartite network, social network, and component interaction matrix heatmap
  plot_bipartite_system_from_mat = function(bipartite_matrix, RSIENA_ITERATION, plot_save = FALSE) {
    # #
    # RSIENA_ITERATION <- length(self$rsiena_model$sims)
    
    N <- ncol(bipartite_matrix)
    M <- nrow(bipartite_matrix)
    
    TS <- round(as.numeric(Sys.time()) * 100)
    
    # Generate labels for components in letter-number sequence
    generate_component_labels <- function(n) {
      letters <- LETTERS  # Uppercase alphabet letters
      labels <- c()
      repeat_count <- ceiling(n / length(letters))
      for (i in 1:repeat_count) {
        labels <- c(labels, paste0(letters, i))
      }
      return(labels[1:n])  # Return only as many as needed
    }
    
    component_labels <- generate_component_labels(self$N)
    
   
    
    # 1. Bipartite network plot using ggraph with vertex labels
    ig_bipartite <- igraph::graph_from_biadjacency_matrix(bipartite_matrix, directed = T, weighted = T)
    
    # E(ig_bipartite)$weight <- 1
    projs <- igraph::bipartite_projection(ig_bipartite, multiplicity = T, which = 'both')
    ig_social <- projs$proj1
    ig_component <- projs$proj2
    
    # Set node attributes for shape and label
    # V(ig_bipartite)$type <- bipartite_mapping(ig_bipartite)$type
    V(ig_bipartite)$shape <- ifelse(V(ig_bipartite)$type, "square", "circle")
    V(ig_bipartite)$color <- ifelse(V(ig_bipartite)$type, "lightblue", "darkorange")
    V(ig_bipartite)$label <- c(as.character(1:self$M), as.character(component_labels))  # Numeric labels for actors
    
    
    # stop('DEBUG')
    
    V(ig_bipartite)$node_text_color <- ifelse(V(ig_bipartite)$type, 'black', 'lightblue') 
    
    
    bipartite_plot <- ggraph(ig_bipartite, layout = "fr") +
      geom_edge_link(color = "gray") +
      geom_node_point(aes(shape = shape, color = color), size = 5) +
      geom_node_text(aes(label = label, color=node_text_color), vjust = 0.5, hjust = 0.5, size = 3) +  # Center vertex labels over vertices
      scale_shape_manual(values = c("circle" = 16, "square" = 15)) +
      scale_color_identity() +
      labs(title = "[DGP] Bipartite Environment\n(Actors and Components)") +  # Split title into two lines
      theme_minimal() +
      theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
      )
    
    # 2. Social network plot using ggraph with different statistics for node color and size
    # ig_social <- self$get_social_igraph()
    
    # Calculate network statistics
    node_size <- igraph::degree(ig_social)  # Degree centrality for node size
    node_color <- eigen_centrality(ig_social)$vector  # Eigenvector centrality for node color
    node_text <- 1:self$M
    
    # V(ig_social)$label <- ifelse(V(ig_social)$type, 1:self$M, 1:self$M)
    
    social_plot <- ggraph(ig_social, layout = "fr") +
      geom_edge_link(color = "gray") +
      geom_node_point(aes(size = node_size, color = node_color)) +
      geom_node_text(aes(label = node_text), vjust = 0.5, hjust = 0.5, size = 3, color='white') +  # Center vertex labels over vertices
      scale_color_gradient(low = "green", high = "red") +
      labs(title = "[Proj1] Actor Social Network\n(Component Overlaps)", 
           color = "Eigenvector\nCentrality", size = "Degree\nCentrality") +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",  # Stack legends vertically
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
      ) +
      guides(
        color = guide_legend(order = 1, nrow = 2),  # Stacks color legend on two rows if needed
        size = guide_legend(order = 2, nrow = 2)    # Stacks size legend on two rows if needed
      )
    
    # 3. Component interaction matrix heatmap
    component_matrix <- igraph::as_adjacency_matrix(ig_component, type = 'both', sparse = F)
    
    
    component_df <- melt(component_matrix)
    colnames(component_df) <- c("Component1", "Component2", "Interaction")
    
    
    # Replace component numbers with labels in the heatmap data frame
    component_df$Component1 <- factor(component_df$Component1, labels = component_labels)
    component_df$Component2 <- factor(component_df$Component2, labels = component_labels)
    
    heatmap_plot <- ggplot(component_df, aes(x = Component1, y = Component2, fill = Interaction)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      labs(title = "[Proj2] Component Interaction\n(Epistasis) Heatmap", x = "Component 1", y = "Component 2") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Calculate average degree (K) for the social space and component interaction space
    avg_degree_social <- mean(igraph::degree(ig_social))
    avg_degree_component <- mean(igraph::degree(ig_component))
    
    # Create a main title using sprintf with simulation parameters
    main_title <- sprintf("Simulation Parameters:\niter = %d, N = %d, M = %d, BI_PROB = %.2f, K_S = %.2f, K_C = %.2f", 
                          RSIENA_ITERATION, N, M, self$BI_PROB, avg_degree_social, avg_degree_component)
    
    # Arrange plots with a main title
    (plt <- grid.arrange(
      social_plot, bipartite_plot, heatmap_plot, ncol = 3,
      top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold"))
    ))
    if(plot_save){
      keystring <- sprintf("%s_sim%.0f_iter%.0f_N%d_M%d_BI_PROB_%.2f_K_S_%.2f_K_C_%.2f", 
                           self$SIM_NAME, TS, RSIENA_ITERATION, N, M, self$BI_PROB, avg_degree_social, avg_degree_component)
      ggsave(filename = sprintf('3plot_SAOM_NK_networks_%s_%s.png', keystring, TS), 
             plot = plt, height = 4.5, width = 9, units = 'in', dpi = 300)
    }
  },
  
  
  
  
    
    
    # Interactive visualization using visNetwork
    visNetwork_social = function() {
      ig_social <- graph_from_adjacency_matrix(as.matrix.network(self$social_network), mode = 'directed')
      vis_data <- toVisNetworkData(ig_social)
      
      visNetwork(nodes = vis_data$nodes, edges = vis_data$edges) %>%
        visNodes(shape = "dot", size = 10) %>%
        visEdges(arrows = "to") %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visLayout(randomSeed = 123)
    },
    
    # Evaluate interactions between social and search spaces using the bipartite matrix
    evaluate_interspace_interactions = function() {
      bipartite_matrix <- as.matrix.network(self$bipartite_igraph)
      actor_influence_on_components <- bipartite_matrix %*% as.matrix.network(self$search_landscape)
      total_interaction_score <- sum(actor_rounfluence_on_components)
      return(total_interaction_score)
    }
  )
)


# ##-----------------



# This value is calculated using the following logic:
#   In a bipartite network with M actors and N components, the total possible ties is M * N.
# For a target density of 20% (0.2), the expected number of ties is 0.2 * M * N.
# The density parameter in RSiena is expressed as the log-odds of the probability of a tie.
# The calculation can be performed using the formula:
#   density parameter
# =
#   log
# ⁡
# (
#   target ties
#   total possible ties
#   −
#   target ties
# )
# density parameter=log( 
#   total possible ties−target ties
#   target ties
#   ​
# )
# For the given example where M = 10 actors, N = 20 components, and a target density of 20%, the resulting density parameter is -1.3862943611198906.


##-- Variable Covariate
# strat_mat <- matrix( c(
#   c(rep(0, 1000)), ## Actor 1
#   c(rep(0, 1000)), ## Actor 2
#   #
#   c(rep(0, 500), rep(1, 500)),  ## Actor 3
#   c(rep(0, 500), rep(-1, 500)), ## Actor 4
#   #
#   c(rep(1, 500), rep(0, 500)),
#   c(rep(-1, 500), rep(0, 500)),
#   #
#   c(rep(1, 1000)),
#   c(rep(-1, 1000))
# ), nrow=environ_params$M, byrow=TRUE)


#,
# varDyadCovar = list(
#   
# ),
# coCovar = list(
#   
# ),
# coDyadCovar = list(
#   
# ),
# dyadicCov = list(
#   
# )


## log( (0.2 * 10 * 20) / (0.8 * 10 * 20)  )
# (0.2 * m1$M * m1$N) / (1 - (0.2 * m1$M * m1$N) )

######
## log( BI_PROB / (1 - BI_PROB)  ) = Beta_density
## log( 0.2 / 0.8  ) = -1.386294 Beta_density







###############  SIM SETUP ####################################



## 1. Environment determines what types of entities and strategies are permissible
environ_params <- list(
  M = 9,        ## Actors
  N = 18,       ## Components
  BI_PROB = .15, ## Environmental Density (DGP hyperparameter)
  component_matrix_start = 'rand', ##**TODO** Implement: 'rand','modular','semi-modular',...
  rand_seed = 123,
  visualize_init = T,
  name = '_TESTrsienaPayoffs3_'#,
  ## use starting matrix param that is character ('modular',etc) or random (random seed 12345)
)


######### CONSTANT ACTOR STRATEGY ########
actor_strat <- c(-1,0,1,-1,0,1,-1,0,1)
component_payoffs <- runif(environ_params$N)

## 2. Strategies sets the objective function as a linear combination of network stats across DVs
structure_model <- list(
  dv_bipartite = list(
    name = 'self$bipartite_rsienaDV',
    effects = list( ##**STRUCTURAL EFFECTS -- dyadic/network endogeneity sources**
      list(effect='density', parameter= -0.1, fix=T, dv_name='self$bipartite_rsienaDV'), ##interaction1 = NULL
      list(effect='inPop',   parameter= .01,  fix=T, dv_name='self$bipartite_rsienaDV'), #interaction1 = NUL
      list(effect='outAct',  parameter= .01, fix=T, dv_name='self$bipartite_rsienaDV')#, #interaction1 = NULL
      # list(effect='outInAss', parameter=0, fix=F, dv_name='self$bipartite_rsienaDV'), #interaction1 = NULL
      # list(effect='cycle4', parameter=.5, fix=T, dv_name='self$bipartite_rsienaDV')#, #interaction1 = NULL
    ),
    coCovars = list( ##**STRATEGY -- MONADIC CONSTANT COVARIATE EFFECTS **
      list(effect='altX',parameter= 1, fix=T, dv_name='self$bipartite_rsienaDV',interaction1='self$component_1_coCovar',
           x = component_payoffs ),
      list(effect='outActX',parameter= .5, fix=T, dv_name='self$bipartite_rsienaDV',interaction1='self$component_2_coCovar',
           x = component_payoffs ), #interaction1 = NULL
      list(effect='egoX',parameter= .5, fix=T, dv_name='self$bipartite_rsienaDV',interaction1='self$strat_1_coCovar',
           x = actor_strat ), #interaction1 = NULL
      list(effect='inPopX',parameter= .5, fix=T, dv_name='self$bipartite_rsienaDV',interaction1='self$strat_2_coCovar',
           x = actor_strat )#, #interaction1 = NULL
      # list(effect=c('outActX','inPopX'),parameter= -.3, fix=T, dv_name='self$bipartite_rsienaDV',interaction1=c('self$strat_1_coCovar','self$strat_2_coCovar'),
      #      x = c(1,0,-1,1,0,-1,1,0,-1) )#, #interaction1 = NULL
      # list(effect='egoXinPop',parameter=.8, fix=T,dv_name='self$bipartite_rsienaDV',interaction1='self$strat_1_coCovar',
      #      x = c(1, 0, -1, 1, 0, -1, 1, 0) )#, #interaction1 = NULL
    ),
    varCovars = list() ##**MONADIC TIME-VARYING COVARIATE EFFECTS -- DYNAMIC STRATEGY PROGRAMS**
  )
)
# strat  <- c(rep(-1, round(environ_params$M / 2)), rep(1, round(environ_params$M / 2)))
# strat <- c(-1, 0, 1, -1, 0, 1, -1, 0)

################## SIM ANALYSIS ########################################
## 1. INIT SIM object
m1 <- SaomNkRSienaBiEnv$new(environ_params)
#
m1$search_rsiena_multiwave_run(structure_model, waves=1, iterations = 2500, rand_seed = 321)
#
m1$search_rsiena_multiwave_process_results()
#
m1$search_rsiena_multiwave_plot('utility_strategy_summary', thin_factor = 25)
#
m1$search_rsiena_multiwave_plot('K_4panel', thin_factor = 25)
#
#
m1$search_rsiena_multiwave_plot('K_A_strategy_summary', thin_factor = 1)
m1$search_rsiena_multiwave_plot('K_B1_strategy_summary', thin_factor = 1)
m1$search_rsiena_multiwave_plot('K_B2_strategy_summary', thin_factor = 1)
m1$search_rsiena_multiwave_plot('K_C_strategy_summary', thin_factor = 1)
#
m1$search_rsiena_multiwave_plot('K_4panel', thin_factor = 1)
#


#
#
m1$search_rsiena_multiwave_plot('utility_density_by_strategy', thin_wave_factor = 1)
#
#
# View(m1$actor_wave_stats)
#
m1$search_rsiena_multiwave_plot('utility_by_strategy', smooth_method = 'loess', show_utility_points = T, thin_wave_factor = 1)
m1$search_rsiena_multiwave_plot('utility_by_strategy', smooth_method = 'loess', show_utility_points = F, thin_wave_factor = 1)

 
m1$bipartite_matrix_waves[[1]]

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


self = m1

forebear_rate <- ( nrow(self$chain_stats) - length(self$rsiena_model$chain) ) / length(self$rsiena_model$chain)


# getStructModNetStatsFromBipartiteMat <- function(bi_env_mat) {
#   # eff <- m1$rsiena_effects[m1$rsiena_effects$include, ]
#   efflist <- self$config_structure_model$dv_bipartite$effects
#   # set the bipartite environment network matrix
#   # bi_env_mat <- self$bipartite_matrix
#   ### empty matrix to hold actor network statistics 
#   mat <- matrix(rep(0, m1$M * length(efflist) ), nrow=self$M, ncol=length(efflist) )
#   effnames <- sapply(efflist, function(x) x$effect, simplify = T)
#   effparams <- sapply(efflist, function(x) x$parameter, simplify = T)
#   colnames(mat) <- effnames
#   rownames(mat) <- 1:self$M
#   #
#   for (i in 1:length(efflist)) 
#   {
#     eff_name <- efflist[[ i ]]$effect
#     #
#     xActorDegree  <- rowSums(bi_env_mat, na.rm=T)
#     xComponentDegree  <- colSums(bi_env_mat, na.rm=T)
#     ## network statistics dataframe
#     #
#     if (eff_name == 'density' ) 
#     {
#       mat[ , i] <- c( xActorDegree )
#     } 
#     else if (eff_name == 'outAct' )
#     {
#       mat[ , i]  <- c( xActorDegree^2  )
#     } 
#     else if (eff_name == 'inPop' ) 
#     {
#       mat[ , i]  <- c( bi_env_mat %*% (xComponentDegree + 1) )  ## MxN %*% N
#     } 
#     else  
#     {
#       print(sprintf('Effect not yet implemented: `%s`', eff_name))
#     }
#   }
#   return(mat)
# }

# get_struct_mod_net_stats_list_from_bi_mat <- function(bi_env_mat, type='all') {
#   # eff <- m1$rsiena_effects[m1$rsiena_effects$include, ]
#   efflist <- self$config_structure_model$dv_bipartite$effects
#   # set the bipartite environment network matrix
#   # bi_env_mat <- self$bipartite_matrix
#   ### empty matrix to hold actor network statistics
#   df <- data.frame() #nrow=self$M, ncol=length(efflist)
#   effnames <- sapply(efflist, function(x) x$effect, simplify = T)
#   effparams <- sapply(efflist, function(x) x$parameter, simplify = T)
#   mat <- matrix(rep(0, m1$M * length(efflist) ), nrow=self$M, ncol=length(efflist) )
#   colnames(mat) <- effnames
#   rownames(mat) <- 1:self$M
#   #
#   for (i in 1:length(efflist))
#   {
#     eff_name <- efflist[[ i ]]$effect
#     #
#     xActorDegree  <- rowSums(bi_env_mat, na.rm=T)
#     xComponentDegree  <- colSums(bi_env_mat, na.rm=T)
#     ## network statistics dataframe
#     #
#     if (eff_name == 'density' )
#     {
#       stat <- c( xActorDegree )
#     }
#     else if (eff_name == 'outAct' )
#     {
#       stat <- c( xActorDegree^2 )
#     }
#     else if (eff_name == 'inPop' )
#     {
#       stat <- c( bi_env_mat %*% (xComponentDegree + 1) )
#     }
#     else
#     {
#       cat(sprintf('\n\nEffect not yet implemented: `%s\n\n`', eff_name))
#     }
#     #
#     mat[ , i] <- stat
#     #
#     df <- rbind(df, data.frame(
#       statistic = stat, 
#       actor_id = factor(1:self$M), 
#       effect_id = factor(i), 
#       effect_name = effnames[i] 
#     ))
#   }
#   #
#   if (type %in% c('df','data.frame'))   return(df)
#   if (type %in% c('mat','matrix'))      return(mat)
#   if (type %in% c('all','both','list')) return(list(df=df, mat=mat))
#   cat(sprintf('specified return type %s not found', type))
# }


##**TODO** Check if updating preallocated matrix is faster then rbind to data.frame rows ?
# get_chain_stats_list <- function() { ## 'all','statdf','bi_env_arr', 'util', 'util_diff')
#   if (is.null(self$chain_stats))
#     stop('chain_stats missing; run simulation with returnChains=TRUE before computing actor utility')
#   # ## actor Utility vector
#   # au <- c( mat %*% self$rsiena_model$theta )
#   # hist( util )
#   effnames <- sapply(self$config_structure_model$dv_bipartite$effects, function(x) x$effect)
#   ## remove chain entries where no change was made (i.e., keep if !stability )
#   tiechdf <- self$chain_stats[ !self$chain_stats$stability, ]
#   
#   ## get matrix timeseries and network statistics timeseries
#   nchains <- length(self$rsiena_model$chain)
#   theta   <- self$rsiena_model$theta
#   ntheta <- length(theta)
#   #
#   bi_env_arr <- array(NA, dim=c(self$M, self$N, nchains))
#   # bi_env_long <- data.frame()
#   statdf <- data.frame()
#   utildf <- data.frame()
#   util_diff <- data.frame()
#   #
#   # utility = util,
#   # chain_step_id = i, 
#   # actor_id = factor(1:self$M )
#   #
#   nrows_1step_stat <- ntheta * self$M
#   nrows_1step_util <- self$M
#   #
#   statmat_long  <- matrix(NA, nrow= (nchains * nrows_1step_stat), ncol=4) ## chain_step_id, actor_id, stat, value
#   utilmat_long  <- matrix(NA, nrow= (nchains * nrows_1step_util), ncol=3) ## chain_step_id, actor_id, utility
#   util_diff_long <- matrix(NA, nrow=(nchains * nrows_1step_util), ncol=3)
#   #
#   for (i in 1:nrow(tiechdf)) {
#     mstep <- tiechdf[i,]
#     cat(sprintf('\n %.2f%s i=%s, j=%s', 100*i/nrow(tiechdf),'%', mstep$id_from,  mstep$id_to))
#     ## update bipartite environment matrix for one step (toggle one dyad)
#     actor_row <- mstep$id_from
#     component_col <- (mstep$id_to - self$M)
#     bi_env_mat <- toggleBiMat(bi_env_mat, actor_row,  component_col)
#     #
#     # bi_env_mat
#     #
#     # Get New Statistics
#     statsobj <- get_struct_mod_net_stats_list_from_bi_mat( bi_env_mat )
#     #
#     statsdf  <- rbind(statsdf,  statsobj$df )
#     #
#     # stat_step_fill <- matrix(
#     #   c(
#     #     rep(i, nrows_1step_stat), 
#     #     rep(1:self$M, each=ntheta),
#     #     rep(1:ntheta, self$M), 
#     #     rep(NA, nrows_1step_stat)
#     #   ), 
#     #   ncol = 4 ##nrow = nrows_1step_stat
#     # )
#     stat_step_mat <-  expand.grid(chain_step_id=i, actor_id=1:self$M, stat_eff_id=1:ntheta, value=NA)
#     util_step_fill <- c()
#     
#     #
#     c(statsobj$mat)
#     #
#     statrow_ids <- ( 1 + (i-1)*nrows_1step_stat ):( i*nrows_1step_stat )
#     statmat_long[ statrow_ids, ] <- matrix(c(), nrow=nrows_1step_stat, ncol=4)
#     #
#     statsmat <- statsobj$mat  #get_n_by_m_mat_from_long_df
#     #
#     ## Add ministep updated bipartite environment to array
#     bi_env_arr[ , , i]  <- bi_env_mat
#     #
#     # stat_arr[ , , i] <- statsmat
#     ## Add utilities to array
#     util <- c( statsmat %*% theta )
#     utildf <- rbind(utildf, data.frame(
#       utility = util,
#       chain_step_id = i, 
#       actor_id = factor(1:self$M )
#     ))
#     util_diff <- rbind(util_diff, data.frame(
#       utility = if(i == 1){ NA } else {util - util_lag }, ## diff 
#       chain_step_id = i, 
#       actor_id = factor(1:self$M )
#     ))
#     ######
#     ## update utility lag for next period difference
#     util_lag <- util
#     ######
#   }
#   #####
#   return(list(
#     utildf = utildf,
#     util_diff = util_diff,
#     statsdf = statsdf,
#     bi_env_arr = bi_env_arr
#   ))
# }




# get_chain_stats_list <- function() { ## 'all','statdf','bi_env_arr', 'util', 'util_diff')
#   if (is.null(self$chain_stats))
#     stop('chain_stats missing; run simulation with returnChains=TRUE in order to compute actor utility')
#   # ## actor Utility vector
#   # au <- c( mat %*% self$rsiena_model$theta )
#   # hist( util )
#   effnames <- sapply(self$config_structure_model$dv_bipartite$effects, function(x) x$effect)
#   ## remove chain entries where no change was made (keep if !stability )
#   tiechdf <- self$chain_stats[ !self$chain_stats$stability, ]
#   
#   ## get matrix timeseries and network statistics timeseries
#   nchains <- length(self$rsiena_model$chain)
#   theta   <- self$rsiena_model$theta
#   #
#   bi_env_arr <- array(NA, dim=c(self$M, self$N, nchains))
#   #
#   stats_li <- list()
#   util_li  <- list()
#   util_diff_li <- list()
#   # bi_env_long <- data.frame()
#   # statdf <- data.frame()
#   # utildf <- data.frame()
#   # util_diff <- data.frame()
#   for (i in 1:nrow(tiechdf)) {
#     mstep <- tiechdf[i,]
#     cat(sprintf('\n %.2f%s i=%s, j=%s', 100*i/nrow(tiechdf),'%', mstep$id_from,  mstep$id_to))
#     ## update bipartite environment matrix for one step (toggle one dyad)
#     bi_env_mat <- toggleBiMat(bi_env_mat, mstep$id_from,  (mstep$id_to - self$M)  )
#     #
#     # bi_env_mat
#     #
#     # Get New Statistics
#     # statsobj <- get_struct_mod_net_stats_list_from_bi_mat( bi_env_mat )
#     statmat <- get_struct_mod_stats_mat_from_bi_mat( bi_env_mat )
#     #
#     step_statgrid <- expand.grid(chain_step_id=i, actor_id=1:self$M, effect_id=1:ntheta)
#     step_statgrid$value <- c( statmat )
#     #
#     stats_li[[i]] <- step_statgrid   ##**CHECK list**
#     #
#     # statsdf  <- rbind(statsdf,  statsobj$df )
#     #
#     # statsmat <- statsobj$mat  #get_n_by_m_mat_from_long_df
#     #
#     
#     #
#     # tmpstatdf <- as.data.frame( statsmat )
#     # tmpstatdf$chain_step_id <- i
#     # tmpstatdf$actor_id <- factor(1:self$M)
#     # statdf <- rbind(statdf,  tmpstatdf )
#     ## network statistics matrix added to array
#     # stat_long <- rbind(stat_long, statdf)
#     ## Add ministep updated bipartite environment to array
#     bi_env_arr[ , , i]  <- bi_env_mat
#     #
#     # stat_arr[ , , i] <- statsmat
#     ## Add utilities to array
#     util <- c( statsmat %*% theta )
#     # utildf <- rbind(utildf, data.frame(
#     #   utility = util,
#     #   chain_step_id = i, 
#     #   actor_id = factor(1:self$M )
#     # ))
#     ##
#     step_utilgrid <- expand.grid(chain_step_id=i, actor_id=1:self$M)
#     step_utilgrid$utility <- util
#     util_li[[i]] <- step_utilgrid  ##**CHECK list**
#     ##
#     ##
#     # util_diff <- rbind(util_diff, data.frame(
#     #   utility = if(i == 1){ NA } else {util - util_lag }, ## diff 
#     #   chain_step_id = i, 
#     #   actor_id = factor(1:self$M )
#     # ))
#     ##
#     step_util_diffgrid <- expand.grid(chain_step_id=i, actor_id=1:self$M)
#     step_util_diffgrid$utility <- if(i == 1){ NA } else { util - util_lag }
#     util_diff_li[[i]] <- step_util_diffgrid  ##**CHECK list**
#     ##
#     ######
#     ## update utility lag for next period difference
#     util_lag <- util
#     ######
#   }
#   #####
#   return(list(
#     stats_li = stats_li,
#     util_li = util_li,
#     util_diff_li = util_diff_li,
#     bi_env_arr = bi_env_arr
#   ))
# }
# 
# chain_stats_li <- self$get_chain_stats_list()
# 
# ## Actor Network Statistics long dataframe
# stats_long <- data.table::rbindlist( chain_stats_li$stats_li )
# stats_long$actor_id <- as.factor(stats_long$actor_id)
# stats_long$effect_name <- as.factor(effnames[ stats_long$effect_id ])
# ## Actor Utility  long dataframe
# util_long <- data.table::rbindlist( chain_stats_li$util_li ) 
# util_long$actor_id <- as.factor(util_long$actor_id)
# ## Actor Utility Difference long dataframe
# util_diff_long <- data.table::rbindlist( chain_stats_li$util_diff_li ) 
# util_diff_long$actor_id <- as.factor(util_diff_long$actor_id)
# 
# 
# csl <- get_chain_stats_list()


# utildf %>% group_by(actor_id) %>% mutate(rollmean= zoo::rollmean(utility, k=10)) %>%

# avgutil <- utildf %>% group_by(chain_step_id ) %>% summarize(mean=mean(utility)) %>% pull(mean)
# actutil1 <- utildf %>% filter(actor_id==1) %>% pull(utility)
# 
# 
# ff <- fft(avgutil, inverse = F)
# freq <- seq(0, 1, length.out = (length(avgutil)/2) + 1)
# ampl <- Mod(ff[1:((length(avgutil)/2)+1)])
# plot(freq, ampl, type='l', log='y')
# # mvfft() ## multivariate fft 
# 
# ## Compare 2 actors utilty
# util_df <- self$actor_stats
# ggplot(util_df%>%filter(actor_id %in% c(1,2)), aes(x=chain_step_id, y=utility, color=actor_id, shape=actor_id)) + 
#   geom_point(alpha=.2, size=3) + # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
#   theme_bw()
# 
# 
# ## Compare 2 actors utilty
# ggplot(utildf%>%filter(actor_id %in% c(1,2)), aes(x=chain_step_id, y=utility, color=actor_id, shape=actor_id)) + 
#   geom_point(alpha=.2, size=3) + # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
#   theme_bw()

## plot actors utility points and smoothed loess lines
thin_factor <- 1
ggplot(utildf%>%filter(chain_step_id%%thin_factor==0), aes(x=chain_step_id, y=utility, color=actor_id)) + 
  geom_point(alpha=.2, shape=1, size=2) +
  geom_smooth(aes(linetype=actor_id), method = 'loess', linewidth=.8, alpha=.15) + 
  theme_bw()
## thin factor = 10 (takes every 10th chain step)
thin_factor <- 10
ggplot(utildf%>%filter(chain_step_id%%thin_factor==0), aes(x=chain_step_id, y=utility, color=actor_id)) + 
  geom_point(alpha=.5, shape=1, size=3) +
  geom_smooth(aes(linetype=actor_id), method = 'loess', linewidth=.8, alpha=.15) + 
  theme_bw()


## plot actors utility smoothed loess lines (no points)
thin_factor <- 1
ggplot(utildf%>%filter(chain_step_id%%thin_factor==0), aes(x=chain_step_id, y=utility, color=actor_id)) + 
  geom_smooth(aes(linetype=actor_id), method = 'loess', linewidth=1, alpha=.1) + 
  geom_hline(yintercept =  mean(utildf%>%filter(chain_step_id%%thin_factor==0)%>%pull(utility))) +
  theme_bw()
#
thin_factor <- 10
ggplot(utildf%>%filter(chain_step_id%%thin_factor==0), aes(x=chain_step_id, y=utility, color=actor_id)) + 
  geom_smooth(aes(linetype=actor_id), method = 'loess', linewidth=1, alpha=.05) + 
  geom_hline(yintercept =  mean(utildf%>%filter(chain_step_id%%thin_factor==0)%>%pull(utility))) +
  theme_bw()

# ggplot(utildf, aes(x=chain_step_id, y=utility, color=actor_id)) + 
#   geom_line(alpha=.3, shape=1) + 
#   theme_bw()

## Densities of actor utilities
ggplot(utildf, aes(x=utility, color=actor_id, fill=actor_id, linetype=actor_id)) + 
  geom_density(alpha=.1, linewidth=1)  + 
  # facet_wrap(~chain_step_id)+ 
  theme_bw()

## Actor density fact plots comparing H1 to H2 utility distribution
utildf$chain_half <- factor(1 + 1*(utildf$chain_step_id >= median(utildf$chain_step_id)), levels = c(2,1)) ## reverse order for linetype_1 used for Half_2
ggplot(utildf, aes(x=utility, color=actor_id, fill=actor_id, linetype=chain_half)) + 
  geom_density(alpha=.1, linewidth=1)  + 
  facet_wrap(~actor_id)+
  theme_bw()

## Actor utility difference timeseries
ggplot(util_diff, aes(x=chain_step_id, y=utility, color=actor_id)) + 
  geom_line(alpha=.9, shape=1) +  facet_grid(actor_id ~ .) +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw()

## Actor utility difference timeseries
ggplot(util_diff, aes(x=utility, color=actor_id, fill=actor_id)) + 
  geom_histogram(alpha=.3, bins=35) + 
  facet_grid(.~actor_id ) +
  scale_y_log10() +
  # geom_hline(yintercept = 0, linetype=2) +
  theme_bw()


# ggplot(util_diff, aes(x=chain_step_id, y=utility, color=actor_id, fill=actor_id, linetype=actor_id)) + 
#   geom_line()  + 
#   # facet_wrap(~chain_step_id)+ 
#   theme_bw()

util_diff_wide <- tidyr::pivot_wider(util_diff, 
    values_from='utility', 
    names_from = 'chain_step_id',  
    id_cols = 'actor_id' 
  ) %>% 
  select( !c(actor_id) )

util_diff_t    <- t(util_diff_wide )
util_diff_m1_t <- t(util_diff_wide %>% select( !c(`1`) ))  ## drop NA's in first columns representing first chain step

util_diff_cor <- cor(util_diff_m1_t)
diag(util_diff_cor) <- 0
#
heatmap(util_diff_cor)
#
rownames(util_diff_cor) <- paste0('A',1:self$M)
colnames(util_diff_cor) <- paste0('A',1:self$M)
#
cor_long <- as.data.frame(util_diff_cor) %>%
  tibble::rownames_to_column(var = "Actor_1") %>%
  pivot_longer(cols = -Actor_1, names_to = "Actor_2", values_to = "correlation")

ggplot(cor_long, aes(x=Actor_1, y=Actor_2, fill=correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0) +
  geom_text(aes(label = round(correlation, 2)), size = 4, color='white') +  # Text labels
  labs(title = "Actor Utility Correlation Heatmap", x = "Actor  i", y = "Actor  j") +
  theme_bw()

## Get dataframe of nchains (rows) by n_summary_stats (columns) for summarying actor utility changes
util_diff_t_summary <- data.frame( 
  ministep_count_id = 1:nrow(util_diff_t),
  sum = rowSums( util_diff_t ),
  mean = rowMeans( util_diff_t ),
  sd = apply( util_diff_t, 1, sd ),
  # skew = apply( util_diff_t, 1, skewness ),
  # kurt = apply( util_diff_t, 1, kurtosis ),
  min = apply( util_diff_t, 1, min ),
  max = apply( util_diff_t, 1, max ),
  range_magnitude = apply( util_diff_t, 1, function(x) abs(diff(range(x))) ),
  cnt_change = apply( util_diff_t, 1, function(x) sum(x != 0) ),  ##sum of logical counts the TRUEs
  cnt_stable = apply( util_diff_t, 1, function(x) sum(x == 0) ),  ##sum of logical counts the TRUEs
  cnt_pos    = apply( util_diff_t, 1, function(x) sum(x > 0) ),   ##sum of logical counts the TRUEs
  cnt_neg  = apply( util_diff_t, 1, function(x) sum(x < 0) )   ##sum of logical counts the TRUEs
)
 ## convert to long for plotting  
util_diff_t_summary_long <- util_diff_t_summary %>% 
  pivot_longer(names_to = 'utility_stat', 
               values_to = 'utility_stat_value', 
               cols = c('sum', 'mean', 'sd', 'min', 'max', 'range_magnitude', 
                        'cnt_change', 'cnt_stable', 'cnt_pos', 'cnt_neg'))

# ggplot(util_diff_t_summary_long, aes(x=utility_stat_value, color=utility_stat, fill=utility_stat)) + 
#   geom_density(alpha=.2) + facet_wrap(~utility_stat, scales = 'free') + 
#   theme_bw()
### Histograph facet_wrap of statistics of actor_utilty change from one ministep of any other actor
ggplot(util_diff_t_summary_long, aes(x=utility_stat_value, color=utility_stat, fill=utility_stat)) + 
  geom_histogram(alpha=.2) + facet_wrap(~utility_stat, scales = 'free_x') + 
  geom_vline(xintercept = 0) +
  theme_bw()


# statmat_arr_long <- melt(statmat_arr, varnames = effnames, )

statmat_arr_long <- statmat_arr #%>% rename(actor = Var1, statistic = Var2, period = Var3, value = Freq)

# Preview the result
head(long_df_tidyr)

# melt(wavemat, varnames = c("Simulation", "Wave", "NetworkEffect"), value.name = "Mean")


# statdf <- ldply(statlist) 


# # eff <- m1$rsiena_effects[m1$rsiena_effects$include, ]
# efflist <- m1$config_structure_model$dv_bipartite$effects
# # mat <- matrix(NA, nrow=m1$M, ncol=length(efflist$efflist))
# statlist <- list()
# for (i in 1:length(efflist)) {
#   #
#   eff <- efflist[[ i ]]
#   bi_env_mat <- m1$bipartite_matrix
#   
#   xActorDegree  <- rowSums(bi_env_mat, na.rm=T)
#   xComponentDegree  <- colSums(bi_env_mat, na.rm=T)
#   ## network statistics dataframe
#   # xmat <- matrix(NA, )
#   #
#   if (eff$effect == 'density') {
#     statlist[['density']] <- xActorDegree 
#   } else if (eff$effect == 'outAct'){
#     statlist[['outAct']]  <- xActorDegree^2 
#   } else if (eff$effect == 'inPop') {
#     statlist[['inPop']]  <- c( bi_env_mat %*% (xComponentDegree + 1) )  ## MxN %*% N
#   } else  {
#     print(sprintf('Effect not yet implemented: `%s`', eff$effect))
#   }
# }
# statdf <- ldply(statlist) 



# ## 2. Strategies sets the objective function as a linear combination of network stats across DVs
# structure_model <- list(
#   dv_social = list(
#     dv_name = 'self$social_rsienaDV',
#     effects = list(
#       list(effect='density', parameter= 0, fix=T, dv_name='self$social_rsienaDV')#, #interaction1 = NULL
#       # list(effect='transTriads', parameter=1, fix=F, dv_name='self$social_rsienaDV')#, #interaction1 = NULL
#       # list(effect='density', parameter=0, fix=T),
#     )
#   ),
#   dv_search = list(
#     dv_name = 'self$search_rsienaDV',
#     effects = list(
#       list(effect='density', parameter= 0, fix=T, dv_name='self$search_rsienaDV')#, #interaction1 = NULL
#       # list(effect='transTriads', parameter=1, fix=F, dv_name='self$search_rsienaDV')#, #interaction1 = NULL
#       # list(effect='density', parameter=0, fix=T),
#     )
#   ),
#   dv_bipartite = list(
#     name = 'self$bipartite_rsienaDV',
#     effects = list(
#       list(effect='density', parameter= -2, fix=T, dv_name='self$bipartite_rsienaDV'), ##interaction1 = NULL
#       list(effect='inPop', parameter= 1, fix=F, dv_name='self$bipartite_rsienaDV'), #interaction1 = NUL
#       list(effect='outAct', parameter= -1, fix=F, dv_name='self$bipartite_rsienaDV')#, #interaction1 = NULL
#       # list(effect='outInAss', parameter=0, fix=F, dv_name='self$bipartite_rsienaDV'), #interaction1 = NULL
#       # list(effect='cycle4', parameter=.5, fix=T, dv_name='self$bipartite_rsienaDV')#, #interaction1 = NULL
#       # list(effect='density', parameter=0, fix=T),
#     )
#   )
# )






































































# getChangeContributions <- function(algorithm, data, effects)
# {
#   ## Gets the simulated statistics.
#   ## The following initializations data, effects, and model
#   ## for calling "getTargets" in "siena07.setup.h"
#   ## is more or less copied from "getTargets" in "getTargets.r".
#   ## However, some modifications have been necessary to get it to work.
#   f <- unpackData(data,algorithm)
#   
#   effects <- effects[effects$include,]
#   if (!is.null(algorithm$settings))
#   {
#     stop('not implemented: RI together with settings')
#     # effects <- addSettingsEffects(effects, algorithm)
#   }
#   else
#   {
#     effects$setting <- rep("", nrow(effects))
#   }
#   pData <- .Call(C_setupData, PACKAGE=pkgname,
#                  list(as.integer(f$observations)),
#                  list(f$nodeSets))
#   ## register a finalizer
#   ans <- reg.finalizer(pData, clearData, onexit = FALSE)
#   ans<- .Call(C_OneMode, PACKAGE=pkgname,
#               pData, list(f$nets))
#   ans <- .Call(C_Bipartite, PACKAGE=pkgname, # added 1.1-299
#                pData, list(f$bipartites))
#   ans<- .Call(C_Behavior, PACKAGE=pkgname, pData,
#               list(f$behavs))
#   ans<-.Call(C_ConstantCovariates, PACKAGE=pkgname,
#              pData, list(f$cCovars))
#   ans<-.Call(C_ChangingCovariates,PACKAGE=pkgname,
#              pData,list(f$vCovars))
#   ans<-.Call(C_DyadicCovariates,PACKAGE=pkgname,
#              pData,list(f$dycCovars))
#   ans<-.Call(C_ChangingDyadicCovariates,PACKAGE=pkgname,
#              pData, list(f$dyvCovars))
#   
#   storage.mode(effects$parm) <- 'integer'
#   storage.mode(effects$group) <- 'integer'
#   storage.mode(effects$period) <- 'integer'
#   
#   effects$effectPtr <- rep(NA, nrow(effects))
#   depvarnames <- names(data$depvars)
#   tmpeffects <- split(effects, effects$name)
#   myeffectsOrder <- match(depvarnames, names(tmpeffects))
#   ans <- .Call(C_effects, PACKAGE=pkgname, pData, tmpeffects)
#   pModel <- ans[[1]][[1]]
#   for (i in 1:length(ans[[2]]))
#   {
#     effectPtr <- ans[[2]][[i]]
#     tmpeffects[[i]]$effectPtr <- effectPtr
#   }
#   myeffects <- tmpeffects
#   for(i in 1:length(myeffectsOrder)){
#     myeffects[[i]]<-tmpeffects[[myeffectsOrder[i]]]
#   }
#   ans <- .Call(C_getTargets, PACKAGE=pkgname, pData, pModel, myeffects,
#                parallelrun=TRUE, returnActorStatistics=FALSE,
#                returnStaticChangeContributions=TRUE)
#   # See getTargets in siena07setup.cpp; also see rTargets in StatisticsSimulation.cpp
#   ans
# }
# 
# unpackData <- function(data, x)
# {
#   f <- NULL
#   observations<- data$observations
#   types <- sapply(data$depvars, function(x) attr(x, "type"))
#   f$nDepvars <- length(data$depvars)
#   oneModes <- data$depvars[types == "oneMode"]
#   Behaviors <- data$depvars[types == "behavior"]
#   continuousBehaviors <- data$depvars[types == "continuous"]
#   bipartites <- data$depvars[types == "bipartite"]
#   ## add the settings
#   # oneModes <- lapply(oneModes, function(depvar) {
#   #                    name <- attr(depvar, "name")
#   #                    if (name %in% names(x$settings)) {
#   #                      # attr(depvar, "settings") <- c("universal", "primary", x$settings[[name]])
#   #                      attr(depvar, "settings") <- c(x$settings[[name]])
#   #                    }
#   #                    depvar
#   #                    })
#   f$nets <- lapply(oneModes, function(x, n, comp)
#     unpackOneMode(x, n, comp),
#     n = observations, comp=data$compositionChange)
#   names(f$nets) <- names(oneModes)
#   f$bipartites <- lapply(bipartites, function(x, n, comp)
#     unpackBipartite(x, n, comp),
#     n = observations, comp=data$compositionChange)
#   names(f$bipartites) <- names(bipartites)
#   f$behavs <-  lapply(Behaviors, function(x, n) unpackBehavior(x, n),
#                       n = observations)
#   names(f$behavs) <- names(Behaviors)
#   f$contbehavs <- lapply(continuousBehaviors, function(x, n)
#     unpackBehavior(x, n), n = observations)
#   names(f$contbehavs) <- names(continuousBehaviors)
#   f$observations <- observations
#   f$seed<- vector("list", observations - 1)
#   f$depvars <- data$depvars
#   f$nodeSets <- data$nodeSets
#   f$oneModes <- oneModes
#   f$Behaviors <- Behaviors
#   f$continuousBehaviors <- continuousBehaviors
#   f$oneModeUpOnly <- sapply(oneModes, function(x) attr(x, "uponly"))
#   f$oneModeDownOnly <- sapply(oneModes, function(x) attr(x, "downonly"))
#   f$behaviorUpOnly <- sapply(Behaviors, function(x) attr(x, "uponly"))
#   f$behaviorDownOnly <- sapply(Behaviors, function(x) attr(x,
#                                                            "downonly"))
#   f$distances <- sapply(data$depvars, function(x) attr(x, "distance"))
#   f$cCovars <- data$cCovars
#   f$vCovars <- data$vCovars
#   ## dyadic covars need to be edgelists
#   f$dycCovars <- lapply(data$dycCovars, function(x) unpackCDyad(x))
#   f$dyvCovars <- lapply(data$dyvCovars, function(x,n) unpackVDyad(x,n),
#                         n=observations)
#   ## create the composition change event lists
#   f$exog <- lapply(data$compositionChange, function(x)
#     unpackCompositionChange(x))
#   f
# }
# 
# unpackBipartite <- function(depvar, observations, compositionChange)
# {
#   edgeLists <- vector("list", observations)
#   networks <- vector("list", observations)
#   actorSet <- attr(depvar, "nodeSet")
#   compActorSets <- sapply(compositionChange, function(x)attr(x, "nodeSet"))
#   thisComp <- match(actorSet, compActorSets)
#   compChange <- any(!is.na(thisComp))
#   if (compChange)
#   {
#     #  stop("Composition change is not yet implemented for bipartite",
#     #       "networks")
#     action <- attr(compositionChange[[thisComp]], "action")
#     ccOption <- attr(compositionChange[[thisComp]], "ccOption")
#   }
#   else
#   {
#     ccOption <- 0
#     action <- matrix(0, nrow=attr(depvar, "netdims")[1], ncol=observations)
#   }
#   sparse <- attr(depvar, "sparse")
#   allowOnly <- attr(depvar, "allowOnly")
#   if (sparse)
#   {
#     ## require(Matrix)
#     ## have a list of sparse matrices in triplet format
#     ## with missings and structurals embedded and 0 based indices!
#     netmiss <- vector("list", observations)
#     for (i in 1:observations)
#     {
#       ## extract this matrix
#       networks[[i]] <- depvar[[i]]
#       nActors <- nrow(depvar[[i]])
#       nReceivers <- ncol(depvar[[i]])
#       ## stop if any duplicates
#       netmat <- cbind(networks[[i]]@i+1, networks[[i]]@j+1,
#                       networks[[i]]@x)
#       if (any(duplicated(netmat[, 1:2])))
#       {
#         stop("duplicate entries in sparse matrix")
#       }
#       ## extract missing entries
#       netmiss[[i]] <- netmat[is.na(netmat[, 3]), , drop = FALSE]
#       ## carry forward missing values if any
#       if (i == 1) # set missings to zero
#       {
#         netmat <- netmat[!is.na(netmat[,3]), ]
#         networks[[i]] <- spMatrix(nActors, nReceivers, netmat[, 1],
#                                   netmat[, 2], netmat[,3])
#       }
#       else
#       {
#         netmiss1 <- netmiss[[i]][, 1:2]
#         storage.mode(netmiss1) <- "integer"
#         networks[[i]][netmiss1[, 1:2]] <-
#           networks[[i-1]][netmiss1[, 1:2]]
#       }
#     }
#     for (i in 1:observations)
#     {
#       mat1 <- networks[[i]]
#       mat1 <- cbind(mat1@i + 1, mat1@j + 1, mat1@x)
#       ##missing edgelist
#       mat2 <- netmiss[[i]]
#       mat2[, 3] <- 1
#       ## rows of mat1 with structural values
#       struct <- mat1[, 3] %in% c(10, 11)
#       ## reset real data
#       mat1[struct, 3] <- mat1[struct, 3] - 10
#       ## copy reset data to structural edgelist
#       mat3 <- mat1[struct, , drop = FALSE]
#       ## now remove the zeros from reset data
#       mat1 <- mat1[!mat1[, 3] == 0, ]
#       ## do comp change
#       if (compChange)
#       {
#         ## revert to sparse matrices temporarily
#         mat1 <- spMatrix(nrow=nActors, ncol=nReceivers, i = mat1[, 1],
#                          j=mat1[, 2], x=mat1[, 3])
#         mat2 <- spMatrix(nrow=nActors, ncol=nReceivers, i = mat2[, 1],
#                          j=mat2[, 2], x=mat2[, 3])
#         mat3 <- spMatrix(nrow=nActors, ncol=nReceivers, i = mat3[, 1],
#                          j=mat3[, 2], x=mat3[, 3])
#         ones <- which(action[, i] == 1)
#         twos <- which(action[, i] == 2)
#         threes <- which(action[, i] == 3)
#         for (j in ones) ## False data is not preceded by anything real
#         {
#           if (ccOption %in% c(1, 2))
#           {
#             ## find missing values for this actor
#             use <- mat2[j, ] > 0
#             ## remove from real data (i.e. zero)
#             mat1[j, use] <- 0
#             ## remove from missing data
#             mat2[j, use] <- 0
#             ## remove from raw data for distances later
#             depvar[[i]][j, use] <- 0 ## zero
#           }
#           else if (ccOption == 3)
#           {
#             ## add the row  to the missing data
#             mat2[j, ] <- 1
#             ## set to missing in raw data for distances later
#             depvar[[i]][j, ] <- NA
#           }
#         }
#         for (j in threes) ## False data is preceded and followed by real
#         {
#           if (ccOption %in% c(1, 2))
#           {
#             ## find missing values for this actor
#             use <- mat2[j, ] > 0
#             ## remove these from mat2, the missing data
#             mat2[j, use] <- 0
#             ## carry forward
#             if (i == 1)
#             {
#               ## 0 any matches from mat1, the real data
#               mat1[j, use] <- 0
#             }
#             else
#             {
#               mat1[j, use] <- networks[[i-1]][j, use]
#             }
#             depvar[[i]][j, use] <- 0 ##  not missing
#           }
#           else if (ccOption == 3)
#           {
#             ## add the row to the missing data
#             mat2[j, ] <- 1
#             depvar[[i]][j, ] <- NA
#           }
#         }
#         for (j in twos) ## False data is not followed by anything real
#         {
#           if (ccOption == 1)
#           {
#             ## find missing values for this actor
#             use <- mat2[j, ] > 0
#             ## remove these from mat2, the missing data
#             mat2[j, use] <- 0
#             depvar[[i]][j, use] <- 0 ##  not missing
#             ## carry forward
#             if (i == 1)
#             {
#               ## 0 any matches from mat1, the real data
#               mat1[j, use] <- 0
#             }
#             else
#             {
#               mat1[j, use] <- networks[[i-1]][j , use]
#             }
#           }
#           else if (ccOption %in% c(2, 3))
#           {
#             ## add the row  to the missing data
#             mat2[j, ] <- 1
#             depvar[[i]][j, ] <- NA
#           }
#         }
#         
#         ## now revert to triplet matrices, after updating networks
#         networks[[i]] <- mat1
#         mat1 <- cbind(mat1@i + 1, mat1@j + 1, mat1@x)
#         mat2 <- cbind(mat2@i + 1, mat2@j + 1, mat2@x)
#         mat3 <- cbind(mat3@i + 1, mat3@j + 1, mat3@x)
#         if (any (mat1[, 3] == 0) || any (mat2[, 3] == 0) ||
#             any (mat3[, 3] == 0))
#         {
#           stop("zero values in sparse matrices")
#         }
#         if (any (duplicated(mat1[, -3])) ||
#             any (duplicated(mat2[, -3])) ||
#             any (duplicated(mat3[, -3])))
#         {
#           stop("duplicate values in sparse matrices")
#         }
#       }
#       ##fix up storage mode to be integer
#       storage.mode(mat1) <- "integer"
#       storage.mode(mat2) <- "integer"
#       storage.mode(mat3) <- "integer"
#       ## add attribute of size
#       attr(mat1,"nActors") <- c(nActors, nReceivers)
#       attr(mat2,"nActors") <- c(nActors, nReceivers)
#       attr(mat3,"nActors") <- c(nActors, nReceivers)
#       if (i < observations)
#       {
#         ## recreate the distance etc
#         mymat1 <- depvar[[i]]
#         mymat2 <- depvar[[i + 1]]
#         ##remove structural values
#         x1 <- mymat1@x
#         x2 <- mymat2@x
#         x1[x1 %in% c(10, 11)] <- NA
#         x2[x2 %in% c(10, 11)] <- NA
#         mymat1@x <- x1
#         mymat2@x <- x2
#         mydiff <- mymat2 - mymat1
#         attr(depvar, "distance")[i] <- sum(mydiff != 0,
#                                            na.rm = TRUE)
#         if (allowOnly)
#         {
#           if (all(mydiff@x >= 0, na.rm=TRUE))
#           {
#             attr(depvar, "uponly")[i] <- TRUE
#           }
#           if (all(mydiff@x <= 0, na.rm=TRUE))
#           {
#             attr(depvar, "downonly")[i] <- TRUE
#           }
#         }
#       }
#       edgeLists[[i]] <- list(mat1 = t(mat1), mat2 = t(mat2),
#                              mat3 = t(mat3))
#     }
#   }
#   else
#   {
#     for (i in 1:observations) ## carry missings forward  if exist
#     {
#       networks[[i]] <- depvar[, , i]
#       if (i == 1)
#         networks[[i]][is.na(depvar[, , i])] <-0
#       else ##carry missing forward!
#         networks[[i]][is.na(depvar[, , i])] <-
#           networks[[i-1]][is.na(depvar[, , i])]
#     }
#     for (i in 1:observations)
#     {
#       ones <- which(action[, i] == 1)
#       twos <- which(action[, i] == 2)
#       threes <- which(action[, i] == 3)
#       for (j in ones) ## False data is not preceded by anything real
#       {
#         if (ccOption %in% c(1, 2))
#         {
#           use <- is.na(depvar[j, , i])
#           depvar[j, use, i] <- 0 ## not missing
#           networks[[i]][j, use] <- 0 ## zero
#         }
#         else if (ccOption == 3)
#         {
#           depvar[j, , i] <- NA ## missing
#         }
#       }
#       for (j in threes) ## False data is preceded and followed by real
#       {
#         
#         if (ccOption %in% c(1, 2))
#         {
#           use <- is.na(depvar[j, , i])
#           depvar[j, use, i] <- 0 ##  not missing
#           ## carry forward already done
#           if (i == 1)
#           {
#             networks[[i]][j, use] <- 0
#           }
#           else
#           {
#             networks[[i]][j, use] <- networks[[i-1]][j, use]
#           }
#         }
#         else if (ccOption == 3)
#         {
#           depvar[j, , i] <- NA ## missing
#         }
#       }
#       for (j in twos) ## False data is not followed by anything real
#       {
#         if (ccOption == 1)
#         {
#           use <- is.na(depvar[j, , i])
#           depvar[j, use, i] <- 0 ##  not missing
#           ## carry forward already done
#           if (i == 1)
#           {
#             networks[[i]][j, use] <- 0
#           }
#           else
#           {
#             networks[[i]][j, use] <- networks[[i-1]][j, use]
#             
#           }
#         }
#         else if (ccOption %in% c(2, 3))
#         {
#           depvar[j, , i] <- NA ## missing
#         }
#       }
#     }
#     for (i in 1:observations)
#     {
#       if (i < observations)
#       {
#         ## recreate distances, as we have none in c++. (no longer true)
#         mymat1 <- depvar[,,i, drop=FALSE]
#         mymat2 <- depvar[,,i + 1,drop=FALSE]
#         ##remove structural values
#         mymat1[mymat1 %in% c(10,11)] <- NA
#         mymat2[mymat2 %in% c(10,11)] <- NA
#         mydiff <- mymat2 - mymat1
#         attr(depvar, "distance")[i] <- sum(mydiff != 0,
#                                            na.rm = TRUE)
#         if (allowOnly)
#         {
#           if (all(mydiff >= 0, na.rm=TRUE))
#           {
#             attr(depvar, "uponly")[i] <- TRUE
#           }
#           if (all(mydiff <= 0, na.rm=TRUE))
#           {
#             attr(depvar, "downonly")[i] <- TRUE
#           }
#         }
#       }
#       
#       edgeLists[[i]] <- createEdgeLists(networks[[i]], depvar[, , i], TRUE)
#     }
#   }
#   ## add attribute of nodeset
#   attr(edgeLists, "nodeSet") <- attr(depvar, "nodeSet")
#   ## add attribute of name
#   attr(edgeLists, "name") <- attr(depvar, "name")
#   ## add attribute of distance
#   attr(edgeLists, "distance") <- attr(depvar, "distance")
#   ## attr uponly and downonly
#   attr(edgeLists, "uponly") <- attr(depvar, "uponly")
#   attr(edgeLists, "downonly") <- attr(depvar, "downonly")
#   ## attr symmetric
#   attr(edgeLists, "symmetric") <- attr(depvar, "symmetric")
#   ## attr balmean
#   attr(edgeLists, "balmean") <- attr(depvar, "balmean")
#   ## attr structmean
#   attr(edgeLists, "structmean") <- attr(depvar, "structmean")
#   attr(edgeLists, "averageOutDegree") <- attr(depvar, "averageOutDegree")
#   return(edgeLists = edgeLists)
# }
# 
# createEdgeLists<- function(mat, matorig, bipartite)
# {
#   ## mat1 is basic values, with missings and structurals replaced
#   tmp <- lapply(1 : nrow(mat), function(x, y)
#   {
#     mymat <- matrix(0, nrow = sum(y[x, ] > 0), ncol = 3)
#     mymat[, 1] <- x
#     mymat[, 2] <- which(y[x, ] != 0)
#     mymat[, 3] <- y[x, mymat[, 2]]
#     mymat
#   }, y = mat)
#   mat1 <- do.call(rbind, tmp)
#   ## mat2 reverts to matorig to get the missing values
#   tmp <- lapply(1 : nrow(matorig), function(x, y)
#   {
#     mymat <- matrix(0, nrow = sum(is.na(y[x, ])), ncol = 3)
#     mymat[, 1] <- x
#     mymat[, 2] <- which(is.na(y[x, ]))
#     mymat[, 3] <- 1
#     mymat
#   }, y = matorig)
#   mat2 <- do.call(rbind, tmp)
#   ## remove the diagonal if not bipartite
#   if (!bipartite)
#   {
#     mat2 <- mat2[mat2[, 1] != mat2[, 2], , drop=FALSE]
#   }
#   ## mat3 structurals
#   struct <- mat1[,3] %in% c(10, 11)
#   mat1[struct, 3] <- mat1[struct,3] - 10
#   mat3 <- mat1[struct, , drop=FALSE]
#   mat3[, 3] <- 1
#   mat1 <- mat1[!mat1[,3] == 0, , drop=FALSE] ##remove any zeros just created
#   ##fix up storage mode to be integer
#   storage.mode(mat1) <- "integer"
#   storage.mode(mat2) <- "integer"
#   storage.mode(mat3) <- "integer"
#   ## add attribute of size
#   if (bipartite)
#   {
#     attr(mat1, "nActors") <- c(nrow(mat), ncol(mat))
#     attr(mat2, "nActors") <- c(nrow(mat), ncol(mat))
#     attr(mat3, "nActors") <- c(nrow(mat), ncol(mat))
#   }
#   else
#   {
#     attr(mat1, "nActors") <- nrow(mat)
#     attr(mat2, "nActors") <- nrow(mat)
#     attr(mat3, "nActors") <- nrow(mat)
#   }
#   
#   list(mat1 = t(mat1), mat2 = t(mat2), mat3 = t(mat3))
# }
# 
# ##@createCovarEdgeLists siena07 Reformat data for C++
# createCovarEdgeList<- function(mat, matorig)
# {
#   tmp <- lapply(1 : nrow(mat), function(x, y)
#   {
#     mymat <- matrix(0, nrow = sum(y[x, ] != 0), ncol = 3)
#     mymat[, 1] <- x
#     mymat[, 2] <- which(y[x, ] != 0)
#     mymat[, 3] <- y[x, mymat[, 2]]
#     mymat
#   }, y = mat)
#   mat1 <- do.call(rbind, tmp)
#   ##mat2 reverts to matorig to get the missing values
#   tmp <- lapply(1 : nrow(matorig), function(x, y)
#   {
#     mymat <- matrix(0, nrow = sum(is.na(y[x, ])), ncol = 3)
#     mymat[, 1] <- x
#     mymat[, 2] <- which(is.na(y[x, ]))
#     mymat[, 3] <- 1
#     mymat
#   }, y = matorig)
#   mat2 <- do.call(rbind, tmp)
#   ## add attribute of size
#   attr(mat1, "nActors1") <- nrow(mat)
#   attr(mat1, "nActors2") <- ncol(mat)
#   list(mat1=t(mat1), mat2=t(mat2))
# }



#
# getChangeContributions(m1$rsiena_algorithm, m1$rsiena_data, m1$rsiena_effects)

































# ##**TODO** Action count data frames (explore, exploit) 
# ##          rows (actors) by columns (iterations), counts: 0,1,2,3,...
# 
# rolling_window = 20
# 
# self <- m1
# 
# m1$chain_stats %>% group_by(iteration_id) %>% count() %>% plot(type='l')
# 
# m1$chain_stats %>% group_by(X__dv_name, X__id_from) %>% count() %>% print(n=60)
# 
# util_after <- m1$chain_stats$X__utility_after
# plot(util_after, type='l')
# plot(cummean(util_after), type='l')
# 
# 
# ## Get Actor-period count of ministep decisions
# step_summary_df <- m1$chain_stats %>%
#   filter( X__dv_name=='self$bipartite_rsienaDV') %>% 
#   group_by(X__id_from)
# iters_with_change <- unique(step_summary_df$iteration_id)
# cnt_mat <- matrix(0, nrow=self$M, ncol=length(self$rsiena_model$chain) )
# for (iter in iters_with_change) {
#   cat(sprintf(' %s ', iter))
#   iter_df <- step_summary_df %>% filter(iteration_id == iter)
#   actor_cnts <- iter_df %>% group_by(X__id_from) %>% count()
#   id_not_0s <- ( ! actor_cnts$X__id_from %in% c(0,'0') )
#   if ( any( ! id_not_0s ) ) {
#     print('zero in id_from ?:')
#     print(id_not_0s)
#     actor_cnts <- actor_cnts %>% filter( ! X__id_from %in% c(0,'0') )
#   }
#   cnt_mat[ actor_cnts$X__id_from, iter ] <- actor_cnts$n
# }
# # View(.)
# 
# cnt_df <- as.data.frame(cnt_mat)
# colnames(cnt_df) <- paste0("Period_", 1:length(self$rsiena_model$chain) )
# cnt_df$Actor <- 1:self$M
# rownames(cnt_df) <- paste0("Actor_", cnt_df$Actor)
# 
# long_data <- cnt_df %>% 
#     pivot_longer(cols = starts_with("Period_"),
#                  names_to = "period",
#                  values_to = "count") %>% 
#   mutate(
#     period = as.numeric(gsub("Period_", "", period)),
#     actor = as.numeric(gsub("Actor_", "", Actor))
#   ) %>%
#   select(actor, period, count)  # Reorder columns
# 
# ## plot actor event sequences
# ggplot(long_data, aes(y = factor(actor), x=period, fill = count)) + 
#   geom_tile() + 
#   scale_fill_gradient(low = "white", high = "blue", name = "Ministeps Count") +
#   theme_bw() + labs('y' = 'Actor')
# 
# 
# # ggplot(long_data, aes(x=period, y=count, color=factor(actor))) + geom_point() + geom_smooth()
# 
# 
# 
# 
# ## Get timeseries of actor utility function value
# util_step_df <- m1$chain_stats %>% filter( X__dv_name=='self$bipartite_rsienaDV')
# n_util_steps <- nrow(util_step_df)
# util_mat <- matrix(NA, nrow=self$M, ncol=length(self$rsiena_model$chain) )
# for (i in 1:n_util_steps) {
#   cat(sprintf(' %s ', i))
#   x <- util_step_df[i, ]
#   util_mat[ x$X__id_from , x$iteration_id ] <- x$X__utility_after
#   
# }
# # Function to fill forward values along rows
# fill_forward_func <- function(row) {
#   for (i in 2:length(row)) {
#     if (is.na(row[i])) {
#       row[i] <- row[i - 1]
#     }
#   }
#   return(row)
# }
# # Apply the fill-forward function to each row
# util_mat_filled <- t(apply(util_mat, 1, fill_forward_func))
# 
# util_df_fill <- as.data.frame(util_mat_filled)
# colnames(util_df_fill) <- paste0("Period_", 1:length(self$rsiena_model$chain) )
# util_df_fill$Actor <- 1:self$M
# rownames(util_df_fill) <- paste0("Actor_", util_df_fill$Actor)
# 
# long_util <- util_df_fill %>% 
#   pivot_longer(cols = starts_with("Period_"),
#                names_to = "period",
#                values_to = "count") %>% 
#   mutate(
#     period = as.numeric(gsub("Period_", "", period)),
#     actor = as.numeric(gsub("Actor_", "", Actor))
#   ) %>%
#   select(actor, period, count)  # Reorder columns
# 
# ggplot(long_util, aes(x=period, y=count, color=factor(actor))) + 
#   geom_point() + geom_smooth() + theme_bw()
# 
# ggplot(long_util, aes(x=period, y=count, color=factor(actor))) +  # geom_point(shape=1, alpha=.5) + 
#   geom_line(linewidth=.8) + # geom_smooth() + 
#   theme_bw() + labs(y='Utility', legend='Actor')
# 
# # compute utility moving average
# long_util_ma <- long_util %>%
#   group_by(actor) %>%
#   mutate(rolling_avg = zoo::rollmean(count, k = rolling_window, fill = NA, align = "right"))
# #
# util_avg <- long_util_ma %>% group_by(period) %>% summarize(mean=mean(rolling_avg, na.rm=T))
# # plot utility moving average
# ggplot(long_util_ma, aes(x=period, y=rolling_avg, color=factor(actor))) +  
#   geom_point(shape=1, alpha=.1) + 
#   geom_line(linewidth=.9, alpha=.85) + # geom_smooth() + 
#   geom_line(data = util_avg, linewidth=2, color='black', aes(x=period, y=mean, color='Mean')) + 
#   theme_bw() + labs(y='Utility', color='Actor')
# 
# 
# ##
# plot(colMeans(util_mat_filled, na.rm = T), type='l')
# 
# 
# # %>%  mutate(Period = as.numeric(gsub("Period_", "", Period)))  
# 
# ## The chain has the structure chain[[run]][[depvar]][[period]][[ministep]].
# ##ch[[5]] = '[[depvar]][[period]][[ministep]]'
# 
# ch <- m1$rsiena_model$chain
# id_ministep <- 5
# run <- ch[[id_ministep]]
# depvar <- run[[1]]
# period <- depvar[[1]]
# ministep <- period[[1]]
# 
# 
# 
# 
# 
# 
# 
# 
# # 1. INIT SIM
# m2 <- SaomNkRSienaBiEnv$new(environ_params, rand_seed=54321)
# # 2. RUN SIM
# m2$search_rsiena_run(structure_model, payoff_formulas, iterations=1000, 
#                            returnDeps = F,  n_snapshots =1, 
#                            rsiena_phase2_nsub=1, rsiena_n2start_scale = 1, 
#                            get_eff_doc = FALSE)
# # 3. PLOT SIM (EXPLORE RESULTS)
# m2$search_rsiena_plot_stability()
# 
# 
# 
# binet <- network(m1$bipartite_matrix, directed = F)
# ergstats <- ergmMPLE(binet ~ edges + cycle(4), output = 'array')
# bi_edges <- ergpreds[,, 'edges']
# bi_cyc4  <- ergpreds[,, 'cycle4']
# 
# 
# 
# 
# 
#  
# id_wave <- 1
# ### ## m1$rsiena_model$sf2[ <sims> , id_wave, <effects> ]
# 
# m1$rsiena_model$sf2
# 
# mean <-  m1$rsiena_model$sf2 
# rownames(wavemat) <- as.character( 1:nrow(wavemat) )
# colnames(wavemat) <- gsub('\\s', '_' , m1$rsiena_model$effects$effectName )
# # dat[,,3] <- dat[1,1,][ m1$rsiena_model$effects$effectName ]
# 
# df_long <- melt(wavemat, varnames = c("Simulation", "Wave", "NetworkEffect"), value.name = "Mean")
# 
# # Plot using ggplot2
# ggplot(df_long, aes(x = Mean, color=NetworkEffect, fill=NetworkEffect)) +
#   # geom_boxplot(aes(group = Wave), outlier.shape = NA) +
#   geom_density(alpha=.1, bins=11) +
#   facet_wrap(~ NetworkEffect, scales = "free_x") +
#   theme_minimal() +
#   labs(title = "Distribution of Simulations by Effects",
#        x = "X",
#        y = "Y")
# 
# 
# 
# ##$sf2  -> [sims,  waves,  effects]
# ##$sf   -> [sims, effects]
# 
# apply(m1$rsiena_model$sf2[,1,], 2, mean)
# 
# colMeans(m1$rsiena_model$sf)
# 
# 
# 
# ##------------------------------------
# 
# # m1 <-
# 
# # We can also format the predictor matrix into an array:
# mplearray <- ergmMPLE(formula, output="array")
# 
# #
# dvnet <- 
# formula <- as.formula(sprintf('%s ~ edges',dvnet))
# 
# # The resulting matrices are big, so only print the first 8 actors:
# mplearray$response[1:8,1:8]
# mplearray$predictor[1:8,1:8, ]
# mplearray$weights[1:8,1:8]
# 
# ##------------------------------------
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###--------------------------------------------------------------------------
# 
# m1$search_rsiena_run(objective_list, iterations=2000, returnDeps = T, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
# m1$search_rsiena_plot_stability()
# 
# 
# # ##
# # m1$search_rsiena_run(iterations=4000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# # m1$search_rsiena_plot_stability()
# # m1$search_rsiena_run(iterations=400, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# # m1$search_rsiena_plot_stability()
# 
# 
# 
# 
# m1 <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .11, sim_name = '_TESTrsiena_')
# m1$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# m1$search_rsiena_plot_stability()
# #
# m2 <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .11, sim_name = '_TESTrsiena_')
# m2$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
# m2$search_rsiena_plot_stability()
# 
# #
# ms1 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .11, sim_name = '_TESTrsiena_')
# ms1$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# ms1$search_rsiena_plot_stability()
# #
# ms2 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .11, sim_name = '_TESTrsiena_')
# ms2$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
# ms2$search_rsiena_plot_stability()
# 
# 
# 
# md1 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .2, sim_name = '_TESTrsiena_')
# md1$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# md1$search_rsiena_plot_stability()
# #
# md2 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .2, sim_name = '_TESTrsiena_')
# md2$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
# md2$search_rsiena_plot_stability()
# 
# # #
# # mc1 <- SAOM_NK_RSiena$new(M = 12, N = 30, BI_PROB = .05, sim_name = '_TESTrsiena_')
# # mc1$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# # mc1$search_rsiena_plot_stability()
# # #
# # mc2 <- SAOM_NK_RSiena$new(M = 12, N = 30, BI_PROB = .05, sim_name = '_TESTrsiena_')
# # mc2$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
# # mc2$search_rsiena_plot_stability()
# 
# # #
# # ms2 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .15, sim_name = '_TESTrsiena_')
# # ms2$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = .1)
# # ms2$search_rsiena_plot_stability()
# 
# 
# #
# # m3 <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .15, sim_name = '_TESTrsiena_')
# # m3$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=3, rsiena_n2start_scale = 1)
# # m3$search_rsiena_plot_stability()
# #
# # m4 <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .15, sim_name = '_TESTrsiena_')
# # m4$search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=4, rsiena_n2start_scale = 1)
# # m4$search_rsiena_plot_stability()
# 
# 
# saomnkrsiena <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .15, sim_name = '_TESTrsiena_')
# ##
# # saomnkrsiena$search_rsiena_run(iterations=10, rsiena_phase2_nsub=1) 
# # saomnkrsiena$search_rsiena_run(iterations=100, rsiena_phase2_nsub=2)
# # saomnkrsiena$search_rsiena_run(iterations=5000, rsiena_phase2_nsub=3)
# saomnkrsiena$search_rsiena_run(iterations=1000, n_snapshots =2, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# saomnkrsiena$search_rsiena_plot_stability()
# # saomnkrsiena$search_rsiena_run(iterations=1000, rsiena_phase2_nsub=5)
# 
# # saom_nk_enhanced$search_rsiena(1)
# saomnkrsiena$search_rsiena_run(iterations=100) 
# 
# saomnkrsiena$search_rsiena_run(iterations=10000)
# saomnkrsiena$search_rsiena_run(iterations=100000)
# # saomnkrsiena$search_rsiena_run(iterations=20, overwrite = F) 
# 
# saomnkrsiena2 <- SAOM_NK_RSiena$new(M = 15, N = 24, BI_PROB = .15, sim_name = '_TESTrsiena_')
# # saom_nk_enhanced$search_rsiena(1)
# saomnkrsiena2$search_rsiena_run(iterations=1000) 
# 
# 
# 
# 
# ########################################
# self <- saomnkrsiena$clone(deep=T)
# sims = self$rsiena_model$sims
# n <- length(sims)
# outlist <- list()
# difflist <- list()
# jaccardlist <- list()
# K_env_list <- list()
# K_soc_list <- list()
# for(i in 1:length(sims)) {
#   cat(sprintf(' %s ', i))
#   ##
#   el_bi_env <- sims[[ i ]][[1]][[3]]$`1`
#   ## update numbering of second mode (the comonent integer names shift upward by the number of actors)
#   el_bi_env[,2] <- el_bi_env[,2] + self$M
#   ##
#   MplusN <- self$M + self$N
#   #############
#   ## Bipartite matrix space (N+M by N+M)
#   ## Undirected --> Upper right rectangle of full bipartite matrix
#   ##   M N
#   ## M[0,X], for X in [0,1]
#   ## N[0,0]
#   bi_env_mat_sp <- sparseMatrix(i = el_bi_env[,1],
#                                 j = el_bi_env[,2],
#                                 x = el_bi_env[,3],
#                                 dims = c(MplusN, MplusN))
#   ## get the top-right triangle for the bipartite environment
#   bi_env_mat_new  <- as.matrix(bi_env_mat_sp)[ 1:self$M, (self$M+1):(MplusN) ]
#   ## Add ACTOR NAMES on rows
#   rownames(bi_env_mat_new) <- as.character( 1:self$M )
#   ## Add COMPONENT NAMES on columns (N+1 ... N+M)
#   colnames(bi_env_mat_new) <- as.character( (1:self$N) + self$M  )
#   ##
#   # self$plot_bipartite_system_from_mat(bi_env_mat_new, i, plot_save = TRUE)
#   ##
#   outlist[[ sprintf('sim%d',i) ]] <- bi_env_mat_new
#   
#   if (i == 1) {
#     difflist[[ sprintf('diff%d-%d', i-1, i) ]] <- 0
#     jaccardlist[[sprintf('jac%d-%d', i-1, i)]] <- 0
#   } else {
#     difflist[[ sprintf('diff%d-%d', i-1, i) ]] <-  bi_env_mat_new - outlist[[ (i-1) ]]
#     jaccardlist[[sprintf('jac%d-%d', i-1, i)]] <- get_jaccard_index(m0 = outlist[[ (i-1) ]], m1 = bi_env_mat_new )
#   }
#   
#   new_bi_g <- igraph::graph_from_biadjacency_matrix(bi_env_mat_new, 
#                                                       directed = F, mode = 'all', 
#                                                       multiple = T, weighted = T, 
#                                                       add.names = T)
#   projections <- igraph::bipartite_projection(new_bi_g, multiplicity = T, which = 'both')
#   new_g_soc <- projections$proj1
#   new_g_env <- projections$proj2
# 
#   K_soc_list[[i]] <- igraph::degree(new_g_soc)
#   K_env_list[[i]] <- igraph::degree(new_g_env)
#   
#   
#   # print(sim_i[[1]][[3]]`1`)
# }
# 
# saomnkrsiena$plot_degree_progress_from_sims(K_soc_list = K_soc_list, K_env_list = K_env_list, plot_save = T)
# 
# 
# # par(mfrow=c(1,3))
# # stability_vec <- cumsum(plyr::ldply(jaccardlist)[,2])/(1:n)
# # plot(stability_vec , type='l' , main='Change in Stability  [t-1, t]')
# # 
# # 
# degree_summary <- do.call(rbind, lapply(1:length(K_soc_list), function(iter) {
#   data.frame(
#     Iteration = iter,
#     Mean_K_S = mean(K_soc_list[[iter]]),
#     Q25_K_S = quantile(K_soc_list[[iter]], 0.25),
#     Q75_K_S = quantile(K_soc_list[[iter]], 0.75),
#     Mean_K_E = mean(K_env_list[[iter]]),
#     Q25_K_E = quantile(K_env_list[[iter]], 0.25),
#     Q75_K_E = quantile(K_env_list[[iter]], 0.75)
#   )
# }))
# (plt <- ggplot(degree_summary, aes(x = Iteration)) +
#     geom_line(aes(y = Mean_K_S, color = "Mean K_S")) +
#     geom_ribbon(aes(ymin = Q25_K_S, ymax = Q75_K_S, fill = "K_S"), alpha = 0.1) +
#     geom_line(aes(y = Mean_K_E, color = "Mean K_E")) +
#     geom_ribbon(aes(ymin = Q25_K_E, ymax = Q75_K_E, fill = "K_E"), alpha = 0.1) +
#     scale_color_manual(values = c("Mean K_S" = "blue", "Mean K_E" = "red")) +
#     scale_fill_manual(values = c("K_S" = "blue", "K_E" = "red")) +
#     labs(title = "Degree Progress Over Iterations",
#          x = "Iteration",
#          y = "Degree",
#          color = "Mean Degree",
#          fill = "IQR (Mid-50%)") +
#     theme_minimal()
# )
# 
# 
# 
# ########################################
# 
# # sims = self$rsiena_model$sims
# # n <- length(sims)
# # outlist <- list()
# # difflist <- list()
# # jaccardlist <- list()
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
# #   ##
# #   bi_env_mat_new  <- as.matrix(bi_env_mat_sp)[ 1:self$M, (self$M+1):(MplusN) ]
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
# #   # print(sim_i[[1]][[3]]`1`)
# # }
# # 
# # par(mfrow=c(1,3))
# # stability_vec <- cumsum(plyr::ldply(jaccardlist)[,2])/(1:n)
# # plot(stability_vec , type='l' , main='Change in Stability  [t-1, t]')
# 
# # stability_delta <- stability_vec / (1:n)
# # plot(stability_delta, type='o', log='y',
# #      xlab='t', 
# #      # ylim=c( stability_delta[n-1],  1 ),
# #      ylab='Ln Stability Change [t-1, t]', main='Stabilization Rate (Ln Change in Stability)' 
# # ); abline(h = tol, col='pink', lty=2)
# # 
# # burn_prop <- 0.2
# # iter_postburn <- round(c( 1-burn_prop, burn_prop) * n )
# # hist(stability_vec[ iter_postburn[1]:iter_postburn[2] ], main='Stability (post-burn)')
# 
# 
# # ########################################################################
# # sims = saomnkrsiena$rsiena_model$sims
# # outlist <- list()
# # difflist <- list()
# # jaccardlist <- list()
# # for(i in 1:length(sims)) {
# #   cat(sprintf(' %s ', i))
# #   ##
# #   el_bi_env <- sims[[ i ]][[1]][[3]]$`1`
# #   ## update numbering of second mode (the comonent integer names shift upward by the number of actors)
# #   el_bi_env[,2] <- el_bi_env[,2] + saomnkrsiena$M
# #   ##
# #   MplusN <- saomnkrsiena$M + saomnkrsiena$N
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
# #   ##
# #   bi_env_mat_new  <- as.matrix(bi_env_mat_sp)[ 1:saomnkrsiena$M, (saomnkrsiena$M+1):(MplusN) ]
# #   ##
# #   # saomnkrsiena$plot_bipartite_system_from_mat(bi_env_mat_new, i, plot_save = TRUE) 
# #   ##
# #   outlist[[ sprintf('sim%d',i) ]] <- bi_env_mat_new
# #   
# #   if (i > 1) {
# #     difflist[[ sprintf('diff%d-%d', i-1, i) ]] <-  bi_env_mat_new - outlist[[ (i-1) ]] 
# #     
# #     jaccardlist[[sprintf('jac%d-%d', i-1, i)]] <- get_jaccard_index(bi_env_mat_new, outlist[[ (i-1) ]] )
# #   }
# #   
# #   # print(sim_i[[1]][[3]]`1`)
# # }
# # 
# # par(mfrow=c(1,3))
# # stability_vec <- cumsum(plyr::ldply(jaccardlist)[,2])/(1:length(jaccardlist))
# # plot(stability_vec , type='l' )
# # 
# # stability_delta <- stability_vec / 1:length(stability_vec)
# # plot(stability_delta, type='o', log='y',
# #      xlab='t', ylab='Ln Stability Change [t-1, t]', main='Change in Stability' )
# # 
# # hist(stability_delta[100:length(stability_delta)])
# # 
# # ########################################################################
# 
# 
# saom_nk_enhanced$local_search_saom_batchrun(batches = 2, iterations = 3, plot_save = T)
# saom_nk_enhanced$plot_fitness_progress()
# saom_nk_enhanced$plot_degree_progress()
# # saom_nk_enhanced$visualize_networks()
# 
# 
# saom_nk_enhanced$local_search_saom_batchrun(batches = 4, iterations = 10, 
#                                             plot_save=TRUE, filename='__TEST_PLOT2__') 
# # saom_nk_enhanced$visualize_networks(plot_save = T)
# 
# 
# 
# #### SANDBOX ########
# 
# saomnkrsiena$rsiena_model$sims[[1]][[1]][[2]]
# 
# M <- saomnkrsiena$M
# N <- saomnkrsiena$N
# actors <- 1:saomnkrsiena$M 
# components <- 1:saomnkrsiena$N
# el_bi_env <- saomnkrsiena$rsiena_model$sims[[1]][[1]][[3]]$`1`
# el_bi_env[,2] <- el_bi_env[,2] + M
# # bi_env_el <- el_bi_env[,1:2]
# ##
# el_proj1  <- saomnkrsiena$rsiena_model$sims[[1]][[1]][[1]]$`1`
# el_proj2  <- saomnkrsiena$rsiena_model$sims[[1]][[1]][[2]]$`1`
# el_proj2 <- el_proj2 + M
# 
# ## Undirected --> Upper right rectangle of full bipartite matrix
# ##   M N
# ## M[0,X] for X=[0,1]
# ## N[0,0]
# bi_env_mat_sp <- sparseMatrix(i = el_bi_env[,1],
#                               j = el_bi_env[,2],
#                               x = el_bi_env[,3],
#                               dims = c(M+N, M+N))
# bi_env_mat <- as.matrix(bi_env_mat_sp)[1:M,(M+1):(M+N)]
# 
# # self$set_system_from_bipartite_igraph( self$random_bipartite_igraph() )
# 
# 
# #####################
# 
# 
# 
# #########################################
# 
# ## Simple environment 
# env_simple <- SAOM_NK_Enhanced$new(M = 12, N = 8, BI_PROB = .15, sim_name = '_ENV_SIMPLE_')
# env_simple$local_search_saom_batchrun(batches = 10, iterations = 5, plot_save=TRUE) 
# 
# ## Mid environment 
# env_mid <- SAOM_NK_Enhanced$new(M = 12, N = 12, BI_PROB = .15, sim_name = '_ENV_MID_')
# env_mid$local_search_saom_batchrun(batches = 10, iterations = 5, plot_save=TRUE)  
# 
# ## Complex environment 
# env_complex <- SAOM_NK_Enhanced$new(M = 12, N = 18, BI_PROB = .15, sim_name = '_ENV_COMPLEX_')
# env_complex$local_search_saom_batchrun(batches = 10, iterations = 5, plot_save=TRUE) 
# 
# ## Very Complex environment 
# env_verycomplex <- SAOM_NK_Enhanced$new(M = 12, N = 30, BI_PROB = .15, sim_name = '_ENV_VERYCOMPLEX_')
# env_verycomplex$local_search_saom_batchrun(batches = 15, iterations = 5, plot_save=TRUE) 
# 
# 
# #########################################




# library(RSiena) # or RSienaTest
# 
# ###############################################################################
# ###                                                                         ###
# ###         First main function: SimulateNetworks                           ###
# ###                                                                         ###
# ###############################################################################
# 
# 
# SimulateNetworks <- function(n, M, rate, dens, rec, tt, c3,
#                              Vaego, Vaalt, Vasim, Vbego, Vbalt, Vbsim){
#   # Simulates M consecutive network waves, with n actors,
#   # according to a stochastic actor-oriented model
#   # with parameter values rate for rate,
#   # dens for outdegree, rec for reciprocity,
#   # tt for transitive triplets, c3 for 3-cycles,
#   # an actor covariate Va with values alternating between 0 and 1,
#   # with parameter values Vaego, Vaalt, Vasim
#   # for egoX, altX, and simX with respect to Va,
#   # and an actor covariate Vb with a standard normal distribution,
#   # with parameter values Vbego, Vbalt, Vbsim
#   # for egoX, altX, and simX with respect to Vb.
#   ##
#   # Create actor covariates
#   V0 <- rep(0, n)
#   V0[2*(1:(n %/% 2))] <- 1 # equal binary
#   V1 <- rnorm(n, 0, 1)
#   # Create initial 2-wave data to get a suitable data structure.
#   # arbitrarily, this initial network has an expected average degree of 3
#   X0 <- matrix(rbinom(n*n,1,3/(n-1)),n,n)
#   diag(X0) <- 0
#   X1 <- X0
#   # but X0 and X1 should not be identical for use in sienaDependent
#   X0[1,2] <- 0
#   X0[2,1] <- 1
#   X1[1,2] <- 1
#   X1[2,1] <- 0
#   XX <- array(NA,c(n,n,2))
#   XX[,,1] <- X0
#   XX[,,2] <- X1
#   # With this data structure, we now can create the data.
#   Va <- coCovar(V0)
#   Vb <- coCovar(V1)
#   X   <- sienaDependent(XX, allowOnly = FALSE)
#   InitData <- sienaDataCreate(X, Va, Vb)
#   InitEff0 <- getEffects(InitData)
#   # sink to avoid printing to the screen
#   sink("eff.txt")
#   # Specify the parameters.
#   # The rate parameter is first multiplied by 10,
#   # which will be used only to get from the totally random network XX[,,1] = X0
#   # to the network that will be the simulated first wave.
#   InitEff0 <- setEffect(InitEff0, Rate, type="rate", initialValue = 10*rate)
#   InitEff0 <- setEffect(InitEff0, density, initialValue = dens)
#   InitEff0 <- setEffect(InitEff0, recip, initialValue = rec)
#   InitEff0 <- setEffect(InitEff0, transTrip, initialValue = tt)
#   InitEff0 <- setEffect(InitEff0, cycle3, initialValue = c3)
#   InitEff0 <- setEffect(InitEff0, egoX, interaction1="Va", initialValue = Vaego)
#   InitEff0 <- setEffect(InitEff0, altX, interaction1="Va", initialValue = Vaalt)
#   InitEff0 <- setEffect(InitEff0, simX, interaction1="Va", initialValue = Vasim)
#   InitEff0 <- setEffect(InitEff0, egoX, interaction1="Vb", initialValue = Vbego)
#   InitEff0 <- setEffect(InitEff0, altX, interaction1="Vb", initialValue = Vbalt)
#   InitEff0 <- setEffect(InitEff0, simX, interaction1="Vb", initialValue = Vbsim)
#   # The parameter given for n3 should be larger than sum(InitEff0$include)
#   nthree <- sum(InitEff0$include)	+ 5
#   InitAlg <- sienaAlgorithmCreate(projname="Init", useStdInits=FALSE,
#                                   cond=FALSE, nsub=0, n3=nthree, simOnly=TRUE)
#   # Simulate the first wave.
#   InitSim   <- siena07(InitAlg, data=InitData, eff=InitEff0,
#                        returnDeps=TRUE, batch=TRUE, silent=TRUE)
#   # Now prepare for simulating waves 2 to M.
#   # Create empty result network.
#   Xs <- array(0, dim=c(n,n,M))
#   # The rate parameter value from the function call is reinstated in InitEff.
#   InitEff <- InitEff0
#   InitEff <- setEffect(InitEff, Rate, type="rate", initialValue = rate)
#   sink()
#   for (m in 1:M){
#     # Note that we start this loop with a previously simulated network.
#     # Transform the previously simulated network
#     # from edge list into adjacency matrix
#     XXsim <- matrix(0,n,n)
#     nsim  <- InitAlg$n3
#     XXsim[InitSim$sims[[nsim]][[1]]$X[[1]][,1:2]]  <- InitSim$sims[[nsim]][[1]]$X[[1]][,3]
#     # Put simulated network into the result matrix.
#     Xs[,,m] <- XXsim
#     # Put simulated network in desired places for the next simulation
#     XX[,,2] <- XX[,,1] # used only to get the data structure
#     XX[,,1] <- XXsim
#     if (m < M){
#       # The following is only to prevent the error that would occur
#       # in the very unlikely event XX[,,1] == XX[,,2].
#       if (identical(XX[,,1], XX[,,2])){XX[1,2,1] <- 1 - XX[1,2,2]}
#       # Specify the two-wave network data set starting with XX[,,1].
#       X <- sienaDependent(XX, allowOnly = FALSE)
#       # Simulate wave m+1 starting at XX[,,1] which is the previous XXsim
#       InitData  <- sienaDataCreate(X, Va, Vb)
#       InitSim <- siena07(InitAlg, data=InitData, eff=InitEff,
#                          returnDeps=TRUE, batch=TRUE, silent=TRUE)
#     }
#   }
#   # Present the average degrees to facilitate tuning the outdegree parameter
#   # to achieve a desired average value for the average degrees.
#   cat("Average degrees ", round(colSums(Xs,dims=2)/n, digits=2), "\n")
#   # Result: simulated data set; covara and covarb are vectors of length n;
#   # networks is an array of dimension nxnxM
#   list(covara = V0, covarb = V1, networks = Xs)
# }
# 
# 
# ###############################################################################
# ###                                                                         ###
# ###         Examples                                                        ###
# ###                                                                         ###
# ###############################################################################
# 
# 
# # Trial values:
# n <- 20
# M <- 4
# rate <- 2
# dens <- -1.9
# rec <- 2
# tt <- 0.3
# c3 <- -0.3
# Vaego <- 0
# Vaalt <- 0
# Vasim <- 0.6
# Vbego <- 0.5
# Vbalt <- -0.5
# Vbsim <- 0.5
# 
# # Example call:
# SN <- SimulateNetworks(n, M, rate, dens, rec, tt, c3, Vaego, Vaalt, Vasim,
#                        Vbego, Vbalt, Vbsim)
# # You can repeat this call a few times, and then see the varying values
# # reported for the average degrees.
# # You can also experiment this with other values for dens,
# # keeping everything else the same,
# # varying dens by values of +/- 0.05 to +/- 0.2, for example.
# 
# # For larger n, slightly lower values of dens are required
# # to achieve roughly the same average degrees. For example:
# n <- 30
# dens <- -2.0
# SN <- SimulateNetworks(n, M, rate, dens, rec, tt, c3, Vaego, Vaalt, Vasim,
#                        Vbego, Vbalt, Vbsim)
# 
# # Results:
# SN[[1]]
# # the same as
# SN$covara
# 
# SN[[2]]
# # the same as
# SN$covarb
# 
# SN[[3]]
# # the same as
# SN$networks




