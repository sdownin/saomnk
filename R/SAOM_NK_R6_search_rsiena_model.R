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
SaomNkRSienaBiEnv_base <- R6Class( 
  ##
  "SaomNkRSienaBiEnv_base",
  ##
  
  public = list(
    #
    ITERATION = 0,
    TIMESTAMP = NULL,
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
    
    
    ##------ Helper functions ---------------
    
    # Define the bipartite (2-mode) network matrix self$toggle function:
    toggleBiMat = function(m,i,j){
      m[i,j] <-  ( 1 - m[i,j] )
      return(m)
    },
    
    ##
    ##
    ##
    get_jaccard_index = function(m0, m1) {
      diffvec <-  c(m1 - m0)   ## new m1 - old m0
      cnt_maintain <- sum( m1 * m0 )
      cnt_change <- sum( diffvec != 0 )  ## sum = count true cases (added + dropped)
      return( cnt_maintain / (cnt_maintain + cnt_change) )
    },
    
    ##
    exists = function(x){
      return(!is.null(x) && !is.na(x) && !is.nan(x))
    },
    
    # Define the one-mode network matrix self$toggle function:
    toggle = function(m,i,j){
      if (i != j) {
        # print(sprintf('test i %s != j %s',i,j))
        m[i,j] <-  ( 1 - m[i,j] )
      }
      return(m)
    },
    
    
    ##-----/helper-------------------------
    
    
    ##
    include_rsiena_effect_from_eff_list = function(eff) {
      ##---------- 2+ Effects Combination --------------------
      if (length(eff$effect)>1)
      {
        if (all(eff$effect %in% c('egoX', 'inPopX'))) {
          self$rsiena_effects <- includeEffects(self$rsiena_effects,  egoX, inPopX, ## get network statistic function from effect name (character)
                                                name = eff$dv_name, 
                                                interaction1 = eff$interaction1,
                                                fix = eff$fix)
          self$rsiena_effects <- setEffect(self$rsiena_effects,  egoX, inPopX, 
                                           interaction1 = eff$interaction1,
                                           name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
        }
        else 
        {
          stop('Effect combination not yet implemented.')
        }
        
        return(NULL)
      }

      
      ##---------- 1 Efect --------------------------
      
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
          ## sqrt( rowSums( covarComponentMat * bi_env_mat, na.rm = T)^2 ) ## square root of the square of the degree == degree
          stat <- rowSums( covarComponentMat * bi_env_mat, na.rm = T)^2 
        }
        else if (item$effect == 'inPopX') ## interaction1 strat_coCovar
        {
          covar <- item$x ## M-vector of actor strategy covariate
          ## MxN matrix holding actor strategy covariate as columns stacked for each component
          covarActorMat <- matrix(rep(covar, self$N), nrow=self$M, ncol=self$N,  byrow = FALSE)
          ## N-vector of square roots of component weights (sum of actor covariate for the component's connected actors)
          component_weights_from_actor_stats <-  colSums(covarActorMat * bi_env_mat, na.rm=T)
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




# ##-----------------




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



