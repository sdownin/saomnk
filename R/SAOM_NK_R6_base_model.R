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
library(uuid)



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
    UUID = NULL,
    DIR_OUTPUT = NULL,
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
    bipartite_matrix_init = NULL, ## cache the init matrix for future reference and plotting
    #
    bipartite_matrix = NULL,
    social_matrix = NULL,
    search_matrix = NULL,
    #
    bipartite_matrix_waves = list(), ## list of simulated networks for extended multi-wave simulation; implements dynamic strategy choice/changes
    #
    bi_env_arr = array(), ## bipartite matrix array from chain of ministeps
    #
    bipartite_rsienaDV = NULL,
    social_rsienaDV = NULL,
    search_rsienaDV = NULL,
    ##----- COVARIATES -----------------
    strat_1_coCovar = NULL,        ## constant strategy covariate  (M-vector)
    strat_2_coCovar = NULL,
    strat_3_coCovar = NULL,
    strat_4_coCovar = NULL,
    #
    strat_1_varCovar = NULL,        ## time-varying strategy covariate  (MxT matrix) for T periods
    strat_2_varCovar = NULL,
    #
    strat_1_coDyadCovar = NULL,    ## constant social network covariate (MxM matrix)
    strat_2_coDyadCovar = NULL,
    #
    strat_1_varDyadCovar = NULL,    ## time varying social network covariate (MxMxT array) for T periods
    strat_2_varDyadCovar = NULL,
    #
    strat_1_interaction = NULL,
    strat_2_interaction = NULL,
    #
    component_1_coCovar = NULL,     ## constant component covariate  (N-vector)
    component_2_coCovar = NULL,
    #
    component_1_varCovar = NULL,     ## time varying component covariate  (NxT matrix) for T periods
    component_2_varCovar = NULL,
    #
    component_1_coDyadCovar = NULL,  ## constant interaction matrix covariate (NxN matrix)
    component_2_coDyadCovar = NULL,
    #
    component_1_varDyadCovar = NULL,  ## time varying interaction matrix covariate (NxNxT array) for T periods
    component_2_varDyadCovar = NULL,
    #
    component_1_interaction = NULL, 
    component_2_interaction = NULL,
    #
    interaction_1 = NULL,
    interaction_2 = NULL,
    interaction_3 = NULL,
    interaction_4 = NULL,
    #
    # component_3_coCovar = NULL,
    ##------/end covariates ------------
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
      self$UUID <- UUIDgenerate(use.time = T)
      self$DIR_OUTPUT <- ifelse(is.null(config_environ_params[['dir_output']]),
                                getwd(),
                                config_environ_params[['dir_output']])
      # self$P_change <- config_environ_params[['P_change']]
      #
      # self$bipartite_igraph <- self$generate_bipartite_igraph()
      # self$social_network <- self$project_social_space()
      # self$search_landscape <- self$project_search_space()
      
      ##--------- INIT COMPONENT MATRIX STRUCTURE --------------------
      component_mat_init_type <- config_environ_params[['component_matrix_start']]
      if ( component_mat_init_type  == 'rand' ) 
      {
        default_seed <- 123
        self$rsiena_env_seed <- ifelse(is.null(config_environ_params[['rand_seed']]), 
                                       default_seed,
                                       config_environ_params[['rand_seed']])
        ## 
        start_bipartite_matrix <- self$random_bipartite_matrix(self$rsiena_env_seed)
        
        # } else if (component_mat_init_type == 'modular') {
        # } else if (component_mat_init_type == 'triangular') {  ## 
          
      } 
      else 
      {
        stop(sprintf('Component matrix init type not implemented: %s', component_mat_init_type))
      }
      
      ## Keep init matrix
      self$bipartite_matrix_init <- start_bipartite_matrix
      ## SET BIPARTITE NETWORK SYSTEM FROM INIT matrix
      self$set_system_from_bipartite_matrix( start_bipartite_matrix )
      ##
      self$TIMESTAMP <- round( as.numeric(Sys.time())*100 )
      self$SIM_NAME <- config_environ_params[['name']]
    },
    
    #
    set_system_from_bipartite_matrix = function(bipartite_matrix) {
      bipartite_igraph <- igraph::graph_from_biadjacency_matrix(bipartite_matrix, directed = F, mode = 'all')
      self$set_system_from_bipartite_igraph(bipartite_igraph)
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
      self$bipartite_matrix <- igraph::as_biadjacency_matrix(bipartite_igraph,attr = 'weight', sparse = F)
      self$social_matrix <- igraph::as_adjacency_matrix(self$social_igraph, attr = 'weight', sparse = F)
      self$search_matrix <- igraph::as_adjacency_matrix(self$search_igraph, attr = 'weight', sparse = F)
    },
    
    # Generate a random bipartite network matrix
    random_bipartite_matrix = function(rand_seed = 123) {
      ## If BI_PROB==1 or 0, return equivalent matrix (full or empty) 
      if (self$BI_PROB %in% c(0, 1))
        return(matrix(self$BI_PROB, nrow = self$M, ncol = self$N))
      # Else sample random matrix
      set.seed(rand_seed)  # For reproducibility
      probs <- c( 1 - self$BI_PROB, self$BI_PROB )
      bipartite_matrix <- matrix(sample(0:1, self$M * self$N, replace = TRUE, prob = probs ),
                                 nrow = self$M, ncol = self$N)
      return(bipartite_matrix)
      # # bipartite_net <- network(bipartite_matrix, bipartite = TRUE, directed = TRUE)
      # bipartite_igraph <- igraph::graph_from_biadjacency_matrix(bipartite_matrix, directed = F, mode = 'all')
      # return(bipartite_igraph)
    },
    
    # Generate a random bipartite network object
    random_bipartite_igraph = function(rand_seed = 123) {
      bipartite_matrix <- self$random_bipartite_matrix(rand_seed)
      # bipartite_net <- network(bipartite_matrix, bipartite = TRUE, directed = TRUE)
      bipartite_igraph <- igraph::graph_from_biadjacency_matrix(bipartite_matrix, directed = F, mode = 'all')
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
      if ( i >= 1 && i <= dim(m)[1] && j >= 1 && j <= dim(m)[2] ) {
        m[i,j] <-   1 - m[i,j] 
      }
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
      
      else if (eff$effect == 'totInDist2')
      {
        
        # activity_covar <- varCovar(activity_data)
        # effects <- includeEffects(effects, egoX, interaction1 = "activity_covar")
        
        self$rsiena_effects <- includeEffects(self$rsiena_effects,  totInDist2, ## get network statistic function from effect name (character)
                                              name = eff$dv_name, 
                                              interaction1 = eff$interaction1,
                                              fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  totInDist2, 
                                         interaction1 = eff$interaction1,
                                         name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      }
      
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
      else if (eff$effect == 'altXOutAct')
      {
        self$rsiena_effects <- includeEffects(self$rsiena_effects,  altXOutAct, ## get network statistic function from effect name (character)
                                              name = eff$dv_name, 
                                              interaction1 = eff$interaction1,
                                              fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  altXOutAct, 
                                         interaction1 = eff$interaction1,
                                         name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      }
      else if (eff$effect == 'homXOutAct')
      {
        self$rsiena_effects <- includeEffects(self$rsiena_effects,  homXOutAct, ## get network statistic function from effect name (character)
                                              name = eff$dv_name, 
                                              interaction1 = eff$interaction1,
                                              fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  homXOutAct, 
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
      ## DYADIC COVARIATE EFFECT
      else if (eff$effect == 'XWX')
      {
        self$rsiena_effects <- includeEffects(self$rsiena_effects,  XWX, ## get network statistic function from effect name (character)
                                              name = eff$dv_name, 
                                              interaction1 = eff$interaction1,
                                              fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  XWX, 
                                         interaction1 = eff$interaction1,
                                         name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      }
      else if (eff$effect == 'X')
      {
        self$rsiena_effects <- includeEffects(self$rsiena_effects,  X, ## get network statistic function from effect name (character)
                                              name = eff$dv_name, 
                                              interaction1 = eff$interactin1,
                                              fix = eff$fix)
        self$rsiena_effects <- setEffect(self$rsiena_effects,  X, 
                                         interaction1 = eff$interaction1,
                                         name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      }
      # else if (eff$effect == 'XWX1')
      # {
      #   self$rsiena_effects <- includeEffects(self$rsiena_effects,  XWX1, ## get network statistic function from effect name (character)
      #                                         name = eff$dv_name, 
      #                                         interaction1 = eff$interaction1,
      #                                         fix = eff$fix)
      #   self$rsiena_effects <- setEffect(self$rsiena_effects,  XWX1, 
      #                                    interaction1 = eff$interaction1,
      #                                    name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      # }
      # else if (eff$effect == 'XWX2')
      # {
      #   self$rsiena_effects <- includeEffects(self$rsiena_effects,  XWX2, ## get network statistic function from effect name (character)
      #                                         name = eff$dv_name, 
      #                                         interaction1 = eff$interaction1,
      #                                         fix = eff$fix)
      #   self$rsiena_effects <- setEffect(self$rsiena_effects,  XWX2, 
      #                                    interaction1 = eff$interaction1,
      #                                    name = eff$dv_name, parameter = eff$parameter,  fix = eff$fix)
      # }
      else 
      {
        ##
        print(eff)
        stop('effect printed above was not included in rsiena_effects.')
      }
      
    },
    
    
    
    
    include_rsiena_interaction_from_eff_list = function(interact) {
      #
      if ( ! 'manual_interaction' %in% names(self$rsiena_effects) ) {
        self$rsiena_effects$manual_interaction <- NA
      }
      #
      if(all(c('totInDist2','X') %in% interact$effects))  {
        # ##**altX**
        # component_cov <- self$rsiena_data$cCovars[[ interact$interaction1 ]]
        # if (attr(component_cov, 'centered')) {
        #   ## return to uncentered original data
        #   component_cov <-  component_cov + attr(component_cov, 'mean')
        # }
        
        # ##**X**
        # dyad_cov <- self$rsiena_data$dycCovars[[ interact$interaction1 ]]
        # if (attr(component_cov, 'centered')) {
        #   ## return to uncentered original data
        #   dyad_cov <-  dyad_cov + attr(dyad_cov, 'mean')
        # }
        # ##**XWX**
        # component_dyad_cov <- self$rsiena_data$dycCovars[[ interact$interaction2 ]]
        # if ( all(diag(component_dyad_cov)==0)  &  all(c(component_dyad_cov) %in% c(0,1)) ) {
        #   diag(component_dyad_cov) <- 1
        # }
        # ##
        # inter_mat <- outer(component_cov, component_cov, "*") * component_dyad_cov
        #
        #
        self$rsiena_effects <- includeInteraction(self$rsiena_effects, totInDist2, X,  ## get network statistic function from effect name (character)
                                                  name = interact$dv_name, # interaction1 = eff$interaction1,
                                                  # parameter = interact$parameter,
                                                  fix = interact$fix,
                                                  interaction1 = c(interact$interaction1, interact$interaction2) )
        ## This workaround adjusts parameter of interaction without name
        ##**TODO: Check for RSiena official implementation**
        self$rsiena_effects[self$rsiena_effects$include,][sum(self$rsiena_effects$include),'parm'] <- interact$parameter
        ## NEED TO SET INTERACTNG VARIABLE NAMES HERE FOR USE IN LATER COMPUTATIONS (theta_matrix)
        self$rsiena_effects[self$rsiena_effects$include,][sum(self$rsiena_effects$include),'manual_interaction'] <- interact$effect
        # self$rsiena_effects <- setEffect(self$rsiena_effects, unspInt,
        #                                  name = interact$dv_name,
        #                                  parameter = interact$parameter,
        #                                  fix = interact$fix,
        #                                  interaction1 = c(interact$interaction1, interact$interaction2))  #,
        
      } else if(all(c('egoX','XWX') %in% interact$effects))  {
        
        self$rsiena_effects <- includeInteraction(self$rsiena_effects, egoX, XWX,  ## get network statistic function from effect name (character)
                                                  name = interact$dv_name, # interaction1 = eff$interaction1,
                                                  # parameter = interact$parameter,
                                                  fix = interact$fix,
                                                  interaction1 = c(interact$interaction1, interact$interaction2) )
        ## This workaround adjusts parameter of interaction without name
        ##**TODO: Check for RSiena official implementation**
        self$rsiena_effects[self$rsiena_effects$include,][sum(self$rsiena_effects$include),'parm'] <- interact$parameter
        ## NEED TO SET INTERACTNG VARIABLE NAMES HERE FOR USE IN LATER COMPUTATIONS (theta_matrix)
        self$rsiena_effects[self$rsiena_effects$include,][sum(self$rsiena_effects$include),'manual_interaction'] <- interact$effect
        # self$rsiena_effects <- setEffect(self$rsiena_effects, unspInt,
        #                                  name = interact$dv_name,
        #                                  parameter = interact$parameter,
        #                                  fix = interact$fix,
        #                                  interaction1 = c(interact$interaction1, interact$interaction2))  #,
        
      } else if(all(c('inPopX','X') %in% interact$effects))  {
        
        self$rsiena_effects <- includeInteraction(self$rsiena_effects, inPopX, X,  ## get network statistic function from effect name (character)
                                                  name = interact$dv_name, # interaction1 = eff$interaction1,
                                                  # parameter = interact$parameter,
                                                  fix = interact$fix,
                                                  interaction1 = c(interact$interaction1, interact$interaction2) )
        ## This workaround adjusts parameter of interaction without name
        ##**TODO: Check for RSiena official implementation**
        self$rsiena_effects[self$rsiena_effects$include,][sum(self$rsiena_effects$include),'parm'] <- interact$parameter
        ## NEED TO SET INTERACTNG VARIABLE NAMES HERE FOR USE IN LATER COMPUTATIONS (theta_matrix)
        self$rsiena_effects[self$rsiena_effects$include,][sum(self$rsiena_effects$include),'manual_interaction'] <- interact$effect
        # self$rsiena_effects <- setEffect(self$rsiena_effects, unspInt,
        #                                  name = interact$dv_name,
        #                                  parameter = interact$parameter,
        #                                  fix = interact$fix,
        #                                  interaction1 = c(interact$interaction1, interact$interaction2))  #,
        
      } else if(all(c('inPopX','egoX') %in% interact$effects))  {
        
        self$rsiena_effects <- includeInteraction(self$rsiena_effects, inPopX, egoX,  ## get network statistic function from effect name (character)
                                                  name = interact$dv_name, # interaction1 = eff$interaction1,
                                                  # parameter = interact$parameter,
                                                  fix = interact$fix,
                                                  interaction1 = c(interact$interaction1, interact$interaction2) )
        ## This workaround adjusts parameter of interaction without name
        ##**TODO: Check for RSiena official implementation**
        self$rsiena_effects[self$rsiena_effects$include,][sum(self$rsiena_effects$include),'parm'] <- interact$parameter
        ## NEED TO SET INTERACTNG VARIABLE NAMES HERE FOR USE IN LATER COMPUTATIONS (theta_matrix)
        self$rsiena_effects[self$rsiena_effects$include,][sum(self$rsiena_effects$include),'manual_interaction'] <- interact$effect
        # self$rsiena_effects <- setEffect(self$rsiena_effects, unspInt,
        #                                  name = interact$dv_name,
        #                                  parameter = interact$parameter,
        #                                  fix = interact$fix,
        #                                  interaction1 = c(interact$interaction1, interact$interaction2))  #,
        
        
      }  else {
        
        stop(sprintf('Interaction `%s*%s` not yet implemented.', interact$effects[1], interact$effects[2]))
        
      }


    },
    
    
    
    # ##-----------------
    getTotInDist2 = function(bipartite_matrix, actor_covariate, interaction_type = "absdiff") {
      # Ensure input is a matrix
      if (!is.matrix(bipartite_matrix)) {
        stop("Bipartite network must be an M × N matrix.")
      }

      # Ensure actor covariate is a vector of length M
      M <- nrow(bipartite_matrix)
      if (length(actor_covariate) != M) {
        stop("Actor covariate must be a vector of length M (number of actors).")
      }

      # Step 1: Compute the N × N projection (component adjacency matrix)
      A <- t(bipartite_matrix) %*% bipartite_matrix  # Project bipartite network onto components
      diag(A) <- 0  # Remove self-loops

      # Step 2: Compute the dyadic transformation of the actor covariate
      actor_matrix <- outer(actor_covariate, actor_covariate, FUN=switch(
        interaction_type,
        "absdiff" = function(x, y) abs(x - y),
        "product" = function(x, y) x * y,
        "sum" = function(x, y) x + y,
        stop("Invalid interaction type. Choose 'absdiff', 'product', or 'sum'.")
      ))

      # Step 3: Apply the transformed covariate to weight second-degree paths
      A_squared <- A %*% A  # A^2 counts 2-step walks
      weighted_A2 <- A_squared * (t(bipartite_matrix) %*% actor_matrix %*% bipartite_matrix)

      # Step 4: Compute totInDist2 statistic as the sum of incoming 2-step weighted paths
      totInDist2_values <- rowSums(weighted_A2)

      # Return as named vector
      names(totInDist2_values) <- paste0("Component_", seq_along(totInDist2_values))

      return(totInDist2_values)
    },
    
    
    
    ##
    get_struct_mod_stats_mat_from_bi_mat = function(bi_env_mat, type='all') {
      #
      xActorDegree      <- rowSums(bi_env_mat, na.rm=T)
      xComponentDegree  <- colSums(bi_env_mat, na.rm=T)
      #
      # eff <- m1$rsiena_effects[m1$rsiena_effects$include, ]
      efflist <- c(
        self$config_structure_model$dv_bipartite$effects,
        self$config_structure_model$dv_bipartite$coCovars,
        self$config_structure_model$dv_bipartite$varCovars,
        self$config_structure_model$dv_bipartite$coDyadCovars,
        self$config_structure_model$dv_bipartite$varDyadCovars,
        self$config_structure_model$dv_bipartite$interactions
      )
      #
      neffs <- length(efflist)
      ### empty matrix to hold actor network statistics
      effnames <- sapply(efflist, function(x) x$effect, simplify = T)
      effparams <- sapply(efflist, function(x) x$parameter, simplify = T)
      #
      mat <- matrix(rep(0, self$M * neffs ), nrow=self$M, ncol=neffs )
      colnames(mat) <- effnames
      rownames(mat) <- 1:self$M
      #  
      for (i in 1:neffs)
      {
        # print(i)
        item <- efflist[[ i ]]
        # print('DEBUG  get_struct_mod_stats_mat_from_bi_mat() ')
        # print(item)
        
        ## network statistics dataframe
        #
        if (item$effect == 'density' )  {
          
          mat[ , i] <- c( xActorDegree )
          
        } else if (item$effect == 'outAct' )  {
          
          mat[ , i] <- c( xActorDegree^2 )
        } else if (item$effect == 'inPop' ) {
          
          mat[ , i] <- c( bi_env_mat %*% (xComponentDegree + 1) )
          
          # else if (item$effect == 'transTriads' ) {
          #   stat <- 
          # } else if (item$effect == 'cycle4' ) {
          #   stat <- 
          # }
        } else if (item$effect == 'egoX') {
          
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
          mat[ , i] <- c( covar * xActorDegree ) ##**vector element-wise multiplication by rows of covar matrix, or elements of covar array
            
        } else if (item$effect == 'altX') {
          
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
          mat[ , i] <- rowSums( covarComponentMat * bi_env_mat, na.rm=T ) ##**vector element-wise multiplication by rows of covar matrix, or elements of covar array
          
        } else if (item$effect == 'outActX') { ## interaction1 component_coCovar
          
          ## N-vector of component covariate
          covar <- item$x 
          # MxN matrix of row-stacked component covariate (repeated for each actor)
          covarComponentMat <- matrix(rep(covar, self$M), nrow=self$M, ncol=self$N, byrow = TRUE)
          ## M-vector of actor's squared sum of component-covariate-weighted component connections (weighted version of the squared degree)
          mat[ , i] <- xActorDegree * rowSums( covarComponentMat * bi_env_mat, na.rm = T)
        
        } else if (item$effect == 'inPopX') { ## interaction1 strat_coCovar
          covar <- item$x ## M-vector of actor strategy covariate
          ## MxN matrix holding actor strategy covariate as columns stacked for each component
          covarActorMat <- matrix(rep(covar, self$N), nrow=self$M, ncol=self$N,  byrow = FALSE)
          ## N-vector of square roots of component weights (sum of actor covariate for the component's connected actors)
          component_weights_from_actor_stats <-  colSums(covarActorMat * bi_env_mat, na.rm=T)
          ## M-vector of actor sum of it's connected component weights (which are computed as the sum of the connected actor covariates)
          mat[ , i] <- rowSums( bi_env_mat * component_weights_from_actor_stats, na.rm=T ) ##**vector element-wise multiplication by rows of covar matrix, or elements of covar array
        
        } else if (item$effect == 'XWX') { ## interaction1 strat_coCovar
          
          covar <- item$x  ## NxN matrix
          ## MxM matrix of inter-actor connections weighted by component covarite matrix
          interactor_cov_w <- bi_env_mat %*% covar %*% t(bi_env_mat)
          ## covert to M-vector of actor attributes
          mat[ , i] <- rowSums( interactor_cov_w, na.rm=T ) ##**TODO: CHECK**
          
        }  else if (item$effect == 'X') { ## interaction1 strat_coCovar
          
          covar <- item$x  ## NxN matrix
          ## MxM matrix of inter-actor connections weighted by component covarite matrix
          interactor_cov_w <- bi_env_mat * (covar - mean(c(covar, na.rm=T)) )
          ## covert to M-vector of actor attributes
          mat[ , i] <- rowSums( interactor_cov_w, na.rm=T ) ##**TODO: CHECK**
          # stop('implement altX|XWX .')
          
        }   else if (item$effect == 'totInDist2') { ## interaction1 strat_coCovar
          
          ## M-vector
          covar <- item$x  
          # 1xN matrix
          component_sums_w_by_actor_covar <-  covar %*% bi_env_mat 
          #
          compo_w_stacked_mat <- matrix(rep(component_sums_w_by_actor_covar, self$M), byrow=T, ncol=self$N)
          ## covert to M-vector of actor attributes
          mat[ , i ] <-  rowSums( compo_w_stacked_mat * bi_env_mat, na.rm=T )
          # mat[ , i] <- self$getTotInDist2(bi_env_mat, covar, interaction_type = 'absdiff') ##**TODO: CHECK**
          
        }  else if(grepl('[|]', item$effect))  {
          
          next ## Skip interactions to add them as interactions of the stat matrix 
          
        } else {
          
          cat(sprintf('\n\nEffect not yet implemented: `%s\n\n`', item$effect))
          next
          
        }
        
        ####
        # else if (item$effect == 'altX|XWX') ## interaction1 strat_coCovar
        # {
        #   # covar <- item$x  ## NxN matrix
        #   # ## MxM matrix of inter-actor connections weighted by component covarite matrix
        #   # interactor_cov_w <- bi_env_mat %*% covar %*% t(bi_env_mat)
        #   # ## covert to M-vector of actor attributes
        #   # stat <- rowSums( interactor_cov_w, na.rm=T ) ##**TODO: CHECK**
        #   stop('implement altX|XWX .')
        # }
        # else if (item$effect == 'totInDist2|X') ## interaction1 strat_coCovar
        # {
        #   # covar <- item$x  ## NxN matrix
        #   # ## MxM matrix of inter-actor connections weighted by component covarite matrix
        #   # interactor_cov_w <- bi_env_mat %*% covar %*% t(bi_env_mat)
        #   # ## covert to M-vector of actor attributes
        #   # stat <- rowSums( interactor_cov_w, na.rm=T ) ##**TODO: CHECK**
        #   stop('implement altX|XWX .')
        # }
        
        #
        # mat[ , i] <- stat
        
      }
      
      return(mat)
    }
    
    
    
  )
    
    
)




# # ##-----------------
# compute_totInDist2 <- function(bipartite_matrix, actor_covariate, interaction_type = "absdiff") {
#   # Ensure input is a matrix
#   if (!is.matrix(bipartite_matrix)) {
#     stop("Bipartite network must be an M × N matrix.")
#   }
#   
#   # Ensure actor covariate is a vector of length M
#   M <- nrow(bipartite_matrix)
#   if (length(actor_covariate) != M) {
#     stop("Actor covariate must be a vector of length M (number of actors).")
#   }
#   
#   # Step 1: Compute the N × N projection (component adjacency matrix)
#   A <- t(bipartite_matrix) %*% bipartite_matrix  # Project bipartite network onto components
#   diag(A) <- 0  # Remove self-loops
#   
#   # Step 2: Compute the dyadic transformation of the actor covariate
#   actor_matrix <- outer(actor_covariate, actor_covariate, FUN=switch(
#     interaction_type,
#     "absdiff" = function(x, y) abs(x - y),
#     "product" = function(x, y) x * y,
#     "sum" = function(x, y) x + y,
#     stop("Invalid interaction type. Choose 'absdiff', 'product', or 'sum'.")
#   ))
#   
#   # Step 3: Apply the transformed covariate to weight second-degree paths
#   A_squared <- A %*% A  # A^2 counts 2-step walks
#   weighted_A2 <- A_squared * (t(bipartite_matrix) %*% actor_matrix %*% bipartite_matrix)
#   
#   # Step 4: Compute totInDist2 statistic as the sum of incoming 2-step weighted paths
#   totInDist2_values <- rowSums(weighted_A2)
#   
#   # Return as named vector
#   names(totInDist2_values) <- paste0("Component_", seq_along(totInDist2_values))
#   
#   return(totInDist2_values)
# }

# # ===========================
# # 📌 Example Usage
# # ===========================
# 
# # Set parameters
# M <- 10  # Number of actors
# N <- 15  # Number of components
# 
# # Generate a random M x N bipartite adjacency matrix (undirected)
# set.seed(123)
# bipartite_matrix <- matrix(sample(0:1, M * N, replace=TRUE, prob=c(0.7, 0.3)), nrow=M, ncol=N)
# 
# # Generate a random actor covariate (M × 1 vector)
# actor_covariate <- runif(M, 0, 1)
# 
# # Compute the totInDist2 statistic with actor-based interaction
# totInDist2_results <- compute_totInDist2(bipartite_matrix, actor_covariate, interaction_type="absdiff")
# 
# # Print results
# print(totInDist2_results)



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



