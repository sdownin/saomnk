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
library(grid)
library(rvest)

###############   DEPENDENCIES ###############################
# SaomNkRSienaBiEnv_base <- source(file.path(dir_proj, 'SAOM_NK_R6_base_model.R'))$value
# ##    DEPENDENCIES     ######  
# SaomNkRSienaBiEnv_base <- source(file.path(dir_proj, 'SAOM_NK_R6_base_model.R'))$value



###############################################################################
###############################################################################

# Define the SAOM_NK_Enhanced class with local search, complex search strategies, network metrics, and visualization
SaomNkRSienaBiEnv <- R6Class( 
  ##
  "SaomNkRSienaBiEnv",
  ##
  # inherit = SaomNkRSienaBiEnv_base,
  inherit =  source(file.path('C://Users//sdr8y//OneDrive - University of Missouri//Research//Search_networks//SaoMNK//R', 'SAOM_NK_R6_base_model.R'))$value,
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
      plot_init <- ifelse(is.null(config_environ_params[['plot_init']]), 
                               TRUE, 
                               config_environ_params[['plot_init']])
      if (plot_init) {
        self$plot_bipartite_system_from_mat(self$bipartite_matrix, RSIENA_ITERATION=0, return_plot=F)
      }
      # self$visualize_system_bi_env_rsiena(plot_save = TRUE)
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
      
      
      structeffs <- sapply(structure_model$dv_bipartite$effects, function(x)x$parameter)
      names(structeffs) <- sapply(structure_model$dv_bipartite$effects, function(x)x$effect)
      
      
      coCovars     <- sapply(structure_model$dv_bipartite$coCovars, function(x)x$effect)
      varCovars     <- sapply(structure_model$dv_bipartite$varCovars, function(x)x$effect)
      coDyadCovars  <- sapply(structure_model$dv_bipartite$coDyadCovars, function(x)x$effect)
      varDyadCovars <- sapply(structure_model$dv_bipartite$varDyadCovars, function(x)x$effect)

      coCovarTypes     <- sapply(structure_model$dv_bipartite$coCovars, function(x)x$interaction1)
      varCovarTypes     <- sapply(structure_model$dv_bipartite$varCovars, function(x)x$interaction1)
      coDyadCovarTypes  <- sapply(structure_model$dv_bipartite$coDyadCovars, function(x)x$interaction1)
      varDyadCovarTypes <- sapply(structure_model$dv_bipartite$varDyadCovars, function(x)x$interaction1)

      component_coCovar_ids     <- grep('self\\$component.+_coCovar', coCovarTypes) ## ex: "self$component_1_coCovar"
      component_varCovar_ids    <- grep('self\\$component.+_varCovar', varCovarTypes) ## ex: "self$component_1_coCovar"
      component_coDyadCovar_ids  <- grep('self\\$component.+_coDyadCovar', coDyadCovarTypes) ## ex: "self$component_1_coCovar"
      component_varDyadCovar_ids <- grep('self\\$component.+_varDyadCovar', varDyadCovarTypes) ## ex: "self$component_1_coCovar"

      strat_coCovar_ids     <- grep('self\\$strat.+_coCovar', coCovarTypes) ## ex: "self$component_1_coCovar"
      strat_varCovar_ids    <- grep('self\\$strat.+_varCovar', varCovarTypes) ## ex: "self$component_1_coCovar"
      strat_coDyadCovar_ids  <- grep('self\\$strat.+_coDyadCovar', coDyadCovarTypes) ## ex: "self$component_1_coCovar"
      strat_varDyadCovar_ids <- grep('self\\$strat.+_varDyadCovar', varDyadCovarTypes) ## ex: "self$component_1_coCovar"
      # stratDV_ids     <- grep('self\\$strat_\\d{1,2}_coCovar', covDvTypes ) ## ex: "self$strat_1_coCovar"

      ncompo_coCovar     <- length( component_coCovar_ids )
      ncompo_varCovar    <- length( component_varCovar_ids )
      ncompo_coDyadCovar  <- length( component_coDyadCovar_ids )
      ncompo_varDyadCovar <- length( component_varDyadCovar_ids )

      nstrat_coCovar     <- length( strat_coCovar_ids )
      nstrat_varCovar    <- length( strat_varCovar_ids )
      nstrat_coDyadCovar  <- length( strat_coDyadCovar_ids )
      nstrat_varDyadCovar <- length( strat_varDyadCovar_ids )
      
      
      ## input list of variable for RSiena model
      input_varlist <- list(`self$bipartite_rsienaDV`=self$bipartite_rsienaDV)
      ##-------------------------------------------------------------
      ## COMPONENTS ##
      if (ncompo_coCovar) {
        for (i in 1:ncompo_coCovar) {
          property <- sprintf('component_%s_coCovar', i)
          eff <- structure_model$dv_bipartite$coCovars[[ component_coCovar_ids[i] ]]
          self[[property]] <- coCovar(eff$x, nodeSet = c('COMPONENTS'))
          input_varlist[[sprintf('self$%s',property)]] <-  self[[property]]
        }
      }
      if (ncompo_varCovar) {
        for (i in 1:ncompo_varCovar) {
          property <- sprintf('component_%s_coCovar', i)
          eff <- structure_model$dv_bipartite$varCovars[[ component_varCovar_ids[i] ]]
          self[[property]] <- varCovar(eff$x, nodeSet = c('COMPONENTS'))
          input_varlist[[sprintf('self$%s',property)]] <-  self[[property]]
        }
      }
      if (ncompo_coDyadCovar) {
        for (i in 1:ncompo_coDyadCovar) {
          property <- sprintf('component_%s_coDyadCovar', i)
          eff <- structure_model$dv_bipartite$coDyadCovar[[ component_coCovar_ids[i] ]]
          self[[property]] <- coDyadCovar(eff$x, nodeSet = c('COMPONENTS','COMPONENTS'))
          input_varlist[[sprintf('self$%s',property)]] <-  self[[property]]
        }
      }
      if (ncompo_varDyadCovar) {
        for (i in 1:ncompo_varDyadCovar) {
          property <- sprintf('component_%s_varDyadCovar', i)
          eff <- structure_model$dv_bipartite$varDyadCovar[[ component_varDyadCovar_ids[i] ]]
          self[[property]] <- varDyadCovar(eff$x, nodeSet = c('COMPONENTS','COMPONENTS'))
          input_varlist[[sprintf('self$%s',property)]] <-  self[[property]]
        }
      }
      ##-------------------------------------------------------------
      ## STRATEGIES ##
      if (nstrat_coCovar) {
        for (i in 1:nstrat_coCovar) {
          property <- sprintf('strat_%s_coCovar', i)
          eff <- structure_model$dv_bipartite$coCovars[[ strat_coCovar_ids[i] ]]
          self[[property]] <- coCovar(eff$x, nodeSet = c('ACTORS'))
          input_varlist[[sprintf('self$%s',property)]] <-  self[[property]]
        }
      }
      if (nstrat_varCovar) {
        for (i in 1:nstrat_varCovar) {
          property <- sprintf('strat_%s_coCovar', i)
          eff <- structure_model$dv_bipartite$varCovars[[ strat_varCovar_ids[i] ]]
          self[[property]] <- varCovar(eff$x, nodeSet = c('ACTORS'))
          input_varlist[[sprintf('self$%s',property)]] <-  self[[property]]
        }
      }
      if (nstrat_coDyadCovar) {
        for (i in 1:nstrat_coDyadCovar) {
          property <- sprintf('strat_%s_coDyadCovar', i)
          eff <- structure_model$dv_bipartite$coDyadCovar[[ strat_coDyadCovar_ids[i] ]]
          self[[property]] <- coDyadCovar(eff$x, nodeSet = c('ACTORS','ACTORS'))
          input_varlist[[sprintf('self$%s',property)]] <-  self[[property]]
        }
      }
      if (nstrat_varDyadCovar) {
        for (i in 1:nstrat_varDyadCovar) {
          property <- sprintf('strat_%s_varDyadCovar', i)
          eff <- structure_model$dv_bipartite$varDyadCovar[[ strat_varDyadCovar_ids[i] ]]
          self[[property]] <- varDyadCovar(eff$x, nodeSet = c('ACTORS','ACTORS'))
          input_varlist[[sprintf('self$%s',property)]] <-  self[[property]]
        }
      }
      ##-------------------------------------------------------------
      #
      sienaDataCreate_args <- c( list(nodeSets=list(ACTORS, COMPONENTS)), input_varlist ) 
      rsiena_data <- do.call(sienaDataCreate, sienaDataCreate_args)
      
      
      return(rsiena_data)
    },
    
    
    
    
    ##
    get_rsiena_data_from_structure_model_v0 = function(structure_model) {
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
      if (ncomponent==2 & nstrat==1 ) 
      {
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV, self$component_1_coCovar,self$component_2_coCovar, self$strat_1_coCovar), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      if (ncomponent==2 & nstrat==2 ) 
      {
        strat_2_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[2] ]]
        self$strat_2_coCovar <- coCovar(strat_2_eff$x, nodeSet = 'ACTORS')
        #
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV, self$component_1_coCovar,self$component_2_coCovar,  self$strat_1_coCovar,self$strat_2_coCovar), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      if (ncomponent==2 & nstrat==3 ) 
      {
        strat_2_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[2] ]]
        self$strat_2_coCovar <- coCovar(strat_2_eff$x, nodeSet = 'ACTORS')
        strat_3_eff <- structure_model$dv_bipartite$coCovars[[ stratDV_ids[3] ]]
        self$strat_3_coCovar <- coCovar(strat_3_eff$x, nodeSet = 'ACTORS')
        #
        rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV, self$component_1_coCovar,self$component_2_coCovar,  self$strat_1_coCovar,self$strat_2_coCovar,self$strat_3_coCovar), nodeSets = list(ACTORS, COMPONENTS))
        return(rsiena_data)
      }
      if (ncomponent==2 & nstrat==4 ) 
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
      if (ncomponent==2 & nstrat==5 ) 
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
    
    ## SETS $rsiena_data property for RSiena::siena07()
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

      ## Simulation baseline nets should not be same; make one small change (self$toggle one dyad)
      ## @see https://www.stats.ox.ac.uk/~snijders/siena/NetworkSimulation.R
      bipartite_matrix1 <- bipartite_matrix
      bipartite_matrix2 <- bipartite_matrix
      # .i <- sample(1:self$M, 1)
      # .j <- sample(1:self$N, 1)
      # bipartite_matrix2 <-  self$toggleBiMat(bipartite_matrix2, .i, .j )
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
      ## Simulation baseline nets should not be same; make one small change (self$toggle one dyad)
      ## @see https://www.stats.ox.ac.uk/~snijders/siena/NetworkSimulation.R
      #
      # if(identical(bipartite_matrix1, bipartite_matrix2)) {
      #   .i <- sample(1:self$M, 1)
      #   .j <- sample(1:self$N, 1)
      #   bipartite_matrix2 <-  self$toggleBiMat(bipartite_matrix2, .i, .j )
      # }
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
          ##**TODO** Implement variable covariates(?)
        }
        
        for (j in 1:length(dv$coDyadCovars)) {
          cat(sprintf('\n coDyadCovars i=%s, j=%s\n', i, j))
          
          eff <- dv$coDyadCovars[[ j ]]
          
          print(eff)
          # print(eff)
          # stop('DEBUG')
          
          self$include_rsiena_effect_from_eff_list(eff)
        }
        
        for (j in 1:length(dv$varDyadCovars)) {
          ##**TODO** Implement variable covariates(?)
        }
        
      }
      
    },
  
    
    
    
    # #######
    # bipartite_matrix1 <- self$bipartite_matrix
    # bipartite_matrix2 <- self$bipartite_matrix
    # #######
    # # cat('\n\nDEBUG: called init_rsiena_model_from_structure_model_bipartite_matrix()\n\n')
    # #
    # ACTORS     <- sienaNodeSet(self$M, nodeSetName="ACTORS")
    # COMPONENTS <- sienaNodeSet(self$N, nodeSetName="COMPONENTS")
    # #
    # self$config_structure_model <- structure_model
    # structure_model_dvs <- names(structure_model)
    # #
    # # ## Simulation baseline nets should not be same; make one small change (self$toggle one dyad)
    # # ## @see https://www.stats.ox.ac.uk/~snijders/siena/NetworkSimulation.R
    # # #
    # # if(identical(bipartite_matrix1, bipartite_matrix2)) {
    # #   .i <- sample(1:self$M, 1)
    # #   .j <- sample(1:self$N, 1)
    # #   bipartite_matrix2 <-  self$toggleBiMat(bipartite_matrix2, .i, .j )
    # # }
    # ##
    # social_matrix1 <- bipartite_matrix1 %*% t(bipartite_matrix1)
    # search_matrix1 <- t(bipartite_matrix1) %*% bipartite_matrix1
    # ##
    # social_matrix2 <- bipartite_matrix2 %*% t(bipartite_matrix2)
    # search_matrix2 <- t(bipartite_matrix2) %*% bipartite_matrix2
    # ## init networks to duplicate for the init arrays (two network waves)
    # array_bi_net <- array(c(bipartite_matrix1, bipartite_matrix2), dim=c(self$M, self$N, 2) )
    # array_social <- array(c(social_matrix1, social_matrix2), dim=c(self$M, self$M, 2) )
    # array_search <- array(c(search_matrix1, search_matrix2), dim=c(self$N, self$N, 2) )
    # ## Drop information above binary ties for RSiena DVs
    # array_bi_net[ array_bi_net > 1 ] <- 1
    # array_social[ array_social > 1 ] <- 1
    # array_search[ array_search > 1 ] <- 1
    # if ('dv_social' %in% structure_model_dvs ) {
    #   self$social_rsienaDV <- sienaDependent(array_social, type='oneMode', nodeSet = 'ACTORS', allowOnly = F)
    # }
    # if ('dv_search' %in% structure_model_dvs) {
    #   self$search_rsienaDV <- sienaDependent(array_search, type='oneMode', nodeSet = 'COMPONENTS', allowOnly = F)
    # }
    # if ('dv_bipartite' %in% structure_model_dvs) {
    #   self$bipartite_rsienaDV <- sienaDependent(array_bi_net, type='bipartite', nodeSet =c('ACTORS', 'COMPONENTS'), allowOnly = F)
    # }
    
    ##
    preview_effects = function(structure_model, filter=TRUE) {
      #
      self$config_structure_model <- structure_model
      
      ##--------- I. SET BIPARTITE NETWORK DV ARRAY IF NULL -----------------
      # rsiena_model$sims[[sim_id]][[wave_id]][[bi_env_dv_id]] 
      ## Set up bipartite matrix RSiena dependent variable 
      if(is.null(self$bipartite_matrix)) 
          stop('bipartite_matrix is not set and no array_bi_net provided.')
      array_bi_net <- array(c(self$bipartite_matrix, self$bipartite_matrix), 
                            dim = c(self$M, self$N, 2))

      # 
      self$bipartite_rsienaDV <- sienaDependent(array_bi_net, 
                                                type='bipartite',
                                                nodeSet =c('ACTORS', 'COMPONENTS'), 
                                                allowOnly = F)
      
      ##----------- II. SET RSIENA DATA AND EFFECTS ---------------------------
      ##  1. RSiena data object
      self$rsiena_data <- self$get_rsiena_data_from_structure_model(structure_model)
      ##  2. Init effects
      self$rsiena_effects <- getEffects(self$rsiena_data)
      
      eff_filename <- file.path(self$DIR_OUTPUT, '_rsiena_effects_doc_')
      effectsDocumentation(self$rsiena_effects, type = 'html', display = F, filename = eff_filename)
      
      # ##--2. NETWORK: STRUCTURE EVOLUTION (structure_Model)-----
      # ##  2.1. Add effects from model objective function list
      # self$add_rsiena_effects(structure_model)
      
      # Parse the HTML
      html <- read_html(sprintf('%s.html', eff_filename))
      
      # html <- read_html('C:\\Users\\sdr8y\\OneDrive - University of Missouri\\Research\\Search_networks\\SaoMNK\\R\\_rsiena_effects_doc_.html')
      
      # Extract tables (returns a list of data.frames)
      efftab_list <- html %>% html_table(fill = TRUE)
      efftab <- efftab_list[[1]] # effect table is first in list
      
      if ( ! filter )
        return(efftab)
      
      efftable_dt <- datatable(
        efftab,               
        filter = "top",       # Adds filter boxes at the top
        options = list(
          pageLength = 10,   # Number of rows per page
          autoWidth = TRUE,   # Auto-adjust column width
          dom = 'lfrtip',     # Layout controls (search box, filters, etc.)
          scrollX = TRUE      # Horizontal scrolling if needed
        ),
        class = "display"
      )   # Apply default styling
      
      return(efftable_dt)
      
    },
    
    
    ##**TODO: check**
    search_rsiena_v1 = function(structure_model, 
                                 get_eff_doc=FALSE,
                                 rsiena_phase2_nsub=1,
                                 rsiena_n2start_scale=1, 
                                 iterations=1000, 
                                 run_seed=123, 
                                 digits=3,
                                 plot_save=FALSE,
                                 return_plot=FALSE
                                 ) {
      if(is.null(self$bipartite_matrix))
        stop('bipartite_matrix is not set.')

      # print('BEFORE self$init_multiwave_rsiena_model_from_structure_model_bipartite_matrix')
      
      ##**self$rsiena_data property** SET rsiena_data FROM 1st 2 bipartite_matrix and structure_model
      self$init_multiwave_rsiena_model_from_structure_model_bipartite_matrix(
        structure_model,
        self$bipartite_matrix, 
        self$bipartite_matrix, ## one tie will be randomly toggled because RSiena sim instructions suggest not using identical networks
        self$rsiena_env_seed
      )
      
      print('self$rsiena_data : ')
      print(self$rsiena_data)
      
      ##  INIT effects list in simulation model (in RSiena model)
      self$rsiena_effects <- getEffects(self$rsiena_data)
      
      # stop('DEBUG XXXXXXXXXXXXX')
      ##--2. NETWORK: STRUCTURE EVOLUTION (structure_Model)-----
      ##  2.1. Add effects from model objective function list
      self$add_rsiena_effects(structure_model)
    
      # Effects Documentation
      if(get_eff_doc)
        effectsDocumentation(self$rsiena_effects, type = 'html', display = T)
      ##-----------------------------
      
      
      ## set fix to FALSE for all parameters to be simulated
      self$rsiena_effects$fix <- rep(FALSE, length(self$rsiena_effects$fix) )
      # self$rsiena_model$theta
      
      ##-----------------------------
      ## RSiena Algorithm
      self$rsiena_run_seed <- run_seed
      self$rsiena_algorithm <- sienaAlgorithmCreate(projname=sprintf('%s_%s',self$SIM_NAME,self$TIMESTAMP),
                                                    simOnly = T,
                                                    # nsub = rsiena_phase2_nsub * 1,
                                                    nsub = 0,
                                                    # n2start = rsiena_n2start_scale * 2.52 * (7+sum(self$rsiena_effects$include)),
                                                    n3 = iterations,
                                                    seed = run_seed)
      
      # Run RSiena simulation
      self$rsiena_model <- siena07(self$rsiena_algorithm,
                                   data = self$rsiena_data,
                                   effects = self$rsiena_effects,
                                   batch = TRUE,
                                   returnDeps = TRUE,
                                   returnChains = TRUE,
                                   returnThetas = TRUE,
                                   returnDataFrame = TRUE, ##**TODO** CHECK
                                   returnLoglik = TRUE     ##**TODO** CHECK
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
      # plt <- self$plot_bipartite_system_from_mat(self$bipartite_matrix, round(iterations), 
      #                                           plot_save = plot_save, 
      #                                           return_plot=return_plot, 
      #                                           normalize_degree=FALSE)
      
      # if (return_plot)
      #   return(plt)
      # ##
      # search_model <- SaomNkRSienaBiEnv_search_rsiena_model$new(structure_model)
      # 
      # return(search_model)
    },
    
    
    
    
    
    # 
    # 
    # ## 2.b. Component Payoffs vector
    # set.seed(12345)
    # component_payoffs <-  runif(environ_params$N, min = 0, max = 1)
    # ## 2. Strategies sets the objective function as a linear combination of network stats across DVs
    # #
    # actor_strats_list <- lapply(strategies, function(strat) rep(strat,  environ_params$M/length(strat)) )
    # 
    #
    # structure_model <- list(
    #   dv_bipartite = list(
    #     name = 'self$bipartite_rsienaDV',
    #     effects = list( ##**STRUCTURAL EFFECTS -- dyadic/network endogeneity sources**
    #       list(effect='density', parameter= 0, fix=T, dv_name=DV_NAME), ##interaction1 = NULL
    #       list(effect='inPop',   parameter= 0, fix=T, dv_name=DV_NAME), #interaction1 = NUL
    #       list(effect='outAct',  parameter= 0, fix=T, dv_name=DV_NAME)
    #     ),
    #     ## COVARIATE EFFECTS
    #     coCovars = list( 
    #       ##** COMPONENTS : MONADIC CONSTANT COVARIATE EFFECTS **##
    #       list(effect='altX',   parameter= 1, dv_name=DV_NAME, fix=T,
    #            interaction1='self$component_1_coCovar', x = component_payoffs 
    #       ),
    #       ##** STRATEGIES : MONADIC CONSTANT COVARIATE EFFECTS **##
    #       list(effect='egoX',   parameter= .1,  dv_name=DV_NAME, fix=T,
    #            interaction1='self$strat_1_coCovar',   x = actor_strats_list[[1]] 
    #       ), #interaction1 = NULL
    #       list(effect='inPopX', parameter= .5,  dv_name=DV_NAME, fix=T,
    #            interaction1='self$strat_2_coCovar',  x = actor_strats_list[[2]] 
    #       )
    #     ),
    #     varCovars = list() ##**MONADIC TIME-VARYING COVARIATE EFFECTS -- DYNAMIC STRATEGY PROGRAMS**
    #   )
    # )
    
    #
    
    
    get_input_from_structure_model = function(structure_model) {
      # actor_strats <- self$get_actor_strategies()
      if( ! 'dv_bipartite' %in% names(structure_model) ) 
        stop('structure_model does not contain dv_bipartite.')
      #
      structeffs <- structure_model$dv_bipartite$effects %>% ldply(as.data.frame) %>% mutate(interaction1=NA, interaction2=NA) # use previous effects
      coCovars   <- structure_model$dv_bipartite$coCovars %>% ldply(function(x) {
        as.data.frame(x[! names(x)%in%c('x')])
      }) %>% mutate(interaction2=NA)
      coDyadCovars   <- structure_model$dv_bipartite$coDyadCovars %>% ldply(function(x) {
        as.data.frame(x[! names(x)%in%c('x')])
      }) %>% mutate(interaction2=NA) 
      # varCovars   <- structure_model$dv_bipartite$varCovars %>% ldply(function(x) {
      #   as.data.frame(x[! names(x)%in%c('x')])
      # })  %>% mutate(interaction2=NA) 
      # varDyadCovars   <- structure_model$dv_bipartite$varDyadCovars %>% ldply(function(x) {
      #   as.data.frame(x[! names(x)%in%c('x')])
      # }) %>% mutate(interaction2=NA) 
      interactions   <- structure_model$dv_bipartite$interactions %>% ldply(function(x) {
        as.data.frame( x )
      }) 
      #
      effects <- structeffs %>%  
        bind_rows(coCovars) %>%  
        bind_rows(coDyadCovars) %>%
        bind_rows(interactions) ## %>%
        # bind_rows(varCovars) %>%
        # bind_rows(varDyadCovars) %>%
      #
      return(effects)
    },
    
    
    ############################################################################
    ##**TODO: check**
    search_rsiena = function(structure_model, 
                             array_bi_net=NULL, ## starting matrices for the simulation
                             theta_matrix=NULL, ## variable theta matrix replaces parameter values in structure_model
                             iterations=1000,   ## replaced by nrow(theta_matrix) if theta_matrix is not null
                             run_seed=123, 
                             process_chain=TRUE,
                             get_eff_doc=FALSE,
                             plot_save=FALSE,
                             return_plot=FALSE,
                             verbose=TRUE,
                             digits=3
    ) {
      
      self$config_structure_model <- structure_model
      
      input_effs <- self$get_input_from_structure_model(structure_model)

      ##--------- I. SET BIPARTITE NETWORK DV ARRAY IF NULL -----------------
      # rsiena_model$sims[[sim_id]][[wave_id]][[bi_env_dv_id]] 
      ## Set up bipartite matrix RSiena dependent variable 
      if (any(is.null(array_bi_net))) 
      {
        if(is.null(self$bipartite_matrix)) 
          stop('bipartite_matrix is not set and no array_bi_net provided.')
        #
        array_bi_net <- array(c(self$bipartite_matrix, self$bipartite_matrix), 
                              dim = c(self$M, self$N, 2))
      }
      # 
      self$bipartite_rsienaDV <- sienaDependent(array_bi_net, 
                                                type='bipartite',
                                                nodeSet =c('ACTORS', 'COMPONENTS'), 
                                                allowOnly = F)
      
      ##----------- II. SET RSIENA DATA AND EFFECTS ---------------------------
      ##  1. RSiena data object
      self$rsiena_data <- self$get_rsiena_data_from_structure_model(structure_model)
      if(verbose) {
        cat('\n\nself$rsiena_data : \n\n')
        print(self$rsiena_data)
      }
      ##  2. Init effects
      self$rsiena_effects <- getEffects(self$rsiena_data)
      ## UPDATE PARAMETERS TO FIXED=FALSE so we can estimate them
      ##  3. Add effects from model objective function list
      self$add_rsiena_effects(structure_model)
      ## set fix to FALSE for all parameters to be simulated
      # self$rsiena_effects$fix <- rep(FALSE, length(self$rsiena_effects$fix) )
      
      
      ##---------- III. SET THETA MATRIX IF NULL  -----------------------------
      if (any(is.null(theta_matrix))) 
      {
        effs <- as.data.frame(self$rsiena_effects[self$rsiena_effects$include, ])
        ## add short name for convenience as effect name
        effs$effect <- effs$shortName
        ## EXCLUDE RATES PARAMETERS NOT INCLUDED IN INPUTE MODEL
        effs <- effs[ effs$effect %in% input_effs$effect , ]
        # effs <- self$get_rsiena_effects_df_from_structure_model(structure_model)
        # names_theta_in <- effs$effectName[effs$include][params_norates]
        theta_in        <- effs$parm
        names(theta_in) <- effs$effect
        # interaction1s   <- effs$interaction1
        # names(interaction1s) <- effs$effect
        nthetas <- length(theta_in)
        ## Number of decision chain steps to simulate
        theta_matrix <- matrix(NA, nrow = iterations, ncol = nthetas)
        for (j in 1:length(theta_in)) {
          print(names(theta_in)[j])
          theta_matrix[, j] <- rep(theta_in[j], iterations)
        }  
        colnames(theta_matrix) <- names(theta_in)
        #print(theta_matrix)
      }
      #
      if(verbose) {
        cat('\n\n theta_matrix : \n\n')
        print(theta_matrix)
      }
      
      
      
      ##---------- IV. RUN SIMULATION  -----------------------------
      ##  4. RSiena Algorithm
      self$rsiena_run_seed <- run_seed
      self$rsiena_algorithm <- sienaAlgorithmCreate(projname=file.path(self$DIR_OUTPUT, sprintf('%s_%s',self$SIM_NAME,self$TIMESTAMP)),
                                                    simOnly = T,
                                                    # nsub = rsiena_phase2_nsub * 1,
                                                    nsub = 0,
                                                    # n2start = rsiena_n2start_scale * 2.52 * (7+sum(self$rsiena_effects$include)),
                                                    n3 = nrow(theta_matrix),
                                                    seed = run_seed)
      ## 5. Run RSiena simulation
      self$rsiena_model <- siena07(self$rsiena_algorithm,
                                   data = self$rsiena_data, 
                                   effects = self$rsiena_effects,
                                   thetaValues = theta_matrix,
                                   batch = TRUE,
                                   returnDeps = TRUE, 
                                   returnChains = TRUE,
                                   returnThetas = TRUE,
                                   returnDataFrame = TRUE, ##**TODO** CHECK
                                   returnLoglik = TRUE   ##**TODO** CHECK
      )   # returnChains = returnChains
      # 6. Summary of fitted model
      mod_summary <- summary( self$rsiena_model )
      if(verbose &  !is.null(mod_summary))
        print(mod_summary)
      
      
      ##----------- IV. POST-PROCESSING OF SIMULATION RESULTS  ------------------
      if(process_chain)
      {
        ## PRocess the decision ministep chain 
        ##  (C++ output to R data.frame: reindexing C++'s 0-index Actor IDs to R's 1-index)
        self$search_rsiena_process_ministep_chain(verbose)
        ## Compute actor and component statistics 
        self$search_rsiena_process_stats()  
        ##
        new_bi_env_matrix <- self$bi_env_arr[,, dim(self$bi_env_arr)[3] ]
        self$set_system_from_bipartite_matrix( new_bi_env_matrix )
      }
      
      ## Show RSiena model screenreg after chain stats (show users what they expect to see at bottom)
      if(verbose)
        print(screenreg(list(self$rsiena_model), single.row = T, digits = digits))
      
    },
    ############################################################################
    
    
    ##
    ## The output is a 3D array with:
    ##  The first dimension NK[a, :, :] is the number of an iteration (from 0 to n_landscapes)
    ##  The second dimension NK[:, b, :] is the location on a given landscape, where b ranges from 0 to 2^N
    ##  The third dimension NK[:, :, c] captures:
    ##   - the first N columns are for the combinations of N decision variables DV
    ##   - the second N columns are for the contribution values of each DV
    ##   - the next value (index=2*N) is for the total fit (avg of N contributions)
    ##   - the last one (index=2*N+1 is to find out whether it is the local peak (0 or 1)
    compute_fitness_landscape = function(n_landscapes=30, 
                                         component_coCovar=1, ## if not exists, then just use random
                                         normalize_int_mat=TRUE,
                                         project_int_mat=FALSE,
                                         component_value_sd=0.1,
                                         verbose=TRUE) {
      # SYSTEM INPUT
      N <- self$N
      # i <- n_landscapes
      
      has_coDyadCovar <- !is.null(self$component_1_coDyadCovar) & all(dim(self$component_1_coDyadCovar) == self$N)
      
      ##**TODO - Use projected search matrix or exogenous covariate matrix? **
      ## use projected matrix
      if (project_int_mat || !has_coDyadCovar) {
        Int_matrix <- self$search_matrix
      } else if ( has_coDyadCovar ) {
        Int_matrix <- as.matrix( self$component_1_coDyadCovar )
      } else {
        stop('project_int_mat=FALSE but no component_x_coDyadCovar variable included in RSiena data.')
      }
      #
      if (normalize_int_mat) {
        Int_matrix <- sapply(1:nrow(Int_matrix), function(i){
          if(max(Int_matrix[i, ])==0) 
            return(Int_matrix[i, ])
          return( Int_matrix[i, ] / max(Int_matrix) )
        })
      }
      
      # K <- self$K_C_df %>% 
      #   filter(chain_step_id==max(chain_step_id)) %>% 
      #   summarize(mean=mean(value - 1)) %>%
      #   pull(mean)
      
      #________internal_functions________________
      powerkey <- function(N) {
        # Compute the powers of 2 for index location
        return(2^((N-1):0))
      }
      nkland <- function(N, component_cov=NULL, component_value_sd=0.1) {
        ## Fully random landscape
        if(!self$exists(component_cov)) {
          runif_mat <- matrix(runif(2^N * N), nrow=2^N, ncol=N)
          return(runif_mat)
        }
        ##**Use component covariates in fitness computation**
        vals <- self[[sprintf('component_%s_coCovar', component_coCovar)]]
        ## Gaussian noise centered at the component value ~ U(0,1), with sd=component_value_sd
        rnorm_mat <- sapply(1:N, function(i) rnorm(n = 2^N, mean = vals[i], sd=component_value_sd))
        return( rnorm_mat )
      }
      calc_fit <- function(N, NK_land, inter_m, Current_position, Power_key) {
        # Compute fitness contributions for each decision variable
        Fit_vector <- numeric(N)
        for (ad1 in 1:N) {
          Fit_vector[ad1] <- NK_land[sum(Current_position * inter_m[ad1, ] * Power_key) + 1, ad1]
        }
        return(Fit_vector)
      }
      comb_and_values <- function(N, NK_land, Power_key, inter_m) {
        # Compute fitness for all combinations
        Comb_and_value <- matrix(0, nrow=2^N, ncol=2*N+2)
        combinations <- expand.grid(replicate(N, 0:1, simplify = FALSE))
        
        for (c1 in 1:(2^N)) {
          Combination1 <- as.integer(combinations[c1,])
          fit_1 <- calc_fit(N, NK_land, inter_m, Combination1, Power_key)
          
          Comb_and_value[c1, 1:N] <- Combination1
          Comb_and_value[c1, (N+1):(2*N)] <- fit_1
          Comb_and_value[c1, (2*N+1)] <- mean(fit_1)
        }
        
        # Check for local peaks
        for (c3 in 1:(2^N)) {
          loc_p <- 1  # Assume it is a peak
          for (c4 in 1:N) {
            new_comb <- Comb_and_value[c3, 1:N]
            new_comb[c4] <- abs(new_comb[c4] - 1)
            index <- sum(new_comb * Power_key) + 1
            if (Comb_and_value[c3, 2*N+1] < Comb_and_value[index, 2*N+1]) {
              loc_p <- 0  # Not a peak
            }
          }
          Comb_and_value[c3, 2*N+2] <- loc_p
        }
        
        return(Comb_and_value)
      }
      ##___________________________
      
      # *** GENERAL VARIABLES AND OBJECTS ***
      Power_key <- powerkey(N)
      NK <- array(0, dim=c(n_landscapes, 2^N, 2*N+2))
      
      start_time <- Sys.time()
      for (i_1 in 1:n_landscapes) {
        nkland_i1 <- nkland(N, component_coCovar, component_value_sd)
        NK[i_1,,] <- comb_and_values(N, nkland_i1, Power_key, Int_matrix)
      }
      end_time <- Sys.time()
      run_time <-  end_time - start_time
      
      self$fitness_landscape <- NK
      
      # Print average number of peaks per landscape
      # average_peaks <- sum(NK[, , (2*N+2)]) / n_landscapes
      # cat("\nThe average number of peaks per landscape is:", average_peaks, "\n")
      npeaks_vec <- sapply(1:n_landscapes, function(i) sum(NK[i,,(2*N+2)]) )
      peaks <- psych::describe(npeaks_vec)
      cat("\nThe mean number of peaks per landscape is:", peaks$mean, "\n")
      cat("\nThe std. deviation of the number of peaks per landscape is:", peaks$sd, "\n")
      cat("\nThe skewness of the number of peaks per landscape is:", peaks$skew, "\n")
      cat("\nThe kurtosis of the number of peaks per landscape is:", peaks$kurtosis, "\n")
      
      # Print elapsed time
      cat("\nElapsed time:", round(run_time, 2), "sec\n")
      
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
      sim_iterations <- self$rsiena_model$n3 ##**TODO** CHECK
      #
      # self$plot_bipartite_system_from_mat(self$bipartite_matrix, iterations, plot_save = T)
      self$plot_bipartite_system_from_mat(self$bipartite_matrix, sim_iterations, 
                                          plot_save = FALSE, return_plot=TRUE)

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
                                 plot_save = TRUE,
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
                                           rand_seed=123, 
                                           dir_output=NA,
                                           file_output=NA) {
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
      .projdir <- ifelse(is.na(dir_output), getwd(), dir_output)
      .projname <-  ifelse(!is.na(file_output), file_output, as.character(as.numeric(Sys.time())))
      self$rsiena_algorithm <- sienaAlgorithmCreate(projname=file.path( .projdir,  sprintf('%s.log', .projname) ), ## rsiena project log filename
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
      sim_iteration <- ifelse( is.null(sim_iteration), rsiena_model$n3 , sim_iteration)
      new_bi_env_mat <- self$get_bipartite_matrix_from_rsiena_model(rsiena_model, sim_iteration)
      new_bi_env_igraph <- igraph::graph_from_biadjacency_matrix(new_bi_env_mat, directed = F, weighted = T, mode = 'out') ##**TODO: CHECK** all vs. out
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
                       length(rsiena_model$chain), ## default current state is the last simulation in sims list
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
      bi_env_mat_sp <- sparseMatrix(i = el_bi_env[,1],
                                    j = el_bi_env[,2],
                                    x = el_bi_env[,3],
                                    dims = c(MplusN, MplusN))
      ## Undirected network --> only uses 'Upper right' rectangle of full bipartite matrix
      ##   M N
      ## M[0,X], for X in {0,1}
      ## N[0,0]
      new_bi_env_mat <- as.matrix(bi_env_mat_sp)[ 1:self$M, (self$M+1):(MplusN) ]
      #
      return( new_bi_env_mat )
    },
    
    ##
    get_chain_stats_list = function() { ## 'all','statdf','bi_env_arr', 'util', 'util_diff')
      if (is.null(self$chain_stats))
        stop('chain_stats missing; run simulation with returnChains=TRUE in order to compute actor utility')
      # ## actor Utility vector
      # au <- c( mat %*% self$rsiena_model$theta )
      # hist( util )
      struct_effnames <- sapply(self$config_structure_model$dv_bipartite$effects, function(x) x$effect)
      #
      coCovar_effnames  <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x) x$effect)
      varCovar_effnames <- sapply(self$config_structure_model$dv_bipartite$varCovars, function(x) x$effect)
      coDyadCovar_effnames <- sapply(self$config_structure_model$dv_bipartite$coDyadCovars, function(x) x$effect)
      varDyadCovar_effnames <- sapply(self$config_structure_model$dv_bipartite$varDyadCovars, function(x) x$effect)
      ## ALL effect names
      effnames <- unlist(c(struct_effnames, coCovar_effnames, varCovar_effnames, coDyadCovar_effnames, varDyadCovar_effnames))
      
      #
      config_param_vals <- c(
        unlist(sapply(self$config_structure_model$dv_bipartite$effects, function(x) x$parameter)),
        unlist(sapply(self$config_structure_model$dv_bipartite$coCovars, function(x) x$parameter)),
        unlist(sapply(self$config_structure_model$dv_bipartite$varCovars, function(x) x$parameter)),
        unlist(sapply(self$config_structure_model$dv_bipartite$coDyadCovars, function(x) x$parameter)),
        unlist(sapply(self$config_structure_model$dv_bipartite$varDyadCovars, function(x) x$parameter))
      )
      names(config_param_vals) <- effnames
      fixed_params <- c(
        unlist(sapply(self$config_structure_model$dv_bipartite$effects, function(x) x$fix)),
        unlist(sapply(self$config_structure_model$dv_bipartite$coCovars, function(x) x$fix)),
        unlist(sapply(self$config_structure_model$dv_bipartite$varCovars, function(x) x$fix)),
        unlist(sapply(self$config_structure_model$dv_bipartite$coDyadCovars, function(x) x$fix)),
        unlist(sapply(self$config_structure_model$dv_bipartite$varDyadCovars, function(x) x$fix))
      )
      
      ## Use system bipartite matrix as start
      bi_env_mat <- self$bipartite_matrix
      
      ## remove chain entries where no change was made (keep if !stability )
      # tiechdf <- self$chain_stats[ !self$chain_stats$stability, ]
      ## Keep all steps (including no-change for forbearance measures)
      tiechdf <- self$chain_stats
      
      ## get matrix timeseries and network statistics timeseries
      # nchains <- length(self$rsiena_model$chain)
      nchains <- nrow(tiechdf)
      
      ## paramters
      theta   <- self$rsiena_model$theta
      ## replace fixed theta param with given param values (instead of zero default when effect is fixed)
      if (any(fixed_params)) 
        theta[ fixed_params ] <- config_param_vals[ fixed_params ]
      #
      ntheta <- length(theta)
      
      
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
      for (i in 1:nchains) {
        mstep <- tiechdf[i,]
        # cat(sprintf('\n %.2f%s i=%s, j=%s', 100*i/nchains,'%', mstep$id_from,  mstep$id_to))
        if(i %% 100 == 0) cat(sprintf('\n %.2f%s', 100*i/nchains,'%'))
        ## update bipartite environment matrix for one step (self$toggle one dyad)
        # bi_env_mat <- self$toggleBiMat(bi_env_mat, mstep$id_from,  (mstep$id_to - self$M)  )
        
        if ( ! mstep$stability ) {
          bi_env_mat <- self$toggleBiMat(bi_env_mat, mstep$id_from,  mstep$id_to)
        }
        #
        # bi_env_mat
        #
        # Get New Statistics
        # statsobj <- self$get_struct_mod_net_stats_list_from_bi_mat( bi_env_mat )
        statmat <- self$get_struct_mod_stats_mat_from_bi_mat( bi_env_mat ) 
        #
        step_statgrid <- expand.grid(chain_step_id=i, actor_id=1:self$M, effect_id=1:ntheta)
        step_statgrid$value <- c( statmat )
        step_statgrid$value_contributions <- c(sweep(statmat, 2, theta, "*")) ## multiply each actor row (margin=2) of statistics by the theta vector
        step_statgrid$stability <- mstep$stability
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
        step_utilgrid$stability <- mstep$stability
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
        step_util_diffgrid$stability <- mstep$stability
        util_diff_li[[i]] <- step_util_diffgrid  ##**CHECK list**
        ##
        ######
        ## update utility lag for next period difference (applies to i>1)
        util_lag <- util
        ######
        
        #
        K_A_grid <- expand.grid(chain_step_id=i, actor_id=1:self$M)
        # K_A_grid$value <- rowSums( bi_env_mat %*% t(bi_env_mat) , na.rm=T) ## WEIGHTED
        K_A_grid$value <- apply( bi_env_mat %*% t(bi_env_mat), 1, function(x)sum(x>0, na.rm=T) )
        K_A_grid$stability <- mstep$stability
        K_A_li[[i]] <- K_A_grid
        #
        K_B1_grid <- expand.grid(chain_step_id=i, actor_id=1:self$M)
        # K_B1_grid$value <- rowSums( bi_env_mat, na.rm=T)  ## WEIGHTED
        K_B1_grid$value <- apply( bi_env_mat, 1, function(x)sum(x>0, na.rm=T))
        K_B1_grid$stability <- mstep$stability
        K_B1_li[[i]] <- K_B1_grid
        #
        K_B2_grid <- expand.grid(chain_step_id=i, component_id=1:self$N)
        # K_B2_grid$value <- colSums( bi_env_mat, na.rm=T)  ## WEIGHTED
        K_B2_grid$value <- apply( bi_env_mat, 2, function(x)sum(x>0, na.rm=T))
        K_B2_grid$stability <- mstep$stability
        K_B2_li[[i]] <- K_B2_grid
        #
        K_C_grid <- expand.grid(chain_step_id=i, component_id=1:self$N)
        # K_C_grid$value <- colSums( t(bi_env_mat) %*% bi_env_mat  , na.rm=T)  ## WEIGHTED
        K_C_grid$value <- apply( t(bi_env_mat) %*% bi_env_mat, 2, function(x)sum(x>0, na.rm=T))  ## WEIGHTED
        K_C_grid$stability <- mstep$stability
        K_C_li[[i]] <- K_C_grid
      }
      actor_strats <- self$get_actor_strategies()
      ##-----------------
      ## Actor Network Statistics long dataframe
      stats_df <- data.table::rbindlist( stats_li )
      stats_df$actor_id <- as.factor(stats_df$actor_id)
      stats_df$strategy <- as.factor( actor_strats[ stats_df$actor_id ] )
      stats_df$effect_name <- as.factor(effnames[ stats_df$effect_id ])
      ## Actor Utility  long dataframe
      util_df <- data.table::rbindlist( util_li ) 
      util_df$actor_id <- as.factor(util_df$actor_id)
      util_df$strategy <- as.factor( actor_strats[ util_df$actor_id ] )
      ## Actor Utility Difference long dataframe
      util_diff_df <- data.table::rbindlist( util_diff_li ) 
      util_diff_df$actor_id <- as.factor(util_diff_df$actor_id)
      util_diff_df$strategy <- as.factor( actor_strats[ util_diff_df$actor_id ] )
      ##---
      
      K_A_df <- data.table::rbindlist( K_A_li ) 
      K_A_df$actor_id <- as.factor(K_A_df$actor_id)
      K_A_df$strategy <- as.factor( actor_strats[ K_A_df$actor_id ] )
      #
      K_B1_df <- data.table::rbindlist( K_B1_li ) 
      K_B1_df$actor_id <- as.factor(K_B1_df$actor_id)
      K_B1_df$strategy <- as.factor( actor_strats[ K_B1_df$actor_id ] )
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
    # from ego to alter is self$toggled (1 -> 0; 0 -> 1).
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
    search_rsiena_process_ministep_chain = function(verbose=TRUE) {
      
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
      chainDat$chain_step_id <- sapply(1:nrow(chainDat), function(.id) ifelse(chainDat$stability[.id], NA, .id ) )
      
      chainDat <- chainDat %>% tidyr::fill(chain_step_id, .direction = "down") ## a no-change takes previous chain_step_id
      
      ## set to simulation self
      self$chain_stats <- chainDat
        
      if (verbose) {
        cat('\n\n\nSimulated Decision Chain Summary:\n\n')
        print(summary(self$chain_stats))
        print(dim(self$chain_stats))
        print(self$chain_stats)
      }

    },
    
    #
    search_rsiena_process_stats = function() {
      ##
      # actor_stats_df = NULL,
      # actor_util_df = NULL,
      # actor_util_diff_df = NULL,
      # ##
      chain_stats_list <- self$get_chain_stats_list()
      # if (is.null(chain_stats_list$bi_env_arr)) 
      #   stop('search_rsiena_process_stats() bi_env_array is NULL')
      ##
      self$bi_env_arr <- chain_stats_list$bi_env_arr  ## chain array of bipartite matrix from chain of decision steps
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
    
    
    # #
    # process_sims = function() {
    #  
    #   #   'K_A'=list(), ##**TODO** Actor social space [degree = central position in social network]
    #   #   'K_B'=list(), #  Bipartite space  [ degree = resource/affiliation density ]
    #   #   'K_C'=list()  #  Component space  [ degree = interdependence/structuration : epistatic interactions ]
    #   #   ), 
    #   
    #   self$sims_stats <- NULL
    # },
    
    
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
          jaccardlist[[sprintf('jac%d-%d', i-1, i)]] <- self$get_jaccard_index(m0 = outlist[[ (i-1) ]], m1 = bi_env_mat_new )
        }
        
        new_bi_g <- igraph::graph_from_biadjacency_matrix(bi_env_mat_new, 
                                                          directed = F, mode = 'all', 
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
    search_rsiena_multiwave_plot = function(type=c(), 
                                            rolling_window = 10, 
                                            actor_ids=c(),
                                            component_ids=c(),
                                            wave_ids=c(),
                                            thin_factor=1, 
                                            thin_wave_factor=1,
                                            smooth_method='loess',
                                            show_utility_points=TRUE,
                                            show_strategy_means=TRUE,
                                            append_plot=FALSE,
                                            histogram_position='identity',
                                            scale_utility=TRUE,
                                            return_plot=TRUE,
                                            plot_file=NA, plot_dir=NA,
                                            loess_span=0.4
    ) {
      # if (all(is.na(type)) |  'ministep_count' %in% type)
      #   return(self$search_rsiena_plot_actor_ministep_count(rolling_window))
      # if (all(is.na(type)) | type == 'utility')
      #   return(self$search_rsiena_plot_actor_utility(rolling_window))
      plist <- list()
      # if (all(is.na(type)) |  'utility' %in% type)
      #   plist[['utility']] <- self$search_rsiena_plot_actor_utility(actor_ids, thin_factor, smooth_method, show_utility_points, return_plot=T)
      # 
      if (length(type)==0 |  'K_4panel' %in% type)
        plist[['K_4panel']] <- self$search_rsiena_multiwave_plot_K_4panel(actor_ids, component_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, return_plot=T, plot_file=plot_file )
      
      if (length(type)==0 |  'K_A_strategy_summary' %in% type)
        plist[['K_A_strategy_summary']] <- self$search_rsiena_multiwave_plot_K_A_strategy_summary(actor_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, return_plot=T, plot_file=plot_file )
      
      if (length(type)==0 |  'K_B1_strategy_summary' %in% type)
        plist[['K_B1_strategy_summary']] <- self$search_rsiena_multiwave_plot_K_B1_strategy_summary(actor_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, return_plot=T, plot_file=plot_file )
      
      if (length(type)==0 |  'K_B2_strategy_summary' %in% type)
        plist[['K_B2_strategy_summary']] <- self$search_rsiena_multiwave_plot_K_B2_strategy_summary(component_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, return_plot=T, plot_file=plot_file )
      
      if (length(type)==0 |  'K_C_strategy_summary' %in% type)
        plist[['K_C_strategy_summary']] <- self$search_rsiena_multiwave_plot_K_C_strategy_summary(component_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, return_plot=T, plot_file=plot_file )
      
      
      if (length(type)==0 |  'utility_strategy_summary' %in% type)
        plist[['utility_strategy_summary']] <- self$search_rsiena_multiwave_plot_actor_utility_strategy_summary(actor_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, scale_utility, return_plot=T, plot_file=plot_file, loess_span=loess_span )
      
      if (length(type)==0 |  'utility_by_strategy' %in% type)
        plist[['utility_by_strategy']] <- self$search_rsiena_multiwave_plot_actor_utility_by_strategy(actor_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, return_plot=T, plot_file=plot_file )
      # 
      # if (length(type)==0 |  'utility_density' %in% type)
      #   plist[['utility_density']] <- self$search_rsiena_plot_actor_utility_density(return_plot=T, plot_file=plot_file )
      # 
      if (length(type)==0 |  'utility_density_by_strategy' %in% type)
        plist[['utility_density_by_strategy']] <- self$search_rsiena_multiwave_plot_actor_utility_density_by_strategy(thin_wave_factor, return_plot=T, plot_file=plot_file )
      # 
      # if (length(type)==0 |  'utility_histogram_by_strategy' %in% type)
      #   plist[['utility_histogram_by_strategy']] <- self$search_rsiena_plot_actor_utility_histogram_by_strategy(histogram_position, return_plot=T, plot_file=plot_file )
      # 
      # if (length(type)==0 |  'stability' %in% type)
      #   plist[['stability']] <- self$search_rsiena_plot_stability() ##**TODO** Fix return plot
      #
      if (length(type)==0 |  'utility_ridge_density_by_strategy' %in% type)
        plist[['utility_ridge_density_by_strategy']] <- self$search_rsiena_multiwave_plot_utility_ridge_density_by_strategy(actor_ids, wave_ids, thin_factor, thin_wave_factor, show_utility_points, show_strategy_means, return_plot=T, plot_file=plot_file )
      #
      #  SET plots 
      self$multiwave_plots <- if(append_plot) { append(self$multiwave_plots, plist) } else { plist }
      
      # stop(sprintf('Plot type not implemented: %s', type))
      if(return_plot)
        return(plist)
    },
    
    
    
    #
    search_rsiena_multiwave_plot_utility_ridge_density_by_strategy = function(actor_ids=c(),
                                                                              wave_ids=c(),
                                                                              thin_factor=1,
                                                                              thin_wave_factor=1,
                                                                              show_utility_points=T,
                                                                              show_strategy_means=T,
                                                                              scale_utility=TRUE,
                                                                              return_plot=TRUE,
                                                                              plot_file=NA,
                                                                              plot_dir=NA,
                                                                              plot_periods=4) {
      #
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
      ## transition_pds <- plot_periods - 1
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
          stabilization_summary_period = ifelse(actor_component_period <= (plot_periods - 1), 
                                                actor_component_period, 
                                                sprintf('%s+\n(%s-%s)',plot_periods,plot_periods,max(actor_component_period)))
        )
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
      # 
      group_dens_means <- dat_dens_rigde %>% ungroup() %>% group_by(strategy) %>% 
        dplyr::summarize(n=n(),mean=mean(utility,na.rm=T))
      ##---------------------
      ## Start Plot
      plt.dr <- ggplot(dat_dens_rigde, aes(y = PeriodFct, x = utility, color=strategy, fill=strategy)) +
        stat_density_ridges(aes(point_color = strategy, point_fill = strategy, point_shape = strategy),
                            quantile_lines = TRUE, alpha = .3, rel_min_height = density_ridges_rel_min_height,
                            point_size=.4,
                            jittered_points = T, 
                            position = position_raincloud(adjust_vlines = FALSE, ygap = -.1, height = .15),# "raincloud",
                            # position = position_points_sina(rel_min = 0.1, rel_max = 0.9, seed = NULL),
                            quantiles = c(0.5), linewidth=.75 ) +
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
        coord_cartesian(clip = "off") +
        theme_ridges(grid = T, center=T) + 
        theme(legend.position = 'bottom')
      #
      if(show_strategy_means) {
        plt.dr <- plt.dr +  
          geom_vline(data = group_dens_means,  aes(xintercept = mean, color=strategy), 
                     linetype=3, linewidth=1.3) +
          geom_text(data = group_dens_means,
                    aes(x = mean, y = 1+length(unique(dat_dens_rigde$stabilization_summary_period)), label = round(mean, 2), color=strategy),
                    inherit.aes = FALSE, size = 5, nudge_x=-.1, nudge_y=.35  ) 
      }

      #
      if(return_plot)
        return(plt.dr)
    },
    
    
    
    #
    search_rsiena_multiwave_plot_K_4panel = function(actor_ids=c(), 
                                                     component_ids=c(), 
                                                     wave_ids=c(), 
                                                     thin_factor=1, 
                                                     thin_wave_factor=1, 
                                                     smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                     show_utility_points=T, 
                                                     return_plot=TRUE,
                                                     plot_file=NA, plot_dir=NA
                                                     ) {
      K_A  <- self$search_rsiena_multiwave_plot_K_A_strategy_summary(actor_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, show_legend=T, show_title=F, return_plot=T)
      K_B1 <- self$search_rsiena_multiwave_plot_K_B1_strategy_summary(actor_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, show_legend=F, show_title=F, return_plot=T)
      K_B2 <- self$search_rsiena_multiwave_plot_K_B2_strategy_summary(component_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, show_legend=F, show_title=F, return_plot=T)
      K_C  <- self$search_rsiena_multiwave_plot_K_C_strategy_summary(component_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, show_legend=T, show_title=F, return_plot=T)
      # K_B1_leg <- self$search_rsiena_multiwave_plot_K_B1_strategy_summary(actor_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, show_legend=T, show_title=F, return_plot=T)
      # K_B2_leg <- self$search_rsiena_multiwave_plot_K_B2_strategy_summary(component_ids, wave_ids, thin_factor, thin_wave_factor, smooth_method, show_utility_points, show_legend=T, show_title=F, return_plot=T)
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
      
      maintitle <- sprintf('Environment: Actors (M) = %s, Components (N) = %s, Init.Prob. = %.2f\nActor Strategy:  %s\nComponent Payoff:  %s\nStructure:  %s', 
                           self$M, self$N, self$BI_PROB,
                           paste( paste(paste(strateffs[stratDV_ids], stratparams[stratDV_ids], sep='= '), stratfixs[stratDV_ids], sep='' ), collapse = ';  '),
                           paste( paste(paste(strateffs[componentDV_ids], stratparams[componentDV_ids], sep='= '), stratfixs[componentDV_ids], sep=''), collapse = ';  '),
                           paste( paste(paste(structeffs, structparams, sep='= '), structfixs, sep=''), collapse = ';  ')
      )
      
      combined_plot_notitle <- ggarrange(
        K_B1, K_B2, 
        K_A, K_C,
        nrow = 2, ncol = 2, 
        # widths = c(7, 7), # Adjust column widths
        # heights = c(7, 7),
        common.legend = TRUE, # Share a common legend if needed
        legend = "bottom"#,     # Place legend at the bottom
        # align = 'v'
      ) 
      
      combined_plot <- annotate_figure(
        combined_plot_notitle,
        top = text_grob(maintitle, color = "black", size = 14) ## face = "bold", 
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
      
      if(!is.na(plot_file))
        ggsave(file = file.path(ifelse(is.na(plot_dir),getwd(),plot_dir), sprintf("%s_%s.png", self$config_environ_params$name, plot_file)), 
               combined_plot, 
               width = 10, height = 8, units = 'in', dpi = 600)
      
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
                                                                  return_plot=TRUE,
                                                                  plot_file=NA, plot_dir=NA
                                                                 ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar are not  set.")
      if ( attr(self$component_1_coCovar, 'nodeSet') != 'COMPONENTS' )
        stop("Component payoff values in self$component_1_coCovar are not set.")
      # actor_strat <- as.factor( self$strat_coCovar )
      # component_types <- as.factor( self$component_coCovar )
      # component_types <- as.factor( sample(0:1, self$N, replace = T) )
      range_midpoint <- min(self$component_1_coCovar, na.rm=T) + ( abs(diff(range(self$component_1_coCovar, na.rm = T))) / 2 )
      component_types <- as.factor( ifelse(self$component_1_coCovar > range_midpoint, 'High', 'Low') )
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
      #
      density_rng <- range(dat$value, na.rm=T)
      density_absdiff_scale <- abs(diff(density_rng)) * 0.15
      y_lim <- c(density_rng[1] - density_absdiff_scale,  density_rng[2] + density_absdiff_scale)
      #
      point_size <- 10 / log( nstep )
      point_alpha <- min( 1,  .5/log10( nstep ) )
      #
      if(length(component_ids))
        dat <- dat %>% filter(component_id %in% component_ids)
      if(length(wave_ids))
        dat <- dat %>% filter(wave_id %in% wave_ids)
      dat_wave_means <- dat %>% group_by(wave_id) %>% 
        dplyr::summarize(mean=mean(value, na.rm=T))
      plt <- ggplot(dat, aes(x=chain_step_id, y=value)) + 
        geom_hline(data=dat_wave_means, aes(yintercept=mean), linetype=3, col='black' ) +
        facet_grid(wave_id ~ .) 
      if(show_utility_points)
        plt <- plt + geom_point(aes(color=component_id), alpha=point_alpha, shape=1, size=point_size, show.legend = F)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(self$exists(smooth_method))
        plt <- plt + geom_smooth(aes(color=component_id, fill=component_id), method = smooth_method, linewidth=.5, alpha=.09, show.legend = F, se=F)
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
      stratmeans <- dat %>% group_by(component_id, wave_id) %>% 
        dplyr::summarize(n=n(), mean=mean(value, na.rm=T), sd=sd(value, na.rm=T))
      ## Actor density fact plots comparing H1 to H2 utility distribution
      plt2 <- ggplot(dat, aes(x=value, color=component_id, fill=component_id)) + ##linetype=chain_half
        geom_density(alpha=.01, linewidth=.5, show.legend = F)  +
        # geom_histogram(alpha=.1, position = 'dodge') +
        geom_vline(data = stratmeans, aes(xintercept = mean, color=component_id), linetype=2, linewidth=.5, show.legend = F) +
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
        plt2 <- plt2 + ggtitle('\n\n\n')
      
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
      if(!is.na(plot_file))
        ggsave(file = file.path(ifelse(is.na(plot_dir),getwd(),plot_dir), sprintf("%s_%s.png", self$config_environ_params$name, plot_file)), 
               combined_plot, 
               width = 10, height = 8, units = 'in', dpi = 600)
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
                                                                  return_plot=TRUE,
                                                                  plot_file=NA, plot_dir=NA
    ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar are not set.")
      if ( attr(self$component_1_coCovar, 'nodeSet') != 'COMPONENTS' )
        stop("Component payoff values in self$component_1_coCovar are not set.")
      # actor_strat <- as.factor( self$strat_coCovar )
      # component_types <- as.factor( self$component_type_coCovar )
      # component_types <- as.factor( rep(1, self$N ) )
      range_midpoint <- min(self$component_1_coCovar, na.rm=T) + ( abs(diff(range(self$component_1_coCovar, na.rm = T))) / 2 )
      component_types <- as.factor( ifelse(self$component_1_coCovar > range_midpoint, 'High', 'Low') )
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
      density_rng <- range(dat$value, na.rm=T)
      density_absdiff_scale <- abs(diff(density_rng)) * 0.15
      y_lim <- c(density_rng[1] - density_absdiff_scale,  density_rng[2] + density_absdiff_scale)
      #
      point_size <- 10 / log( nstep )
      point_alpha <- min( 1,  .5/log10( nstep ) )
      #
      if(length(component_ids))
        dat <- dat %>% filter(component_id %in% component_ids)
      if(length(wave_ids))
        dat <- dat %>% filter(wave_id %in% wave_ids)
      dat_wave_means <- dat %>% group_by(wave_id) %>% 
        dplyr::summarize(mean=mean(value, na.rm=T))
      plt <- ggplot(dat, aes(x=chain_step_id, y=value)) + 
        geom_hline(data=dat_wave_means, aes(yintercept=mean), linetype=3, col='black' ) +
        facet_grid(wave_id ~ .) 
      if(show_utility_points)
        plt <- plt + geom_point(aes(color=component_id), alpha=point_alpha, shape=1, size=point_size, show.legend = F)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(self$exists(smooth_method))
        plt <- plt + geom_smooth(aes(color=component_id, fill=component_id), method = smooth_method, linewidth=.5, alpha=.09, se=F, show.legend = F)
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
      stratmeans <- dat %>% group_by(component_id, wave_id) %>% 
        dplyr::summarize(n=n(), mean=mean(value, na.rm=T), sd=sd(value, na.rm=T))
      ## Actor density fact plots comparing H1 to H2 utility distribution
      plt2 <- ggplot(dat, aes(x=value, color=component_id, fill=component_id)) + ##linetype=chain_half
        geom_density(alpha=.01, linewidth=.5, show.legend = F)  +
        # geom_histogram(alpha=.1, position = 'dodge') +
        geom_vline(data = stratmeans, aes(xintercept = mean, color=component_id), linetype=2, linewidth=.5, show.legend = F) +
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
        plt2 <- plt2 + ggtitle('\n\n\n')
      
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
      if(!is.na(plot_file))
        ggsave(file = file.path(ifelse(is.na(plot_dir),getwd(),plot_dir), sprintf("%s_%s.png", self$config_environ_params$name, plot_file)), 
               combined_plot, 
               width = 10, height = 8, units = 'in', dpi = 600)
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
                                                                 return_plot=TRUE,
                                                                 plot_file=NA, plot_dir=NA
    ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <-  as.factor( self$get_actor_strategies() )
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
      #
      density_rng <- range(dat$value, na.rm=T)
      density_absdiff_scale <- abs(diff(density_rng)) * 0.15
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
        dplyr::summarize(mean=mean(value, na.rm=T))
      plt <- ggplot(dat, aes(x=chain_step_id, y=value)) + 
        geom_hline(data=dat_wave_means, aes(yintercept=mean), linetype=3, col='black' ) +
        facet_grid(wave_id ~ .) 
      if(show_utility_points)
        plt <- plt + geom_point(aes(color=strategy), alpha=point_alpha, shape=1, size=point_size)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(self$exists(smooth_method))
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
        dplyr::summarize(n=n(), mean=mean(value, na.rm=T), sd=sd(value, na.rm=T))
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
        plt2 <- plt2 + ggtitle('\n\n\n')
      
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
      if(!is.na(plot_file))
        ggsave(file = file.path(ifelse(is.na(plot_dir),getwd(),plot_dir), sprintf("%s_%s.png", self$config_environ_params$name, plot_file)), 
               combined_plot, 
               width = 10, height = 8, units = 'in', dpi = 600)
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
                                                                 return_plot=TRUE,
                                                                 plot_file=NA, plot_dir=NA
    ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <-  as.factor( self$get_actor_strategies() )
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
      density_absdiff_scale <- abs(diff(density_rng)) * 0.15
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
        dplyr::summarize(mean=mean(value, na.rm=T))
      plt <- ggplot(dat, aes(x=chain_step_id, y=value)) + 
        geom_hline(data=dat_wave_means, aes(yintercept=mean), linetype=3, col='black' ) +
        facet_grid(wave_id ~ .) 
      if(show_utility_points)
        plt <- plt + geom_point(aes(color=strategy), alpha=point_alpha, shape=1, size=point_size)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(self$exists(smooth_method))
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
        dplyr::summarize(n=n(), mean=mean(value, na.rm=T), sd=sd(value, na.rm=T))
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
        plt2 <- plt2 + ggtitle('\n\n\n')
      
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
      if(!is.na(plot_file))
        ggsave(file = file.path(ifelse(is.na(plot_dir),getwd(),plot_dir), sprintf("%s_%s.png", self$config_environ_params$name, plot_file)), 
               combined_plot, 
               width = 10, height = 8, units = 'in', dpi = 600)
      #
      if(return_plot)
        return(combined_plot)
    },
    
    get_actor_strategies = function() {
      coveffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
      covDvTypes <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$interaction1)
      stratDV_ids     <- grep('self\\$strat_\\d{1,2}_coCovar', covDvTypes ) ## ex: "self$strat_1_coCovar" 
      strat_effs <- coveffs[ stratDV_ids ]
      nstrat <- length(strat_effs)
      if (nstrat == 1) return(self$strat_1_coCovar)
      if (nstrat == 2) return(paste(self$strat_1_coCovar, 
                                    self$strat_2_coCovar, sep='_'))
      if (nstrat == 3) return(paste(self$strat_1_coCovar, 
                                    self$strat_2_coCovar, 
                                    self$strat_3_coCovar, sep='_'))
      if (nstrat == 4) return(paste(self$strat_1_coCovar, 
                                    self$strat_2_coCovar, 
                                    self$strat_3_coCovar, 
                                    self$strat_4_coCovar, sep='_'))
      if (nstrat == 5) return(paste(self$strat_1_coCovar, 
                                    self$strat_2_coCovar, 
                                    self$strat_3_coCovar, 
                                    self$strat_4_coCovar, 
                                    self$strat_5_coCovar, sep='_'))
      # if (nstrat == 1) return(paste(strat_effs, self$strat_1_coCovar, sep='= '))
      # if (nstrat == 2) return(paste(strat_effs, paste(self$strat_1_coCovar, 
      #                                                 self$strat_2_coCovar, sep='_'), sep='= '))
      # if (nstrat == 3) return(paste(strat_effs, paste(self$strat_1_coCovar, 
      #                                                 self$strat_2_coCovar, 
      #                                                 self$strat_3_coCovar, sep='_'), sep='= '))
      # if (nstrat == 4) return(paste(strat_effs, paste(self$strat_1_coCovar, 
      #                                                 self$strat_2_coCovar, 
      #                                                 self$strat_3_coCovar, 
      #                                                 self$strat_4_coCovar, sep='_'), sep='= '))
      # if (nstrat == 5) return(paste(strat_effs, paste(self$strat_1_coCovar, 
      #                                                 self$strat_2_coCovar, 
      #                                                 self$strat_3_coCovar, 
      #                                                 self$strat_4_coCovar, 
      #                                                 self$strat_5_coCovar, sep='_'), sep='= '))
    },
    
    # actor_ids=c()
    # wave_ids=c()
    # thin_factor=1 
    # thin_wave_factor=1
    # smooth_method='loess'  ##"lm", "glm", "gam", "loess","auto"
    # show_utility_points=T
    # scale_utility=TRUE
    # return_plot=TRUE
    # plot_file=NA
    # plot_dir=NA
    
    
    ##
    search_rsiena_multiwave_plot_actor_utility_strategy_summary = function(actor_ids=c(), 
                                                                           wave_ids=c(),
                                                                            thin_factor=1, 
                                                                            thin_wave_factor=1,
                                                                            smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                                            show_utility_points=T,
                                                                            scale_utility=TRUE,
                                                                            return_plot=TRUE,
                                                                            plot_file=NA, plot_dir=NA,
                                                                           loess_span=0.4
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
      plt <- ggplot(dat, aes(x=chain_step_id, y=utility)) + 
        geom_hline(data=dat_wave_means, aes(yintercept=mean), linetype=3, col='black' ) +
        facet_grid(wave_id ~ .) 
      if(show_utility_points)
        plt <- plt + geom_point(aes(color=strategy), alpha=point_alpha, shape=1, size=point_size)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(self$exists(smooth_method))
        plt <- plt + geom_smooth(aes(color=strategy, fill=strategy, shape=actor_id), 
                                 method = smooth_method, span=loess_span, linewidth=1, alpha=.09)
      #
      # plt <- plt + geom_text(aes(x=chain_step_id, y=utility, linetype=factor(actor_id)), size = 3, vjust = -1) 
      # Population Mean in black
      plt <- plt + geom_smooth(aes(x=chain_step_id, y=utility_mean), method = smooth_method, span=loess_span, 
                               data=dat %>% group_by(chain_step_id)%>%summarize(utility_mean=mean(utility)),
                               color='black', linetype=1, alpha=.1) 
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
    },
    
    ##
    search_rsiena_multiwave_plot_actor_utility_by_strategy = function(actor_ids=c(), 
                                                                      thin_factor=1, 
                                                                      thin_wave_factor=1,
                                                                      smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                                      show_utility_points=TRUE,
                                                                      return_plot=TRUE,
                                                                      plot_file=NA, plot_dir=NA
    ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <- as.factor( self$get_actor_strategies() )
      ## Compare 2 actors utilty
      dat <- self$actor_wave_util %>% 
        filter(chain_step_id %% thin_factor == 0) %>% 
        filter(wave_id %% thin_wave_factor == 0 ) %>% 
        mutate(strategy = actor_strat[ actor_id ] )
      if(length(actor_ids))
        dat <- dat %>% filter(actor_id %in% actor_ids)
      plt <- ggplot(dat, aes(x=chain_step_id, y=utility)) + 
        geom_hline(data=dat%>%group_by(wave_id)%>%dplyr::summarize(mean=mean(utility, na.rm=T)), aes(yintercept=mean), linetype=2, col='black' ) +
        facet_wrap( ~ wave_id)
      if(show_utility_points)
        plt <- plt + geom_point(aes(color=strategy), alpha=.25, shape=1, size=2)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
      if(self$exists(smooth_method))
        plt <- plt + geom_smooth(aes(linetype=actor_id, color=strategy), method = smooth_method, linewidth=1, alpha=.15)
      #
      plt <- plt + theme_bw()
      #
      # Add marginal density plots
      plt <- ggMarginal(plt, type = "density", margins = "y")
      
      #
      if(!is.na(plot_file))
        ggsave(file = file.path(ifelse(is.na(plot_dir),getwd(),plot_dir), sprintf("%s_%s.png", self$config_environ_params$name, plot_file)), 
               plt, 
               width = 10, height = 8, units = 'in', dpi = 600)
      #
      if(return_plot)
        return(plt)
    },
    
    
    ##
    search_rsiena_multiwave_plot_actor_utility_density_by_strategy = function(thin_wave_factor=1, 
                                                                              return_plot=TRUE, 
                                                                              plot_file=NA, 
                                                                              plot_dir=NA) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <-  as.factor( self$get_actor_strategies() )
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
        dplyr::summarize(n=n(), mean=mean(utility, na.rm=T), sd=sd(utility, na.rm=T))
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
      
      if(!is.na(plot_file))
        ggsave(file = file.path(ifelse(is.na(plot_dir),getwd(),plot_dir), sprintf("%s_%s.png", self$config_environ_params$name, plot_file)), 
               plt, 
               width = 10, height = 8, units = 'in', dpi = 600)
      
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
                                  return_plot=TRUE,
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
      if(self$exists(smooth_method))
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
                                                        return_plot=TRUE
                                                        ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <-  as.factor( self$get_actor_strategies() )
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
      if(self$exists(smooth_method))
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
    search_rsiena_plot_actor_utility_density = function(return_plot=TRUE) {
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
    search_rsiena_plot_actor_utility_density_by_strategy = function(return_plot=TRUE) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <-  as.factor( self$get_actor_strategies() )
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
        dplyr::summarize(n=n(), mean=mean(utility, na.rm=T), sd=sd(utility, na.rm=T))
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
                                                                      return_plot=TRUE
                                                                      ) {
      ## actor strategy
      if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
        stop("Actor Strategy self$strat_1_coCovar not set.")
      actor_strat <-  as.factor( self$get_actor_strategies() )
      ## Compare 2 actors utilty
      dat <- self$actor_util_df %>%  
        mutate(
          strategy = actor_strat[ actor_id ] ,
          chain_below_med =  chain_step_id < median(chain_step_id)
        )
      dat$chain_half <- factor(ifelse(dat$chain_below_med, '1st Half', '2nd Half'))
      ##
      stratmeans <- dat %>% group_by(strategy, chain_half) %>% 
        dplyr::summarize(n=n(), mean=mean(utility, na.rm=T), sd=sd(utility, na.rm=T))
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
    
    
    
    # # Visualization using ggplot2 and ggraph: bipartite network, social network, and component interaction matrix heatmap
    # visualize_system_bi_env_rsiena = function(plot_save = FALSE) {
    #   #
    #   RSIENA_ITERATION <- length(self$rsiena_model$sims)
    #   
    #   # Generate labels for components in letter-number sequence
    #   generate_component_labels <- function(n) {
    #     letters <- LETTERS  # Uppercase alphabet letters
    #     if (n <= 26)
    #       return(letters[1:n])
    #     labels <- c()
    #     repeat_count <- ceiling(n / length(letters))
    #     for (i in 1:repeat_count) {
    #       labels <- c(labels, paste0(letters, i))
    #     }
    #     return(labels[1:n])  # Return only as many as needed
    #   }
    #   
    #   component_labels <- generate_component_labels(self$N)
    #   
    #   # 1. Bipartite network plot using ggraph with vertex labels
    #   ig_bipartite <- self$bipartite_igraph
    #   
    #   # Set node attributes for shape and label
    #   # V(ig_bipartite)$type <- V(ig_bipartite)$type
    #   V(ig_bipartite)$shape <- ifelse(V(ig_bipartite)$type, "square", "circle")
    #   V(ig_bipartite)$color <- ifelse(V(ig_bipartite)$type, "lightblue", "darkorange")
    #   V(ig_bipartite)$label <- c( 1:self$M, component_labels)  # Numeric labels for actors
    #   V(ig_bipartite)$node_text_color <- ifelse(V(ig_bipartite)$type, 'black', 'lightblue') 
    #   
    #   bipartite_plot <- ggraph(ig_bipartite, layout = "fr") +
    #     geom_edge_link(color = "gray") +
    #     geom_node_point(aes(shape = shape, color = color), size = 5) +
    #     geom_node_text(aes(label = label, color=node_text_color), vjust = 0.5, hjust = 0.5, size = 3) +  # Center vertex labels over vertices
    #     scale_shape_manual(values = c("circle" = 16, "square" = 15)) +
    #     scale_color_identity() +
    #     labs(title = "[DGP] Bipartite Environment\n(Actors and Components)") +  # Split title into two lines
    #     theme_minimal() +
    #     theme(
    #       legend.position = "none",
    #       panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(),
    #       axis.title.x = element_blank(),
    #       axis.title.y = element_blank(),
    #       axis.text.x = element_blank(),
    #       axis.text.y = element_blank(),
    #       axis.ticks = element_blank()
    #     )
    #   
    #   # 2. Social network plot using ggraph with different statistics for node color and size
    #   ig_social <- self$social_igraph
    #   
    #   # Calculate network statistics
    #   node_size <- degree(ig_social)  # Degree centrality for node size
    #   node_color <- eigen_centrality(ig_social)$vector  # Eigenvector centrality for node color
    #   node_text <- 1:vcount(ig_social)
    #   
    #   # V(ig_social)$label <- ifelse(V(ig_social)$type, 1:self$M, 1:self$M)
    #   
    #   social_plot <- ggraph(ig_social, layout = "fr") +
    #     geom_edge_link(color = "gray") +
    #     geom_node_point(aes(size = node_size, color = node_color)) +
    #     geom_node_text(aes(label = node_text), vjust = 0.5, hjust = 0.5, size = 3, color='white') +  # Center vertex labels over vertices
    #     scale_color_gradient(low = "green", high = "red") +
    #     labs(title = "[Proj1] Actor Social Network\n(Component Overlaps)", 
    #          color = "Eigenvector\nCentrality", size = "Degree\nCentrality") +
    #     theme_minimal() +
    #     theme(
    #       legend.position = "bottom",
    #       legend.box = "vertical",  # Stack legends vertically
    #       panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(),
    #       axis.title.x = element_blank(),
    #       axis.title.y = element_blank(),
    #       axis.text.x = element_blank(),
    #       axis.text.y = element_blank()
    #     ) +
    #     guides(
    #       color = guide_legend(order = 1, nrow = 2),  # Stacks color legend on two rows if needed
    #       size = guide_legend(order = 2, nrow = 2)    # Stacks size legend on two rows if needed
    #     )
    #   
    #   
    #   # 3. Component interaction matrix heatmap
    #   component_matrix <- self$search_matrix
    #   component_df <- melt(component_matrix)
    #   colnames(component_df) <- c("Component1", "Component2", "Interaction")
    #   
    #   ##**DEBUG**
    #   # stop('DEBUG print')
    # 
    #   # Replace component numbers with labels in the heatmap data frame
    #   component_df$Component1 <- factor(component_df$Component1, labels = component_labels)
    #   component_df$Component2 <- factor(component_df$Component2, labels = component_labels)
    #   
    #   heatmap_plot <- ggplot(component_df, aes(x = Component1, y = Component2, fill = Interaction)) +
    #     geom_tile() +
    #     scale_fill_gradient(low = "white", high = "red") +
    #     labs(title = "[Proj2] Component Interaction\n(Epistasis) Heatmap", x = "Component 1", y = "Component 2") +
    #     theme_minimal() +
    #     theme(legend.position = "bottom")
    #   
    #   # Calculate average degree (K) for the social space and component interaction space
    #   avg_degree_social <- mean(degree(ig_social))
    #   avg_degree_component <- mean(degree(graph_from_adjacency_matrix(component_matrix)))
    #   
    #   # Create a main title using sprintf with simulation parameters
    #   main_title <- sprintf("Simulation Parameters:\niter = %d, N = %d, M = %d, BI_PROB = %.2f, K_S = %.2f, K_C = %.2f", 
    #                         RSIENA_ITERATION, self$N, self$M, self$BI_PROB, avg_degree_social, avg_degree_component)
    #   
    #   # Arrange plots with a main title
    #   (plt <- grid.arrange(
    #     social_plot, bipartite_plot, heatmap_plot, ncol = 3,
    #     top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold"))
    #   ))
    #   if(plot_save){
    #     keystring <- sprintf("%s_sim%.0f_iter%.0f_N%d_M%d_BI_PROB_%.2f_K_S_%.2f_K_C_%.2f", 
    #                          self$SIM_NAME, self$TIMESTAMP, RSIENA_ITERATION, self$N, self$M, self$BI_PROB, avg_degree_social, avg_degree_component)
    #     ggsave(filename = sprintf('3plot_SAOM-NK_networks_%s_%s.png', keystring, self$TIMESTAMP), 
    #            plot = plt, height = 4.5, width = 9, units = 'in', dpi = 300)
    #   }
    # },
    
  
  
  
  
  # Visualization using ggplot2 and ggraph: bipartite network, social network, and component interaction matrix heatmap
  plot_bipartite_system_from_mat = function(bipartite_matrix, RSIENA_ITERATION, 
                                            plot_save = FALSE, return_plot=TRUE,
                                            normalize_degree = FALSE) {
    # #
    # RSIENA_ITERATION <- length(self$rsiena_model$sims)
    
    N <- ncol(bipartite_matrix)
    M <- nrow(bipartite_matrix)
    
    TS <- round(as.numeric(Sys.time()) * 100)
    
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
    ig_bipartite <- igraph::graph_from_biadjacency_matrix(bipartite_matrix, directed = F, weighted = T)
    
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
    component_matrix <- igraph::as_adjacency_matrix(ig_component, type = 'both', sparse = F, attr = 'weight')
    
    component_df <- melt(component_matrix)
    colnames(component_df) <- c("Component1", "Component2", "Interaction")
    
    
    # Replace component numbers with labels in the heatmap data frame
    component_df$Component1 <- factor(component_df$Component1, labels = component_labels)
    component_df$Component2 <- factor(component_df$Component2, labels = component_labels)
    
    heatmap_plot <- ggplot(component_df, aes(x = Component1, y = Component2, fill = Interaction)) +
      geom_tile() +
      # scale_fill_gradient(low = "white", high = "red") +
      labs(title = "[Proj2] Component Heatmap\n(Epistatic Interactions)", x = "Component 1", y = "Component 2") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    if ( sum(component_df$Interaction) == 0 ) {
      heatmap_plot <- heatmap_plot + scale_fill_gradient(low = "white", high = "white")
    } else {
      heatmap_plot <- heatmap_plot + scale_fill_gradient(low = "white", high = "red")
    }
    
    # Calculate average degree (K) for the social space and component interaction space
    avg_degree_social <- mean(igraph::degree(ig_social, mode = 'all', loops = F, normalized = normalize_degree))
    avg_degree_component <- mean(igraph::degree(ig_component, mode = 'all', loops = F, normalized = normalize_degree))
    
    # Create a main title using sprintf with simulation parameters
    main_title <- sprintf("Simulation Parameters:\niter = %s, N = %s, M = %s, BI_PROB = %.2f, K_A = %.2f, K_C = %.2f", 
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
    if(return_plot)
      return(plt)
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
  },


  ## 
  plot_actor_utility = function(xints=c(), thin_factor=1, loess_span=0.5, return_plot=TRUE) {
    ## Individual actors
    # xints <- (1:shock.splits) * floor(max(self$actor_util_df$chain_step_id)/shock.splits)
    # xints <- c( length(beta_seq_1), length(beta_seq_1) + length(beta_seq_2))
    # thin_factor <- 1
    actthin <- self$actor_util_df %>%  filter(chain_step_id %% thin_factor == 0)
    print(dim(actthin))
    tmpdf <-  actthin %>% mutate(actor_id=as.factor(actor_id))
    #
    # nstep <- length(unique(actthin$chain_step_id))
    npoints <- nrow(self$chain_stats) * self$M
    #
    point_size <- 6 / log10( npoints )
    point_alpha <- min( 1,  .55/log10( npoints ) )
    #
    suppressMessages({
      plt.act <-  ggplot(aes(x=chain_step_id, y=utility, color=strategy), 
                         data=tmpdf) + # linetype=actor_id, color=strategy, point=actor_id, shape=actor_id
        geom_point(alpha=point_alpha, shape=1, size= point_size ) +
        geom_smooth(aes(fill=actor_id, linetype=actor_id), method='loess', span=loess_span, alpha=.1) +
        geom_smooth(aes(x=chain_step_id, y=mean), 
                    data=actthin %>% group_by(chain_step_id,actor_id) %>% summarize(mean=mean(utility, na.rm=T)),
                    method='loess', color='black', span=loess_span, alpha=.05, linewidth=1.1) +
        geom_hline(yintercept = 0, linetype=4, color='black') +
        theme_bw()
    })
    
    if (length(xints)) 
      plt.act <- plt.act + geom_vline(xintercept=xints, linetype=1, color='black') 
    
    if (return_plot)
      return(plt.act)
  },
  
  ##
  plot_strategy_utility = function(xints=c(), thin_factor=1, loess_span=0.5, return_plot=TRUE) {
    #
    # xints <- (1:shock.splits) * floor(max(self$actor_util_df$chain_step_id)/shock.splits)
    # xints <- c( length(beta_seq_1), length(beta_seq_1) + length(beta_seq_2))
    # thin_factor <- 1
    actthin <- self$actor_util_df %>%  filter(chain_step_id %% thin_factor == 0)
    print(dim(actthin))
    tmpdf <-  actthin %>% mutate(actor_id=actor_id)
    #
    # nstep <- length(unique(actthin$chain_step_id))
    npoints <- nrow(self$chain_stats) * self$M
    #
    point_size <- 6 / log10( npoints )
    point_alpha <- min( 1,  .55/log10( npoints ) )
    #
    suppressMessages({
      plt.act <- tmpdf %>%       #ungroup() %>%
        ggplot(aes(x=chain_step_id, y=utility, color=strategy)) +  #
        geom_point(aes(shape=actor_id), alpha=point_size, shape=1, size= point_size) +
        geom_smooth(aes(fill=strategy), 
                    method='loess', span=loess_span, alpha=.1) +
        geom_smooth(aes(x=chain_step_id, y=mean), 
                    method='loess', color='black', span=loess_span, alpha=.05, linewidth=1.1,
                    data=actthin %>% group_by(chain_step_id) %>% summarize(mean=mean(utility, na.rm=T)) %>% 
                      mutate(strategy=NA)) +
        geom_hline(yintercept = 0, linetype=4, color='black') +
        theme_bw()
    })
    
    if (length(xints)) 
      plt.act <- plt.act + geom_vline(xintercept=xints, linetype=1, color='black') 
    
    if (return_plot)
      return(plt.act)
    
  },


  ##
  plot_utility_components = function(use_contributions=TRUE, loess_span=0.5, return_plot=TRUE, save_plot=FALSE) {
    #
    structure_model <- self$config_structure_model
    ##==================================================
    #
    # eff <- m1$rsiena_effects[m1$rsiena_effects$include, ]
    efflist <- c(
      structure_model$dv_bipartite$effects,
      structure_model$dv_bipartite$coCovars,
      structure_model$dv_bipartite$varCovars,
      structure_model$dv_bipartite$coDyadCovars,
      structure_model$dv_bipartite$varDyadCovars
    )
    #
    neffs <- length(efflist)
    ### empty matrix to hold actor network statistics
    effnames <- sapply(efflist, function(x) x$effect, simplify = T)
    effparams <- sapply(efflist, function(x) x$parameter, simplify = T)
    
    names(effparams) <- effnames
    
    
    
    
    
    ACTORS     <- sienaNodeSet(self$M, nodeSetName="ACTORS")
    COMPONENTS <- sienaNodeSet(self$N, nodeSetName="COMPONENTS")
    # strat_mat_varCovar <- varCovar(strat_eff$x, nodeSet = 'ACTORS')
    if (! 'coCovars' %in% names(structure_model$dv_bipartite) ) {
      rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV), nodeSets = list(ACTORS, COMPONENTS))
      return(rsiena_data)
    }
    dv_name <- structure_model$dv_bipartite$name
    
    
    structeffs <- sapply(structure_model$dv_bipartite$effects, function(x)x$parameter)
    names(structeffs) <- sapply(structure_model$dv_bipartite$effects, function(x)x$effect)
    
    # coCovarParams     <- sapply(structure_model$dv_bipartite$coCovars, function(x)x$effect)
    # varCovarParams     <- sapply(structure_model$dv_bipartite$varCovars, function(x)x$effect)
    # coDyadCovarParams  <- sapply(structure_model$dv_bipartite$coDyadCovars, function(x)x$effect)
    # varDyadCovarParams <- sapply(structure_model$dv_bipartite$varDyadCovars, function(x)x$effect)
    
    coCovars     <- sapply(structure_model$dv_bipartite$coCovars, function(x)x$effect)
    varCovars     <- sapply(structure_model$dv_bipartite$varCovars, function(x)x$effect)
    coDyadCovars  <- sapply(structure_model$dv_bipartite$coDyadCovars, function(x)x$effect)
    varDyadCovars <- sapply(structure_model$dv_bipartite$varDyadCovars, function(x)x$effect)
    
    coCovars_param     <- sapply(structure_model$dv_bipartite$coCovars, function(x)x$parameter)
    varCovars_param     <- sapply(structure_model$dv_bipartite$varCovars, function(x)x$parameter)
    coDyadCovars_param  <- sapply(structure_model$dv_bipartite$coDyadCovars, function(x)x$parameter)
    varDyadCovars_param <- sapply(structure_model$dv_bipartite$varDyadCovars, function(x)x$parameter)
    
    covs <- unlist(c(coCovars_param, varCovars_param, coDyadCovars_param,  varDyadCovars_param))
    names(covs) <- unlist(c(coCovars, varCovars, coDyadCovars, varDyadCovars))
    
    coCovarTypes     <- sapply(structure_model$dv_bipartite$coCovars, function(x)x$interaction1)
    varCovarTypes     <- sapply(structure_model$dv_bipartite$varCovars, function(x)x$interaction1)
    coDyadCovarTypes  <- sapply(structure_model$dv_bipartite$coDyadCovars, function(x)x$interaction1)
    varDyadCovarTypes <- sapply(structure_model$dv_bipartite$varDyadCovars, function(x)x$interaction1)
    
    vartypes <- unlist(c(coCovarTypes, varCovarTypes, coDyadCovarTypes, varDyadCovarTypes))
    component_type_ids <-  grep('self\\$component.+var', vartypes)
    strat_type_ids <-  grep('self\\$strat.+var', vartypes)
    
    # component_coCovar_ids     <- grep('self\\$component.+_coCovar', coCovarTypes) ## ex: "self$component_1_coCovar"
    # component_varCovar_ids    <- grep('self\\$component.+_varCovar', varCovarTypes) ## ex: "self$component_1_coCovar"
    # component_coDyadCovar_ids  <- grep('self\\$component.+_coDyadCovar', coDyadCovarTypes) ## ex: "self$component_1_coCovar"
    # component_varDyadCovar_ids <- grep('self\\$component.+_varDyadCovar', varDyadCovarTypes) ## ex: "self$component_1_coCovar"
    
    # strat_coCovar_ids     <- grep('self\\$strat.+_coCovar', coCovarTypes) ## ex: "self$component_1_coCovar"
    # strat_varCovar_ids    <- grep('self\\$strat.+_varCovar', varCovarTypes) ## ex: "self$component_1_coCovar"
    # strat_coDyadCovar_ids  <- grep('self\\$strat.+_coDyadCovar', coDyadCovarTypes) ## ex: "self$component_1_coCovar"
    # strat_varDyadCovar_ids <- grep('self\\$strat.+_varDyadCovar', varDyadCovarTypes) ## ex: "self$component_1_coCovar"
    # stratDV_ids     <- grep('self\\$strat_\\d{1,2}_coCovar', covDvTypes ) ## ex: "self$strat_1_coCovar"
    
    # ncompo_coCovar     <- length( component_coCovar_ids )
    # ncompo_varCovar    <- length( component_varCovar_ids )
    # ncompo_coDyadCovar  <- length( component_coDyadCovar_ids )
    # ncompo_varDyadCovar <- length( component_varDyadCovar_ids )
    # 
    # nstrat_coCovar     <- length( strat_coCovar_ids )
    # nstrat_varCovar    <- length( strat_varCovar_ids )
    # nstrat_coDyadCovar  <- length( strat_coDyadCovar_ids )
    # nstrat_varDyadCovar <- length( strat_varDyadCovar_ids )
    
    
    # compo_vars <- c()
    
  
    ## actor strategy
    # if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
    #   stop("Actor Strategy self$strat_1_coCovar not set.")
    # actor_strat_lvls <- as.factor( self$get_actor_strategies() )
    # actor_strats <- rep(actor_strat_lvls, times = round(self$M/length(actor_strat_lvls)))
    #
    nstep <- sum(!self$chain_stats$stability)
    dv_bipartite <- self$config_structure_model$dv_bipartite
    # ## levels order
    # coveffs   <- sapply(dv_bipartite$coCovars, function(x)x$effect)
    # covparams <- sapply(dv_bipartite$coCovars, function(x)x$parameter)
    # covfixs   <- sapply(dv_bipartite$coCovars, function(x)ifelse(x$fix,'','(var)'))
    # #
    # strateffs   <- sapply(dv_bipartite$coCovars, function(x)x$effect)
    # stratparams <- sapply(dv_bipartite$coCovars, function(x)x$parameter)
    # stratfixs   <- sapply(dv_bipartite$coCovars, function(x)ifelse(x$fix,'','(var)'))
    # # stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)substr(as.character(x$fix),1,1))
    # #
    # structeffs   <- sapply(dv_bipartite$effects, function(x)x$effect)
    # structparams <- sapply(dv_bipartite$effects, function(x)x$parameter)
    # structfixs   <- sapply(dv_bipartite$effects, function(x)ifelse(x$fix,'','(var)'))
    # # structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)substr(as.character(x$fix),1,1))
    
    # covDvTypes <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$interaction1)
    # componentDV_ids <- grep('self\\$component.+var', covDvTypes) ## ex: "self$component_1_coCovar"
    # stratDV_ids     <- grep('self\\$strat.+var', covDvTypes ) ## ex: "self$strat_1_coCovar"  
    
    
    
    #
    sim_title_str <- sprintf(
      'Environment: Actors (M) = %s, Components (N) = %s, Init.Prob. = %.2f\nActor Strategy:  %s\nComponent Payoff:  %s\nStructure:  %s', 
      self$M, self$N, self$BI_PROB,
      paste( paste(paste(names(covs)[strat_type_ids], covs[strat_type_ids], sep='= '), sep='' ), collapse = ';  '),
      paste( paste(paste(names(covs)[component_type_ids], covs[component_type_ids], sep='= '), sep=''), collapse = ';  '),
      paste( paste(paste(names(structeffs), structeffs, sep='= '), sep=''), collapse = ';  ')
    )
    # #
    # compoeffs   <- coveffs[ componentDV_ids ]
    # compoparams <- covparams[ componentDV_ids ]
    # compofixs   <- covfixs[ componentDV_ids ]
    # #
    # strateffs   <- coveffs[ stratDV_ids ]
    # stratparams <- covparams[ stratDV_ids ]
    # stratfixs   <- covfixs[ stratDV_ids ]
    efflvls <- c('utility', 
                 names(covs)[strat_type_ids], 
                 names(covs)[component_type_ids],
                 names(structeffs)
                 )
    ##=================================================
    ##--------------------------------------------------------------
    ## 2. Strategy Signals
    ##--------------------------------------------------------------
    actor_strats <- self$get_actor_strategies()
    #
    actor_stats_df <- self$actor_stats_df
    actor_stats_df$strategy <- actor_strats[ actor_stats_df$actor_id ]
    #
    actor_util_df <- self$actor_util_df  %>% mutate(
      effect_id=NA, 
      value=utility, 
      effect_name='utility', 
      strategy= actor_strats[ actor_id ], 
      utility=NULL
    )
    
    plt_title <- sprintf('Actor Utility Statistics Decomposition:\n%s', sim_title_str)
    
    ## use actor stats contributions to utility instead of original stats
    if (use_contributions) {
      actor_stats_df <- actor_stats_df %>% mutate(value=value_contributions)
      plt_title <- sprintf('Actor Utility Contributions (statistics * betas):\n%s', sim_title_str)
    } 
    
    ## Add utility as extra 'effect' 
    act_effs <- actor_stats_df %>% bind_rows( actor_util_df ) 
    act_effs$effect_name <- factor(act_effs$effect_name, levels=efflvls)
    
    #plot signals
    thinplt2 <- 1
    
    act_effs2 <- act_effs %>% filter(chain_step_id %% thinplt2 == 0)
    
    # nstep <- length(unique(act_effs2$chain_step_id))
    npoints <- nrow(self$chain_stats) * self$M
    #
    point_size <- 4 / log10( npoints )
    point_alpha <- min( 1,  .55/log10( npoints ) )
    
    suppressMessages({
      plt2 <- act_effs2  %>% 
        ggplot(aes(x=chain_step_id, y=value, color=strategy,fill=strategy, linetype=strategy)) + 
        geom_point(alpha=point_alpha, shape=1, size= point_size )  + 
        geom_smooth(method='loess', alpha=.1, span=loess_span) +
        facet_grid(effect_name ~ ., scales='free_y') +
        theme_bw() + theme(legend.position = 'bottom') +
        ggtitle(plt_title)
    })
    
    plotname2 <- sprintf('plot_actor_utility_components_thin%s_%s.png',thinplt2, round(as.numeric(Sys.time())*100) )
    
    if (save_plot) {
      plot_dir <- getwd()
      ggsave(filename=file.path(plot_dir, plotname2), plt2, 
             height = 12, width = 8, dpi = 600, units = 'in')
    }
    
    if (return_plot)
      return(plt2)
    
  },
  
  
  ##
  plot_component_degrees = function(loess_span=0.5, return_plot=TRUE) {
    
    Kdf2 <- self$K_B2_df %>% mutate(effect='K_B2') %>% 
      bind_rows(
        self$K_C_df %>% mutate(effect='K_C')
      ) 
    
    # nstep <- length(unique(Kdf2$chain_step_id))
    npoints <- nrow(self$chain_stats) * self$M
    #
    point_size <- 6 / log10( npoints )
    point_alpha <- min( 1,  .55/log10( npoints ) )
    
    suppressMessages({
      plt <- Kdf2 %>% 
        ggplot(aes(x=chain_step_id, y=value)) +
        geom_point(aes(fill=component_id, color=component_id), 
                   pch=1, alpha=point_alpha, size=point_size) + 
        geom_smooth(aes(fill=component_id, color=component_id, linetype=component_id), 
                    method='loess', alpha=.1, span=loess_span) + 
        geom_smooth(aes(x=chain_step_id, y=mean), span=loess_span, 
                    data=Kdf2 %>% group_by(chain_step_id, effect) %>% summarize(mean=mean(value)),
                    method='loess', se=F, color='black', size=1) +
        geom_hline(yintercept = 0, linetype=2) +
        facet_grid( effect ~ . ) +
        theme_bw() +
        ggtitle('Component Degrees: K_B2, K_C')
    })
    
    if (return_plot)
      return(plt)
    
  }, 
  
  plot_actor_degrees = function(loess_span=0.5, return_plot=TRUE) {
    
    Kdf1 <- self$K_A_df %>% mutate(effect='K_A', node_id=actor_id) %>% 
      bind_rows(
        self$K_B1_df %>% mutate(effect='K_B1', node_id=actor_id)
      )
    
    # nstep <- length(unique(Kdf1$chain_step_id))
    npoints <- nrow(self$chain_stats) * self$M
    #
    point_size <- 6 / log10( npoints )
    point_alpha <- min( 1,  .55/log10( npoints ) )
    
    suppressMessages({
      plt <- Kdf1 %>% 
        ggplot(aes(x=chain_step_id, y=value)) +
        # stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1) +
        geom_point(aes(color=strategy, fill=strategy), 
                   pch=1, alpha=point_alpha, size=point_size) + 
        geom_smooth(aes(color=strategy, fill=strategy, linetype=actor_id),
                    method='loess', alpha=.1, span=loess_span) + 
        geom_smooth(aes(x=chain_step_id, y=mean), span=loess_span, 
                    data=Kdf1 %>% group_by(chain_step_id, effect) %>% summarize(mean=mean(value)),
                    method='loess', se=F, color='black', size=1) +
        geom_hline(yintercept = 0, linetype=2) +
        facet_grid( effect ~ . ) +
        theme_bw() + 
        ggtitle('Actor Degrees: K_A, K_B1')
    })

    
    
    if (return_plot)
      return(plt)
  }, 
  
  
  ##
  plot_actor_utility_strategy_summary = function(actor_ids=c(), 
                                                 thin_factor=1, 
                                                 thin_wave_factor=1,
                                                 smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                 show_utility_points=T,
                                                 scale_utility=TRUE,
                                                 return_plot=TRUE,
                                                 plot_file=NA, 
                                                 plot_dir=NA,
                                                 loess_span=0.5
  ) {
    ## actor strategy
    if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
      stop("Actor Strategy self$strat_1_coCovar not set.")
    actor_strat <- as.factor( self$get_actor_strategies() )
    #
    npoints <- nrow(self$chain_stats) * self$M
    # nstep <- nrow(self$chain_stats) * self$M
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
    ## Compare 2 actors utilty
    dat <- self$actor_util_df %>% 
      filter(chain_step_id %% thin_factor == 0) %>% 
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
      if (!all(dat$utility == 0))
        dat <- dat %>% mutate(utility = c(scale(utility)))
    }
    #
    density_rng <- range(dat$utility, na.rm=T)
    density_absdiff_scale <- abs(diff(density_rng)) * 0.15
    util_lim <- c(density_rng[1] - density_absdiff_scale,  density_rng[2] + density_absdiff_scale)
    #
    point_size <- 6 / log10( npoints )
    point_alpha <- min( 1,  .55/log10( npoints ) )
    #
    if(length(actor_ids))
      dat <- dat %>% filter(actor_id %in% actor_ids)
    plt <- ggplot(dat, aes(x=chain_step_id, y=utility)) + 
      geom_hline(yintercept = dat$utility%>%mean(), linetype=3, col='black' ) 
    if(show_utility_points)
      plt <- plt + geom_point(aes(color=strategy), alpha=point_alpha, shape=1, size=point_size)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
    if(self$exists(smooth_method))
      plt <- plt + geom_smooth(aes(color=strategy, fill=strategy, linetype=actor_id), 
                               method = smooth_method, span=loess_span, linewidth=1, alpha=.09)
    #
    # plt <- plt + geom_text(aes(x=chain_step_id, y=utility, linetype=factor(actor_id)), size = 3, vjust = -1) 
    # Population Mean in black
    plt <- plt + geom_smooth(aes(x=chain_step_id, y=utility_mean), method = smooth_method, span=loess_span, 
                             data=dat %>% group_by(chain_step_id)%>%summarize(utility_mean=mean(utility)),
                             color='black', linetype=1, alpha=.1) 
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
    if ( ! 'coCovars' %in% names(self$config_structure_model$dv_bipartite) ) {
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
    stratmeans <- dat %>% group_by(strategy) %>% 
      dplyr::summarize(n=n(), mean=mean(utility, na.rm=T), sd=sd(utility, na.rm=T))
    ## Actor density fact plots comparing H1 to H2 utility distribution
    plt2 <- ggplot(dat, aes(x=utility, color=strategy, fill=strategy)) + ##linetype=chain_half
      geom_density(alpha=.1, linewidth=1)  +
      # geom_histogram(alpha=.1, position = 'dodge') +
      geom_vline(data = stratmeans, aes(xintercept = mean, color=strategy), linetype=2, linewidth=.9) +
      labs(y='', x='') +
      xlim(util_lim) +
      coord_flip() +
      # facet_grid(wave_id ~ .) +
      ylab('Actor Utility Density') +
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
    
    suppressMessages({
      combined_plot <- ggarrange(
        plt, plt2, 
        ncol = 2, 
        widths = c(4.1,0.9), # Adjust column widths
        common.legend = TRUE, # Share a common legend if needed
        legend = "bottom"#,     # Place legend at the bottom
        # align = 'h'
      )
    })
 
    # print(combined_plot)
    #
    # print(plt)
    #
    
    
    # if(!is.na(plot_file))
    #   ggsave(file = file.path(ifelse(is.na(plot_dir),getwd(),plot_dir), sprintf("%s_%s.png", self$config_environ_params$name, plot_file)), 
    #          combined_plot, 
    #          width = 10, height = 8, units = 'in', dpi = 600)
    
    if(return_plot)
      return(combined_plot)
  }
  
  
  
  #  ##
  # plot_utility_components = function(use_contributions=FALSE, loess_span=0.5, return_plot=TRUE, save_plot=FALSE) {
  #   ##==================================================
  #   ## actor strategy
  #   # if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
  #   #   stop("Actor Strategy self$strat_1_coCovar not set.")
  #   # actor_strat_lvls <- as.factor( self$get_actor_strategies() )
  #   # actor_strats <- rep(actor_strat_lvls, times = round(self$M/length(actor_strat_lvls)))
  #   #
  #   nstep <- sum(!self$chain_stats$stability)
  #   dv_bipartite <- self$config_structure_model$dv_bipartite
  #   ## levels order
  #   coveffs   <- sapply(dv_bipartite$coCovars, function(x)x$effect)
  #   covparams <- sapply(dv_bipartite$coCovars, function(x)x$parameter)
  #   covfixs   <- sapply(dv_bipartite$coCovars, function(x)ifelse(x$fix,'','(var)'))
  #   #
  #   strateffs   <- sapply(dv_bipartite$coCovars, function(x)x$effect)
  #   stratparams <- sapply(dv_bipartite$coCovars, function(x)x$parameter)
  #   stratfixs   <- sapply(dv_bipartite$coCovars, function(x)ifelse(x$fix,'','(var)'))
  #   # stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)substr(as.character(x$fix),1,1))
  #   #
  #   structeffs   <- sapply(dv_bipartite$effects, function(x)x$effect)
  #   structparams <- sapply(dv_bipartite$effects, function(x)x$parameter)
  #   structfixs   <- sapply(dv_bipartite$effects, function(x)ifelse(x$fix,'','(var)'))
  #   # structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)substr(as.character(x$fix),1,1))
  #   covDvTypes <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$interaction1)
  #   componentDV_ids <- grep('self\\$component_\\d{1,2}_coCovar', covDvTypes) ## ex: "self$component_1_coCovar"
  #   stratDV_ids     <- grep('self\\$strat_\\d{1,2}_coCovar', covDvTypes ) ## ex: "self$strat_1_coCovar"     
  #   #
  #   sim_title_str <- sprintf(
  #     'Environment: Actors (M) = %s, Components (N) = %s, Init.Prob. = %.2f\nActor Strategy:  %s\nComponent Payoff:  %s\nStructure:  %s', 
  #     self$M, self$N, self$BI_PROB,
  #     paste( paste(paste(strateffs[stratDV_ids], stratparams[stratDV_ids], sep='= '), stratfixs[stratDV_ids], sep='' ), collapse = ';  '),
  #     paste( paste(paste(strateffs[componentDV_ids], stratparams[componentDV_ids], sep='= '), stratfixs[componentDV_ids], sep=''), collapse = ';  '),
  #     paste( paste(paste(structeffs, structparams, sep='= '), structfixs, sep=''), collapse = ';  ')
  #   )
  #   #
  #   compoeffs   <- coveffs[ componentDV_ids ]
  #   compoparams <- covparams[ componentDV_ids ]
  #   compofixs   <- covfixs[ componentDV_ids ]
  #   #
  #   strateffs   <- coveffs[ stratDV_ids ]
  #   stratparams <- covparams[ stratDV_ids ]
  #   stratfixs   <- covfixs[ stratDV_ids ]
  #   efflvls <- c('utility', strateffs, compoeffs, structeffs)
  #   ##=================================================
  #   ##--------------------------------------------------------------
  #   ## 2. Strategy Signals
  #   ##--------------------------------------------------------------
  #   actor_strats <- self$get_actor_strategies()
  #   #
  #   actor_stats_df <- self$actor_stats_df
  #   actor_stats_df$strategy <- actor_strats[ actor_stats_df$actor_id ]
  #   #
  #   actor_util_df <- self$actor_util_df  %>% mutate(
  #     effect_id=NA, 
  #     value=utility, 
  #     effect_name='utility', 
  #     strategy= actor_strats[ actor_id ], 
  #     utility=NULL
  #   )
  #   
  #   ## Add utility as extra 'effect' 
  #   act_effs <- actor_stats_df %>% bind_rows( actor_util_df ) 
  #   act_effs$effect_name <- factor(act_effs$effect_name, levels=efflvls)
  #   
  #   #plot signals
  #   thinplt2 <- 1
  #   
  #   act_effs2 <- act_effs %>% filter(chain_step_id %% thinplt2 == 0)
  #   
  #   if (use_contributions) {
  #     plt2 <- act_effs2 %>% ggplot(aes(x=chain_step_id, y=value_contributions, color=strategy,fill=strategy, linetype=strategy))
  #   } else {
  #     plt2 <- act_effs2  %>% ggplot(aes(x=chain_step_id, y=value, color=strategy,fill=strategy, linetype=strategy))
  #   }
  #   
  #   plt2 <- plt2 + geom_point(alpha=.25, shape=16)  + 
  #     geom_smooth(method='loess', alpha=.1, span=loess_span) +
  #     facet_grid(effect_name ~ ., scales='free_y') +
  #     theme_bw() + theme(legend.position = 'left') +
  #     ggtitle(sprintf('Actor Utility Signals Decomposition by Effects:\n%s',sim_title_str))
  #   
  #   plotname2 <- sprintf('plot_actor_utility_components_thin%s_%s.png',thinplt2, round(as.numeric(Sys.time())*100) )
  #   
  #   if (save_plot) {
  #     plot_dir <- getwd()
  #     ggsave(filename=file.path(plot_dir, plotname2), plt2, 
  #            height = 12, width = 8, dpi = 600, units = 'in')
  #   }
  #   
  #   if (return_plot)
  #     return(plt2)
  #   
  # }



  )  ##/end pubic list
  
)  ##/end class


# ##-----------------





