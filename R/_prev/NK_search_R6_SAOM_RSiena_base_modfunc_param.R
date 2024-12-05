
# Load necessary libraries
library(R6)
library(Matrix)
library(network)
library(igraph)
library(visNetwork)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggraph)
library(RSiena)
library(texreg)
# Load the necessary library for grid.text
library(grid)


##
##
##
get_jaccard_index <- function(m0, m1) {
  mdiff <-  m1 - m0   ## new m1 - old m0
  vecdiff <- c(mdiff)
  # sumdiff <- sum( vecdiff )
  cnt_maintain <- sum( m1 * m0 )
  cnt_add <- length(which(vecdiff > 0))
  cnt_drop <- length(which(vecdiff < 0))
  return( cnt_maintain / (cnt_maintain + cnt_add + cnt_drop) )
}

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
    ## K = NULL,             # Avg. degree of component-[actor]-component interactions
    bipartite_igraph = NULL, # Bipartite network object 
    social_igraph = NULL,   # Social space projection (network object)
    search_igraph = NULL, # Search space projection (network object)
    #
    bipartite_matrix = NULL,
    social_matrix = NULL,
    search_matrix = NULL,
    #
    bipartite_rsienaDV = NULL,
    social_rsienaDV = NULL,
    search_rsienaDV = NULL,
    #
    fitness_landscape = NULL, # Fitness landscape matrix
    progress_scores = c(),  # Track scores over iterations
    P_change = NULL,          # Probability of changing the component interaction matrix
    #
    fitness_history = list(),
    degree_history_env = list(),  ##**TODO** Add
    degree_history_K_S = list(),
    degree_history_K_E = list(),
    #
    rsiena_model = NULL,   
    rsiena_data = NULL,
    rsiena_effects = NULL,
    rsiena_algorithm = NULL,
    
    # Constructor to initialize the SAOM-NK model
    initialize = function(params) {
      cat('\nCALLED _BASE_ INIT\n')
      ## ----- prevent clashes with sna package-----------
      sna_err_check <-tryCatch(expr = { detach('package:sna') }, error=function(e)e )
      ## -------------------------------------------------
      ##**TODO:  LOAD DEPENDENCY FUNCTIONS ETC**
      ##
      self$M <- params[['M']]
      self$N <- params[['N']]
      self$BI_PROB <- params[['BI_PROB']]
      self$P_change <- params[['P_change']]
      #
      # self$bipartite_igraph <- self$generate_bipartite_igraph()
      # self$social_network <- self$project_social_space()
      # self$search_landscape <- self$project_search_space()
      self$set_system_from_bipartite_igraph( self$random_bipartite_igraph() )
      #
      self$TIMESTAMP <- round( as.numeric(Sys.time())*100 )
      self$SIM_NAME <- params[['sim_name']]
      #
      visualize_init <- ifelse(is.null(params[['visualize_init']]), TRUE, params[['visualize_init']])
      if (visualize_init)
        self$visualize_bi_env_rsiena(plot_save = TRUE)
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
      # bipartite_net <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
      bipartite_igraph <- igraph::graph_from_biadjacency_matrix(bipartite_matrix, directed = F, mode = 'all')
      return(bipartite_igraph)
    },
    
    # get_bipartite_projections_from_matrix = function(bipartite_matrix) {
    #   # 1. Bipartite network plot using ggraph with vertex labels
    #   ig_bipartite <- igraph::graph_from_biadjacency_matrix(bipartite_matrix, directed = F, mode = 'all')
    #   E(ig_bipartite)$weight <- 1
    #   projs <- igraph::bipartite.projection(ig_bipartite, multiplicity = T, which = 'both')
    #   # ig_social <- projs$proj1
    #   # ig_component <- projs$proj2
    #   return(projs)
    # },
    
    get_bipartite_projections = function(ig_bipartite) {
      projs <- igraph::bipartite.projection(ig_bipartite, multiplicity = T, which = 'both')
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
    
    # # Project the social space (actor network) from the bipartite network
    # project_social_space = function() {
    #   bipartite_matrix <- as.matrix.network(self$bipartite_igraph)
    #   actor_matrix <- bipartite_matrix %*% t(bipartite_matrix)
    #   diag(actor_matrix) <- 0  # Remove self-loops
    #   social_net <- network(actor_matrix, directed = FALSE)
    #   return(social_net)
    # },
    # 
    # # Project the search landscape space (component interaction network)
    # project_search_space = function() {
    #   bipartite_matrix <- as.matrix.network(self$bipartite_igraph)
    #   component_matrix <- t(bipartite_matrix) %*% bipartite_matrix
    #   diag(component_matrix) <- 0  # Remove self-loops
    #   search_net <- network(component_matrix, directed = FALSE)
    #   return(search_net)
    # },
    
    # get_social_igraph = function() {
    #   return(graph_from_adjacency_matrix(as.matrix.network(self$social_network), mode = "undirected"))
    # },
    # 
    # get_component_igraph = function() {
    #   return(graph_from_adjacency_matrix(as.matrix.network(self$search_landscape), mode = "undirected"))
    # },
    
    # Update Iteration progress
    increment_sim_iter = function(val = 1) {
      self$ITERATION <- self$ITERATION + val
    }, 
    
    
        ##
    include_effect_from_eff_list = function(eff) {
      
      if (eff$effect == 'density') {
        self$rsiena_effects <- includeEffects(self$rsiena_effects,  density, ## get network statistic function from effect name (character)
                                             name = eff$dv_name, parameter = eff$parameter,
                                             # interaction1 = eff$interaction1,
                                             fix = eff$fix)
      }
      if (eff$effect == 'transTriads') {
        self$rsiena_effects <- includeEffects(self$rsiena_effects,  transTriads, ## get network statistic function from effect name (character)
                                             name = eff$dv_name, parameter = eff$parameter,
                                             # interaction1 = eff$interaction1,
                                             fix = eff$fix)
      }
      if (eff$effect == 'inPop') {
        self$rsiena_effects <- includeEffects(self$rsiena_effects, inPop, ## get network statistic function from effect name (character)
                                             name = eff$dv_name, parameter = eff$parameter,
                                             # interaction1 = eff$interaction1,
                                             fix = eff$fix)
      }
      if (eff$effect == 'outAct') {
        self$rsiena_effects <- includeEffects(self$rsiena_effects, outAct, ## get network statistic function from effect name (character)
                                             name = eff$dv_name, parameter = eff$parameter,
                                             # interaction1 = eff$interaction1,
                                             fix = eff$fix)
      }
      if (eff$effect == 'cycle4') {
        self$rsiena_effects <- includeEffects(self$rsiena_effects, cycle4, ## get network statistic function from effect name (character)
                                             name = eff$dv_name,
                                             parameter = eff$parameter,
                                             # interaction1 = eff$interaction1,
                                             fix = eff$fix)
      }
      
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
    initialize = function(params) {
      ##
      cat('\nTEST FROM CALLED CLASS INIT *BEFORE* BASE INIT\n')
      ##
      super$initialize(params)
      ##
      cat('\nTEST FROM CALLED CLASS INIT *AFTER* BASE INIT\n')
      ##
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
    
    #
    init_rsiena_model_from_bipartite_matrix = function(bipartite_matrix) {
      #
      ACTORS     <- sienaNodeSet(self$M, nodeSetName="ACTORS")
      COMPONENTS <- sienaNodeSet(self$N, nodeSetName="COMPONENTS")
      
      ## init networks to duplicate for the init arrays (two network waves)
      array_bi_net <- array(c(bipartite_matrix, bipartite_matrix), 
                            dim=c(self$M, self$N, 2) )
      array_social <- array(c(self$social_matrix, self$social_matrix), 
                            dim=c(self$M, self$M, 2) )
      array_search <- array(c(self$search_matrix, self$search_matrix), 
                            dim=c(self$N, self$N, 2) )
      # Drop information above binary ties for RSiena DVs
      array_bi_net[ array_bi_net > 1 ] <- 1
      array_social[ array_social > 1 ] <- 1
      array_search[ array_search > 1 ] <- 1
      #
      self$social_rsienaDV    <- sienaDependent(array_social, type='oneMode', nodeSet = 'ACTORS', allowOnly = F)
      self$search_rsienaDV    <- sienaDependent(array_search, type='oneMode', nodeSet = 'COMPONENTS', allowOnly = F)
      self$bipartite_rsienaDV <- sienaDependent(array_bi_net, type='bipartite', nodeSet =c('ACTORS', 'COMPONENTS'), allowOnly = F)
      #

      ###
      self$rsiena_data <- sienaDataCreate(list(self$social_rsienaDV,
                                                self$search_rsienaDV,
                                                self$bipartite_rsienaDV),
                                         nodeSets = list(ACTORS, COMPONENTS))
      ###
      # self$rsiena_data <- rsienaDataCreate(list(self$bipartite_rsienaDV), 
      #                                    nodeSets = list(ACTORS, COMPONENTS))
      
    },
    
    ##
    ##**TODO**
    ##**Create custom RSiena interaction functions for only bipartite DV, **
    ##**but objective function includes statistics of the projections (social net, search landscape)**
    ##
    add_rsiena_effects = function(objective_list) {
      
      if (is.null(self$rsiena_effects))
        stop('initiate self$rsiena_effects before adding effects.')

      
      for (i in 1:length(objective_list)) {
        
        dv <- objective_list[[ i ]]
        
        for (j in 1:length(dv$effects)) {
          
          cat(sprintf('\n i=%s, j=%s\n', i, j))
          
          eff <- dv$effects[[ j ]]
          
          # print(eff)
          # stop('DEBUG')
          
          self$include_effect_from_eff_list(eff)
          # self$rsiena_effects <- setEffect(self$rsiena_effects, density, 
          #                                 name = 'self$bipartite_rsienaDV', parameter = 0, fix=T)  ## effect parameter 
        }
        
      }
      
    },
    
  
    
    # RSIENA
    local_search_rsiena_init = function(objective_list,  
                                        returnDeps = TRUE, ## TRUE = simulation only
                                        get_eff_doc = FALSE) {
      
      self$init_rsiena_model_from_bipartite_matrix(self$bipartite_matrix)

      # Step 3: Define the effects for the model
      self$rsiena_effects <- getEffects(self$rsiena_data)
      
      # Effects Documentation
      if(get_eff_doc)
        effectsDocumentation(self$rsiena_effects)
      ##-----------------------------
      
      
      # stop('DEBUG XXXXXXXXXXXXX')
      
      ## Add effects from model objective function list
      self$add_rsiena_effects(objective_list)
      
    },
    
    
    #
    local_search_rsiena_execute_sim = function(iterations, returnDeps=T, rsiena_phase2_nsub=1, rsiena_n2start_scale=1) {
      ##-----------------------------
      ## RSiena Algorithm
      self$rsiena_algorithm <- sienaAlgorithmCreate(projname=sprintf('%s_%s',self$SIM_NAME,self$TIMESTAMP),
                                                   simOnly = T,
                                                   nsub = rsiena_phase2_nsub,
                                                   n2start = rsiena_n2start_scale * 2.52 * (7+sum(self$rsiena_effects$include)),
                                                   n3 = iterations )
      
      
      
      # Step 4: Run RSiena simulation
      self$rsiena_model <- siena07(self$rsiena_algorithm, data = self$rsiena_data, effects = self$rsiena_effects,
                                    batch = TRUE, returnDeps = returnDeps)
      
      # Step 5: Summarize and plot results
      mod_summary <- summary(self$rsiena_model)
      if(!is.null(mod_summary)) 
        print(mod_summary)
      # plot(self$rsiena_model)
      
      
      
      print(screenreg(list(self$rsiena_model), single.row = T, digits = 3))
      
      ## update simulation object environment from 
      new_bi_env_igraph <- self$get_bipartite_igraph_from_rsiena_sim()
      

      self$set_system_from_bipartite_igraph( new_bi_env_igraph )

      #
      self$plot_bipartite_system_from_mat(self$bipartite_matrix, iterations, plot_save = T)
      
      # stop('OK_DEBG')
      
    
    },  
  
    #
    get_bipartite_matrix_from_rsiena_sim = function(sim_iteration = NULL) {
      ######### UPDATE SYSTEM ENVIRONMENT SNAPSHOT FROM EVOLVED Bipartite Network DV #####################
      # sim_id_last <- length( self$rsiena_model$sims ) ##
      sim_id <- ifelse(is.null(sim_iteration), 
                       length(self$rsiena_model$sims), ## default current state is the last simulation in sims list
                       sim_iteration)
      ##
      # actors <- 1:self$M 
      # components <- 1:self$N
      MplusN <- self$M + self$N
      ## Get 3rd DV (bi-partite network) from simulation iteration=sim_id
      el_bi_env <- self$rsiena_model$sims[[ sim_id ]][[1]][[3]]$`1`
      ## update numbering of second mode (the comonent integer names shift upward by the number of actors)
      el_bi_env[,2] <- el_bi_env[,2] + self$M
      ### get networks from other prjected space ties 
      # el_proj1  <- self$rsiena_model$sims[[1]][[1]][[1]]$`1` ## Social network
      # el_proj2  <- self$rsiena_model$sims[[1]][[1]][[2]]$`1` ## 
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
    
    #
    get_bipartite_igraph_from_rsiena_sim = function(sim_iteration = NULL) {
      new_bi_env_mat <- self$get_bipartite_matrix_from_rsiena_sim(sim_iteration)
      new_bi_env_igraph <- igraph::graph_from_biadjacency_matrix(new_bi_env_mat, directed = F, weighted = T, mode = 'all')
      return(new_bi_env_igraph)
    },
    
    #
    get_rsiena_model_snapshots = function(n=4, plot_save = TRUE, return_list=TRUE) {
      #
      snapshot_sim_ids <- seq(0,1,length.out= n) * length(self$rsiena_model$sims)
      snapshot_sim_ids[1] <- 1
      ##
      MplusN <- self$M + self$N
      ##
      snapshots <- list()
      for (sim_id in snapshot_sim_ids) {
        ##
        el_bi_env <- self$rsiena_model$sims[[ sim_id ]][[1]][[3]]$`1`
        ## update numbering of second mode (the comonent integer names shift upward by the number of actors)
        el_bi_env[,2] <- el_bi_env[,2] + self$M
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
        
          
        bi_env_mat_new  <- as.matrix(bi_env_mat_sp)[ 1:self$M, (self$M+1):(MplusN) ]
        
        self$plot_bipartite_system_from_mat(bi_env_mat_new, sim_id, plot_save = TRUE) 
        
        snapshots[[ sprintf('sim%d',sim_id) ]] <- bi_env_mat_new
      }
     
      ##
      if(return_list) return( snapshots )
      
    },
  
    #
    plot_rsiena_sim_stability = function(tol=1e-5) {
      
      sims = self$rsiena_model$sims
      n <- length(sims)
      outlist <- list()
      difflist <- list()
      jaccardlist <- list()
      K_soc_list <- list()
      K_env_list <- list()
      for(i in 1:n) {
        cat(sprintf(' %s ', i))
        ##
        el_bi_env <- sims[[ i ]][[1]][[3]]$`1`
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
                                                          directed = F, mode = 'all', 
                                                          multiple = T, weighted = T, 
                                                          add.names = T)
        projections <- igraph::bipartite.projection(new_bi_g, multiplicity = T, which = 'both')
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
      jaccard_vec <- plyr::ldply(jaccardlist)[2:n, 2] ## skip first period (no change yet)
      n_changes <- length(jaccard_vec)
      
      stability_vec <- cumsum(jaccard_vec) / (1:n_changes)
      plot(stability_vec , type='l' , ylab='Jaccard Index [t-1, t]', main='Jaccard Index Moving Average')
      
      # stability_delta <- cumsum(stability_vec) / (1:n_changes)
      stability_delta <- c(0, abs(diff(stability_vec))) / abs(stability_vec)
      plot(stability_delta, type='l', log='y',
           xlab='t', 
           # ylim=c( stability_delta[n-1],  1 ),
           ylab='Ln Stability Change [t-1, t]', main='Stabilization Rate\n(Change of Stability Moving Average)' 
      ); abline(h = tol, col='pink', lty=2)
      
      burn_prop <- 0.2
      iter_postburn <- round(c( 1-burn_prop, 1) * n_changes )
      idx <- iter_postburn[1]:iter_postburn[2]
      if (any(idx))
        hist(stability_vec[ idx ], main='Stability (post-burn)', breaks=31)
      
    },
    
    # Extend existing simulation environment
    local_search_rsiena_extend = function(iterations, returnDeps = TRUE) {
      self$rsiena_algorithm$n3 <- iterations
      self$rsiena_model <- rsiena07(self$rsiena_algorithm, data = self$rsiena_data, effects = self$rsiena_effects,
                                  prevAns = self$rsiena_model,
                                  batch = TRUE, returnDeps = returnDeps)
      ##
      mod_summary <- summary(self$rsiena_model)
      if(!is.null(mod_summary)) 
        print(mod_summary)
      ##
      print(screenreg(list(self$rsiena_model), single.row = T, digits = 3))
    },
  

    ##**TODO: SIENA MODEL SIMULATION SNAPSHOTS**
    ## Single Simulation Run
    local_search_rsiena_run = function(objective_list, iterations, n_snapshots=1, 
                                       returnDeps=TRUE, ## TRUE=simulation only
                                       plot_save = TRUE, overwrite=TRUE, 
                                       rsiena_phase2_nsub=1, rsiena_n2start_scale=1,
                                       get_eff_doc = FALSE
                                       ) {
      if( overwrite | is.null(self$rsiena_model) ) {
        self$local_search_rsiena_init(objective_list, returnDeps=returnDeps, get_eff_doc=get_eff_doc)
        self$local_search_rsiena_execute_sim(iterations, 
                                             rsiena_phase2_nsub=rsiena_phase2_nsub, 
                                             rsiena_n2start_scale=rsiena_n2start_scale)
      } else {
        # self$local_search_rsiena_extend(objective_list, iterations)
        stop('_extend() method not yet implemented.')
      }
      # # self$visualize_bi_env_rsiena(plot_save = plot_save)
      # self$get_rsiena_model_snapshots(n=n_snapshots)
    },

    # ## Batch of Multiple Simulation Extensions to 
    # local_search_rsiena_batchrun = function(batches=10, iterations=100, plot_save = TRUE, overwrite=FALSE, rsiena_phase2_nsub=1) {
    #   
    #   for (i in 1:batches) {
    #   
    #     overwrite_by_batch <- ifelse(i == 1, TRUE, overwrite)
    #     if(overwrite_by_batch) {
    #       self$local_search_rsiena_extend(iterations)
    #     } else {
    #       self$local_search_rsiena_init(iterations, rsiena_phase2_nsub=rsiena_phase2_nsub)
    #     }
    # 
    #    
    #     #rsiena_run(iterations=iterations, plot_save=plot_save, overwrite=overwrite_by_batch)
    #     
    #     ## Take snapshot of evolving system at this iteration of progress
    #     self$visualize_bi_env_rsiena(plot_save = plot_save)
    #     
    #   }
    #   
    # },
    
    
    ##
    ##
    ##
    ##**TODO**
    ##**CREATE VISUALIZE_NETWORKS plots for the RSiena approach local_search_rsiena_batchrun **
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
    visualize_bi_env_rsiena = function(plot_save = FALSE) {
      #
      RSIENA_ITERATION <- length(self$rsiena_model$sims)
      
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
    ig_bipartite <- igraph::graph_from_biadjacency_matrix(bipartite_matrix, directed = F, weighted = T)
    
    # E(ig_bipartite)$weight <- 1
    projs <- igraph::bipartite.projection(ig_bipartite, multiplicity = T, which = 'both')
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
      ig_social <- graph_from_adjacency_matrix(as.matrix.network(self$social_network), mode = "undirected")
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
# ##-----------------------------
# # bipartite effects
# self$rsiena_effects <- includeEffects(self$rsiena_effects, 
#                                      density, 
#                                      name = 'self$bipartite_rsienaDV')
# self$rsiena_effects <- setEffect(self$rsiena_effects, density, 
#                                 name = 'self$bipartite_rsienaDV', parameter = 0, fix=T)  ## effect parameter vs. initial value ?
# 
# # # 
# # self$rsiena_effects <- includeEffects(self$rsiena_effects, 
# #                                      inPop, 
# #                                      name = 'self$bipartite_rsienaDV')
# # self$rsiena_effects <- setEffect(self$rsiena_effects, inPop, 
# #                                 name = 'self$bipartite_rsienaDV', parameter = .5)  ## effect parameter vs. initial value ?
# # # 
# # self$rsiena_effects <- includeEffects(self$rsiena_effects, 
# #                                      outAct, 
# #                                      name = 'self$bipartite_rsienaDV')
# # self$rsiena_effects <- setEffect(self$rsiena_effects, outAct, 
# #                                 name = 'self$bipartite_rsienaDV', parameter = .5)  ## effect parameter vs. initial value ?
# # # 
# # self$rsiena_effects <- includeEffects(self$rsiena_effects, 
# #                                      cycle4, 
# #                                      name = 'self$bipartite_rsienaDV')
# # self$rsiena_effects <- setEffect(self$rsiena_effects, cycle4, 
# #                                 name = 'self$bipartite_rsienaDV', parameter = -.2, fix=T)  ## effect parameter vs. initial value ?
# # #-------------------------------
# ##-----------------------------
# ## SOCIAL SPACE
# self$rsiena_effects <- includeEffects(self$rsiena_effects, 
#                                      inPop, outAct, include = FALSE, 
#                                      name = 'self$social_rsienaDV')
# ## SOcial density
# # self$rsiena_effects <- includeEffects(self$rsiena_effects,
# #                                      density,
# #                                      name = 'self$social_rsienaDV')
# self$rsiena_effects <- setEffect(self$rsiena_effects, density,
#                                 name = 'self$social_rsienaDV', parameter = 0, fix=T)  ## effect parameter vs. initial value ?
# # ## Social transitivity
# # self$rsiena_effects <- includeEffects(self$rsiena_effects,
# #                                      transTriads,
# #                                      name = 'self$social_rsienaDV')
# # self$rsiena_effects <- setEffect(self$rsiena_effects, transTriads,
# #                                 name = 'self$social_rsienaDV', parameter = 0.4, fix = T)  ## effect parameter vs. initial value ?
# # 
# # ##
# ##-----------------------------
# # # Search landscape DV interactions with bipartite     
# # self$rsiena_effects <- includeEffects(self$rsiena_effects, 
# #                                      from.w.ind, 
# #                                      name = 'self$search_rsienaDV',
# #                                      interaction1 = "self$bipartite_rsienaDV",
# #                                      parameter = .5)
# # self$rsiena_effects <- includeEffects(self$rsiena_effects, 
# #                                      density, 
# #                                      name = 'self$search_rsienaDV')
# self$rsiena_effects <- setEffect(self$rsiena_effects, density,
#                                 name = 'self$search_rsienaDV', parameter = 0, fix=T)  ## effect parameter vs. initial value ?
# ##-----------------------------
# },
# ##-----------------



##**TODO**
## 1. [list] environ_params     Environmental parameters list
## 2. [list] structure_model    SAOM for environmental structural evolution (network ties: social, search, bipartite)
## 3. [list] payoff_formulas    [1] structure_model if searching=shaping; [2] arbitrary payoffs if search=/=shaping
## 4. [list] < other search specific attributes, like environ shock from regulatory change ? >
##**/**


## Environment determines what types of entities and strategies are permissible
params <- list(M = 10, N = 20, BI_PROB = .2, sim_name = '_TESTrsiena_')
## Strategies sets the objective function as a linear combination of network stats across DVs
objective_list <- list(
  dv_social = list(
    dv_name = 'self$social_rsienaDV',
    effects = list(
      list(effect='density', parameter=0, fix=F, dv_name='self$social_rsienaDV'), #interaction1 = NULL
      list(effect='transTriads', parameter=0, fix=F, dv_name='self$social_rsienaDV')#, #interaction1 = NULL
      # list(effect='density', parameter=0, fix=T),
    )
  ),
  dv_search = list(
    dv_name = 'self$search_rsienaDV',
    effects = list(
      list(effect='density', parameter=0, fix=F, dv_name='self$search_rsienaDV'), #interaction1 = NULL
      list(effect='transTriads', parameter=0, fix=F, dv_name='self$search_rsienaDV')#, #interaction1 = NULL
      # list(effect='density', parameter=0, fix=T),
    )
  ),
  dv_bipartite = list(
    name = 'self$bipartite_rsienaDV',
    effects = list(
      list(effect='density', parameter=0, fix=F, dv_name='self$bipartite_rsienaDV'), #interaction1 = NULL
      # list(effect='inPop', parameter=0, fix=F, dv_name='self$bipartite_rsienaDV'), #interaction1 = NUL
      list(effect='outAct', parameter=0, fix=F, dv_name='self$bipartite_rsienaDV'), #interaction1 = NULL
      list(effect='outInAss', parameter=0, fix=F, dv_name='self$bipartite_rsienaDV'), #interaction1 = NULL
      list(effect='cycle4', parameter=0, fix=F, dv_name='self$bipartite_rsienaDV')#, #interaction1 = NULL
      # list(effect='density', parameter=0, fix=T),
    )
  )
)
####
# objective_list <- list(
#   dv_social = list(
#     dv_name = 'self$social_rsienaDV',
#     effects = list(
#       list(effect='density', parameter=1, fix=T, dv_name='self$social_rsienaDV'), #interaction1 = NULL
#       list(effect='transTriads', parameter=1, fix=T, dv_name='self$social_rsienaDV')#, #interaction1 = NULL
#       # list(effect='density', parameter=0, fix=T),
#     )
#   ),
#   dv_search = list(
#     dv_name = 'self$search_rsienaDV',
#     effects = list(
#       list(effect='density', parameter=1, fix=T, dv_name='self$search_rsienaDV'), #interaction1 = NULL
#       list(effect='transTriads', parameter=1, fix=T, dv_name='self$search_rsienaDV')#, #interaction1 = NULL
#       # list(effect='density', parameter=0, fix=T),
#     )
#   ),
#   dv_bipartite = list(
#     name = 'self$bipartite_rsienaDV',
#     effects = list(
#       list(effect='density', parameter=1, fix=T, dv_name='self$bipartite_rsienaDV'), #interaction1 = NULL
#       # list(effect='inPop', parameter=0, fix=F, dv_name='self$bipartite_rsienaDV'), #interaction1 = NUL
#       list(effect='outAct', parameter=1, fix=T, dv_name='self$bipartite_rsienaDV'), #interaction1 = NULL
#       list(effect='outInAss', parameter=1, fix=T, dv_name='self$bipartite_rsienaDV'), #interaction1 = NULL
#       list(effect='cycle4', parameter=1, fix=T, dv_name='self$bipartite_rsienaDV')#, #interaction1 = NULL
#       # list(effect='density', parameter=0, fix=T),
#     )
#   )
# )

#
m1 <- SaomNkRSienaBiEnv$new(params)
# m1$get_effects_doc()
m1$local_search_rsiena_run(objective_list, iterations=1000, returnDeps = F,  n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
m1$plot_rsiena_sim_stability()

 
id_wave <- 1
### ## m1$rsiena_model$sf2[ <sims> , id_wave, <effects> ]

m1$rsiena_model$sf2

mean <-  m1$rsiena_model$sf2 
rownames(wavemat) <- as.character( 1:nrow(wavemat) )
colnames(wavemat) <- gsub('\\s', '_' , m1$rsiena_model$effects$effectName )
# dat[,,3] <- dat[1,1,][ m1$rsiena_model$effects$effectName ]

df_long <- melt(wavemat, varnames = c("Simulation", "Wave", "NetworkEffect"), value.name = "Mean")

# Plot using ggplot2
ggplot(df_long, aes(x = Mean, color=NetworkEffect, fill=NetworkEffect)) +
  # geom_boxplot(aes(group = Wave), outlier.shape = NA) +
  geom_density(alpha=.1, bins=11) +
  facet_wrap(~ NetworkEffect, scales = "free_x") +
  theme_minimal() +
  labs(title = "Distribution of Simulations by Effects",
       x = "X",
       y = "Y")



##$sf2  -> [sims,  waves,  effects]
##$sf   -> [sims, effects]

apply(m1$rsiena_model$sf2[,1,], 2, mean)

colMeans(m1$rsiena_model$sf)



##------------------------------------

# m1 <-

# We can also format the predictor matrix into an array:
mplearray <- ergmMPLE(formula, output="array")

#
dvnet <- 
formula <- as.formula(sprintf('%s ~ edges',dvnet))

# The resulting matrices are big, so only print the first 8 actors:
mplearray$response[1:8,1:8]
mplearray$predictor[1:8,1:8, ]
mplearray$weights[1:8,1:8]

##------------------------------------









###--------------------------------------------------------------------------

m1$local_search_rsiena_run(objective_list, iterations=2000, returnDeps = T, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
m1$plot_rsiena_sim_stability()


# ##
# m1$local_search_rsiena_run(iterations=4000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# m1$plot_rsiena_sim_stability()
# m1$local_search_rsiena_run(iterations=400, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# m1$plot_rsiena_sim_stability()




m1 <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .11, sim_name = '_TESTrsiena_')
m1$local_search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
m1$plot_rsiena_sim_stability()
#
m2 <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .11, sim_name = '_TESTrsiena_')
m2$local_search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
m2$plot_rsiena_sim_stability()

#
ms1 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .11, sim_name = '_TESTrsiena_')
ms1$local_search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
ms1$plot_rsiena_sim_stability()
#
ms2 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .11, sim_name = '_TESTrsiena_')
ms2$local_search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
ms2$plot_rsiena_sim_stability()



md1 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .2, sim_name = '_TESTrsiena_')
md1$local_search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
md1$plot_rsiena_sim_stability()
#
md2 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .2, sim_name = '_TESTrsiena_')
md2$local_search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
md2$plot_rsiena_sim_stability()

# #
# mc1 <- SAOM_NK_RSiena$new(M = 12, N = 30, BI_PROB = .05, sim_name = '_TESTrsiena_')
# mc1$local_search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
# mc1$plot_rsiena_sim_stability()
# #
# mc2 <- SAOM_NK_RSiena$new(M = 12, N = 30, BI_PROB = .05, sim_name = '_TESTrsiena_')
# mc2$local_search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=2, rsiena_n2start_scale = 1)
# mc2$plot_rsiena_sim_stability()

# #
# ms2 <- SAOM_NK_RSiena$new(M = 20, N = 10, BI_PROB = .15, sim_name = '_TESTrsiena_')
# ms2$local_search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=1, rsiena_n2start_scale = .1)
# ms2$plot_rsiena_sim_stability()


#
# m3 <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .15, sim_name = '_TESTrsiena_')
# m3$local_search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=3, rsiena_n2start_scale = 1)
# m3$plot_rsiena_sim_stability()
#
# m4 <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .15, sim_name = '_TESTrsiena_')
# m4$local_search_rsiena_run(iterations=2000, n_snapshots =1, rsiena_phase2_nsub=4, rsiena_n2start_scale = 1)
# m4$plot_rsiena_sim_stability()


saomnkrsiena <- SAOM_NK_RSiena$new(M = 10, N = 20, BI_PROB = .15, sim_name = '_TESTrsiena_')
##
# saomnkrsiena$local_search_rsiena_run(iterations=10, rsiena_phase2_nsub=1) 
# saomnkrsiena$local_search_rsiena_run(iterations=100, rsiena_phase2_nsub=2)
# saomnkrsiena$local_search_rsiena_run(iterations=5000, rsiena_phase2_nsub=3)
saomnkrsiena$local_search_rsiena_run(iterations=1000, n_snapshots =2, rsiena_phase2_nsub=1, rsiena_n2start_scale = 1)
saomnkrsiena$plot_rsiena_sim_stability()
# saomnkrsiena$local_search_rsiena_run(iterations=1000, rsiena_phase2_nsub=5)

# saom_nk_enhanced$local_search_rsiena(1)
saomnkrsiena$local_search_rsiena_run(iterations=100) 

saomnkrsiena$local_search_rsiena_run(iterations=10000)
saomnkrsiena$local_search_rsiena_run(iterations=100000)
# saomnkrsiena$local_search_rsiena_run(iterations=20, overwrite = F) 

saomnkrsiena2 <- SAOM_NK_RSiena$new(M = 15, N = 24, BI_PROB = .15, sim_name = '_TESTrsiena_')
# saom_nk_enhanced$local_search_rsiena(1)
saomnkrsiena2$local_search_rsiena_run(iterations=1000) 




########################################
self <- saomnkrsiena$clone(deep=T)
sims = self$rsiena_model$sims
n <- length(sims)
outlist <- list()
difflist <- list()
jaccardlist <- list()
K_env_list <- list()
K_soc_list <- list()
for(i in 1:length(sims)) {
  cat(sprintf(' %s ', i))
  ##
  el_bi_env <- sims[[ i ]][[1]][[3]]$`1`
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
                                                      directed = F, mode = 'all', 
                                                      multiple = T, weighted = T, 
                                                      add.names = T)
  projections <- igraph::bipartite.projection(new_bi_g, multiplicity = T, which = 'both')
  new_g_soc <- projections$proj1
  new_g_env <- projections$proj2

  K_soc_list[[i]] <- igraph::degree(new_g_soc)
  K_env_list[[i]] <- igraph::degree(new_g_env)
  
  
  # print(sim_i[[1]][[3]]`1`)
}

saomnkrsiena$plot_degree_progress_from_sims(K_soc_list = K_soc_list, K_env_list = K_env_list, plot_save = T)


# par(mfrow=c(1,3))
# stability_vec <- cumsum(plyr::ldply(jaccardlist)[,2])/(1:n)
# plot(stability_vec , type='l' , main='Change in Stability  [t-1, t]')
# 
# 
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



########################################

# sims = self$rsiena_model$sims
# n <- length(sims)
# outlist <- list()
# difflist <- list()
# jaccardlist <- list()
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
#   ##
#   bi_env_mat_new  <- as.matrix(bi_env_mat_sp)[ 1:self$M, (self$M+1):(MplusN) ]
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
#   # print(sim_i[[1]][[3]]`1`)
# }
# 
# par(mfrow=c(1,3))
# stability_vec <- cumsum(plyr::ldply(jaccardlist)[,2])/(1:n)
# plot(stability_vec , type='l' , main='Change in Stability  [t-1, t]')

# stability_delta <- stability_vec / (1:n)
# plot(stability_delta, type='o', log='y',
#      xlab='t', 
#      # ylim=c( stability_delta[n-1],  1 ),
#      ylab='Ln Stability Change [t-1, t]', main='Stabilization Rate (Ln Change in Stability)' 
# ); abline(h = tol, col='pink', lty=2)
# 
# burn_prop <- 0.2
# iter_postburn <- round(c( 1-burn_prop, burn_prop) * n )
# hist(stability_vec[ iter_postburn[1]:iter_postburn[2] ], main='Stability (post-burn)')


# ########################################################################
# sims = saomnkrsiena$rsiena_model$sims
# outlist <- list()
# difflist <- list()
# jaccardlist <- list()
# for(i in 1:length(sims)) {
#   cat(sprintf(' %s ', i))
#   ##
#   el_bi_env <- sims[[ i ]][[1]][[3]]$`1`
#   ## update numbering of second mode (the comonent integer names shift upward by the number of actors)
#   el_bi_env[,2] <- el_bi_env[,2] + saomnkrsiena$M
#   ##
#   MplusN <- saomnkrsiena$M + saomnkrsiena$N
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
#   ##
#   bi_env_mat_new  <- as.matrix(bi_env_mat_sp)[ 1:saomnkrsiena$M, (saomnkrsiena$M+1):(MplusN) ]
#   ##
#   # saomnkrsiena$plot_bipartite_system_from_mat(bi_env_mat_new, i, plot_save = TRUE) 
#   ##
#   outlist[[ sprintf('sim%d',i) ]] <- bi_env_mat_new
#   
#   if (i > 1) {
#     difflist[[ sprintf('diff%d-%d', i-1, i) ]] <-  bi_env_mat_new - outlist[[ (i-1) ]] 
#     
#     jaccardlist[[sprintf('jac%d-%d', i-1, i)]] <- get_jaccard_index(bi_env_mat_new, outlist[[ (i-1) ]] )
#   }
#   
#   # print(sim_i[[1]][[3]]`1`)
# }
# 
# par(mfrow=c(1,3))
# stability_vec <- cumsum(plyr::ldply(jaccardlist)[,2])/(1:length(jaccardlist))
# plot(stability_vec , type='l' )
# 
# stability_delta <- stability_vec / 1:length(stability_vec)
# plot(stability_delta, type='o', log='y',
#      xlab='t', ylab='Ln Stability Change [t-1, t]', main='Change in Stability' )
# 
# hist(stability_delta[100:length(stability_delta)])
# 
# ########################################################################


saom_nk_enhanced$local_search_saom_batchrun(batches = 2, iterations = 3, plot_save = T)
saom_nk_enhanced$plot_fitness_progress()
saom_nk_enhanced$plot_degree_progress()
# saom_nk_enhanced$visualize_networks()


saom_nk_enhanced$local_search_saom_batchrun(batches = 4, iterations = 10, 
                                            plot_save=TRUE, filename='__TEST_PLOT2__') 
# saom_nk_enhanced$visualize_networks(plot_save = T)



#### SANDBOX ########

saomnkrsiena$rsiena_model$sims[[1]][[1]][[2]]

M <- saomnkrsiena$M
N <- saomnkrsiena$N
actors <- 1:saomnkrsiena$M 
components <- 1:saomnkrsiena$N
el_bi_env <- saomnkrsiena$rsiena_model$sims[[1]][[1]][[3]]$`1`
el_bi_env[,2] <- el_bi_env[,2] + M
# bi_env_el <- el_bi_env[,1:2]
##
el_proj1  <- saomnkrsiena$rsiena_model$sims[[1]][[1]][[1]]$`1`
el_proj2  <- saomnkrsiena$rsiena_model$sims[[1]][[1]][[2]]$`1`
el_proj2 <- el_proj2 + M

## Undirected --> Upper right rectangle of full bipartite matrix
##   M N
## M[0,X] for X=[0,1]
## N[0,0]
bi_env_mat_sp <- sparseMatrix(i = el_bi_env[,1],
                              j = el_bi_env[,2],
                              x = el_bi_env[,3],
                              dims = c(M+N, M+N))
bi_env_mat <- as.matrix(bi_env_mat_sp)[1:M,(M+1):(M+N)]

# self$set_system_from_bipartite_igraph( self$random_bipartite_igraph() )


#####################



#########################################

## Simple environment 
env_simple <- SAOM_NK_Enhanced$new(M = 12, N = 8, BI_PROB = .15, sim_name = '_ENV_SIMPLE_')
env_simple$local_search_saom_batchrun(batches = 10, iterations = 5, plot_save=TRUE) 

## Mid environment 
env_mid <- SAOM_NK_Enhanced$new(M = 12, N = 12, BI_PROB = .15, sim_name = '_ENV_MID_')
env_mid$local_search_saom_batchrun(batches = 10, iterations = 5, plot_save=TRUE)  

## Complex environment 
env_complex <- SAOM_NK_Enhanced$new(M = 12, N = 18, BI_PROB = .15, sim_name = '_ENV_COMPLEX_')
env_complex$local_search_saom_batchrun(batches = 10, iterations = 5, plot_save=TRUE) 

## Very Complex environment 
env_verycomplex <- SAOM_NK_Enhanced$new(M = 12, N = 30, BI_PROB = .15, sim_name = '_ENV_VERYCOMPLEX_')
env_verycomplex$local_search_saom_batchrun(batches = 15, iterations = 5, plot_save=TRUE) 


#########################################







