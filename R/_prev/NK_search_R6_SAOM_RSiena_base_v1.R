
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




  



# Define the SAOM_NK_Enhanced class with local search, complex search strategies, network metrics, and visualization
SAOM_NK_RSiena <- R6Class(
  "SAOM_NK_Rsiena",
  
  public = list(
    #
    ITERATION = NULL,
    TIMESTAMP = NULL,
    SIM_NAME = NULL,
    #
    M = NULL,                # Number of actors
    N = NULL,                # Number of components
    BI_PROB = NULL,          # probability of tie in bipartite network dyads (determines bipartite density)
    ## K = NULL,             # Avg. degree of component-[actor]-component interactions
    bipartite_network = NULL, # Bipartite network object 
    social_network = NULL,   # Social space projection (network object)
    search_landscape = NULL, # Search space projection (network object)
    #
    bipartite_sienaDV = NULL,
    social_sienaDV = NULL,
    search_sienaDV = NULL,
    #
    fitness_landscape = NULL, # Fitness landscape matrix
    progress_scores = NULL,  # Track scores over iterations
    P_change = NULL,          # Probability of changing the component interaction matrix
    #
    fitness_history = list(),
    degree_history_K_S = list(),
    degree_history_K_E = list(),
    #
    siena_model = NULL,   
    siena_data = NULL,
    siena_effects = NULL,
    siena_algorithm = NULL,

    # Constructor to initialize the SAOM-NK model
    initialize = function(M, N, BI_PROB=0.2, sim_name = '_', P_change = 0.1, visualize_init=T) {
      ## ----- prevent clashes with sna package-----------
      sna_err_check <-tryCatch(expr = { detach('package:sna') }, error=function(e)e )
      ## -------------------------------------------------
      self$M <- M 
      self$N <- N
      self$BI_PROB <- BI_PROB
      self$P_change <- P_change
      #
      # self$bipartite_network <- self$generate_bipartite_network()
      # self$social_network <- self$project_social_space()
      # self$search_landscape <- self$project_search_space()
      self$set_system_from_bipartite_network( self$random_bipartite_network() )
      #
      self$progress_scores <- c()  # Initialize empty vector to track scores
      self$ITERATION <- 0
      self$TIMESTAMP <- round( as.numeric(Sys.time())*100  )
      self$SIM_NAME <- sim_name
      #
      if (visualize_init)
        self$visualize_networks_rsiena(plot_save = TRUE)
    },
    
    #
    set_system_from_bipartite_network = function(bipartite_network) {
      if( ! 'network' %in% class(bipartite_network)) 
        stop(sprintf('\nbipartite_network is not a network object; class %s\n', class(bipartite_network)))
      self$bipartite_network <- bipartite_network
      #
      self$social_network <- self$project_social_space()
      self$search_landscape <- self$project_search_space()
    },
    
    # Generate a random bipartite network object
    random_bipartite_network = function(rand_seed = 123) {
      set.seed(rand_seed)  # For reproducibility
      probs <- c( 1 - self$BI_PROB, self$BI_PROB )
      bipartite_matrix <- matrix(sample(0:1, self$M * self$N, replace = TRUE, prob = probs ),
                                 nrow = self$M, ncol = self$N)
      bipartite_net <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
      return(bipartite_net)
    },
    
    # Project the social space (actor network) from the bipartite network
    project_social_space = function() {
      bipartite_matrix <- as.matrix.network(self$bipartite_network)
      actor_matrix <- bipartite_matrix %*% t(bipartite_matrix)
      diag(actor_matrix) <- 0  # Remove self-loops
      social_net <- network(actor_matrix, directed = FALSE)
      return(social_net)
    },
    
    # Project the search landscape space (component interaction network)
    project_search_space = function() {
      bipartite_matrix <- as.matrix.network(self$bipartite_network)
      component_matrix <- t(bipartite_matrix) %*% bipartite_matrix
      diag(component_matrix) <- 0  # Remove self-loops
      search_net <- network(component_matrix, directed = FALSE)
      return(search_net)
    },
    
    get_social_igraph = function() {
      return(graph_from_adjacency_matrix(as.matrix.network(self$social_network), mode = "undirected"))
    },
    
    get_component_igraph = function() {
      return(graph_from_adjacency_matrix(as.matrix.network(self$search_landscape), mode = "undirected"))
    },
    
    # Update Iteration progress
    increment_sim_iter = function(val = 1) {
      self$ITERATION <- self$ITERATION + val
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
    init_rsiena_model_from_bipartite_matrix = function() {
      #
      bipartite_matrix <- as.matrix.network(self$bipartite_network)
      
      
      ACTORS     <- sienaNodeSet(self$M, nodeSetName="ACTORS")
      COMPONENTS <- sienaNodeSet(self$N, nodeSetName="COMPONENTS")
      
      ## init networks to duplicate for the init arrays (two network waves)
      array_bi_net <- array(c(bipartite_matrix, bipartite_matrix), 
                            dim=c(self$M, self$N, 2) )
      array_social <- array(c(self$social_network[,], self$social_network[,]), 
                            dim=c(self$M, self$M, 2) )
      array_search <- array(c(self$search_landscape[,], self$search_landscape[,]), 
                            dim=c(self$N, self$N, 2) )
      #
      self$bipartite_sienaDV <- sienaDependent(array_bi_net, type='bipartite', nodeSet =c('ACTORS', 'COMPONENTS'), allowOnly = F)
      self$social_sienaDV    <- sienaDependent(array_social, type='oneMode', nodeSet = 'ACTORS', allowOnly = F)
      self$search_sienaDV    <- sienaDependent(array_search, type='oneMode', nodeSet = 'COMPONENTS', allowOnly = F)
      #

      ###
      self$siena_data <- sienaDataCreate(list(self$social_sienaDV,
                                              self$search_sienaDV,
                                              self$bipartite_sienaDV),
                                         nodeSets = list(ACTORS, COMPONENTS))
      ###
      # self$siena_data <- sienaDataCreate(list(self$bipartite_sienaDV), 
      #                                    nodeSets = list(ACTORS, COMPONENTS))
      
    },
    
  ##
  ##**TODO**
  ##**Create custom RSiena interaction functions for only bipartite DV, **
  ##**but objective function includes statistics of the projections (social net, search landscape)**
  ##
    
    # RSIENA
    local_search_rsiena_init = function(iterations, show_effects_doc=FALSE) {
      
      self$init_rsiena_model_from_bipartite_matrix()

      # Step 3: Define the effects for the model
      self$siena_effects <- getEffects(self$siena_data)
      
      # Effects Documentation
      if(show_effects_doc)
        effectsDocumentation(self$siena_effects)
      ##-----------------------------
      ##-----------------------------
      # bipartite effects
      self$siena_effects <- includeEffects(self$siena_effects, 
                                           density, 
                                           name = 'self$bipartite_sienaDV')
      self$siena_effects <- setEffect(self$siena_effects, density, 
                                      name = 'self$bipartite_sienaDV', initialValue = -1, fix=T)  ## effect parameter vs. initial value ?
      
      # 
      self$siena_effects <- includeEffects(self$siena_effects, 
                                           inPop, 
                                           name = 'self$bipartite_sienaDV')
      self$siena_effects <- setEffect(self$siena_effects, inPop, 
                                      name = 'self$bipartite_sienaDV', initialValue = .1)  ## effect parameter vs. initial value ?
      # 
      self$siena_effects <- includeEffects(self$siena_effects, 
                                           outAct, 
                                           name = 'self$bipartite_sienaDV')
      self$siena_effects <- setEffect(self$siena_effects, outAct, 
                                      name = 'self$bipartite_sienaDV', initialValue = .1)  ## effect parameter vs. initial value ?
      # 
      self$siena_effects <- includeEffects(self$siena_effects, 
                                           cycle4, 
                                           name = 'self$bipartite_sienaDV')
      self$siena_effects <- setEffect(self$siena_effects, cycle4, 
                                      name = 'self$bipartite_sienaDV', initialValue = -.1, fix=T)  ## effect parameter vs. initial value ?
      #-------------------------------
      ##-----------------------------
      ## SOCIAL SPACE
      self$siena_effects <- includeEffects(self$siena_effects, 
                                           inPop, outAct, include = FALSE, 
                                           name = 'self$social_sienaDV')
      ## SOcial density
      # self$siena_effects <- includeEffects(self$siena_effects,
      #                                      density,
      #                                      name = 'self$social_sienaDV')
      self$siena_effects <- setEffect(self$siena_effects, density,
                                      name = 'self$social_sienaDV', initialValue = -0.2, fix = T)  ## effect parameter vs. initial value ?
      # ## Social transitivity
      # self$siena_effects <- includeEffects(self$siena_effects,
      #                                      transTriads,
      #                                      name = 'self$social_sienaDV')
      # self$siena_effects <- setEffect(self$siena_effects, transTriads,
      #                                 name = 'self$social_sienaDV', initialValue = 0.3, fix = T)  ## effect parameter vs. initial value ?
      
      ##
      ##-----------------------------
      # # Search landscape DV interactions with bipartite     
      # self$siena_effects <- includeEffects(self$siena_effects, 
      #                                      from.w.ind, 
      #                                      name = 'self$search_sienaDV',
      #                                      interaction1 = "self$bipartite_sienaDV",
      #                                      initialValue = .5)
      # self$siena_effects <- includeEffects(self$siena_effects, 
      #                                      density, 
      #                                      name = 'self$search_sienaDV')
      self$siena_effects <- setEffect(self$siena_effects, density,
                                      name = 'self$search_sienaDV', initialValue = 0.1, fix = T)  ## effect parameter vs. initial value ?
      ##-----------------------------
      ##-----------------------------
      
      ## RSiena Algorithm
      self$siena_algorithm <- sienaAlgorithmCreate(projname=sprintf('%s_%s',self$SIM_NAME,self$TIMESTAMP),
                                                   simOnly = T,
                                                   n3 = iterations)
      
      # Step 4: Run RSiena simulation
      self$siena_model <- siena07(self$siena_algorithm, data = self$siena_data, effects = self$siena_effects,
                                  batch = TRUE, returnDeps = TRUE)
      
      # Step 5: Summarize and plot results
      mod_summary <- summary(self$siena_model)
      if(!is.null(mod_summary)) 
        print(mod_summary)
      # plot(self$siena_model)
      
      print(screenreg(list(self$siena_model), single.row = T, digits = 3))
      
      ## update simulation object environment from 
      new_bi_env_mat <- self$get_bipartite_network_from_siena_sim()
      new_bi_env_net <- network(new_bi_env_mat, bipartite = TRUE, directed = FALSE)
      self$set_system_from_bipartite_network( new_bi_env_net )
    
    },  
  
    #
    get_bipartite_network_from_siena_sim = function(sim_iteration = NULL) {
      ######### UPDATE SYSTEM ENVIRONMENT SNAPSHOT FROM EVOLVED Bipartite Network DV #####################
      # sim_id_last <- length( self$siena_model$sims ) ##
      sim_id <- ifelse(is.null(sim_iteration), 
                       length(self$siena_model$sims), ## default current state is the last simulation in sims list
                       sim_iteration)
      ##
      # actors <- 1:self$M 
      # components <- 1:self$N
      MplusN <- self$M + self$N
      el_bi_env <- self$siena_model$sims[[ sim_id ]][[1]][[3]]$`1`
      ## update numbering of second mode (the comonent integer names shift upward by the number of actors)
      el_bi_env[,2] <- el_bi_env[,2] + M
      ### get networks from other prjected space ties 
      # el_proj1  <- self$siena_model$sims[[1]][[1]][[1]]$`1` ## Social network
      # el_proj2  <- self$siena_model$sims[[1]][[1]][[2]]$`1` ## 
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
    get_siena_model_snapshots = function(n=4, return_list=FALSE) {
      #
      snapshot_sim_ids <- seq(0,1,length.out= n) * length(self$siena_model$sims)
      snapshot_sim_ids[1] <- 1
      ##
      MplusN <- self$M + self$N
      ##
      snapshots <- list()
      for (sim_id in snapshot_sim_ids) {
        ##
        el_bi_env <- self$siena_model$sims[[ sim_id ]][[1]][[3]]$`1`
        ## update numbering of second mode (the comonent integer names shift upward by the number of actors)
        el_bi_env[,2] <- el_bi_env[,2] + M
        #############
        ## Undirected --> Upper right rectangle of full bipartite matrix
        ##   M N
        ## M[0,X], for X in [0,1]
        ## N[0,0]
        bi_env_mat_sp <- sparseMatrix(i = el_bi_env[,1],
                                      j = el_bi_env[,2],
                                      x = el_bi_env[,3],
                                      dims = c(MplusN, MplusN))
        snapshots[[ sprintf('sim%d',sim_id) ]] <- as.matrix(bi_env_mat_sp)[ 1:self$M, (self$M+1):(MplusN) ]
        
      }
     
      ##
      if(return_list) 
        return( snapshots )
      
    },
    
    # Extend existing simulation environment
    local_search_rsiena_extend = function(iterations) {
      self$siena_algorithm$n3 <- iterations
      self$siena_model <- siena07(self$siena_algorithm, data = self$siena_data, effects = self$siena_effects,
                                  prevAns = self$siena_model,
                                  batch = TRUE, returnDeps = TRUE)
      ##
      mod_summary <- summary(self$siena_model)
      if(!is.null(mod_summary)) 
        print(mod_summary)
      ##
      print(screenreg(list(self$siena_model), single.row = T, digits = 3))
    },
  

    ## Single Simulation Run
    local_search_rsiena_run = function(iterations=100, plot_save = TRUE, overwrite=TRUE) {
      if( overwrite | is.null(self$siena_model) ) {
        self$local_search_rsiena_init(iterations)
      } else {}
      # self$visualize_networks_rsiena(plot_save = plot_save)
      self$get_siena_model_snapshots(n=5)
    },

    ## Batch of Multiple Simulation Extensions to 
    local_search_rsiena_batchrun = function(batches=10, iterations=100, plot_save = TRUE, overwrite=FALSE) {
      
      for (i in 1:batches) {
      
        overwrite_by_batch <- ifelse(i == 1, TRUE, overwrite)
        if(overwrite_by_batch) {
          self$local_search_rsiena_extend(iterations)
        } else {
          self$local_search_rsiena_init(iterations)
        }

       
        #siena_run(iterations=iterations, plot_save=plot_save, overwrite=overwrite_by_batch)
        
        ## Take snapshot of evolving system at this iteration of progress
        self$visualize_networks_rsiena(plot_save = plot_save)
        
      }
      
    },
    
    
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
    #   basic rate parameter self$social_sienaDV          0.000 (0.000) ***
    #   self$social_sienaDV: degree (density)            -2.000 (0.000) *** [suppressed K_S]
    #   basic rate parameter self$search_sienaDV          0.045 (0.000) ***
    #   self$search_sienaDV: degree (density)            -0.557 (0.000) *** [suppressed K_E]
    #   basic rate parameter self$bipartite_sienaDV       0.000 (0.000) ***
    #   self$bipartite_sienaDV: outdegree (density)      57.398 (0.000) ***
    #   self$bipartite_sienaDV: indegree - popularity     8.086 (0.000) ***
    #   self$bipartite_sienaDV: outdegree - activity     -2.389 (0.000) ***
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
    visualize_networks_rsiena = function(plot_save = FALSE) {
      #
      RSIENA_ITERATION <- length(self$siena_model$sims)
      
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
      ig_bipartite <- graph_from_biadjacency_matrix(as.matrix.network(self$bipartite_network))
      
      # Set node attributes for shape and label
      V(ig_bipartite)$type <- bipartite_mapping(ig_bipartite)$type
      V(ig_bipartite)$shape <- ifelse(V(ig_bipartite)$type, "square", "circle")
      V(ig_bipartite)$color <- ifelse(V(ig_bipartite)$type, "lightblue", "darkorange")
      V(ig_bipartite)$label <- ifelse(V(ig_bipartite)$type, component_labels, 1:self$M)  # Numeric labels for actors
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
      ig_social <- self$get_social_igraph()
      
      # Calculate network statistics
      node_size <- degree(ig_social)  # Degree centrality for node size
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
      component_matrix <- as.matrix.network(self$search_landscape)
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
      bipartite_matrix <- as.matrix.network(self$bipartite_network)
      actor_influence_on_components <- bipartite_matrix %*% as.matrix.network(self$search_landscape)
      total_interaction_score <- sum(actor_influence_on_components)
      return(total_interaction_score)
    }
  )
)


saomnkrsiena <- SAOM_NK_RSiena$new(M = 8, N = 30, BI_PROB = .15, sim_name = '_TESTrsiena_')
# saom_nk_enhanced$local_search_rsiena(1)
saomnkrsiena$local_search_rsiena_run(iterations=100) 
saomnkrsiena$local_search_rsiena_run(iterations=1000) 
saomnkrsiena$local_search_rsiena_run(iterations=10000)
saomnkrsiena$local_search_rsiena_run(iterations=100000)
# saomnkrsiena$local_search_rsiena_run(iterations=20, overwrite = F) 

saomnkrsiena2 <- SAOM_NK_RSiena$new(M = 15, N = 24, BI_PROB = .15, sim_name = '_TESTrsiena_')
# saom_nk_enhanced$local_search_rsiena(1)
saomnkrsiena2$local_search_rsiena_run(iterations=1000) 


###


saom_nk_enhanced$local_search_saom_batchrun(batches = 2, iterations = 3, plot_save = T)
saom_nk_enhanced$plot_fitness_progress()
saom_nk_enhanced$plot_degree_progress()
# saom_nk_enhanced$visualize_networks()


saom_nk_enhanced$local_search_saom_batchrun(batches = 4, iterations = 10, 
                                            plot_save=TRUE, filename='__TEST_PLOT2__') 
# saom_nk_enhanced$visualize_networks(plot_save = T)



#### SANDBOX ########

saomnkrsiena$siena_model$sims[[1]][[1]][[2]]

M <- saomnkrsiena$M
N <- saomnkrsiena$N
actors <- 1:saomnkrsiena$M 
components <- 1:saomnkrsiena$N
el_bi_env <- saomnkrsiena$siena_model$sims[[1]][[1]][[3]]$`1`
el_bi_env[,2] <- el_bi_env[,2] + M
# bi_env_el <- el_bi_env[,1:2]
##
el_proj1  <- saomnkrsiena$siena_model$sims[[1]][[1]][[1]]$`1`
el_proj2  <- saomnkrsiena$siena_model$sims[[1]][[1]][[2]]$`1`
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

# self$set_system_from_bipartite_network( self$random_bipartite_network() )


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







