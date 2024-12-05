
# Load necessary libraries
library(R6)
library(network)
library(igraph)
library(visNetwork)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggraph)
# Load the necessary library for grid.text
library(grid)






# Define the SAOM_NK_Enhanced class with local search, complex search strategies, network metrics, and visualization
SAOM_NK_Enhanced <- R6Class(
  "SAOM_NK_Enhanced",
  
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
    bipartite_network = NULL, # Bipartite network object (actors x components)
    social_network = NULL,   # Social space projection (network object)
    search_landscape = NULL, # Search space projection (network object)
    fitness_landscape = NULL, # Fitness landscape matrix
    progress_scores = NULL,  # Track scores over iterations
    P_change = NULL,          # Probability of changing the component interaction matrix
    fitness_history = list(),
    degree_history_K_S = list(),
    degree_history_K_E = list(),
    
    
    # Constructor to initialize the SAOM-NK model
    initialize = function(M, N, BI_PROB=0.2, sim_name = '_', P_change = 0.1, visualize_init=T) {
      self$M <- M 
      self$N <- N
      self$BI_PROB <- BI_PROB
      self$P_change <- P_change
      self$bipartite_network <- self$generate_bipartite_network()
      self$social_network <- self$project_social_space()
      self$search_landscape <- self$project_search_space()
      self$progress_scores <- c()  # Initialize empty vector to track scores
      self$ITERATION <- 0
      self$TIMESTAMP <- round( as.numeric(Sys.time())*100  )
      self$SIM_NAME <- sim_name
      #
      if (visualize_init)
        self$visualize_networks(plot_save = TRUE)
    },
    
    # Generate a random bipartite network object
    generate_bipartite_network = function() {
      set.seed(123)  # For reproducibility
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
    
    
    
    # Objective function that computes utility based on network statistics
    objective_func = function(actor_id, component_id) {
      # Calculate network statistics from bipartite, social, and component spaces
      
      # Bipartite network degree (e.g., number of connections for the actor)
      bipartite_matrix <- as.matrix.network(self$bipartite_network)
      actor_degree_bipartite <- sum(bipartite_matrix[actor_id, ])
      
      # Social network statistics (e.g., degree centrality)
      social_igraph <- graph_from_adjacency_matrix(as.matrix.network(self$social_network), mode = "undirected")
      degree_centrality_social <- degree(social_igraph)[actor_id]
      transitivity_social <- transitivity(social_igraph, type = "global")
      # cat(sprintf('\n trans = %s', transitivity_social))
      
      # Component interaction network statistics (e.g., degree of the component)
      component_igraph <- graph_from_adjacency_matrix(as.matrix.network(self$search_landscape), mode = "undirected")
      degree_centrality_component <- degree(component_igraph)[component_id]
      
      # Define utility function (example combining all three statistics)
      # utility <- (
      #     1.0 * ( round(self$N / 4) - actor_degree_bipartite )
      #   + 0.1 * ( round(self$M / 5) - degree_centrality_social )
      #   + 0.1 * ( round(self$N / 3) - degree_centrality_component )
      #   ## + 0.1 * transitivity_social
      # )
      utility <- (
        # 2 * ( round(self$N / 3) - actor_degree_bipartite )
           .1   *  actor_degree_bipartite 
        -  .5   *  degree_centrality_social 
        +  .2   *  degree_centrality_component 
        +  .1   *  transitivity_social
      )
      return(utility)
    },
    
    # Perform local search using SAOM objective function for tie changes with normalization over all configurations
    local_search_saom = function(iterations) {
      for (iter in 1:iterations) {
        
        for (actor_id in 1:self$M) {
          # Compute the sum of exponential utilities for all single-component changes
          bipartite_matrix <- as.matrix.network(self$bipartite_network)
          exp_utilities <- numeric(self$N)
          
          for (component_id in 1:self$N) {
            # Temporarily flip the current tie for evaluation
            test_matrix <- bipartite_matrix
            test_matrix[actor_id, component_id] <- 1 - test_matrix[actor_id, component_id]  # Flip the tie
            
            # Temporarily update the bipartite network for utility calculation
            self$bipartite_network <- network(test_matrix, bipartite = TRUE, directed = FALSE)
            self$social_network <- self$project_social_space()  # Update social network projection
            self$search_landscape <- self$project_search_space()  # Update component interaction projection
            
            # Calculate the utility for the current flipped state using objective_func
            exp_utilities[component_id] <- exp(self$objective_func(actor_id, component_id))
          }
          
          # Sum the exponential utilities for normalization
          total_exp_utilities <- sum(exp_utilities)
          
          # Iterate through component changes for the current actor to decide tie changes
          for (component_id in 1:self$N) {
            
            ## update simulation progress iteration count
            self$increment_sim_iter()
            
            # Flip the tie to evaluate the new state
            new_matrix <- bipartite_matrix
            current_tie <- new_matrix[actor_id, component_id]
            new_matrix[actor_id, component_id] <- 1 - current_tie  # Flip the tie
            
            # Update the bipartite network with the modified matrix
            self$bipartite_network <- network(new_matrix, bipartite = TRUE, directed = FALSE)
            self$social_network <- self$project_social_space()  # Update social network projection
            self$search_landscape <- self$project_search_space()  # Update component interaction projection
            
            # Calculate the new utility using objective_func
            new_payoff <- self$objective_func(actor_id, component_id)
            
            # Normalize prob_accept over all possible configurations
            prob_accept <- exp(new_payoff) / total_exp_utilities
            
            cat(sprintf('\nactor %s, component %s, prob_acct = %s', actor_id, component_id, prob_accept))
            
            if (runif(1) > prob_accept) {
              # Revert the change if not accepted
              new_matrix[actor_id, component_id] <- current_tie # return to current tie value
              self$bipartite_network <- network(new_matrix, bipartite = TRUE, directed = FALSE)
              self$social_network <- self$project_social_space()  # Revert social network projection
              self$search_landscape <- self$project_search_space()  # Revert component interaction projection
            } else {
              cat(sprintf("\nIteration %d, Actor %d, Component %d: Tie change accepted with new payoff = %.3f\n",
                          iter, actor_id, component_id, new_payoff))
            }
            
            ## UPDATE SIM PROGRESS
            self$record_fitness()
            self$record_Ks()
          }
        }
      }
    },
    
    # Batch runs with intermittent snapshots (network plots) progress 
    local_search_saom_batchrun = function(batches, batch_iters, plot_save=TRUE) {
      
      for (batch_iter in 1:batches) {
        cat(sprintf('\n batch %s\n', batch_iter))
        self$local_search_saom(iterations = batch_iters)
        #
        self$visualize_networks(plot_save = plot_save)
        ##
        self$plot_fitness_progress(plot_save = plot_save)
        #
        self$plot_degree_progress(plot_save = plot_save)
      }
      
    },
    
    
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
    
    
    
    # Complex search strategies: hill climbing
    hill_climbing_search = function(actor_id, component_id, iterations) {
      current_score <- self$evaluate_interspace_interactions()
      self$progress_scores <- c(current_score)  # Track the initial score
      
      for (i in 1:iterations) {
        bipartite_matrix <- as.matrix.network(self$bipartite_network)
        
        current_tie <- bipartite_matrix[actor_id, component_id]
        new_tie <- ifelse(current_tie == 1, 0, 1)
        bipartite_matrix[actor_id, component_id] <- new_tie
        
        self$bipartite_network <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
        self$social_network <- self$project_social_space()
        self$search_landscape <- self$project_search_space()
        
        new_score <- self$evaluate_interspace_interactions()
        
        if (new_score > current_score) {
          current_score <- new_score
          cat("Iteration", i, ": Improved interaction score to", current_score, "\n")
        } else {
          bipartite_matrix[actor_id, component_id] <- current_tie
          self$bipartite_network <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
        }
        
        self$progress_scores <- c(self$progress_scores, current_score)
      }
    },
    
    # Calculate network metrics using igraph
    calculate_network_metrics = function() {
      ig_social <- graph_from_adjacency_matrix(as.matrix.network(self$social_network), mode = "undirected")
      degree_centrality <- degree(ig_social)
      betweenness_centrality <- betweenness(ig_social)
      eigenvector_centrality <- eigen_centrality(ig_social)$vector
      
      cat("Network Metrics:\n")
      cat("Degree Centrality:\n")
      print(degree_centrality)
      cat("Betweenness Centrality:\n")
      print(betweenness_centrality)
      cat("Eigenvector Centrality:\n")
      print(eigenvector_centrality)
      
      par(mfrow=c(2,2))
      hist(degree_centrality, main = 'Degree Centrality')
      hist(betweenness_centrality, main = 'Betweenness Centrality')
      hist(eigenvector_centrality, main = 'Eigenvector Centrality')
      
      return(list(
        degree = degree_centrality,
        betweenness = betweenness_centrality,
        eigenvector = eigenvector_centrality
      ))
    },
    
    
    
    # Visualization using ggplot2 and ggraph: bipartite network, social network, and component interaction matrix heatmap
    visualize_networks = function(plot_save = FALSE) {
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
      ig_bipartite <- graph_from_incidence_matrix(as.matrix.network(self$bipartite_network))
      
      # Set node attributes for shape and label
      V(ig_bipartite)$type <- bipartite_mapping(ig_bipartite)$type
      V(ig_bipartite)$shape <- ifelse(V(ig_bipartite)$type, "square", "circle")
      V(ig_bipartite)$color <- ifelse(V(ig_bipartite)$type, "lightblue", "orange")
      V(ig_bipartite)$label <- ifelse(V(ig_bipartite)$type, component_labels, 1:self$M)  # Numeric labels for actors
      
      bipartite_plot <- ggraph(ig_bipartite, layout = "bipartite") +
        geom_edge_link(color = "gray") +
        geom_node_point(aes(shape = shape, color = color), size = 5) +
        geom_node_text(aes(label = label), vjust = 0.5, hjust = 0.5, size = 3) +  # Center vertex labels over vertices
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
      
      social_plot <- ggraph(ig_social, layout = "fr") +
        geom_edge_link(color = "gray") +
        geom_node_point(aes(size = node_size, color = node_color)) +
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
                            self$ITERATION, self$N, self$M, self$BI_PROB, avg_degree_social, avg_degree_component)
      
      # Arrange plots with a main title
      (plt <- grid.arrange(
        social_plot, bipartite_plot, heatmap_plot, ncol = 3,
        top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold"))
      ))
      if(plot_save){
        keystring <- sprintf("%s_sim%.0f_iter%.0f_N%d_M%d_BI_PROB_%.2f_K_S_%.2f_K_C_%.2f", 
                             self$SIM_NAME, self$TIMESTAMP, self$ITERATION, self$N, self$M, self$BI_PROB, avg_degree_social, avg_degree_component)
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


saom_nk_enhanced <- SAOM_NK_Enhanced$new(M = 12, N = 12, BI_PROB = .15, sim_name = '_TEST_')
saom_nk_enhanced$local_search_saom_batchrun(batches = 2, batch_iters = 3, plot_save = T)
saom_nk_enhanced$plot_fitness_progress()
saom_nk_enhanced$plot_degree_progress()
# saom_nk_enhanced$visualize_networks()


saom_nk_enhanced$local_search_saom_batchrun(batches = 4, batch_iters = 10, 
                                            plot_save=TRUE, filename='__TEST_PLOT2__') 
# saom_nk_enhanced$visualize_networks(plot_save = T)





#########################################

## Simple environment 
env_simple <- SAOM_NK_Enhanced$new(M = 12, N = 8, BI_PROB = .15, sim_name = '_ENV_SIMPLE_')
env_simple$local_search_saom_batchrun(batches = 10, batch_iters = 5, plot_save=TRUE) 

## Mid environment 
env_mid <- SAOM_NK_Enhanced$new(M = 12, N = 12, BI_PROB = .15, sim_name = '_ENV_MID_')
env_mid$local_search_saom_batchrun(batches = 10, batch_iters = 5, plot_save=TRUE)  

## Complex environment 
env_complex <- SAOM_NK_Enhanced$new(M = 12, N = 18, BI_PROB = .15, sim_name = '_ENV_COMPLEX_')
env_complex$local_search_saom_batchrun(batches = 10, batch_iters = 5, plot_save=TRUE) 

## Very Complex environment 
env_verycomplex <- SAOM_NK_Enhanced$new(M = 12, N = 30, BI_PROB = .15, sim_name = '_ENV_VERYCOMPLEX_')
env_verycomplex$local_search_saom_batchrun(batches = 15, batch_iters = 5, plot_save=TRUE) 


#########################################


# 
# 
# # Run the local search method
# saom_nk_enhanced$local_search_saom(iterations = 1)
# saom_nk_enhanced$visualize_networks()
# 
# saom_nk_enhanced$local_search_saom(iterations = 1)
# saom_nk_enhanced$visualize_networks()
# 
# saom_nk_enhanced$local_search_saom(iterations = 1)
# saom_nk_enhanced$visualize_networks()
# 
# # Calculate and print network metrics
# network_metrics <- saom_nk_enhanced$calculate_network_metrics()
# 
# # Plot the search progress
# # saom_nk_enhanced$plot_search_progress()
# saom_nk_enhanced$visualize_networks()
# 
# 
# 
# saom_nk_enhanced$local_search_saom(iterations = 5)
# saom_nk_enhanced$visualize_networks()
# 
# 
# 
# saom_nk_enhanced$local_search_saom(iterations = 5)
# saom_nk_enhanced$visualize_networks()



















# # Example usage of the SAOM-NK model with enhanced local search
# saom_nk_enhanced <- SAOM_NK_Enhanced$new(M = 10, N = 15, BI_PROB = .15)
# saom_nk_enhanced$visualize_networks()
# 
# 
# # Run the local search method
# saom_nk_enhanced$local_search(iterations = 10, mutation_rate = 0.01)
# 
# # Calculate and print network metrics
# network_metrics <- saom_nk_enhanced$calculate_network_metrics()
# # Plot the search progress
# saom_nk_enhanced$plot_search_progress()
# 
# saom_nk_enhanced$visualize_networks()








