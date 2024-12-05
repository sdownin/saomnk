
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
    M = NULL,                # Number of actors
    N = NULL,                # Number of components
    BI_PROB = NULL,          # probability of tie in bipartite network dyads (determines bipartite density)
    ## K = NULL,                # Avg. degree of component-[actor]-component interactions
    bipartite_network = NULL, # Bipartite network object (actors x components)
    social_network = NULL,   # Social space projection (network object)
    search_landscape = NULL, # Search space projection (network object)
    fitness_landscape = NULL, # Fitness landscape matrix
    progress_scores = NULL,  # Track scores over iterations
    P_change = NULL,          # Probability of changing the component interaction matrix
    
    # Constructor to initialize the SAOM-NK model
    initialize = function(M, N, BI_PROB=0.2, P_change = 0.1) {
      self$M <- M
      self$N <- N
      self$BI_PROB <- BI_PROB
      self$P_change <- P_change
      self$bipartite_network <- self$generate_bipartite_network()
      self$social_network <- self$project_social_space()
      self$search_landscape <- self$project_search_space()
      self$progress_scores <- c()  # Initialize empty vector to track scores
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
           1 *  actor_degree_bipartite 
        +  1 *  degree_centrality_social 
        +  1 *  degree_centrality_component 
        # -  0.1 *  transitivity_social
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
              new_matrix[actor_id, component_id] <- current_tie
              self$bipartite_network <- network(new_matrix, bipartite = TRUE, directed = FALSE)
              self$social_network <- self$project_social_space()  # Revert social network projection
              self$search_landscape <- self$project_search_space()  # Revert component interaction projection
            } else {
              cat(sprintf("Iteration %d, Actor %d, Component %d: Tie change accepted with new payoff = %.3f\n",
                          iter, actor_id, component_id, new_payoff))
            }
          }
        }
      }
    },
    # 
    # local_search_saom_runbatch = function(batch_iters, run_iters, plot, file) {
    # },
    
  
    # # Perform local search using SAOM objective function for tie changes with probability-based matrix updates
    # local_search_saom = function(iterations) {
    #   for (iter in 1:iterations) {
    #     for (actor_id in 1:self$M) {
    #       for (component_id in 1:self$N) {
    #         # Calculate current utility
    #         utility_current <- self$objective_func(actor_id, component_id)
    # 
    #         # Flip the tie to evaluate the new state
    #         bipartite_matrix <- as.matrix.network(self$bipartite_network)
    #         current_tie <- bipartite_matrix[actor_id, component_id]
    #         bipartite_matrix[actor_id, component_id] <- ifelse(current_tie == 1, 0, 1)
    # 
    #         # Update the bipartite network with the modified matrix
    #         self$bipartite_network <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
    # 
    #         # With probability P_change, update the component interaction matrix
    #         if (runif(1) < self$P_change) {
    #           self$search_landscape <- self$project_search_space()  # Update the component interaction projection
    #         }
    # 
    #         self$social_network <- self$project_social_space()  # Always update the social network projection
    # 
    #         # Calculate new utility
    #         utility_new <- self$objective_func(actor_id, component_id)
    # 
    #         # Calculate probability of accepting the new tie state
    #         prob_accept <- exp(utility_new) / (exp(utility_current) + exp(utility_new))
    # 
    #         if (runif(1) > prob_accept) {
    #           # Revert the change if not accepted
    #           bipartite_matrix[actor_id, component_id] <- current_tie
    #           self$bipartite_network <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
    #           self$social_network <- self$project_social_space()  # Revert social network projection
    # 
    #           # Revert the component interaction matrix if it was updated
    #           if (runif(1) < self$P_change) {
    #             self$search_landscape <- self$project_search_space()  # Revert component interaction projection
    #           }
    #         } else {
    #           cat(sprintf("Iteration %d, Actor %d, Component %d: Tie change accepted with utility = %.3f\n",
    #                       iter, actor_id, component_id, utility_new))
    #         }
    #       }
    #     }
    #   }
    # } ,
    
    
    # # Perform local search using SAOM objective function for tie changes with full normalization
    # local_search_saom = function(iterations) {
    #   for (iter in 1:iterations) {
    #     for (actor_id in 1:self$M) {
    #       # Enumerate all possible configurations for the actor
    #       num_configs <- 2^self$N
    #       all_configs <- expand.grid(rep(list(c(0, 1)), self$N))  # Generate all possible combinations
    #       utility_values <- numeric(num_configs)
    #       
    #       # Calculate utility for each configuration
    #       for (config_idx in 1:num_configs) {
    #         # Apply configuration to the actor's row in the bipartite matrix
    #         temp_matrix <- as.matrix.network(self$bipartite_network)
    #         temp_matrix[actor_id, ] <- as.numeric(all_configs[config_idx, ])
    #         temp_network <- network(temp_matrix, bipartite = TRUE, directed = FALSE)
    #         
    #         # Temporarily update the bipartite network for utility calculation
    #         self$bipartite_network <- temp_network
    #         self$social_network <- self$project_social_space()  # Update social network projection
    #         self$search_landscape <- self$project_search_space()  # Update component interaction projection
    #         
    #         # Store utility value for this configuration
    #         utility_values[config_idx] <- self$objective_func(actor_id, config_idx)
    #       }
    #       
    #       # Sum the exponential utilities for normalization
    #       exp_utilities <- exp(utility_values)
    #       total_exp_utilities <- sum(exp_utilities)
    #       
    #       # Iterate through component changes for the current actor
    #       for (component_id in 1:self$N) {
    #         # Calculate current utility
    #         current_utility <- self$objective_func(actor_id, component_id)
    #         
    #         # Flip the tie to evaluate the new state
    #         bipartite_matrix <- as.matrix.network(self$bipartite_network)
    #         current_tie <- bipartite_matrix[actor_id, component_id]
    #         bipartite_matrix[actor_id, component_id] <- ifelse(current_tie == 1, 0, 1)
    #         
    #         # Update the bipartite network with the modified matrix
    #         self$bipartite_network <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
    #         self$social_network <- self$project_social_space()  # Update social network projection
    #         self$search_landscape <- self$project_search_space()  # Update component interaction projection
    #         
    #         # Calculate new utility
    #         new_utility <- self$objective_func(actor_id, component_id)
    #         
    #         # Normalize prob_accept over all possible configurations
    #         prob_accept <- exp(new_utility) / total_exp_utilities
    #         
    #         prob_exists <- ( !is.null(prob_accept) & !is.nan(prob_accept) & !is.na(prob_accept) )
    #         
    #         if ( prob_exists &  runif(1) > prob_accept ) {
    #           # Revert the change if not accepted
    #           bipartite_matrix[actor_id, component_id] <- current_tie
    #           self$bipartite_network <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
    #           self$social_network <- self$project_social_space()  # Revert social network projection
    #           self$search_landscape <- self$project_search_space()  # Revert component interaction projection
    #         } else {
    #           cat(sprintf("Iteration %d, Actor %d, Component %d: Tie change accepted with utility = %.3f\n",
    #                       iter, actor_id, component_id, new_utility))
    #         }
    #       }
    #     }
    #   }
    # }    ,
    
    # Visualization of the search progress
    plot_search_progress = function() {
      if (length(self$progress_scores) == 0) {
        stop("No search progress data available. Run a search method first.")
      }
      plot.new()
      plot(self$progress_scores, type = "l", col = "blue", lwd = 2,
           main = "Search Progress",
           xlab = "Iteration", ylab = "Interaction Score",
           ylim = range(self$progress_scores))
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
    visualize_networks = function() {
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
        labs(title = "Bipartite Network\n(Actors and Components)") +  # Split title into two lines
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
      ig_social <- graph_from_adjacency_matrix(as.matrix.network(self$social_network), mode = "undirected")
      
      # Calculate network statistics
      node_size <- degree(ig_social)  # Degree centrality for node size
      node_color <- eigen_centrality(ig_social)$vector  # Eigenvector centrality for node color
      
      social_plot <- ggraph(ig_social, layout = "fr") +
        geom_edge_link(color = "gray") +
        geom_node_point(aes(size = node_size, color = node_color)) +
        scale_color_gradient(low = "green", high = "red") +
        labs(title = "Social Network", color = "Eigenvector\nCentrality", size = "Degree\nCentrality") +
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
        labs(title = "Component Interaction Matrix Heatmap", x = "Component 1", y = "Component 2") +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      # Calculate average degree (K) for the social space and component interaction space
      avg_degree_social <- mean(degree(ig_social))
      avg_degree_component <- mean(degree(graph_from_adjacency_matrix(component_matrix)))
      
      # Create a main title using sprintf with simulation parameters
      main_title <- sprintf("Simulation Parameters: N = %d, M = %d, BI_PROB = %.2f, K_S = %.2f, K_C = %.2f", 
                            self$N, self$M, self$BI_PROB, avg_degree_social, avg_degree_component)
      
      # Arrange plots with a main title
      grid.arrange(
        social_plot, bipartite_plot, heatmap_plot, ncol = 3,
        top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold"))
      )
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


saom_nk_enhanced <- SAOM_NK_Enhanced$new(M = 10, N = 10, BI_PROB = .15)
saom_nk_enhanced$visualize_networks()

# Run the local search method
saom_nk_enhanced$local_search_saom(iterations = 1)
saom_nk_enhanced$visualize_networks()

saom_nk_enhanced$local_search_saom(iterations = 1)
saom_nk_enhanced$visualize_networks()

saom_nk_enhanced$local_search_saom(iterations = 1)
saom_nk_enhanced$visualize_networks()

# Calculate and print network metrics
network_metrics <- saom_nk_enhanced$calculate_network_metrics()

# Plot the search progress
# saom_nk_enhanced$plot_search_progress()
saom_nk_enhanced$visualize_networks()



saom_nk_enhanced$local_search_saom(iterations = 5)
saom_nk_enhanced$visualize_networks()



saom_nk_enhanced$local_search_saom(iterations = 5)
saom_nk_enhanced$visualize_networks()


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








