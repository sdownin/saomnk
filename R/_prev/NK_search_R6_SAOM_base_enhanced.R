
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
    progress_scores = NULL,  # Track scores over iterations
    
    # Constructor to initialize the SAOM-NK model
    initialize = function(M, N, BI_PROB=0.2) {
      self$M <- M
      self$N <- N
      self$BI_PROB <- BI_PROB
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
    
    # Local search method for broader exploration
    local_search = function(iterations, mutation_rate = 0.1) {
      current_score <- self$evaluate_interspace_interactions()
      self$progress_scores <- c(current_score)  # Track the initial score
      
      for (i in 1:iterations) {
        # Convert the bipartite network to an adjacency matrix for direct modification
        bipartite_matrix <- as.matrix.network(self$bipartite_network)
        
        # Mutate a random set of actor-component ties based on the mutation rate
        for (actor_id in 1:self$M) {
          for (component_id in 1:self$N) {
            if (runif(1) < mutation_rate) {
              current_tie <- bipartite_matrix[actor_id, component_id]
              bipartite_matrix[actor_id, component_id] <- ifelse(current_tie == 1, 0, 1)
            }
          }
        }
        
        # Update the bipartite network and projections
        self$bipartite_network <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
        self$social_network <- self$project_social_space()
        self$search_landscape <- self$project_search_space()
        
        # Evaluate new interaction score
        new_score <- self$evaluate_interspace_interactions()
        
        # Accept the new configuration if it improves the score
        if (new_score > current_score) {
          current_score <- new_score
          cat("Iteration", i, ": Improved interaction score to", current_score, "\n")
        } else {
          # Optionally, keep some non-improving moves for exploration
          if (runif(1) < 0.1) {  # Example exploration rate
            current_score <- new_score
            cat("Iteration", i, ": Accepted non-improving score for exploration:", current_score, "\n")
          } else {
            # Revert to the previous configuration if not accepted
            bipartite_matrix <- as.matrix.network(self$bipartite_network)
          }
        }
        
        # Track the score at this iteration
        self$progress_scores <- c(self$progress_scores, current_score)
      }
    },
    
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
    visualize_networks = function(BI_PROB = 0.3) {
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
                            self$N, self$M, BI_PROB, avg_degree_social, avg_degree_component)
      
      # Arrange plots with a main title
      grid.arrange(
        social_plot, bipartite_plot, heatmap_plot, ncol = 3,
        top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold"))
      )
    }
    
    
    ,
    
    
    # # Visualization using ggplot2 and ggraph: bipartite network, social network, and component interaction matrix heatmap
    # visualize_networks = function() {
    #   # 1. Bipartite network plot
    #   bipartite_matrix <- as.matrix.network(self$bipartite_network)
    #   bipartite_df <- data.frame(
    #     Actor = rep(1:self$M, each = self$N),
    #     Component = rep(1:self$N, times = self$M),
    #     Tie = as.vector(bipartite_matrix)
    #   )
    #   bipartite_df <- bipartite_df[bipartite_df$Tie == 1, ]  # Only keep edges with ties
    #   
    #   # bipartite_plot <- ggplot(bipartite_df, aes(x = as.factor(Actor), y = as.factor(Component))) +
    #   #   geom_point(color = "blue", size = 3) +
    #   #   labs(title = "Bipartite Network\n(Actors vs. Components)", x = "Actors", y = "Components") +
    #   #   theme_minimal()
    #   #####
    #   ig_bipartite <- graph_from_biadjacency_matrix(as.matrix.network(self$bipartite_network))
    #   # Set node attributes for shape
    #   V(ig_bipartite)$type <- bipartite_mapping(ig_bipartite)$type
    #   V(ig_bipartite)$shape <- ifelse(V(ig_bipartite)$type, "square", "circle")
    #   V(ig_bipartite)$color <- ifelse(V(ig_bipartite)$type, "lightblue", "orange")
    #   # plot bipartite style
    #   bipartite_plot <- ggraph(ig_bipartite, layout = "bipartite") +
    #     geom_edge_link(color = "gray") +
    #     geom_node_point(aes(shape = shape, color = color), size = 5) +
    #     scale_shape_manual(values = c("circle" = 16, "square" = 15)) +
    #     scale_color_identity() +
    #     labs(title = "Bipartite Network\n(Actors and Components)") +
    #     theme_minimal() +
    #     theme(
    #       legend.position = "none",
    #       panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(),
    #       axis.title.x = element_blank(),
    #       axis.title.y = element_blank(),
    #       axis.text.x = element_blank(),       
    #       axis.text.y = element_blank()      
    #     )
    #   
    #   # 2. Social network plot using ggraph with different statistics for node color and size
    #   ig_social <- graph_from_adjacency_matrix(as.matrix.network(self$social_network), mode = "undirected")
    #   
    #   # Calculate network statistics
    #   node_size <- degree(ig_social)  # degree centrality for node size
    #   node_color <- eigen_centrality(ig_social)$vector  # Eigenvector centrality for node color
    #   
    #   social_plot <- ggraph(ig_social, layout = "fr") +
    #     geom_edge_link(color = "gray") +
    #     geom_node_point(aes(size = node_size, color = node_color)) +
    #     scale_color_gradient(low = "green", high = "red") +
    #     labs(title = "Social Network", color = "Eigenvector\nCentrality", size = "Degree\nCentrality") +
    #     theme_minimal() +
    #     theme(
    #       legend.position = "bottom",
    #       legend.box = "vertical",  # Stack legends vertically
    #       panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(),
    #       axis.title.x = element_blank(),
    #       axis.title.y = element_blank(),
    #       axis.text.x = element_blank(),       # Remove x-axis numbering
    #       axis.text.y = element_blank()        # Remove y-axis numbering
    #     ) +
    #     guides(
    #       color = guide_legend(order = 1, nrow = 2),
    #       size = guide_legend(order = 2, nrow = 2)
    #     )
    #   
    #   # 3. Component interaction matrix heatmap
    #   component_matrix <- as.matrix.network(self$search_landscape)
    #   component_df <- melt(component_matrix)
    #   colnames(component_df) <- c("Component1", "Component2", "Interaction")
    #   
    #   heatmap_plot <- ggplot(component_df, aes(x = Component1, y = Component2, fill = Interaction)) +
    #     geom_tile() +
    #     scale_fill_gradient(low = "white", high = "red") +
    #     labs(title = "Component Interaction Matrix Heatmap", x = "Component 1", y = "Component 2") +
    #     theme_minimal() +
    #     theme(legend.position = "bottom")
    #   
    #   # Combine plots into a panel
    #   grid.arrange(social_plot, bipartite_plot, heatmap_plot, ncol = 3)
    # },
    
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

# Example usage of the SAOM-NK model with enhanced local search
saom_nk_enhanced <- SAOM_NK_Enhanced$new(M = 12, N = 25, BI_PROB = .05)
saom_nk_enhanced$visualize_networks()

###
# saom_nk_enhanced$local_search(iterations = 5, mutation_rate = 0.02)
# saom_nk_enhanced$visualize_networks()
# saom_nk_enhanced$local_search(iterations = 5, mutation_rate = 0.02)
# saom_nk_enhanced$visualize_networks()

# Run the local search method
saom_nk_enhanced$local_search(iterations = 10, mutation_rate = 0.01)

# Plot the search progress
saom_nk_enhanced$plot_search_progress()

# Calculate and print network metrics
network_metrics <- saom_nk_enhanced$calculate_network_metrics()

saom_nk_enhanced$visualize_networks()



# # Visualize the social network
# saom_nk_enhanced$visNetwork_social()












# # Load necessary libraries
# library(R6)
# library(network)
# library(igraph)
# library(visNetwork)
# 
# # Define the SAOM_NK_Enhanced class with enhancements and visualization
# SAOM_NK_Enhanced <- R6Class(
#   "SAOM_NK_Enhanced",
#   
#   public = list(
#     M = NULL,                # Number of actors
#     N = NULL,                # Number of components
#     bipartite_network = NULL, # Bipartite network object (actors x components)
#     social_network = NULL,   # Social space projection (network object)
#     search_landscape = NULL, # Search space projection (network object)
#     progress_scores = NULL,  # Track scores over iterations
#     
#     # Constructor to initialize the SAOM-NK model
#     initialize = function(M, N) {
#       self$M <- M
#       self$N <- N
#       self$bipartite_network <- self$generate_bipartite_network()
#       self$social_network <- self$project_social_space()
#       self$search_landscape <- self$project_search_space()
#       self$progress_scores <- c()  # Initialize empty vector to track scores
#     },
#     
#     # Generate a random bipartite network object
#     generate_bipartite_network = function() {
#       set.seed(123)  # For reproducibility
#       bipartite_matrix <- matrix(sample(0:1, self$M * self$N, replace = TRUE, prob = c(0.7, 0.3)),
#                                  nrow = self$M, ncol = self$N)
#       bipartite_net <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
#       return(bipartite_net)
#     },
#     
#     # Project the social space (actor network) from the bipartite network
#     project_social_space = function() {
#       bipartite_matrix <- as.matrix.network(self$bipartite_network)
#       actor_matrix <- bipartite_matrix %*% t(bipartite_matrix)
#       diag(actor_matrix) <- 0  # Remove self-loops
#       social_net <- network(actor_matrix, directed = FALSE)
#       return(social_net)
#     },
#     
#     # Project the search landscape space (component interaction network)
#     project_search_space = function() {
#       bipartite_matrix <- as.matrix.network(self$bipartite_network)
#       component_matrix <- t(bipartite_matrix) %*% bipartite_matrix
#       diag(component_matrix) <- 0  # Remove self-loops
#       search_net <- network(component_matrix, directed = FALSE)
#       return(search_net)
#     },
#     
#     # Modified hill climbing search to handle early peaks
#     hill_climbing_search = function(actor_id, component_id, iterations) {
#       current_score <- self$evaluate_interspace_interactions()
#       self$progress_scores <- c(current_score)  # Track the initial score
#       peak_found <- FALSE
#       
#       for (i in 1:iterations) {
#         if (peak_found) break  # Stop if no improvements are found
#         
#         # Convert the bipartite network to an adjacency matrix for direct modification
#         bipartite_matrix <- as.matrix.network(self$bipartite_network)
#         
#         # Flip the tie between the specified actor and component
#         current_tie <- bipartite_matrix[actor_id, component_id]
#         new_tie <- ifelse(current_tie == 1, 0, 1)
#         bipartite_matrix[actor_id, component_id] <- new_tie
#         
#         # Update the bipartite network and projections
#         self$bipartite_network <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
#         self$social_network <- self$project_social_space()
#         self$search_landscape <- self$project_search_space()
#         
#         # Evaluate new interaction score
#         new_score <- self$evaluate_interspace_interactions()
#         
#         # Accept or revert change based on new score and track progress
#         if (new_score > current_score) {
#           current_score <- new_score
#           cat("Iteration", i, ": Improved interaction score to", current_score, "\n")
#         } else {
#           # Revert the change
#           bipartite_matrix[actor_id, component_id] <- current_tie
#           self$bipartite_network <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
#           peak_found <- TRUE  # Mark that no improvement was found in this iteration
#           cat("Iteration", i, ": No improvement found, potential peak reached.\n")
#         }
#         
#         # Track the score at this iteration
#         self$progress_scores <- c(self$progress_scores, current_score)
#       }
#     },
#     
#     # Visualization of the search progress
#     plot_search_progress = function() {
#       if (length(self$progress_scores) == 0) {
#         stop("No search progress data available. Run a search method first.")
#       }
#       
#       plot(self$progress_scores, type = "l", col = "blue", lwd = 2,
#            main = "Hill Climbing Search Progress",
#            xlab = "Iteration", ylab = "Interaction Score",
#            ylim = range(self$progress_scores))
#     },
#     
#     # Evaluate interactions between social and search spaces using the bipartite matrix
#     evaluate_interspace_interactions = function() {
#       bipartite_matrix <- as.matrix.network(self$bipartite_network)
#       actor_influence_on_components <- bipartite_matrix %*% as.matrix.network(self$search_landscape)
#       total_interaction_score <- sum(actor_influence_on_components)
#       return(total_interaction_score)
#     }
#   )
# )
# 
# # Example usage of the SAOM-NK model with visualization
# saom_nk_enhanced <- SAOM_NK_Enhanced$new(M = 5, N = 10)
# 
# # Perform hill climbing search
# saom_nk_enhanced$hill_climbing_search(actor_id = 3, component_id = 1, iterations = 10)
# 
# # Plot the search progress
# saom_nk_enhanced$plot_search_progress()
# 
# 
# 
# 
# # Perform hill climbing search
# saom_nk_enhanced$hill_climbing_search(actor_id = 2, component_id = 5, iterations = 5)
# saom_nk_enhanced$plot_search_progress()




# # Load necessary libraries
# library(R6)
# library(network)
# library(igraph)
# library(visNetwork)
# 
# # Define the SAOM_NK_Enhanced class with all enhancements
# SAOM_NK_Enhanced <- R6Class(
#   "SAOM_NK_Enhanced",
#   
#   public = list(
#     M = NULL,                # Number of actors
#     N = NULL,                # Number of components
#     bipartite_network = NULL, # Bipartite network object (actors x components)
#     social_network = NULL,   # Social space projection (network object)
#     search_landscape = NULL, # Search space projection (network object)
#     
#     # Constructor to initialize the SAOM-NK model
#     initialize = function(M, N) {
#       self$M <- M
#       self$N <- N
#       self$bipartite_network <- self$generate_bipartite_network()
#       self$social_network <- self$project_social_space()
#       self$search_landscape <- self$project_search_space()
#     },
#     
#     # Generate a random bipartite network object
#     generate_bipartite_network = function() {
#       set.seed(123)  # For reproducibility
#       bipartite_matrix <- matrix(sample(0:1, self$M * self$N, replace = TRUE, prob = c(0.7, 0.3)),
#                                  nrow = self$M, ncol = self$N)
#       bipartite_net <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
#       return(bipartite_net)
#     },
#     
#     # Project the social space (actor network) from the bipartite network
#     project_social_space = function() {
#       bipartite_matrix <- as.matrix.network(self$bipartite_network)
#       actor_matrix <- bipartite_matrix %*% t(bipartite_matrix)
#       diag(actor_matrix) <- 0  # Remove self-loops
#       social_net <- network(actor_matrix, directed = FALSE)
#       return(social_net)
#     },
#     
#     # Project the search landscape space (component interaction network)
#     project_search_space = function() {
#       bipartite_matrix <- as.matrix.network(self$bipartite_network)
#       component_matrix <- t(bipartite_matrix) %*% bipartite_matrix
#       diag(component_matrix) <- 0  # Remove self-loops
#       search_net <- network(component_matrix, directed = FALSE)
#       return(search_net)
#     },
#     
#     # Complex search strategies: hill climbing
#     hill_climbing_search = function(actor_id, component_id, iterations) {
#       current_score <- self$evaluate_interspace_interactions()
#       for (i in 1:iterations) {
#         # Convert the bipartite network to an adjacency matrix for direct modification
#         bipartite_matrix <- as.matrix.network(self$bipartite_network)
#         
#         # Flip the tie between the specified actor and component
#         current_tie <- bipartite_matrix[actor_id, component_id]
#         new_tie <- ifelse(current_tie == 1, 0, 1)
#         bipartite_matrix[actor_id, component_id] <- new_tie
#         
#         # Update the bipartite network and projections
#         self$bipartite_network <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
#         self$social_network <- self$project_social_space()
#         self$search_landscape <- self$project_search_space()
#         
#         # Evaluate new interaction score
#         new_score <- self$evaluate_interspace_interactions()
#         
#         # Accept or revert change based on new score
#         if (new_score > current_score) {
#           current_score <- new_score
#           cat("Iteration", i, ": Improved interaction score to", current_score, "\n")
#         } else {
#           # Revert the change
#           bipartite_matrix[actor_id, component_id] <- current_tie
#           self$bipartite_network <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
#         }
#       }
#     },
#     
#     # Complex search strategies: simulated annealing
#     simulated_annealing_search = function(actor_id, component_id, iterations, initial_temp = 1, cooling_rate = 0.95) {
#       current_score <- self$evaluate_interspace_interactions()
#       temperature <- initial_temp
#       
#       for (i in 1:iterations) {
#         # Convert the bipartite network to an adjacency matrix for direct modification
#         bipartite_matrix <- as.matrix.network(self$bipartite_network)
#         
#         # Flip the tie between the specified actor and component
#         current_tie <- bipartite_matrix[actor_id, component_id]
#         new_tie <- ifelse(current_tie == 1, 0, 1)
#         bipartite_matrix[actor_id, component_id] <- new_tie
#         
#         # Update the bipartite network and projections
#         self$bipartite_network <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
#         self$social_network <- self$project_social_space()
#         self$search_landscape <- self$project_search_space()
#         
#         # Evaluate new interaction score
#         new_score <- self$evaluate_interspace_interactions()
#         
#         # Accept based on annealing probability
#         if (new_score > current_score || runif(1) < exp((new_score - current_score) / temperature)) {
#           current_score <- new_score
#           cat("Iteration", i, ": Accepted new score", current_score, "\n")
#         } else {
#           # Revert the change
#           bipartite_matrix[actor_id, component_id] <- current_tie
#           self$bipartite_network <- network(bipartite_matrix, bipartite = TRUE, directed = FALSE)
#         }
#         
#         # Cool down temperature
#         temperature <- temperature * cooling_rate
#       }
#     },
#     
#     # Calculate network metrics using igraph
#     calculate_network_metrics = function() {
#       ig_social <- graph_from_adjacency_matrix(as.matrix.network(self$social_network), mode = "undirected")
#       degree_centrality <- degree(ig_social)
#       betweenness_centrality <- betweenness(ig_social)
#       eigenvector_centrality <- eigen_centrality(ig_social)$vector
#       
#       cat("Network Metrics:\n")
#       cat("Degree Centrality:\n")
#       print(degree_centrality)
#       cat("Betweenness Centrality:\n")
#       print(betweenness_centrality)
#       cat("Eigenvector Centrality:\n")
#       print(eigenvector_centrality)
#       
#       return(list(
#         degree = degree_centrality,
#         betweenness = betweenness_centrality,
#         eigenvector = eigenvector_centrality
#       ))
#     },
#     
#     # Interactive visualization using visNetwork
#     visualize_network = function() {
#       ig_social <- graph_from_adjacency_matrix(as.matrix.network(self$social_network), mode = "undirected")
#       vis_data <- toVisNetworkData(ig_social)
#       
#       visNetwork(nodes = vis_data$nodes, edges = vis_data$edges) %>%
#         visNodes(shape = "dot", size = 10) %>%
#         visEdges(arrows = "to") %>%
#         visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
#         visLayout(randomSeed = 123)
#     },
#     
#     # Method to evaluate interactions between social and search spaces using the bipartite matrix
#     evaluate_interspace_interactions = function() {
#       # Convert the bipartite network to an adjacency matrix
#       bipartite_matrix <- as.matrix.network(self$bipartite_network)
#       
#       # Check the dimensions of the matrices for debugging
#       cat("Dimensions of bipartite_matrix:", dim(bipartite_matrix), "\n")
#       cat("Dimensions of social network (M x M):", dim(as.matrix.network(self$social_network)), "\n")
#       cat("Dimensions of search landscape (N x N):", dim(as.matrix.network(self$search_landscape)), "\n")
#       
#       # Calculate influence from actors to components using the bipartite matrix
#       actor_influence_on_components <- bipartite_matrix %*% as.matrix.network(self$search_landscape)
#       
#       # Sum the influence to get a total interaction score
#       total_interaction_score <- sum(actor_influence_on_components)
#       
#       return(total_interaction_score)
#     }
#     
#   )
# )
# 
# # Example usage of the SAOM-NK model with enhancements
# saom_nk_enhanced <- SAOM_NK_Enhanced$new(M = 5, N = 10)
# 
# # Perform hill climbing search
# saom_nk_enhanced$hill_climbing_search(actor_id = 2, component_id = 3, iterations = 10)
# 
# # Perform simulated annealing search
# saom_nk_enhanced$simulated_annealing_search(actor_id = 2, component_id = 4, iterations = 10)
# 
# # Calculate and print network metrics
# network_metrics <- saom_nk_enhanced$calculate_network_metrics()
# 
# # Visualize the social network
# saom_nk_enhanced$visualize_network()
