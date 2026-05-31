# Load necessary libraries
library(R6)
library(network)
library(igraph)

# Define the SAOM_NK_Base class with local search methods
SAOM_NK_Base <- R6Class(
  "SAOM_NK_Base",
  
  public = list(
    M = NULL,                # Number of actors
    N = NULL,                # Number of components
    bipartite_network = NULL, # Bipartite network object (actors x components)
    social_network = NULL,   # Social space projection (actors linked by common components)
    search_landscape = NULL, # Search space projection (components linked by common actors)
    
    # Constructor to initialize the SAOM-NK model
    initialize = function(M, N) {
      self$M <- M
      self$N <- N
      self$bipartite_network <- self$generate_bipartite_network()
      self$social_network <- self$project_social_space()
      self$search_landscape <- self$project_search_space()
    },
    
    # Generate a random bipartite network object
    generate_bipartite_network = function() {
      set.seed(123)  # For reproducibility
      bipartite_matrix <- matrix(sample(0:1, self$M * self$N, replace = TRUE, prob = c(0.7, 0.3)),
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
    
    # Method for performing a local search by modifying bipartite ties
    local_search = function(actor_id, component_id, iterations) {
      for (i in 1:iterations) {
        # Flip a random tie between the actor and the component
        current_tie <- get.edge.value(self$bipartite_network, actor_id, component_id)
        new_tie <- ifelse(current_tie == 1, 0, 1)
        set.edge.value(self$bipartite_network, actor_id, component_id, new_tie)
        
        # Update the projections after modifying the bipartite network
        self$social_network <- self$project_social_space()
        self$search_landscape <- self$project_search_space()
        
        # Print the update status
        cat("Iteration", i, ": Updated bipartite tie between Actor", actor_id, "and Component", component_id, "\n")
      }
    },
    
    # Print current state of networks
    print_networks = function() {
      cat("Bipartite Network:\n")
      print(as.matrix.network(self$bipartite_network))
      cat("\nSocial Network:\n")
      print(as.matrix.network(self$social_network))
      cat("\nSearch Landscape Network:\n")
      print(as.matrix.network(self$search_landscape))
    }
  )
)

# Example usage of the SAOM-NK model with local search
saom_nk <- SAOM_NK_Base$new(M = 5, N = 10)
saom_nk$print_networks()

# Perform a local search modifying ties in the bipartite network
actor_id <- 2
component_id <- 3
saom_nk$local_search(actor_id, component_id, iterations = 5)

# Print networks after performing local search
saom_nk$print_networks()
