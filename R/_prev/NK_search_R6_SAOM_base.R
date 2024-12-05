
# Load necessary libraries
library(R6)
library(network)
library(RSiena)

# Define the SAOM_NK_Base class for handling bipartite networks and projections
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
    
    # Generate a random bipartite matrix and convert to a network object
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
    
    # Compute network statistics using RSiena
    compute_network_statistics = function() {
      social_data <- sienaDependent(as.sociomatrix(self$social_network))
      siena_effects <- getEffects(social_data)
      siena_model <- sienaAlgorithmCreate(projname = "SAOM_NK_Model")
      effects <- includeEffects(siena_effects, recip, transTrip, outdegreePop)
      
      # Run a simple RSiena model for statistics
      model_result <- siena07(siena_model, data = sienaDataCreate(social_data), effects = effects)
      summary(model_result)
    },
    
    # Method to evaluate interactions between social and search spaces
    evaluate_interspace_interactions = function() {
      # Convert the network objects to adjacency matrices
      social_adj <- as.matrix.network(self$social_network)
      search_adj <- as.matrix.network(self$search_landscape)
      
      # Ensure conformable arguments for interaction analysis
      interaction_matrix <- social_adj %*% t(search_adj)
      total_interaction_score <- sum(interaction_matrix)
      
      return(total_interaction_score)
    }
  )
)

# Example usage of the SAOM-NK model with network and RSiena integration
saom_nk <- SAOM_NK_Base$new(M = 5, N = 10)
print("Bipartite Network:")
plot(saom_nk$bipartite_network, main = "Bipartite Network of Actors and Components")

print("Social Space Network:")
plot(saom_nk$social_network, main = "Social Space Projection")

print("Search Landscape Network:")
heatmap(saom_nk$search_landscape[,], main = "Search Landscape Projection")

# Run network statistics using RSiena
cat("Network statistics using RSiena:\n")
saom_nk$compute_network_statistics()

# Evaluate interspace interactions
interspace_score <- saom_nk$evaluate_interspace_interactions()
cat("Interspace interaction score:", interspace_score, "\n")

