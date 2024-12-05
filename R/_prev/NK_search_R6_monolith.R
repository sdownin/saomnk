# Load the R6 package
library(R6)

# Define the NKModel base class
NKModel <- R6Class(
  "NKModel",
  
  public = list(
    N = NULL,                # Number of components
    K = NULL,                # Number of interacting components
    interaction_matrix = NULL,
    power_key = NULL,
    landscape = NULL,
    
    # Constructor: initializes the NK model properties
    initialize = function(N, K) {
      self$N <- N
      self$K <- K
      self$interaction_matrix <- self$generate_interaction_matrix()
      self$power_key <- self$generate_power_key()
      self$landscape <- self$generate_nk_landscape()
    },
    
    # Method to generate a random interaction matrix
    generate_interaction_matrix = function() {
      interaction_matrix <- matrix(0, nrow = self$N, ncol = self$N)
      for (i in 1:self$N) {
        indices <- setdiff(1:self$N, i)
        chosen_indices <- sample(indices, self$K, replace = FALSE)
        interaction_matrix[i, chosen_indices] <- 1
        interaction_matrix[i, i] <- 1  # Include self-loop
      }
      return(interaction_matrix)
    },
    
    # Method to generate the power key for binary indexing
    generate_power_key = function() {
      return(2 ^ (seq(self$N - 1, 0)))
    },
    
    # Method to generate the NK landscape with random fitness values
    generate_nk_landscape = function() {
      num_combinations <- 2^self$N
      landscape <- matrix(runif(num_combinations * self$N), nrow = num_combinations, ncol = self$N)
      return(landscape)
    },
    
    # Method to calculate fitness for a given configuration
    calculate_fitness = function(current_position) {
      fitness_vector <- numeric(self$N)
      for (i in 1:self$N) {
        sub_bit_string <- current_position[which(self$interaction_matrix[i, ] == 1)]
        index <- sum(sub_bit_string * self$power_key[1:length(sub_bit_string)]) + 1  # 1-based index
        
        if (index > nrow(self$landscape)) {
          stop("Error: Index out of bounds. Check landscape dimensions.")
        }
        fitness_vector[i] <- self$landscape[index, i]
      }
      return(mean(fitness_vector))
    },
    
    # Scenario 1: Basic NK Simulation
    basic_simulation = function() {
      # Evaluate all configurations and identify local peaks
      num_combinations <- 2^self$N
      result_matrix <- matrix(0, nrow = num_combinations, ncol = (2 * self$N) + 2)
      
      for (comb_index in 0:(num_combinations - 1)) {
        current_combination <- as.integer(intToBits(comb_index)[1:self$N])
        current_fitness <- self$calculate_fitness(current_combination)
        result_matrix[comb_index + 1, 1:self$N] <- current_combination
        result_matrix[comb_index + 1, (self$N + 1):(2 * self$N)] <- current_fitness
        result_matrix[comb_index + 1, (2 * self$N) + 1] <- mean(current_fitness)
        
        is_local_peak <- TRUE
        for (j in 1:self$N) {
          neighbor_combination <- current_combination
          neighbor_combination[j] <- 1 - neighbor_combination[j]
          neighbor_index <- sum(neighbor_combination * self$power_key) + 1
          
          if (neighbor_index <= num_combinations) {
            if (current_fitness < result_matrix[neighbor_index, (2 * self$N) + 1]) {
              is_local_peak <- FALSE
              break
            }
          }
        }
        result_matrix[comb_index + 1, (2 * self$N) + 2] <- as.integer(is_local_peak)
      }
      
      return(result_matrix)
    },
    
    # Scenario 2: Local Search with Long Jumps
    local_search_with_jumps = function(initial_position, iterations, p_jump = 0.1) {
      current_position <- initial_position
      current_fitness <- self$calculate_fitness(current_position)
      fitness_over_time <- numeric(iterations)
      
      for (t in 1:iterations) {
        fitness_over_time[t] <- current_fitness
        
        if (runif(1) < p_jump) {
          # Perform a long jump to a new random position
          new_position <- sample(0:1, self$N, replace = TRUE)
        } else {
          # Flip one random bit for local search
          new_position <- current_position
          flip_index <- sample(1:self$N, 1)
          new_position[flip_index] <- 1 - new_position[flip_index]
        }
        
        new_fitness <- self$calculate_fitness(new_position)
        if (new_fitness > current_fitness) {
          current_position <- new_position
          current_fitness <- new_fitness
        }
      }
      
      return(fitness_over_time)
    },
    
    # Scenario 3: Decentralized Local Search
    decentralized_local_search = function(initial_position, iterations, reorg_point) {
      current_position <- initial_position
      fitness_over_time <- numeric(iterations)
      
      for (t in 1:iterations) {
        if (t < reorg_point) {
          # Split the vector for decentralized search
          part1 <- 1:(self$N / 2)
          part2 <- ((self$N / 2) + 1):self$N
          new_combA <- current_position[part1]
          new_combB <- current_position[part2]
          
          # Mutate one bit in each part
          choice_varA <- sample(part1, 1)
          choice_varB <- sample(part2, 1)
          new_combA[choice_varA] <- 1 - new_combA[choice_varA]
          new_combB[choice_varB] <- 1 - new_combB[choice_varB]
          
          # Combine parts and evaluate fitness
          new_position <- c(new_combA, new_combB)
          new_fitness <- self$calculate_fitness(new_position)
          
          # Accept the change if it improves fitness
          current_fitness <- max(current_fitness, new_fitness)
        } else {
          # Centralized search after reorg_point
          new_position <- current_position
          choice_var <- sample(1:self$N, 1)
          new_position[choice_var] <- 1 - new_position[choice_var]
          new_fitness <- self$calculate_fitness(new_position)
          
          if (new_fitness > current_fitness) {
            current_position <- new_position
            current_fitness <- new_fitness
          }
        }
        
        fitness_over_time[t] <- current_fitness
      }
      
      return(fitness_over_time)
    }
  )
)

# Example usage
nk_model <- NKModel$new(N = 6, K = 2)
initial_position <- sample(0:1, nk_model$N, replace = TRUE)

# Run basic simulation
basic_result <- nk_model$basic_simulation()
cat("Basic simulation result:\n")
print(basic_result)

# Run local search with long jumps
local_search_result <- nk_model$local_search_with_jumps(initial_position, iterations = 50)

# Run decentralized local search
decentralized_result <- nk_model$decentralized_local_search(initial_position, iterations = 50, reorg_point = 25)

# Plot the fitness over time for different scenarios
plot(local_search_result, type = "l", col = "blue", lwd = 2, main = "Local Search with Long Jumps", xlab = "Iteration", ylab = "Fitness")
