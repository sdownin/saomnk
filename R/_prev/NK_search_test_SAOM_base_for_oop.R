

#################################
##
##  1.  NK Landscape Creation
##
##################################


# Load necessary libraries
library(Matrix)

# Function to generate a random interaction matrix
generate_interaction_matrix <- function(N, K) {
  interaction_matrix <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    indices <- setdiff(1:N, i)  # Exclude self
    chosen_indices <- sample(indices, K, replace = FALSE)
    interaction_matrix[i, chosen_indices] <- 1
    interaction_matrix[i, i] <- 1  # Ensure the diagonal is set
  }
  return(interaction_matrix)
}

# Function to generate the power key (used for indexing combinations)
generate_power_key <- function(N) {
  return(2 ^ (seq(N - 1, 0)))
}

# Function to create an NK landscape with random fitness values
generate_nk_landscape <- function(N) {
  num_combinations <- 2^N
  landscape <- matrix(runif(num_combinations * N), nrow = num_combinations, ncol = N)
  return(landscape)
}

# Function to calculate fitness using the NK landscape and interaction matrix
evaluate_fitness <- function(N, landscape, interaction_matrix, current_position, power_key) {
  fitness_vector <- numeric(N)
  for (i in 1:N) {
    # Extract the sub-bit string using the interaction matrix
    interactions <- which(interaction_matrix[i, ] == 1)
    sub_position <- current_position[interactions]
    
    # Convert the sub-position to an index
    index <- sum(sub_position * power_key[1:length(sub_position)]) + 1  # 1-based index
    
    # Ensure the index is within the landscape bounds
    if (index > nrow(landscape)) {
      stop("Error: Calculated index is out of bounds. Check landscape dimensions.")
    }
    
    # Assign fitness value from the landscape
    fitness_vector[i] <- landscape[index, i]
  }
  return(fitness_vector)
}

# Example usage
N <- 6  # Number of components
K <- 2  # Number of interacting components

# Generate the interaction matrix and power key
interaction_matrix <- generate_interaction_matrix(N, K)
power_key <- generate_power_key(N)

# Create an NK landscape
landscape <- generate_nk_landscape(N)

# Define a current position (bit string)
current_position <- sample(0:1, N, replace = TRUE)

# Calculate fitness
fitness_vector <- evaluate_fitness(N, landscape, interaction_matrix, current_position, power_key)
cat("Fitness vector for the given position:\n")
print(fitness_vector)
print(mean(fitness_vector))










#################################
##
##  2.  Local Search
##
##################################

# Function for local search with long jumps in the NK model
local_search_nk <- function(N, K, landscape, interaction_matrix, power_key, iterations, p_jump = 0.1) {
  fitness_over_time <- numeric(iterations)
  
  # Generate initial random configuration (bit string)
  current_position <- sample(0:1, N, replace = TRUE)
  
  # Evaluate the initial fitness
  current_fitness_vector <- evaluate_fitness(N, landscape, interaction_matrix, current_position, power_key)
  current_fitness <- mean(current_fitness_vector)
  
  # Iterate for the specified number of time periods
  for (t in 1:iterations) {
    # Store the normalized fitness value
    fitness_min <- min(landscape)
    fitness_max <- max(landscape)
    fitness_over_time[t] <- (current_fitness - fitness_min) / (fitness_max - fitness_min)
    
    # Check for long jump probability
    if (runif(1) < p_jump) {
      # Perform a long jump: generate a completely new random position
      new_position <- sample(0:1, N, replace = TRUE)
    } else {
      # Perform a local search: flip one random bit in the current position
      new_position <- current_position
      flip_index <- sample(1:N, 1)
      new_position[flip_index] <- 1 - new_position[flip_index]  # Flip the bit
    }
    
    # Evaluate the new fitness
    new_fitness_vector <- evaluate_fitness(N, landscape, interaction_matrix, new_position, power_key)
    new_fitness <- mean(new_fitness_vector)
    
    # Accept the new position if it improves the fitness
    if (new_fitness > current_fitness) {
      current_position <- new_position
      current_fitness <- new_fitness
      cat("Iteration", t, ": Found a better configuration with fitness", current_fitness, "\n")
    }
  }
  
  return(fitness_over_time)
}

# Function to generate the power key (used for indexing combinations)
generate_power_key <- function(N) {
  return(2 ^ (seq(N - 1, 0)))
}

# Example usage
N <- 6  # Number of components
K <- 2  # Number of interacting components
iterations <- 50  # Number of time periods for local search
p_jump <- 0.1  # Probability of performing a long jump

# Generate the interaction matrix, power key, and fitness landscape
interaction_matrix <- generate_interaction_matrix(N, K)
power_key <- generate_power_key(N)
landscape <- generate_nk_landscape(N)

# Run the local search and plot fitness over time
fitness_over_time <- local_search_nk(N, K, landscape, interaction_matrix, power_key, iterations, p_jump)

# Plotting fitness over time
plot(fitness_over_time, type = "l", col = "green", lwd = 2, 
     main = "Results of Local Search with Long Jumps",
     xlab = "Time Periods", ylab = "Normalized Fitness")






##########
##########
## Mutation Rates
##########
# Function for local search with mutation rates and long jumps in the NK model
local_search_nk_with_mutation <- function(N, K, landscape, interaction_matrix, power_key, iterations, p_jump = 0.1, mutation_rate = 0.05) {
  fitness_over_time <- numeric(iterations)
  
  # Generate initial random configuration (bit string)
  current_position <- sample(0:1, N, replace = TRUE)
  
  # Evaluate the initial fitness
  current_fitness_vector <- evaluate_fitness(N, landscape, interaction_matrix, current_position, power_key)
  current_fitness <- mean(current_fitness_vector)
  
  # Iterate for the specified number of time periods
  for (t in 1:iterations) {
    # Store the normalized fitness value
    fitness_min <- min(landscape)
    fitness_max <- max(landscape)
    fitness_over_time[t] <- (current_fitness - fitness_min) / (fitness_max - fitness_min)
    
    # Check for long jump probability
    if (runif(1) < p_jump) {
      # Perform a long jump: generate a completely new random position
      new_position <- sample(0:1, N, replace = TRUE)
    } else {
      # Apply mutation based on the mutation rate to create a new position
      new_position <- current_position
      for (i in 1:N) {
        if (runif(1) < mutation_rate) {
          new_position[i] <- 1 - new_position[i]  # Flip the bit with mutation probability
        }
      }
    }
    
    # Evaluate the new fitness
    new_fitness_vector <- evaluate_fitness(N, landscape, interaction_matrix, new_position, power_key)
    new_fitness <- mean(new_fitness_vector)
    
    # Accept the new position if it improves the fitness
    if (new_fitness > current_fitness) {
      current_position <- new_position
      current_fitness <- new_fitness
      cat("Iteration", t, ": Found a better configuration with fitness", current_fitness, "\n")
    }
  }
  
  return(fitness_over_time)
}

# Example usage
N <- 6  # Number of components
K <- 2  # Number of interacting components
iterations <- 50  # Number of time periods for local search
p_jump <- 0.1  # Probability of performing a long jump
mutation_rate <- 0.05  # Mutation rate for bit flipping

# Generate the interaction matrix, power key, and fitness landscape
interaction_matrix <- generate_interaction_matrix(N, K)
power_key <- generate_power_key(N)
landscape <- generate_nk_landscape(N)

# Run the local search with mutation
fitness_over_time <- local_search_nk_with_mutation(N, K, landscape, interaction_matrix, power_key, iterations, p_jump, mutation_rate)

# Plotting fitness over time
plot(fitness_over_time, type = "l", col = "blue", lwd = 2, 
     main = "Results of Local Search with Mutation and Long Jumps",
     xlab = "Time Periods", ylab = "Normalized Fitness")


##########
### visualize
##########
# Function to visualize the impact of mutation rates on fitness over time
visualize_mutation_impact <- function(N, K, landscape, interaction_matrix, power_key, iterations, p_jump, mutation_rates) {
  fitness_results <- list()
  
  # Run local search for each mutation rate and store the results
  for (mutation_rate in mutation_rates) {
    fitness_over_time <- local_search_nk_with_mutation(N, K, landscape, interaction_matrix, power_key, iterations, p_jump, mutation_rate)
    fitness_results[[paste("Mutation Rate =", mutation_rate)]] <- fitness_over_time
  }
  
  # Plot fitness results for each mutation rate
  plot(1:iterations, fitness_results[[1]], type = "n", ylim = c(0, 1), 
       main = "Impact of Mutation Rates on Fitness Over Time",
       xlab = "Time Periods", ylab = "Normalized Fitness")
  
  colors <- rainbow(length(mutation_rates))
  for (i in 1:length(mutation_rates)) {
    lines(1:iterations, fitness_results[[i]], col = colors[i], lwd = 2)
  }
  
  # Add legend
  legend("bottomright", legend = names(fitness_results), col = colors, lwd = 2, cex = 0.8)
}

# Example parameters and setup
N <- 6  # Number of components
K <- 2  # Number of interacting components
iterations <- 50  # Number of time periods for local search
p_jump <- 0.1  # Probability of performing a long jump
mutation_rates <- c(0.01, 0.05, 0.1, 0.2)  # Different mutation rates to test

# Generate the interaction matrix, power key, and fitness landscape
interaction_matrix <- generate_interaction_matrix(N, K)
power_key <- generate_power_key(N)
landscape <- generate_nk_landscape(N)

# Visualize the impact of mutation rates
visualize_mutation_impact(N, K, landscape, interaction_matrix, power_key, iterations, p_jump, mutation_rates)





####
####
## Comparing Local Search with Different Jumps
####
# Function to compare the impact of different jump probabilities on fitness over time
compare_jump_probabilities <- function(N, K, landscape, interaction_matrix, power_key, iterations, mutation_rate, jump_probs) {
  fitness_results <- list()
  
  # Run local search for each jump probability and store the results
  for (p_jump in jump_probs) {
    fitness_over_time <- local_search_nk_with_mutation(N, K, landscape, interaction_matrix, power_key, iterations, p_jump, mutation_rate)
    fitness_results[[paste("Jump Prob =", p_jump)]] <- fitness_over_time
  }
  
  # Plot fitness results for each jump probability
  plot(1:iterations, fitness_results[[1]], type = "n", ylim = c(0, 1), 
       main = "Impact of Jump Probabilities on Fitness Over Time",
       xlab = "Time Periods", ylab = "Normalized Fitness")
  
  colors <- rainbow(length(jump_probs))
  for (i in 1:length(jump_probs)) {
    lines(1:iterations, fitness_results[[i]], col = colors[i], lwd = 2)
  }
  
  # Add legend
  legend("bottomright", legend = names(fitness_results), col = colors, lwd = 2, cex = 0.8)
}

# Example parameters and setup
N <- 6  # Number of components
K <- 2  # Number of interacting components
iterations <- 50  # Number of time periods for local search
mutation_rate <- 0.05  # Fixed mutation rate for all runs
jump_probs <- c(0.0, 0.1, 0.2, 0.5)  # Different jump probabilities to test

# Generate the interaction matrix, power key, and fitness landscape
interaction_matrix <- generate_interaction_matrix(N, K)
power_key <- generate_power_key(N)
landscape <- generate_nk_landscape(N)

# Visualize the impact of jump probabilities
compare_jump_probabilities(N, K, landscape, interaction_matrix, power_key, iterations, mutation_rate, jump_probs)

##







#################################
##
##  3. Decentralized search
##
##################################

# Function for decentralized local search in the NK model with debugging
decentralized_search_nk <- function(N, K, landscape, interaction_matrix, power_key, iterations, reorg) {
  fitness_over_time <- matrix(0, nrow = iterations, ncol = 2)
  
  # Generate initial random configuration (bit string)
  current_position <- sample(0:1, N, replace = TRUE)
  
  # Evaluate the initial fitness
  current_fitness_vector <- evaluate_fitness(N, landscape, interaction_matrix, current_position, power_key)
  current_fitness <- mean(current_fitness_vector)
  
  # Normalize the initial fitness
  max_fit <- max(landscape)
  min_fit <- min(landscape)
  fitness_norm <- (current_fitness - min_fit) / (max_fit - min_fit)
  
  for (t in 1:iterations) {
    fitness_over_time[t, ] <- c(t, fitness_norm)
    
    if (t < reorg) {
      # Decentralized search: split the vector into two halves
      new_position <- current_position
      part1 <- 1:(N / 2)
      part2 <- (N / 2 + 1):N
      
      # Randomly flip one bit in each half
      choice_varA <- sample(part1, 1)
      choice_varB <- sample(part2, 1)
      new_position[choice_varA] <- 1 - new_position[choice_varA]
      new_position[choice_varB] <- 1 - new_position[choice_varB]
      
      # Check that the row index for new_position is valid
      row <- sum(new_position * power_key) + 1  # 1-based index for row lookup
      if (row > nrow(landscape)) {
        stop(paste("Error: Row index", row, "out of bounds for landscape with", nrow(landscape), "rows."))
      }
      
      # Evaluate the new fitness for each half
      new_fitA <- mean(landscape[row, 1:(N / 2)])
      new_fitB <- mean(landscape[row, ((N / 2) + 1):N])
      
      # Accept new half if it improves fitness
      current_row <- sum(current_position * power_key) + 1  # Current row index
      current_fitA <- mean(landscape[current_row, 1:(N / 2)])
      current_fitB <- mean(landscape[current_row, ((N / 2) + 1):N])
      
      if (new_fitA > current_fitA) {
        current_position[part1] <- new_position[part1]
      }
      if (new_fitB > current_fitB) {
        current_position[part2] <- new_position[part2]
      }
      
    } else {
      # Centralized search: flip one random bit in the full vector
      new_position <- current_position
      choice_var <- sample(1:N, 1)
      new_position[choice_var] <- 1 - new_position[choice_var]
      
      # Check that the row index for new_position is valid
      row <- sum(new_position * power_key) + 1  # 1-based index for row lookup
      if (row > nrow(landscape)) {
        stop(paste("Error: Row index", row, "out of bounds for landscape with", nrow(landscape), "rows."))
      }
      
      # Evaluate the new fitness
      new_fitness_vector <- evaluate_fitness(N, landscape, interaction_matrix, new_position, power_key)
      new_fitness <- mean(new_fitness_vector)
      
      # Accept the new configuration if it improves fitness
      if (new_fitness > current_fitness) {
        current_position <- new_position
        current_fitness <- new_fitness
      }
    }
    
    # Normalize fitness for plotting
    fitness_norm <- (current_fitness - min_fit) / (max_fit - min_fit)
  }
  
  return(fitness_over_time)
}


# Example usage
N <- 6  # Number of components
K <- 2  # Number of interacting components
iterations <- 50  # Number of time periods for decentralized search
reorg <- 25  # Iteration at which to switch from decentralized to centralized search

# Generate the interaction matrix, power key, and fitness landscape
interaction_matrix <- generate_interaction_matrix(N, K)
power_key <- generate_power_key(N)
landscape <- generate_nk_landscape(N)

# Run the decentralized search
fitness_over_time <- decentralized_search_nk(N, K, landscape, interaction_matrix, power_key, iterations, reorg)

# Plotting the results
plot(fitness_over_time[, 1], fitness_over_time[, 2], type = "l", col = "blue", lwd = 2,
     main = "Decentralized vs Centralized Search",
     xlab = "Time Periods", ylab = "Normalized Fitness")







