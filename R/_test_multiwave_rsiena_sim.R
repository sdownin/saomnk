

###############################################################################
###                                                                         ###
###       SimulateNetworks for Two-Mode Bipartite Networks                  ###
###                                                                         ###
###############################################################################

library(RSiena) # or RSienaTest

# SimulateBipartiteNetworks <- function(M, N, W, rate, dens, alt, ego, sim, groupMembership, strategy){
#   # Simulates W consecutive network waves for a two-mode bipartite network
#   # with M actors and N components (bipartite structure).
#   #
#   # Arguments:
#   # - M: Number of actors (rows in the bipartite matrix).
#   # - N: Number of components (columns in the bipartite matrix).
#   # - W: Number of waves to simulate.
#   # - rate: Rate parameter for simulation.
#   # - dens: Outdegree effect for actors.
#   # - alt, ego, sim: altX, egoX, and simX effect parameters for covariates.
#   # - groupMembership: Constant covariate indicating actor group membership.
#   # - strategy: Time-varying covariate for actor strategies.
#   #
#   # Returns:
#   # - List containing covariates, simulated networks, and summary statistics.
#   
#   ## Create Initial Network
#   set.seed(123)
#   network_array <- array(NA, dim = c(M, N, W)) # Start with 2 waves
#   network_array[, , 1] <- matrix(rbinom(M * N, 1, prob = 0.2), nrow = M, ncol = N)
#   for (i in 2:W) {
#     network_array[, , W] <- network_array[, , W-1] # Start with same network
#     ##
#     ##**TODO** check changing 1 dyad vs. identical networks; check sensitivity 
#   }
#   
#   ACTORS <- sienaNodeSet(M, nodeSetName = "Actors")
#   COMPONENTS <- sienaNodeSet(N, nodeSetName = "Components")
#   
#   # Ensure it is a bipartite structure
#   siena_network <- sienaDependent(network_array, type = "bipartite", nodeSet = c("Actors", "Components"), allowOnly = FALSE)
#   
#   ## Create Covariates
#   group_covar <- coCovar(groupMembership) # Constant covariate
#   strategy_covar <- varCovar(strategy)    # Time-varying covariate
#   
#   ## RSiena Data Object
#   siena_data <- sienaDataCreate(siena_network, group_covar, strategy_covar, nodeSets = list(ACTORS,COMPONENTS))
#   
#   ## Define Effects
#   effects <- getEffects(siena_data)
#   effectsDocumentation(effects)
#   effects <- includeEffects(effects, density, parameter = dens)
#   effects <- includeEffects(effects, egoX, interaction1 = "strategy_covar") ## parameter = ego
#   effects <- includeEffects(effects, sameEgoInDist2, interaction1 = "group_covar")
#   
#   ## Algorithm Settings
#   siena_algorithm <- sienaAlgorithmCreate(projname = "bipartite_sim", useStdInits = FALSE, 
#                                           cond = FALSE, nsub = 0, simOnly = TRUE, n3 = 100)
#   
#   ## Simulate the First Wave
#   first_sim <- siena07(siena_algorithm, data = siena_data, effects = effects, returnDeps = TRUE)
#   
#   ## Simulate Consecutive Waves
#   network_array_full <- array(0, dim = c(M, N, W)) # Full wave array
#   network_array_full[, , 1:W] <- network_array      # First two waves
#   
#   for (w in 2:W) {
#     # Extract simulated network from the previous wave
#     last_wave <- matrix(0, nrow = M, ncol = N)
#     sim_data <- first_sim$sims[[ length(first_sim$sims) ]][[ 1 ]][[ 1 ]]$`1`
#     for (i in 1:nrow(sim_data)) {
#       last_wave[sim_data[i, 1], sim_data[i, 2]] <- sim_data[i, 3]
#     }
#     
#     network_array_full[, , w] <- last_wave  # Add the wave to the array
#     
#     # last_waves_w <- c(last_wave, last_wave)
#     
#     # Prepare the next wave
#     if (w < W) { 
#       cat(sprintf('\n w=%s',w))
#       siena_network <- sienaDependent(array(c(last_wave, last_wave), dim = c(M, N, W)), 
#                                       type = "bipartite", nodeSet = c("Actors", "Components"), allowOnly = FALSE)
#       siena_data <- sienaDataCreate(siena_network, group_covar, strategy_covar, nodeSets = list(ACTORS,COMPONENTS))
#       
#       first_sim <- siena07(siena_algorithm, data = siena_data, effects = effects, returnDeps = TRUE)
#     }
#   }
#   
#   ## Summary Statistics
#   avg_degree <- apply(network_array_full, 3, function(mat) mean(rowSums(mat)))
#   cat("Average degrees by wave:", round(avg_degree, 2), "\n")
#   
#   ## Return Results
#   return(list(
#     groupMembership = groupMembership,
#     strategy = strategy,
#     networks = network_array_full,
#     avg_degree = avg_degree
#   ))
# }





siena_algorithm
siena_data
siena_effects

# SimulateBipartiteNetworks <- function(M, N, W, rate, dens, alt, ego, sim, groupMembership, strategy){
  # Simulates W consecutive network waves for a two-mode bipartite network
  # with M actors and N components (bipartite structure).
  #
  # Arguments:
  # - M: Number of actors (rows in the bipartite matrix).
  # - N: Number of components (columns in the bipartite matrix).
  # - W: Number of waves to simulate.
  # - rate: Rate parameter for simulation.
  # - dens: Outdegree effect for actors.
  # - alt, ego, sim: altX, egoX, and simX effect parameters for covariates.
  # - groupMembership: Constant covariate indicating actor group membership.
  # - strategy: Time-varying covariate for actor strategies.
  #
  # Returns:
  # - List containing covariates, simulated networks, and summary statistics.

  ## Create Initial Network
  set.seed(123)
  network_array <- array(NA, dim = c(M, N, W)) # Start with 2 waves
  network_array[, , 1] <- matrix(rbinom(M * N, 1, prob = 0.2), nrow = M, ncol = N)
  for (i in 2:W) {
    network_array[, , W] <- network_array[, , W-1] # Start with same network
    ##
    ##**TODO** check changing 1 dyad vs. identical networks; check sensitivity
  }

  ACTORS <- sienaNodeSet(M, nodeSetName = "Actors")
  COMPONENTS <- sienaNodeSet(N, nodeSetName = "Components")

  # Ensure it is a bipartite structure
  siena_network <- sienaDependent(network_array, type = "bipartite", nodeSet = c("Actors", "Components"), allowOnly = FALSE)

  ## Create Covariates
  group_covar <- coCovar(groupMembership) # Constant covariate
  strategy_covar <- varCovar(strategy)    # Time-varying covariate

  ## RSiena Data Object
  siena_data <- sienaDataCreate(siena_network, group_covar, strategy_covar, nodeSets = list(ACTORS,COMPONENTS))

  ## Define Effects
  effects <- getEffects(siena_data)
  effectsDocumentation(effects)
  effects <- includeEffects(effects, density, parameter = dens)
  effects <- includeEffects(effects, egoX, interaction1 = "strategy_covar") ## parameter = ego
  effects <- includeEffects(effects, sameEgoInDist2, interaction1 = "group_covar")

  
  
  
  ##**--------BEGIN MULTIWAVE NETWORK SIMULATION FUNCTION HERE----------**

  ## Algorithm Settings ## input argument or default saved algorithm (?)
  siena_algorithm <- sienaAlgorithmCreate(projname = "bipartite_sim", useStdInits = FALSE,
                                          cond = FALSE, nsub = 0, simOnly = TRUE, n3 = 100)

  ## Simulate the First Wave
  first_sim <- siena07(siena_algorithm, data = siena_data, effects = effects, returnDeps = TRUE)

  ## Simulate Consecutive Waves
  network_array_full <- array(0, dim = c(M, N, W)) # Full wave array
  network_array_full[, , 1:2] <- network_array      # First two waves

  for (w in 2:W) {
    # Extract simulated network from the previous wave
    last_wave <- matrix(0, nrow = M, ncol = N)
    sim_data <- first_sim$sims[[ length(first_sim$sims) ]][[ 1 ]][[ 1 ]]$`1`
    for (i in 1:nrow(sim_data)) {
      last_wave[sim_data[i, 1], sim_data[i, 2]] <- sim_data[i, 3]
    }

    network_array_full[, , w] <- last_wave  # Add the wave to the array

    # last_waves_w <- c(last_wave, last_wave)

    # Prepare the next wave
    if (w < W) {
      cat(sprintf('\n w=%s',w))
      siena_network <- sienaDependent(array(c(last_wave, last_wave), dim = c(M, N, W)),
                                      type = "bipartite", nodeSet = c("Actors", "Components"), allowOnly = FALSE)
      siena_data <- sienaDataCreate(siena_network, group_covar, strategy_covar, nodeSets = list(ACTORS,COMPONENTS))

      first_sim <- siena07(siena_algorithm, data = siena_data, effects = effects, returnDeps = TRUE)
    }
  }

  ## Summary Statistics
  avg_degree <- apply(network_array_full, 3, function(mat) mean(rowSums(mat)))
  cat("Average degrees by wave:", round(avg_degree, 2), "\n")

  ## Return Results
  return(list(
    groupMembership = groupMembership,
    strategy = strategy,
    networks = network_array_full,
    avg_degree = avg_degree
  ))
# }







###############################################################################
###                                                                         ###
###       Example Usage                                                     ###
###                                                                         ###
###############################################################################

# Define Parameters
M <- 10  # Number of actors
N <- 15  # Number of components
W <- 5   # Number of waves

rate <- 3
dens <- -1.5
alt <- 0.7
ego <- 0.5
sim <- 0.4

# Create Covariates
groupMembership <- sample(1:3, M, replace = TRUE) # Constant group membership
strategy <- matrix(sample(0:1, M * W, replace = TRUE), nrow = M, ncol = W) # Time-varying covariate

# Simulate Networks
bipartite_sim <- SimulateBipartiteNetworks(M, N, W, rate, dens, alt, ego, sim, groupMembership, strategy)

# View Results
print(bipartite_sim$avg_degree)































###############################################################################
###                                                                         ###
###         First main function: SimulateNetworks                           ###
###                                                                         ###
###############################################################################


SimulateNetworks <- function(n, M, rate, dens, rec, tt, c3,
                             Vaego, Vaalt, Vasim, Vbego, Vbalt, Vbsim){
  # Simulates M consecutive network waves, with n actors,
  # according to a stochastic actor-oriented model
  # with parameter values rate for rate,
  # dens for outdegree, rec for reciprocity,
  # tt for transitive triplets, c3 for 3-cycles,
  # an actor covariate Va with values alternating between 0 and 1,
  # with parameter values Vaego, Vaalt, Vasim
  # for egoX, altX, and simX with respect to Va,
  # and an actor covariate Vb with a standard normal distribution,
  # with parameter values Vbego, Vbalt, Vbsim
  # for egoX, altX, and simX with respect to Vb.
  ##
  # Create actor covariates
  V0 <- rep(0, n)
  V0[2*(1:(n %/% 2))] <- 1 # equal binary
  V1 <- rnorm(n, 0, 1)
  
  
  # Create initial 2-wave data to get a suitable data structure.
  # arbitrarily, this initial network has an expected average degree of 3
  X0 <- matrix(rbinom(n*n,1,3/(n-1)),n,n)
  diag(X0) <- 0
  X1 <- X0
  # but X0 and X1 should not be identical for use in sienaDependent
  X0[1,2] <- 0
  X0[2,1] <- 1
  X1[1,2] <- 1
  X1[2,1] <- 0
  XX <- array(NA,c(n,n,2))
  XX[,,1] <- X0
  XX[,,2] <- X1
  
  
  # With this data structure, we now can create the data.
  Va <- coCovar(V0)
  Vb <- coCovar(V1)
  X   <- sienaDependent(XX, allowOnly = FALSE)
  InitData <- sienaDataCreate(X, Va, Vb)
  InitEff0 <- getEffects(InitData)
  
  
  # sink to avoid printing to the screen
  sink("eff.txt")
  
  
  # Specify the parameters.
  # The rate parameter is first multiplied by 10,
  # which will be used only to get from the totally random network XX[,,1] = X0
  # to the network that will be the simulated first wave.
  InitEff0 <- setEffect(InitEff0, Rate, type="rate", initialValue = 10*rate)
  InitEff0 <- setEffect(InitEff0, density, initialValue = dens)
  InitEff0 <- setEffect(InitEff0, recip, initialValue = rec)
  InitEff0 <- setEffect(InitEff0, transTrip, initialValue = tt)
  InitEff0 <- setEffect(InitEff0, cycle3, initialValue = c3)
  InitEff0 <- setEffect(InitEff0, egoX, interaction1="Va", initialValue = Vaego)
  InitEff0 <- setEffect(InitEff0, altX, interaction1="Va", initialValue = Vaalt)
  InitEff0 <- setEffect(InitEff0, simX, interaction1="Va", initialValue = Vasim)
  InitEff0 <- setEffect(InitEff0, egoX, interaction1="Vb", initialValue = Vbego)
  InitEff0 <- setEffect(InitEff0, altX, interaction1="Vb", initialValue = Vbalt)
  InitEff0 <- setEffect(InitEff0, simX, interaction1="Vb", initialValue = Vbsim)
  
  
  # The parameter given for n3 should be larger than sum(InitEff0$include)
  nthree <- sum(InitEff0$include)	+ 5
  InitAlg <- sienaAlgorithmCreate(projname="Init", useStdInits=FALSE,
                                  cond=FALSE, nsub=0, n3=nthree, simOnly=TRUE)
  
  
  # Simulate the first wave.
  InitSim   <- siena07(InitAlg, data=InitData, eff=InitEff0,
                       returnDeps=TRUE, batch=TRUE, silent=TRUE)
  
  
  # Now prepare for simulating waves 2 to M.
  # Create empty result network.
  Xs <- array(0, dim=c(n,n,M))
  # The rate parameter value from the function call is reinstated in InitEff.
  InitEff <- InitEff0
  InitEff <- setEffect(InitEff, Rate, type="rate", initialValue = rate)
  
  
  sink()
  for (m in 1:M){
    # Note that we start this loop with a previously simulated network.
    # Transform the previously simulated network
    # from edge list into adjacency matrix
    XXsim <- matrix(0,n,n)
    nsim  <- InitAlg$n3
    XXsim[ InitSim$sims[[nsim]][[1]]$X[[1]][,1:2] ]  <- InitSim$sims[[nsim]][[1]]$X[[1]][,3]
    # Put simulated network into the result matrix.
    Xs[,,m] <- XXsim
    # Put simulated network in desired places for the next simulation
    XX[,,2] <- XX[,,1] # used only to get the data structure
    XX[,,1] <- XXsim
    
    
    if (m < M){
      # The following is only to prevent the error that would occur
      # in the very unlikely event XX[,,1] == XX[,,2].
      if (identical(XX[,,1], XX[,,2])){XX[1,2,1] <- 1 - XX[1,2,2]}
      # Specify the two-wave network data set starting with XX[,,1].
      X <- sienaDependent(XX, allowOnly = FALSE)
      # Simulate wave m+1 starting at XX[,,1] which is the previous XXsim
      InitData  <- sienaDataCreate(X, Va, Vb)
      InitSim <- siena07(InitAlg, data=InitData, eff=InitEff,
                         returnDeps=TRUE, batch=TRUE, silent=TRUE)
    }
    
    
  }
  
  
  # Present the average degrees to facilitate tuning the outdegree parameter
  # to achieve a desired average value for the average degrees.
  cat("Average degrees ", round(colSums(Xs,dims=2)/n, digits=2), "\n")
  
  
  # Result: simulated data set; covara and covarb are vectors of length n;
  # networks is an array of dimension nxnxM
  list(covara = V0, covarb = V1, networks = Xs)
}



