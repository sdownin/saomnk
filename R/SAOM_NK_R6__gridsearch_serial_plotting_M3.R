rm(list=ls())  ## uncomment to clear environ before run
TIME_START <- Sys.time()
########################
##
##   SAOM-NK Runscript:
##
##   Serial Plotting
##
##
#######################
library(readxl)
library(ggrain)
library(ggplot2)
library(plyr); library(dplyr)


##****************************************************##
## GRID SEARCH DIRECTORY
grid.dir <- '__SAOMNK_gridsearch_input_start0_M9__'
## INPUT FILE
# grid.input.file <- 'SAOMNK_gridsearch_input_start0_M3.xlsx'
## ITERATIONS
ITERATIONS <- 6000
##****************************************************##


###############  Load R DEPENDENCIES ############################
## Directories
dir_ext <- 'D:\\Search_networks'
dir_r <- 'C:\\Users\\sdr8y\\OneDrive - University of Missouri\\Research\\Search_networks\\SaoMNK\\R'
dir_proj <- file.path(dir_ext, grid.dir )
dir_proj_results <- file.path(dir_proj, 'grid_results')


## Biparite Environment Search Simulation Class
SaomNkRSienaBiEnv      <- source(file.path(dir_r, 'SAOM_NK_R6_model.R'))$value
# SaomNkRSienaBiEnv_base <- source(file.path(dir_proj, 'SAOM_NK_R6_base_model.R'))$value


## Working director
setwd(dir_proj)



## Simulation files
simfiles <- dir(dir_proj_results, pattern = '\\.rds$')



#######################################
## FUNCTIONS
#######################################
##
## Plot the bipartite system (social network, bipartite network, component matrix)
##
plot_bipartite_system_from_mat <- function(bipartite_matrix, RSIENA_ITERATION) {
  # #
  # RSIENA_ITERATION <- length(self$rsiena_model$sims)
  N <- ncol(bipartite_matrix)
  M <- nrow(bipartite_matrix)
  #
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
  #
  component_labels <- generate_component_labels(self$N)
  # 1. Bipartite network plot using ggraph with vertex labels
  ig_bipartite <- igraph::graph_from_biadjacency_matrix(bipartite_matrix, directed = T, weighted = T)
  # E(ig_bipartite)$weight <- 1
  projs <- igraph::bipartite_projection(ig_bipartite, multiplicity = T, which = 'both')
  ig_social <- projs$proj1
  ig_component <- projs$proj2
  # Set node attributes for shape and label
  # V(ig_bipartite)$type <- bipartite_mapping(ig_bipartite)$type
  V(ig_bipartite)$shape <- ifelse(V(ig_bipartite)$type, "square", "circle")
  V(ig_bipartite)$color <- ifelse(V(ig_bipartite)$type, "lightblue", "darkorange")
  V(ig_bipartite)$label <- c(as.character(1:self$M), as.character(component_labels))  # Numeric labels for actors
  # stop('DEBUG')
  V(ig_bipartite)$node_text_color <- ifelse(V(ig_bipartite)$type, 'black', 'lightblue') 
  #
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
  #
  component_df <- melt(component_matrix)
  colnames(component_df) <- c("Component1", "Component2", "Interaction")
  # Replace component numbers with labels in the heatmap data frame
  component_df$Component1 <- factor(component_df$Component1, labels = component_labels)
  component_df$Component2 <- factor(component_df$Component2, labels = component_labels)
  #
  heatmap_plot <- ggplot(component_df, aes(x = Component1, y = Component2, fill = Interaction)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    labs(title = "[Proj2] Component Heatmap\n(Epistatic Interactions)", x = "Component 1", y = "Component 2") +
    theme_minimal() +
    theme(legend.position = "bottom")
  # Calculate average degree (K) for the social space and component interaction space
  avg_degree_social <- mean(igraph::degree(ig_social))
  avg_degree_component <- mean(igraph::degree(ig_component))
  # Create a main title using sprintf with simulation parameters
  main_title <- sprintf("Simulation Parameters:\niter = %d, N = %d, M = %d, BI_PROB = %.2f, K_A = %.2f, K_C = %.2f", 
                        RSIENA_ITERATION, N, M, self$BI_PROB, avg_degree_social, avg_degree_component)
  return(grid.arrange(
    social_plot, bipartite_plot, heatmap_plot, ncol = 3,
    top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold"))
  ))
}

#########################################################
## FILE MANAGEMENT
#########################################################
## RDS files saved in grid_results directory within the gridsearch job directory
dir_proj_plotting <- file.path(dir_proj, 'img')
if (!dir.exists(dir_proj_plotting)) {
  dir.create(dir_proj_plotting)
}



#########################################################
## GRIDSEARCH INDEX
#########################################################
# full index 
paramgrid_all <- readRDS(file.path(dir_proj, 'paramgrid_index_ALL.rds'))
# written files to index for completed scenarios
paramgrid <- read.csv(file.path(dir_proj, 'paramgrid_index.csv'))


#########################################################
## BEGIN OUTPUT TXT FILE
#########################################################
sink(file = file.path(dir_proj_plotting, sprintf('output_plotting_%s.txt',grid.dir) ) , type='output' )

cat(sprintf('\nJob start time: %s\n', as.character(TIME_START)))
cat(sprintf('\nMain loop start time: %s\n', as.character(Sys.time())))


cat(sprintf('\ndir_proj: %s\n', dir_proj))



cat('\nparamgrid dimensions:\n')
print(dim(paramgrid))

cat('\nMissing paramgrid scenarios from paramgrid:\n')
print( paramgrid$uuid[ ! paramgrid$uuid %in% paramgrid_all$uuid ])


#########################################################
##
## MAIN LOOP
##
#########################################################
# output_foreach <- foreach(i = 1:nrow(paramgrid), .combine = 'rbind', .packages = "base", .errorhandling = "pass", .verbose=T) %dopar% {
#
for (i in 1:length(simfiles)) {

  ## Simulation scenario i (inputs list for next simulation) 
  simfile <- simfiles[i]
  
  parts <- strsplit(gsub('\\.rds','',simfile), '_')[[1]]
  uuid <- parts[length(parts)]
  
  ## Load the Simulation 
  sim <- readRDS(file.path(dir_proj_results, simfile ))
  
  ## alias during debugging before finalizing package
  self <- sim
  
  ## get output file name base of this scenario for log and data output, etc.
  # sim_filebase <-  sprintf('%s_%s',grid.dir, scenario$uuid)
  sim_filebase <-  sprintf('%s', uuid )
  
  

  
  ##==================================================
  ## actor strategy
  # if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
  #   stop("Actor Strategy self$strat_1_coCovar not set.")
  actor_strat_lvls <- as.factor( self$get_actor_strategies() )
  actor_strats <- rep(actor_strat_lvls, times = round(self$M/length(actor_strat_lvls)))
  #
  nstep <- sum(!self$chain_stats$stability)
  ## levels order
  coveffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
  covparams <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$parameter)
  covfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)ifelse(x$fix,'','(var)'))
  #
  strateffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
  stratparams <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$parameter)
  stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)ifelse(x$fix,'','(var)'))
  # stratfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)substr(as.character(x$fix),1,1))
  #
  structeffs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$effect)
  structparams <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$parameter)
  structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)ifelse(x$fix,'','(var)'))
  # structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)substr(as.character(x$fix),1,1))
  covDvTypes <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$interaction1)
  componentDV_ids <- grep('self\\$component_\\d{1,2}_coCovar', covDvTypes) ## ex: "self$component_1_coCovar"
  stratDV_ids     <- grep('self\\$strat_\\d{1,2}_coCovar', covDvTypes ) ## ex: "self$strat_1_coCovar"     
  #
  sim_title_str <- sprintf(
    'Environment: Actors (M) = %s, Components (N) = %s, Init.Prob. = %.2f\nActor Strategy:  %s\nComponent Payoff:  %s\nStructure:  %s', 
     self$M, self$N, self$BI_PROB,
     paste( paste(paste(strateffs[stratDV_ids], stratparams[stratDV_ids], sep='= '), stratfixs[stratDV_ids], sep='' ), collapse = ';  '),
     paste( paste(paste(strateffs[componentDV_ids], stratparams[componentDV_ids], sep='= '), stratfixs[componentDV_ids], sep=''), collapse = ';  '),
     paste( paste(paste(structeffs, structparams, sep='= '), structfixs, sep=''), collapse = ';  ')
  )
  #
  compoeffs   <- coveffs[ componentDV_ids ]
  compoparams <- covparams[ componentDV_ids ]
  compofixs   <- covfixs[ componentDV_ids ]
  #
  strateffs   <- coveffs[ stratDV_ids ]
  stratparams <- covparams[ stratDV_ids ]
  stratfixs   <- covfixs[ stratDV_ids ]
  efflvls <- c('utility', strateffs, compoeffs, structeffs)
  ##=================================================

  
  
  ##****************************************##
  ## II. PLOTTING STARTS
  ##****************************************##
  
  ##--------------------------------------------------------------
  ## 1.  Utility by Strategy Chain Overview
  ##--------------------------------------------------------------
  ## 1.1 Utility Strategy Summary timeseries and density of whole chain
  plt1.1 <- self$search_rsiena_multiwave_plot('utility_strategy_summary',
    thin_factor = 1,
    return_plot=TRUE#,
    # plot_dir = dir_proj_plotting,  ## Not yet implemented when gridsearch was run
    # plot_file = sprintf('%s_strategy_summary_thin1.png',sim_filebase)   ## Not yet implemented when gridsearch was run
  )[[1]] ## get first item of list of plots
  plotname1.1 <- sprintf('%s_plot1-1_strat_overview.png',sim_filebase) 
  ggsave(filename=file.path(dir_proj_plotting, plotname1.1), plt1.1, 
         height = 7, width = 9, dpi = 600, units = 'in')
  
  ## 1.2 thinned / 10
  plt1.2 <- self$search_rsiena_multiwave_plot('utility_strategy_summary',
    thin_factor = 10,
    return_plot=TRUE#,
    # plot_dir = dir_proj_plotting,  ## Not yet implemented when gridsearch was run
    # plot_file = sprintf('%s_strategy_summary_thin10.png',sim_filebase)   ## Not yet implemented when gridsearch was run
  )[[1]] ## get first item of list of plots
  plotname1.2 <- sprintf('%s_plot1-2_strat_overview.png',sim_filebase) 
  ggsave(filename=file.path(dir_proj_plotting, plotname1.2), plt1.2, 
         height = 7, width = 9, dpi = 600, units = 'in')
  
  
  ##--------------------------------------------------------------
  ## 2. Strategy Signals
  ##--------------------------------------------------------------
  #
  actor_stats_df <- self$actor_stats_df
  actor_stats_df$strategy <- actor_strats[ actor_stats_df$actor_id ]
  #
  actor_util_df <- self$actor_util_df  %>% mutate(
    effect_id=NA, 
    value=utility, 
    effect_name='utility', 
    strategy= actor_strats[ actor_id ], 
    utility=NULL
  )
  
  ## Add utility as extra 'effect' 
  act_effs <- actor_stats_df %>% bind_rows( actor_util_df ) 
  act_effs$effect_name <- factor(act_effs$effect_name, levels=efflvls)
  
  #plot signals
  thinplt2 <- 3
  plt2 <- act_effs %>% filter(chain_step_id %% thinplt2 == 0) %>%
    ggplot(aes(x=chain_step_id, y=value, color=strategy,fill=strategy, linetype=strategy)) + 
    geom_point(alpha=.04, shape=1)  + 
    geom_smooth(method='loess', alpha=.1) +
    facet_grid(effect_name ~ ., scales='free_y') +
    theme_bw() + theme(legend.position = 'left') +
    ggtitle(sprintf('Actor Utility Signals Decomposition by Effects:\n%s',sim_title_str))
  plotname2 <- sprintf('%s_plot2_strat_signals_thin%s.png',sim_filebase, thinplt2) 
  ggsave(filename=file.path(dir_proj_plotting, plotname2), plt2, 
         height = 12, width = 8, dpi = 600, units = 'in')
  
  
  
  ##--------------------------------------------------------------
  ## 3.  K4
  ##--------------------------------------------------------------
  # 3. K4 Panel Structural Distribution
  plt3 <- self$search_rsiena_multiwave_plot('K_4panel',
      thin_factor = 3,
      return_plot=TRUE,
      # plot_dir = dir_proj_plotting,
      # plot_file = sprintf('%s_K4panel_thin1.png',sim_filebase) 
  )[[1]] ## get first item of list of plots
  plotname3 <- sprintf('%s_plot2_K4panel_thin3.png',sim_filebase) 
  ggsave(filename=file.path(dir_proj_plotting, plotname3), plt3, 
         height = 8, width = 8, dpi = 600, units = 'in')
  
  
  
  
  
  #######################################################################
  
  #
  config_param_vals <- c(
    unlist(sapply(self$config_structure_model$dv_bipartite$effects, function(x) x$parameter)),
    unlist(sapply(self$config_structure_model$dv_bipartite$coCovars, function(x) x$parameter)),
    unlist(sapply(self$config_structure_model$dv_bipartite$varCovars, function(x) x$parameter))
  )
  fixed_params <- c(
    unlist(sapply(self$config_structure_model$dv_bipartite$effects, function(x) x$fix)),
    unlist(sapply(self$config_structure_model$dv_bipartite$coCovars, function(x) x$fix)),
    unlist(sapply(self$config_structure_model$dv_bipartite$varCovars, function(x) x$fix))
  )
  #
  #UTILITY
  theta   <- self$rsiena_model$theta
  ## replace fixed theta param with given param values (instead of zero default when effect is fixed)
  if (any(fixed_params)) theta[ fixed_params ] <- config_param_vals[ fixed_params ]
  
  
  
  ## PLOT ACTOR UTILITY LANDSCAPE
  utilist <- list()
  act_counterfacts <- list()
  pltlist <- list()
  sys_snapshots <- list()
  
  # step_ids <- round(c(1,8,64,512, self$rsiena_model$n3/2, self$rsiena_model$n3))\
  pd_nsteps <- self$M*self$N
  acpds <- seq(1, self$rsiena_model$n3 - 1, by= pd_nsteps )
  plot_npds <- 6
  step_ids <- c(acpds[1:plot_npds])
  for (step_id in step_ids) {
    # step_id <- 1
    cat(sprintf('\nstep %s: actors ', step_id))
    
    ## get network for this step_id in the decision chain
    bi_env_mat_step <- self$bi_env_arr[, , step_id]
    
    ## Save snapshot plot
    sys_snapshots[[sprintf('step_%s',step_id)]] <- plot_bipartite_system_from_mat(bi_env_mat_step, step_id)
    
    ## Actor loop
    for (ii in 1:self$M) {
      ##**ACTOR ii DECISION PERSPECTIVE in scenario i** 0000000000000000000000000000000000
      cat(sprintf(' %s ', ii))
      
      ##**TODO**
      ## ALL actor-component counterfactual configuations for Actor i  (2^N rows by N cols)
      iland = expand.grid(lapply(1:self$N, function(x) 0:1 ))
      iland_config_step_row_id <- which(apply(iland, 1, function(x) all(x == bi_env_mat_step[ii,]) ))
      
      tmpmat <- bi_env_mat_step
      ## ifit dimensions [ M, 2^N ]
      ifit <- apply(iland, 1, function(x){
        tmpmat[ii,] <- x  ## set counterfactual actor-component configuration
        return( self$get_struct_mod_stats_mat_from_bi_mat( tmpmat ) %*% theta ) ## compute utility vector
      }) 
      
      #
      ids.max <- which(ifit == max(ifit), arr.ind = TRUE) 
      fits.max <- ifit[ ids.max ]
      nmax <- length(fits.max)
      
      act_counterfacts[[ ii ]] <- list(
        ifit = ifit, #Given all other ties, Actor i's configurations applied to utility func for all actors
        ids.max = ids.max,
        fits.max = fits.max,
        nmax = nmax
      )
      
      ## distances of each counterfactual configuration 
      dist_counterfac <- iland - bi_env_mat_step[ ii, ]
      #
      z <- apply( dist_counterfac, 1, function(x) {
        c(dist=sum(x!=0),`drop`=sum(x==-1), `nochange`=sum(x==0), `add`=sum(x==1), 
          counterfac=paste(x, collapse = '|'), start=paste(bi_env_mat_step[ ii, ], collapse = '|') )
      })
      # z <- z[,order(z['d',], decreasing = F)]
      ##**TODO**
      ##**INFORMATION LOSS**
      ##  This uses only ego's own counterfactual fits (landscape)
      ##  but the corresponding other actor's affected fits (landscapes are not currently used )
      z <- rbind(z, utility_ego=ifit[ ii , ] )
      z <- rbind(z, utility_alter_mean= colMeans(ifit[ -ii , ], na.rm = T) )
      z <- rbind(z, utility_alter_sd  = apply(ifit[ -ii , ], 2, function(x) sd(x, na.rm = T)) )
      
      wdf <- as.data.frame( t(z) )
      
      actfit_long <- wdf %>% 
        mutate(
          dist = as.numeric(dist),
          drop = as.numeric(drop),
          nochange = as.numeric(nochange),
          add = as.numeric(add) ,
          utility_ego=as.numeric(utility_ego),
          utility_alter_mean = as.numeric(utility_alter_mean),
          utility_alter_sd = as.numeric(utility_alter_sd)
        ) %>%
        pivot_longer(cols = c( drop:add, utility_ego:utility_alter_sd )) 
      
      
      plt <- actfit_long %>% ggplot(aes(x=value))+ 
        # geom_density() + 
        geom_histogram() +
        facet_grid(dist ~ name, scales='free_x') + theme_bw() + 
        ggtitle(sprintf('Utility Transition Paths: Step %s, Actor %s', step_id, ii ))
      
      step_actor_key <- sprintf('%s|%s',step_id, ii)
      
      pltlist[[step_actor_key]] <- plt

      
      utilist[[step_actor_key]] <- actfit_long %>% 
        filter(name %in% c('utility_ego','utility_alter_mean')) %>%
        mutate(actor_id=ii, chain_step_id = step_id, strategy=actor_strats[ii])
      
    } ##/end actor loop ii
    
  } ##/end decision chain step loop step_id

  #
  utildf <- data.table::rbindlist(utilist, use.names = T, idcol = 'step_actor_key')
  
  
  
  
  ##--------------------------------------------------------------
  ## 4.  Actor Utility Transition Paths
  ##--------------------------------------------------------------
  # Arrange plots across multiple pages
  plt4 <- marrangeGrob(pltlist, nrow = 1, ncol = 1)
  
  # Save to a PDF
  plotname4pdf <- sprintf('%s_actor_pd_dist2peaks.pdf', sim_filebase)
  ggsave(filename = file.path(dir_proj_plotting, plotname4pdf), plt4, width = 8, height = 8)

  
    
  
  ##--------------------------------------------------------------
  ## 5. 
  ##--------------------------------------------------------------
  varname <- 'utility_ego'
  ndists_max <- 9
  utildfvar <- utildf%>%filter(name==varname, dist <= ndists_max )
  plt5 <- ggplot(utildfvar, aes(x=value, color=strategy, fill=strategy, linetype=strategy)) + 
    # geom_density(alpha=.05, size=.9) +
    geom_histogram(alpha=.15, position = 'identity')+
    geom_vline(aes(xintercept = mean, color=strategy, linetype=strategy),
               data = utildfvar %>%
                 group_by(dist,strategy,chain_step_id)%>%
                 dplyr::summarize(mean=mean(value))
    ) +
    facet_grid(chain_step_id ~ dist)+ theme_cowplot() + 
    # ggtitle('Fitness Landscapes by Strategy over Time (Actor Counterfactual Utility)')
    ggtitle(sprintf('How far are the peaks over time? ( var = %s )\nCounterfactual Utility by search distance (cols) and step count (rows)', varname))
  plotname5 <- sprintf('%s_strategy_dist2peaks.png',sim_filebase) 
  ggsave(filename=file.path(dir_proj_plotting, plotname5), plt5, 
         height = 8, width = 10, dpi = 600, units = 'in')
  
  
  ##--------------------------------------------------------------
  ## 6. 
  ##--------------------------------------------------------------
  ## ACTOR POINT SERIES BY PERIOD (lines by actor; colors by strategy)
  actutilstrat <- actor_util_df %>% filter(effect_name=='utility') %>%
    mutate(
      period = ceiling(chain_step_id / (self$M*self$N)),
      utility = value
    )
  plt6 <- ggplot(actutilstrat, aes(x=chain_step_id, y=utility, color=strategy, linetype = actor_id)) + 
    # geom_density(alpha=.05, size=.9) +
    geom_point(alpha=.4, shape=1)+
    # facet_grid(actor_id ~ ., scales=)+
    geom_vline(xintercept = step_ids) +
    geom_line() +
    # geom_smooth(method='loess') +
    xlim(c(0,self$M * self$N * plot_npds )) +
    theme_bw() + 
    ggtitle('ACTOR UTILITY POINT-LINE SERIES colored by STRATEGY vline by PERIOD')
  plotname6 <- sprintf('%s_util_pointline_strat_pd.png',sim_filebase) 
  ggsave(filename=file.path(dir_proj_plotting, plotname6), plt6, 
         height = 7, width = 9, dpi = 600, units = 'in')
  
  
  
  
  
  ##--------------------------------------------------------------
  ## 7.  Raincloud plots by period
  ##--------------------------------------------------------------
  NtimesM <- self$M * self$N
  nstrats <- length(levels(actor_strats))
  
  plt7 <- ggplot(actutilstrat %>% filter(period <= plot_npds) %>%
                           mutate(period=factor(period, levels=as.character(1:plot_npds))), #%>%   mutate(strategy_pd=unique(paste(self$strat_1_coCovar[actor_id],period,collapse = '_'))),
                         aes(x=period, y=utility, fill = strategy, color=strategy)) + #shape=actor_id, 
    geom_rain(alpha = .15, rain.side = 'r',  # id.long.var = 'actor_pd',
              boxplot.args.pos = list(width = .1,
                                      position = ggpp::position_dodgenudge(width = .1,
                                                                           x = c(.13, .13, # t1 old, t1 young
                                                                                 .13, .13,
                                                                                 .13, .13,
                                                                                 .13, .13,
                                                                                 .13, .13,
                                                                                 .13, .13)
                                      )
              ),
              violin.args.pos = list(width = 1.1,
                                     position = position_nudge(
                                       x = c(rep(.2, 256*2), rep(.2, 256*2),# t1
                                             rep(.2, 256*2), rep(.2, 256*2), # t2
                                             rep(.2, 256*2), rep(.2, 256*2),
                                             rep(.2, 256*2), rep(.2, 256*2),
                                             rep(.2, 256*2), rep(.2, 256*2),
                                             rep(.2, 256*2), rep(.2, 256*2)
                                       )
                                     )
              )
    ) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = strategy, color = strategy)) +
    stat_summary(fun = mean, geom = "point",
                 aes(group = strategy, color = strategy)) +
    geom_hline(yintercept = 0, linetype=2) +
    guides(fill = 'none', shape='none') + 
    theme_classic() + 
    ggtitle(sprintf('Strategy Transition Paths: Utility Distributions by Period (%s steps)\n%s',
                    self$M*self$N, sim_title_str))
  #
  plotname7 <- sprintf('%s_util_raincloud_pd.png',sim_filebase) 
  ggsave(filename=file.path(dir_proj_plotting, plotname7), plt7, 
         height = 7, width = 9, dpi = 600, units = 'in')
  
  
  
  
  
  
  ##--------------------------------------------------------------
  ## 8.  Visualize network
  ##--------------------------------------------------------------
  
  step_end <- dim(self$bi_env_arr)[3] 
  sys_snapshots[[sprintf('step_%s',step_end)]] <- plot_bipartite_system_from_mat(self$bi_env_arr[,,step_end], step_end)
  
  # Arrange plots across multiple pages
  plt8 <- marrangeGrob(sys_snapshots, nrow = 1, ncol = 1)
  
  # Save to a PDF
  plotname8pdf <- sprintf('%s_system_snapshots.pdf', sim_filebase)
  ggsave(filename = file.path(dir_proj_plotting, plotname8pdf), plt8, width = 13, height = 7)
  
  
  
  
  
  
  ##---------------
  plt <- NULL
  plt1.1 <- NULL
  plt1.2 <- NULL
  plt2 <- NULL
  plt3 <- NULL
  plt4 <- NULL
  plt5 <- NULL
  plt6 <- NULL
  plt7 <- NULL
  plt8 <- NULL
  sys_snapshots <- NULL
  gc();
  ##---------------
} ##**end main loop**



## Ending time
TIME_END <- Sys.time()
cat(sprintf('\nEnd time: %s\n', TIME_END))
cat(sprintf('\nElapsed time: %s\n', TIME_END - TIME_START ))
print( TIME_END - TIME_START )


cat('\n\nSIMULATION COMPLETED SUCCESSFULLY.\n\n')


## return output to console
sink()
##**END**




