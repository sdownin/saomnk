rm(list=ls())  ## uncomment to clear environ before run
########################
##
##   SAOM-NK Runscript:
##
##   1. Random Sensitivity
##
##
#######################
library(plyr)
library(dplyr)
library(uuid)
library(ggpubr)
library(data.table)
library(dunn.test)
library(forcats)



##****************************************************##
Ms <- c(3, 18)
Ns <- c(128,512)
BI_PROBs <- c(0, 1)
## ITERATIONS
MAX_ITERATIONS <- 5000
## Replications
nreps <- 30
##****************************************************##
#
.ts <- gsub('\\.','',as.character(as.numeric(Sys.time())))
#
# .uuid <- UUIDgenerate()
##  DIRECTORY
grid.dir <- sprintf('_SAOMNK_sensitivity_AND_grid_rand%s_%s', nreps , .ts )


###############  Load R DEPENDENCIES ############################
## Directories
dir_ext <- 'D:\\Search_networks'
dir_r <- 'C:\\Users\\sdr8y\\OneDrive - University of Missouri\\Research\\Search_networks\\SaoMNK\\R'
dir_proj <- file.path( dir_ext, grid.dir )
dir_proj_results <- file.path(dir_proj, 'sensitivity_results')

## script run directory
if (!dir.exists(dir_proj)) {
  dir.create(dir_proj)
}

## default settings: Users do not change; TODO: implment within restricted class attributes
DV_NAME <- 'self$bipartite_rsienaDV'


###############  Load R6 Class DEPENDENCIES ############################
## Biparite Environment Search Simulation Class
SaomNkRSienaBiEnv <- source(file.path(dir_r, 'SAOM_NK_R6_model.R'))$value
# ## RSiena search Class
# SaomNkRSienaBiEnv_search_rsiena <- source(file.path(dir_proj, 'SAOM_NK_R6_search_rsiena_model.R'))$value
###########




#########################################################
## FILE MANAGEMENT
#########################################################
## RDS files saved in grid_results directory within the gridsearch job directory
if (!dir.exists(dir_proj_results)) {
  dir.create(dir_proj_results)
} else {
  ## if img dir exists, remove simufiles from this run list by UUIDs in existing plots
  plotted <- dir(dir_proj_results, pattern = '\\.png')
  skipuuids <- unique(sapply(strsplit(plotted, '_'), function(x) x[1] ))
  dropids <- unname(sapply(skipuuids, function(x) grep(x, simfiles)))
  if(length(dropids))
    simfiles <- simfiles[ -dropids ] ## drop files with UUIDs already plotted
}







## Working director
setwd(dir_proj)


## Start output
sink(file = file.path(dir_proj, sprintf('%s_out.txt', .ts)), type = 'output')


# #################### FUNCTIONS ####################################
# getGridRunNameFromScenario <- function(scenario, environ_seed_params) {
#   varnames <- names(scenario)
#   vals <- unname(scenario)
#   outvec <- c()
#   for (i in 1:length(varnames)) {
#     if (varnames[i] %in% environ_seed_params) {
#       outvec <- c(outvec, paste(varnames[i], vals[i], sep='='))
#     }
#   }
#   outstr <- paste(outvec, collapse='_')
#   return(sprintf('_%s_', outstr))
# }






strategies <- list(
  egoX   =  c(-1,0, 1), #c(0),
  inPopX =  c(1,0, -1)  #c(0),
)



##_-----------------------------------------------------

# ITERATIONS <- 4000
set.seed(123)
seeds <- sample(1:999999, nreps, replace = F)

actlist <- list()
act_effs_list <- list()

dunntestlist <- list()

utlist <- list()
utdifflist <- list()
outlist <- list()
statlist <- list()

for (l in 1:length(Ms)) {
  M <- Ms[l]
  for (k in 1:length(Ns)) {
    N <- Ns[k]
    for (j in 1:length(BI_PROBs)) {
      BI_PROB <- BI_PROBs[j]
      
      #
      
      
      
      #########################################################
      ## GRIDSEARCH INDEX
      #########################################################
      # # full index 
      # paramgrid_all <- readRDS(file.path(dir_proj, 'paramgrid_index_ALL.rds'))
      # # written files to index for completed scenarios
      # paramgrid <- read.csv(file.path(dir_proj, 'paramgrid_index.csv'))
      
      # paramgrid$uuid <- sapply(1:nrow(paramgrid), uuid::UUIDgenerate )
      
      
      
      
      for (i in 1:length(seeds)) 
      {
        
        rep_uuid <- UUIDgenerate()
        
        seed_str <- sprintf('%s-%s-%s-%s', seeds[i], M, N, BI_PROB)
        
        cat(sprintf('\n\n---------------------- SEED %s --------------------\n\n', seeds[i]))
        
        
        environ_params <- list(
          M = M,        ## Actors
          N = N,       ## Components
          BI_PROB = BI_PROB, ## Environmental Density (DGP hyperparameter)
          component_matrix_start = 'rand', ##**TODO** Implement: 'rand','modular','semi-modular',...
          rand_seed = 1234,
          visualize_init = F,
          name = sprintf('%s__%s', rep_uuid, seed_str)
        )
        
        ## 2.b. Component Payoffs vector
        component_payoffs <-  runif(environ_params$N, min = 0, max = 1)
        # component_payoffs <-  rpois(n = environ_params$N, lambda = 1.5)  
        ## 2. Strategies sets the objective function as a linear combination of network stats across DVs
        #
        actor_strats_list <- lapply(strategies, function(strat) rep(strat,  environ_params$M/length(strat)) )
        
        #
        structure_model <- list(
          dv_bipartite = list(
            name = 'self$bipartite_rsienaDV',
            effects = list( ##**STRUCTURAL EFFECTS -- dyadic/network endogeneity sources**
              list(effect='density', parameter= -1.5, fix=T, dv_name=DV_NAME), ##interaction1 = NULL
              list(effect='inPop',   parameter=  .1,  fix=T, dv_name=DV_NAME), #interaction1 = NUL
              list(effect='outAct',  parameter=  .1, fix=T, dv_name=DV_NAME)
            ),
            coCovars = list( ##**STRATEGY -- MONADIC CONSTANT COVARIATE EFFECTS **
              list(effect='altX',   parameter= 1, fix=T,dv_name=DV_NAME,interaction1='self$component_1_coCovar', x= component_payoffs ),
              list(effect='egoX',   parameter= 1,fix=T,dv_name=DV_NAME,interaction1='self$strat_1_coCovar', x= actor_strats_list[[1]] ), #interaction1 = NULL
              list(effect='inPopX', parameter= 1, fix=T,dv_name=DV_NAME,interaction1='self$strat_2_coCovar', x= actor_strats_list[[2]] )
            ),
            varCovars = list() ##**MONADIC TIME-VARYING COVARIATE EFFECTS -- DYNAMIC STRATEGY PROGRAMS**
          )
        )
        
        ##****************************************##
        ## II. SIM ANALYSIS 
        ##****************************************##
        # ## INIT SIM ENVIRONMENT: 
        # env1 <- SaomNkRSienaBiEnv$new(environ_params)
        
        ## INIT SIM ENVIRONMENT: 
        env1 <- SaomNkRSienaBiEnv$new(environ_params)
        
        # pd_steps <- env1$M * env1$N
        pd_steps <- 2 * env1$M 
        
        ITERATIONS <- min( 50 * pd_steps, MAX_ITERATIONS)
        
        ## 1.1. Search 1. SHORT Run
        env1$search_rsiena_multiwave_run(
          structure_model, 
          waves=1, ##id='_short_run',
          iterations = ITERATIONS, #4 * env1$M * env1$N, 
          rand_seed = seeds[i]
        )
        # Process results
        env1$search_rsiena_multiwave_process_results()
        #
        # env1$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T) #+  geom_vline(xintercept=seq(0,440,8*9))
        
        
        
        ##-----------------------
        self <- env1
        
        actor_strats <- self$get_actor_strategies()
        
        tmpdf <- self$actor_util_df %>% mutate(
          seed=seeds[i],
          period = round(ceiling(chain_step_id/ pd_steps))
        )
        tmpdfdiff <- self$actor_util_diff_df %>% mutate(
          seed=seeds[i],
          period = round(ceiling(chain_step_id/ pd_steps))
        )
        
        ## loop over types of effects in structure_model to add to dataframe
        for (k in 1:length(structure_model$dv_bipartite)) {
          # print(k)
          efflist <- structure_model$dv_bipartite[[ k  ]]
          if ( ! 'list' %in% class(efflist) || !length(efflist)) 
            next 
          for (h in 1:length(efflist)) {
            item <- efflist[[ h ]]
            # print(h)
            tmpdf[[item$effect]] <- item$parameter
            tmpdfdiff[[item$effect]] <- item$parameter
            if (!is.null(item$interaction1) && grepl('self\\$strat.+', item$interaction1)) {
              tmpdf$interaction1     <- paste(unique(item$interaction1), unique(item$x), collapse = '_')
              tmpdfdiff$interaction1 <- paste(unique(item$interaction1), unique(item$x), collapse = '_')
            }
          }
        }
        
        
        
        
        #
        #
        act <- self$chain_stats %>% 
          mutate(
            actor_id = as.factor(id_from),
            period =  round(ceiling(chain_step_id/ pd_steps)),
            strategy = as.factor(actor_strats[id_from])
          ) %>%
          group_by(actor_id, period, strategy)  %>% 
          dplyr::summarize(
            n_nochange=sum(stability), ## sum logical = count TRUE
            n_decision=n(), 
            forbearance = sum(stability)/n(),
            n_comps_explored =length(unique(id_to)),
            strategy = as.factor( paste(unique(actor_strats[ actor_id ]), collapse="|")),
            # explore = length(unique(id_to)) / n(), ## more components per decision = explore
            exploitation_rate = n()/length(unique(id_to)),  ## more decisions among fewer number of components = exploit
            entropy= -sum( (plyr::count(id_to)/n()) * log2(plyr::count(id_to)/n()) ), 
            tie_seq=paste(id_to, collapse="|"),
            step_seq = paste(chain_step_id, collapse = '|'),
            step_precedence_sum = sum(1 - chain_step_id/nrow(self$chain_stats), na.rm=T),
            step_precedence_mean=mean(1 - chain_step_id/nrow(self$chain_stats), na.rm=T),
            step_precedence_sd=sd(1 - chain_step_id/nrow(self$chain_stats), na.rm=T),
          ) %>% 
          left_join(
            self$actor_util_df %>% 
              mutate(
                actor_id=as.factor(actor_id), 
                period = round(ceiling(chain_step_id/ pd_steps))#,  strategy = as.factor(actor_strats[actor_id])
              ) %>% 
              group_by(actor_id, period, strategy ) %>%
              dplyr::summarize(
                utility_mean=mean(utility, na.rm=T)
              ), 
            by=c('actor_id', 'period','strategy')
          ) %>% 
          mutate(
            seed = seeds[i]
          )
        
        # act
        
        ## levels order
        coveffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
        covparams <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$parameter)
        covfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)ifelse(x$fix,'','(var)'))
        #
        structeffs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$effect)
        structparams <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)x$parameter)
        structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)ifelse(x$fix,'','(var)'))
        # structfixs   <- sapply(self$config_structure_model$dv_bipartite$effects, function(x)substr(as.character(x$fix),1,1))
        covDvTypes <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$interaction1)
        componentDV_ids <- grep('self\\$component_\\d{1,2}_coCovar', covDvTypes) ## ex: "self$component_1_coCovar"
        stratDV_ids     <- grep('self\\$strat_\\d{1,2}_coCovar', covDvTypes ) ## ex: "self$strat_1_coCovar"
        #
        compoeffs   <- coveffs[ componentDV_ids ]
        compoparams <- covparams[ componentDV_ids ]
        compofixs   <- covfixs[ componentDV_ids ]
        #
        strateffs   <- coveffs[ stratDV_ids ]
        stratparams <- covparams[ stratDV_ids ]
        stratfixs   <- covfixs[ stratDV_ids ]
        efflvls <- c('utility', strateffs, compoeffs, structeffs)
        # facet_strip_colors <- c('black',
        #                         rep('darkgray', length(strateffs)), 
        #                         rep('lightgray', length(compoeffs)),
        #                         rep('white', length(structeffs)))
        # facet_strip_text_colors <- c('white',rep('black', length(efflvls)-1))
        
        actor_stats_df <- self$actor_stats_df
        # actor_stats_df$strategy <- as.factor( actor_strats[ actor_stats_df$actor_id ] )
        ## Add utility as extra 'effect' 
        act_effs <- actor_stats_df %>% bind_rows( 
          self$actor_util_df %>% mutate(
            effect_id=NA, 
            value=utility, 
            effect_name='utility', 
            utility=NULL
          )
        ) %>%
          mutate(
            seed = seeds[i]
          )
        act_effs$effect_name <- factor(act_effs$effect_name, levels=efflvls)
        
        
        
        
        ##---- TESTING ------
        utility <- tmpdf %>% 
          group_by(strategy, chain_step_id) %>% 
          summarize(mean=mean(utility))
        #
        (dt0 <- dunn.test(utility$mean, utility$strategy, method='bonferroni' ))
        
        
        u1  <- tmpdf %>% filter(chain_step_id < median(chain_step_id)) %>%
          group_by(strategy, chain_step_id) %>% 
          summarize(mean=mean(utility))
        (dt1 <- dunn.test(u1$mean, u1$strategy, method='bonferroni' ))
        #
        u2  <- tmpdf %>% filter(chain_step_id >= median(chain_step_id)) %>%
          group_by(strategy, chain_step_id) %>% 
          summarize(mean=mean(utility))
        (dt2 <- dunn.test(u2$mean, u2$strategy, method='bonferroni' ))
        
        
        ##--/end------
        
        
        
        
        datlist <- list(  # plot = self$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T, thin_factor=3, loess_span=.5), #+  geom_vline(xintercept=seq(0,440,8*9))
          sim = self,
          seed = seeds[i],
          M = M, 
          N = N, 
          BI_PROB = BI_PROB,
          utdf = tmpdf,
          utdfdiff = tmpdfdiff,
          act = act,
          act_effs = act_effs,
          proj_ts = .ts,
          rep_uuid = rep_uuid,
          dunn.test = dt0,
          dunn.test.h1 = dt1,
          dunn.test.h2 = dt2
        )
        repfile <- sprintf('%s__%s.rds', .ts, rep_uuid)
        saveRDS(datlist,file = file.path(dir_proj_results, repfile))
        
        
        
        
        
      } ##/end i loop -------------------------------
      
      
      
      
    }
  }
}









## LOAD and re-list 
sensi_files <- dir(dir_proj_results, patter=sprintf('%s.+\\.rds',.ts))
#
lplot <- list()
lsim <- list()
lseed <- list()
lutdf <- list()
lutdfdiff <- list()
lact <- list()
lact_effs <- list()
lproj_ts <- list()
lrep_uuid <- list()
dunn.tests <- list()
dunn.tests.h1 <- list()
dunn.tests.h2 <- list()
for (i in 1:length(sensi_files)) {
  cat(sprintf('\n %s %.1f%s', i, 100*i/length(sensi_files),'%' ) )
  run <- readRDS(file.path(dir_proj_results, sensi_files[ i ] ))
  #
  lsim[[i]] <- run$sim 
  lseed[[i]] <- run$seed 
  lutdf[[i]] <- run$utdf %>% mutate(M=run$M, N=run$N, BI_PROB=run$BI_PROB)
  lutdfdiff[[i]] <- run$utdfdiff   %>% mutate(M=run$M, N=run$N, BI_PROB=run$BI_PROB)
  lact[[i]] <- run$act  %>% mutate(M=run$M, N=run$N, BI_PROB=run$BI_PROB)
  lact_effs[[i]] <- run$act_effs   %>% mutate(M=run$M, N=run$N, BI_PROB=run$BI_PROB)
  lproj_ts[[i]] <- run$proj_ts 
  lrep_uuid[[i]] <- run$rep_uuid 
  #
  dunn.tests[[i]]    <- c(run$dunn.test,     c(M=run$M, N=run$N, BI_PROB=run$BI_PROB))
  dunn.tests.h1[[i]] <- c(run$dunn.test.h1 , c(M=run$M, N=run$N, BI_PROB=run$BI_PROB))
  dunn.tests.h2[[i]] <- c(run$dunn.test.h2 , c(M=run$M, N=run$N, BI_PROB=run$BI_PROB))
}


## PLOTTING

# lplot[[1:2]]

combined_plot <- ggarrange(
  lsim[[1]]$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T, thin_factor=1, loess_span=.5), 
  lsim[[1]]$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T, thin_factor=10, loess_span=.5), 
  lsim[[2]]$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T, thin_factor=1, loess_span=.5), 
  lsim[[2]]$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T, thin_factor=10, loess_span=.5), 
  lsim[[3]]$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T, thin_factor=1, loess_span=.5), 
  lsim[[3]]$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T, thin_factor=10, loess_span=.5), 
  lsim[[4]]$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T, thin_factor=1, loess_span=.5), 
  lsim[[4]]$search_rsiena_multiwave_plot_actor_utility_strategy_summary(return_plot = T, thin_factor=10, loess_span=.5), 
  # plotlist = lplot[1:4], # lplot[[1]], lplot[[2]],  lplot[[3]],  lplot[[4]],
  ncol = 2,  nrow = 4,
  # widths = c(4.1,0.9), # Adjust column widths
  common.legend = TRUE, # Share a common legend if needed
  legend = "bottom"#,     # Place legend at the bottom
  # align = 'h'
) 
# print( combined_plot )
ggsave(sprintf('%s_sensitivity-gridsearch_utility_summary_compare4seeds_color_strategy.png', .ts ), 
       combined_plot,  device = 'png',
       height = 14, width = 14, units = 'in', dpi=800)


#
utdf <- data.table::rbindlist(lutdf, idcol = 'run_seed_str')
utdfdiff <- data.table::rbindlist(lutdfdiff, idcol = 'run_seed_str')
actdf <- data.table::rbindlist(lact, idcol = 'run_seed_str')
act_effsdf <- data.table::rbindlist(lact_effs, idcol = 'run_seed_str')


pdcnts <- plyr::count(actdf$period)
pdskeep <- pdcnts$x[pdcnts$freq > .9*pdcnts$freq[1]]
pdskeep <- pdskeep[ c(1:min(length(pdskeep)-1, 25), length(pdskeep)) ]


##====================== UTILITY TRANSITION PROFILE =============================

for (facet_scale in c('free','free_x','free_y'))
{
  
    
  for (ii in 1:length(Ms)) {
    ii_M <- Ms[ ii ]
    
    
    tran_profile_plots <- list()
    
    for (jj in 1:length(BI_PROBs)) {
      jj_BI_PROB <- BI_PROBs[ jj ]
      
      # DENSITY COMPARISONS BY PERIOD
      utdf2 <- utdf %>% filter(period %in% pdskeep, BI_PROB == jj_BI_PROB,  M == ii_M )
      plt.summary <- utdf2 %>% 
        ggplot(aes(x=utility, color=strategy, fill=strategy)) + 
        geom_histogram(alpha=.1, position='identity') +
        # geom_density(alpha=.1) +
        # scale_x_log10() +
        geom_vline(aes(xintercept = mean, color=strategy), 
                   data=utdf2%>%group_by(strategy,period,N,seed)%>%summarize(mean=mean(utility)),
                   linewidth=.8, linetype=3) + 
        geom_vline(aes(xintercept = mean, color=strategy), 
                   data=utdf2%>%group_by(strategy,period,N)%>%summarize(mean=mean(utility)),
                   linewidth=2) + 
        geom_vline(xintercept=0, linetype=2) +
        facet_grid(period ~ N, scales = facet_scale ) +
        theme_bw() + 
        ggtitle(sprintf('Utility\nBI_PROB=%.2f, M=%d (facet scales %s )', utdf2$BI_PROB[1], utdf2$M[1], facet_scale))
      # print(plt.summary)
      ggsave(sprintf('%s_sensitivity-gridsearch_utility_mean_facet_period_color_strategy_BI_PROB%.2f_M%s_%s.png',
                     .ts , utdf2$BI_PROB[1], utdf2$M[1], facet_scale ), 
             plt.summary,  
             height = 10, width = 8, units = 'in', dpi=600)
      
      tran_profile_plots[[length(tran_profile_plots)+1]] <- plt.summary
    }
    
    
    combined_tran_profile_plot <- ggarrange(
      plotlist = tran_profile_plots, # lplot[[1]], lplot[[2]],  lplot[[3]],  lplot[[4]],
      ncol = length(BI_PROBs),  nrow = 1,
      # widths = c(4.1,0.9), # Adjust column widths
      common.legend = TRUE, # Share a common legend if needed
      legend = "bottom"#,     # Place legend at the bottom
      # align = 'h'
    ) 
    # print( combined_tran_profile_plot )
    ggsave(sprintf('%s_sensitivity-gridsearch_utility_summary_color_strategy_COMBINEDFACETGRID_M%s_%s.png', .ts, ii_M, facet_scale ), 
           combined_tran_profile_plot,  device = 'png',
           height = 12, width = 18, units = 'in', dpi=800)
    
  }

  
}



# utdiff2 <- utdfdiff %>%  filter(period %in% pdskeep[1:10])
# plt.diff <- utdiff2 %>%
#     ggplot(aes(x=utility, color=strategy, fill=strategy)) +
#     geom_histogram(alpha=.1, position='identity') +
#     # geom_density(alpha=.1) +
#     geom_vline(aes(xintercept = mean, color=strategy),
#                data=utdiff2 %>%group_by(strategy,period,seed)%>%summarize(mean=mean(utility)),
#                linewidth=.8, linetype=3) +
#     geom_vline(aes(xintercept = mean, color=strategy),
#                data=utdiff2 %>%group_by(strategy,period)%>%summarize(mean=mean(utility)),
#                linewidth=2) +
#     geom_vline(xintercept=0, linetype=2) +
#     facet_grid(period ~ ., scales='free') +
#     theme_bw() + ggtitle('Utility Difference')
# ggsave(sprintf('%s_utility_diff_facet_period_color_strategy.png', .ts ), 
#        plt.diff,  
#        height = 10, width = 6, units = 'in', dpi=600)
# print(plt.diff)



##====================== ENTROPY TRANSITION PROFILE =============================

for (ii in 1:length(Ms)) {
  ii_M <- Ms[ ii ]
  
  
  entro_profile_plots <- list()
  
  for (jj in 1:length(BI_PROBs)) {
    jj_BI_PROB <- BI_PROBs[ jj ]
    
        
    actdf2 <- actdf %>% filter(period %in% pdskeep, BI_PROB == jj_BI_PROB,  M == ii_M )
    plt.ent <- actdf2 %>%
        ggplot(aes(x=entropy, color=strategy, fill=strategy, linetype=strategy)) + 
        geom_histogram(alpha=.05, position='identity') +
        # geom_density(alpha=.1) +
        geom_vline(aes(xintercept = mean, color=strategy), 
                   data=actdf2 %>%group_by(strategy,period,seed,N)%>%summarize(mean=mean(entropy)),
                   linewidth=.8, linetype=3) + 
        geom_vline(aes(xintercept = mean, color=strategy, linetype=strategy), 
                   data=actdf2%>%group_by(strategy,period,N)%>%summarize(mean=mean(entropy)),
                   linewidth=2) + 
        geom_vline(xintercept=0, linetype=2) +
        facet_grid(period ~ N, scales='free') +
        theme_bw() + ggtitle('Entropy')
    ggsave(sprintf('%s_sensitivity-gridsearch_entropy_facet_period_color_strategy_BI_PROB%.2f_M%s.png',
                   .ts , utdf2$BI_PROB[1], utdf2$M[1] ), 
           plt.ent,  
           height = 10, width = 6, units = 'in', dpi=600)
    print(plt.ent)
    
    
    
    
    entro_profile_plots[[length(entro_profile_plots)+1]] <- plt.ent
    
  }
    
  
  combined_entro_profile_plot <- ggarrange(
    plotlist = entro_profile_plots, # lplot[[1]], lplot[[2]],  lplot[[3]],  lplot[[4]],
    ncol = length(BI_PROBs),  nrow = 1,
    # widths = c(4.1,0.9), # Adjust column widths
    common.legend = TRUE, # Share a common legend if needed
    legend = "bottom"#,     # Place legend at the bottom
    # align = 'h'
  ) 
  # print( combined_tran_profile_plot )
  ggsave(sprintf('%s_sensitivity-gridsearch_entropy_summary_color_strategy_COMBINEDFACETGRID_M%s.png', .ts, ii_M ), 
         combined_entro_profile_plot,  device = 'png',
         height = 12, width = 18, units = 'in', dpi=800)
  
}
  

# (plt.ent <- actdf2 %>%
#     ggplot(aes(x=entropy, color=strategy, fill=strategy, linetype=strategy)) + 
#     # geom_histogram(alpha=.05, position='identity') +
#     geom_density(alpha=.1) +
#     geom_vline(aes(xintercept = mean, color=strategy), 
#                data=actdf2 %>%group_by(strategy,period,seed)%>%summarize(mean=mean(entropy)),
#                linewidth=.8, linetype=3) + 
#     geom_vline(aes(xintercept = mean, color=strategy, linetype=strategy), 
#                data=actdf2%>%group_by(strategy)%>%summarize(mean=mean(entropy)),
#                linewidth=2) + 
#     geom_vline(xintercept=0, linetype=2) +
#     # facet_grid(seed ~ ., scales='free_y') +
#     theme_bw() + ggtitle('Entropy')
# )






##====================== Mean Densities PROFILE =============================


for (facet_scale in c('free','free_x','free_y'))
{
  
    
  for (ii in 1:length(Ms)) {
    ii_M <- Ms[ ii ]
    
    globalmean_plots <- list()
    periodmean_plots <- list()
    
    for (jj in 1:length(BI_PROBs)) {
      jj_BI_PROB <- BI_PROBs[ jj ]
      
      
      utdf2 <-  utdf %>% filter(period %in% pdskeep, BI_PROB == jj_BI_PROB,  M == ii_M )
      
      plt_nseeds <- length(unique(utdf2$seed))
      plt_npds <- length(unique(utdf2$period))
      
      ##
      plt.gm <- utdf2 %>%
        group_by(strategy,seed,N) %>% summarize(utility=mean(utility,na.rm=T)) %>%
        mutate(rank=factor(order(utility, decreasing = TRUE))) %>%
        ggplot(aes(x=utility, color=strategy, fill=strategy)) +
        geom_density(alpha=.2) + 
        geom_vline(aes(xintercept = mean, color=strategy), 
                   data=utdf2%>%group_by(strategy,N)%>%summarize(mean=mean(utility))) +
        facet_wrap(~ N  , scales='free') +
        theme_bw() + ggtitle(sprintf('%s-period Replication Means (n = %d):\nBI_PROB=%.2f, M=%d (facet scales %s )', 
                                     plt_npds, plt_nseeds,  jj_BI_PROB, ii_M, facet_scale))  
      print(plt.gm)
      ggsave(sprintf('%s_sensitivity-gridsearch_global_means_by_strategy_COMBINEDFACETGRID_M%s_%s.png', .ts, ii_M, facet_scale ), 
             plt.gm,  
             height = 6, width = 12, units = 'in', dpi=600)
      
      ##
      plt.pdm <- utdf2 %>%
        group_by(strategy,period,seed,N) %>% summarize(utility=mean(utility,na.rm=T)) %>%
        mutate(rank=factor(order(utility, decreasing = TRUE))) %>%
        ggplot(aes(x=utility, color=strategy, fill=strategy)) +
        geom_density(alpha=.2) + 
        geom_vline(aes(xintercept = mean, color=strategy), 
                   data=utdf2%>%group_by(strategy,N)%>%summarize(mean=mean(utility))) +
        facet_wrap(~ N  , scales='free') +
        theme_bw() + ggtitle(sprintf('Replication-Period Means (n = %d : %s seeds %s periods):\nBI_PROB=%.2f, M=%d (facet scales %s )', 
                                     plt_nseeds * plt_npds, plt_nseeds, plt_npds, jj_BI_PROB, ii_M, facet_scale)) 
      print(plt.pdm)
      ggsave(sprintf('%s_sensitivity-gridsearch_periodwise_means_by_strategy_COMBINEDFACETGRID_M%s_%s.png', .ts, ii_M, facet_scale ), 
             plt.gm,  
             height = 6, width = 12, units = 'in', dpi=600)
      
      
      ##
      globalmean_plots[[length(globalmean_plots)+1]] <- plt.gm
      periodmean_plots[[length(periodmean_plots)+1]] <- plt.pdm
        
    }
    
    ##
    combined_globalmean_plot <- ggarrange(
      plotlist = globalmean_plots, # lplot[[1]], lplot[[2]],  lplot[[3]],  lplot[[4]],
      ncol = 1, nrow = length(BI_PROBs),
      # widths = c(4.1,0.9), # Adjust column widths
      common.legend = TRUE, # Share a common legend if needed
      legend = "bottom"#,     # Place legend at the bottom
      # align = 'h'
    ) 
    # print( combined_tran_profile_plot )
    ggsave(sprintf('%s_sensitivity-gridsearch_global_means_summary_color_strategy_COMBINEDFACETGRID_M%s_%s.png', .ts, ii_M, facet_scale ), 
           combined_globalmean_plot,  device = 'png',
           height = 11, width = 12, units = 'in', dpi=800)
    
    
    
    
    ##
    combined_periodmean_plots <- ggarrange(
      plotlist = periodmean_plots, # lplot[[1]], lplot[[2]],  lplot[[3]],  lplot[[4]],
      ncol = 1, nrow = length(BI_PROBs),
      # widths = c(4.1,0.9), # Adjust column widths
      common.legend = TRUE, # Share a common legend if needed
      legend = "bottom"#,     # Place legend at the bottom
      # align = 'h'
    ) 
    # print( combined_tran_profile_plot )
    ggsave(sprintf('%s_sensitivity-gridsearch_periodwise_means_summary_color_strategy_COMBINEDFACETGRID_M%s_%s.png', .ts, ii_M, facet_scale ), 
           combined_periodmean_plots,  device = 'png',
           height = 11, width = 12, units = 'in', dpi=800)
    
    
  }
  

}


# ## Rank by Seeds
# ralllist <- list()
# for (i in 1:length(seeds)) {
#   z <- utdf %>% filter(seed == seeds[i], BI_PROB == jj_BI_PROB,  M == ii_M ) %>%
#     group_by(strategy) %>% 
#     summarize(mean=mean(utility), sem=sd(utility)/sqrt(n()))
#   z$rank <- factor(order(z$mean, decreasing = T))
#   ralllist[[i]] <- z
# }
# ralldf <- rbindlist(ralllist)
# 
# plt.rall <- ralldf %>% ggplot(aes(y=fct_rev(rank), color=strategy, fill=strategy)) + 
#   geom_bar(position='dodge' ) +
#   theme_bw() + 
#   labs(y='Performance (Utility) Ranking') +
#   ggtitle(sprintf('Count of Replication (seed) Utility Rank by Strategy \n(seeds = %s, periods = %s )', nreps, length(unique(utdf$period)))) 
# print(plt.rall)
# ggsave(sprintf('%s_sensitivity-gridsearch_seed_rank_count_by_strategy.png', .ts ),
#        plt.rall, 
#        height = 6, width = 6, units = 'in', dpi=600)
# 
# 
# ## Rank by Period
# rlist <- list()
# for (i in 1:length(unique(utdf$period))) {
#   z <- utdf %>% filter(period == i) %>% group_by(strategy,period) %>% 
#     summarize(mean=mean(utility), sem=sd(utility)/sqrt(n()))
#   z$rank <- factor(order(z$mean, decreasing = T))
#   rlist[[i]] <- z
# }
# rdf <- rbindlist(rlist)
# 
# plt.r <- rdf %>% ggplot(aes(y=fct_rev(rank), color=strategy, fill=strategy)) + 
#   geom_bar(position='dodge' ) +
#   theme_bw() +
#   ggtitle(sprintf('Count of Periodwise Utility Rank by Strategy \n(seeds = %s, periods = %s )', nreps, length(unique(utdf$period)))) 
# print(plt.r)
# ggsave(sprintf('%s_sensitivity-gridsearch_period_rank_count_by_strategy.png', .ts ),
#        plt.r, 
#        height = 6, width = 6, units = 'in', dpi=600)









##====================== RANK ORDERING COUNTS  =============================

ralllist2 <- list()

for (ii in 1:length(Ms)) {
  ii_M <- Ms[ ii ]
  
  rank_h1h2_plots <- list()
  
  for (jj in 1:length(BI_PROBs)) {
    jj_BI_PROB <- BI_PROBs[ jj ]
    
    
    # utdf2 <-  utdf %>% filter(period %in% pdskeep, BI_PROB == jj_BI_PROB,  M == ii_M )
  
    
    ##** SPLIT HALF**
    ## Rank by Seeds
    for (i in 1:length(seeds)) {
      seed_i <- seeds[ i ]
      
      zlist <- list()
      
      for (j in 1:length(Ns)) {
        
        N_j <- Ns[ j ]
        
        z <- utdf %>% filter(period %in% pdskeep, BI_PROB == jj_BI_PROB,  M == ii_M, N == N_j, seed == seed_i ) %>%
          mutate(half=factor(ifelse(period < median(period), 'H1','H2'))) %>% 
          group_by(half, strategy) %>% 
          summarize(
            mean=mean(utility), 
            sem=sd(utility)/sqrt(n())
          ) %>% 
          mutate(
            rank = factor(order(mean, decreasing = T)),
            N = N_j,
            M = ii_M,
            BI_PROB = jj_BI_PROB
          ) 
          
        # z <- z %>% right_join(
        #     expand.grid(rank=factor(1:length(unique(utdf$strategy))), 
        #                 strategy=unique(utdf$strategy),
        #                 half=c('H1','H2')) 
        #   ) 
        
        
        
        plt_nseeds <- length(unique(z$seed))
        plt_npds <- length(unique(z$period))

      
        zlist[[ length(zlist)+1 ]] <-  z
      }
      
      
      ralllist2[[i]] <- rbindlist(zlist)
      # ralllist2[[i]] <- z1 %>% bind_rows(z2)
      
    }
    ralldf2 <- rbindlist(ralllist2) # %>% filter(!is.na(N))
    
    
    plt.rnk.split <-  expand.grid(
        rank=factor(1:length(unique(utdf$strategy))), 
        strategy=factor(unique(utdf$strategy)),
        half=factor(c('H1','H2')),
        N=unique(utdf$N)
      ) %>% 
      left_join( ralldf2 ) %>%
      group_by(strategy, half, rank, N) %>%
      summarize(count=sum(!is.na(mean))) %>%
      # mutate(N = as.factor(levels(N)[!is.na(levels(N))])) %>%
      ggplot(aes(x=count, y=fct_rev(rank), color=strategy, fill=strategy)) + 
      geom_bar(position='dodge', stat = 'identity', width=.7, alpha=.3) +
      facet_grid(half ~ N) +
      labs(y='Rank (Utility)') + 
      theme_bw() +
      ggtitle(sprintf('Count of Replication Utility Rank by Strategy \n(seeds = %s, periods = %s ) BI_PROB = %.2f, M = %d', 
                      nreps, length(unique(utdf$period)), jj_BI_PROB, ii_M)) 
    print(plt.rnk.split)
    ggsave(sprintf('%s_sensitivity-gridsearch_splithalf_seed_rank_count_by_strategy.png', .ts ),
           plt.rnk.split, 
           height = 6, width = 6, units = 'in', dpi=600)
    


    ##
    rank_h1h2_plots[[length(rank_h1h2_plots)+1]] <- plt.rnk.split
    
  }
  
  ##
  combined_rank_h1h2_plots <- ggarrange(
    plotlist = rank_h1h2_plots, # lplot[[1]], lplot[[2]],  lplot[[3]],  lplot[[4]],
    ncol = length(BI_PROBs),  nrow = 1,
    # widths = c(4.1,0.9), # Adjust column widths
    common.legend = TRUE, # Share a common legend if needed
    legend = "bottom"#,     # Place legend at the bottom
    # align = 'h'
  ) 
  # print( combined_tran_profile_plot )
  ggsave(sprintf('%s_sensitivity-gridsearch_splithalf_seed_rank_count_color_strategy_COMBINEDFACETGRID_M%s.png', .ts, ii_M ), 
         combined_rank_h1h2_plots,  device = 'png',
         height = 12, width = 18, units = 'in', dpi=800)
  
}
  
    
    

sink()
##----------------------------
##   TESTING
##----------------------------
sink(file = file.path(dir_proj, sprintf('%s_StrategyPairwiseWilcoxRankTestout.txt', .ts)), type = 'output')


cat('\n\n\n\nStrategy Pairwise Wilcox Rank Sum Test:\n')


cat('\n\nAll Periods (both H1,H2):\n\n')
# pwtest21 <- pairwise.wilcox.test( r2$rank, r2$strategy  )
# print(pwtest21)
pwtesth1h2 <- pairwise.wilcox.test( as.numeric(ralldf2$rank), ralldf2$strategy, p.adjust.method = 'bonf', paired = T  )
print(pwtesth1h2)


cat('\nHalf 1:\n\n')
r1 <- ralldf2 %>% filter(half=='H1')
# pwtest11 <- pairwise.wilcox.test( r1$rank, r1$strategy  )
# print(pwtest11)
pwtest12 <- pairwise.wilcox.test( as.numeric(r1$rank), r1$strategy, p.adjust.method = 'bonf', paired = T  )
print(pwtest12)


cat('\n\nHalf 2:\n\n')
r2 <- ralldf2 %>% filter(half=='H2')
# pwtest21 <- pairwise.wilcox.test( r2$rank, r2$strategy  )
# print(pwtest21)
pwtest22 <- pairwise.wilcox.test( as.numeric(r2$rank), r2$strategy, p.adjust.method = 'bonf', paired = T  )
print(pwtest22)


# pvaldf <- pwtest12$p.value %>% as.data.frame() %>% mutate(strategy1=as.character(rownames(pwtest12$p.value))) %>% 
#   tidyr::pivot_longer(cols=!strategy1, names_to = 'strategy2') %>% filter(!is.na(value))
# print(pvaldf)


##----------- subsample by conditions --------------------

cat('\n\n\n\n SUBSAMPLE BY CONDITIONS BI_PROB, M -- Strategy Pairwise Wilcox Rank Sum Test:\n')






for (ii in 1:length(Ms)) {
  ii_M <- Ms[ ii ]
  cat(sprintf('\n\n M = %s \n\n', ii_M))
  
  # rank_h1h2_plots <- list()
  
  for (jj in 1:length(BI_PROBs)) {
    jj_BI_PROB <- BI_PROBs[ jj ]
    cat(sprintf('\n   BI_PROB = %s \n', jj_BI_PROB))
    
    
    # utdf2 <-  utdf %>% filter(period %in% pdskeep, BI_PROB == jj_BI_PROB,  M == ii_M )
    
    rsub <- ralldf2 %>% filter(M == ii_M, BI_PROB == jj_BI_PROB)
    
      
    cat('\n\nAll Periods (both H1,H2):\n\n')
    # pwtest21 <- pairwise.wilcox.test( r2$rank, r2$strategy  )
    # print(pwtest21)
    pwtesth1h2sub <- pairwise.wilcox.test( as.numeric(rsub$rank), rsub$strategy, p.adjust.method = 'bonf', paired = T  )
    print(pwtesth1h2sub)
    
    
    cat('\nHalf 1:\n\n')
    r1sub <- rsub %>% filter(half=='H1')
    pwtest12sub <- pairwise.wilcox.test( as.numeric(r1sub$rank), r1sub$strategy, p.adjust.method = 'bonf', paired = T  )
    print(pwtest12sub)
    
    
    cat('\n\nHalf 2:\n\n')
    r2sub <- rsub %>% filter(half=='H2')
    pwtest22sub <- pairwise.wilcox.test( as.numeric(r2sub$rank), r2sub$strategy, p.adjust.method = 'bonf', paired = T  )
    print(pwtest22sub)
    
    
    # pvaldf <- pwtest12$p.value %>% as.data.frame() %>% mutate(strategy1=as.character(rownames(pwtest12$p.value))) %>% 
    #   tidyr::pivot_longer(cols=!strategy1, names_to = 'strategy2') %>% filter(!is.na(value))
    # print(pvaldf)


  }
  
}






cat('\n\n\n\nStrategy Pairwise Dunn Test over GRIDSEARCH:\n')

# cat('\nAll Period (both H1,H2):\n\n')
# # Perform Dunn's test  
# dunn_resultsh1h2 <- dunn.test( as.numeric(ralldf2$rank), ralldf2$strategy,  method = "bonferroni")  
# print(dunn_resultsh1h2)
# 
# cat('\nHalf 1:\n\n')
# # Perform Dunn's test  
# dunn_results1 <- dunn.test( as.numeric(r1$rank), r1$strategy,  method = "bonferroni")  
# # Print the results  
# print(dunn_results1) 
# 
# cat('\n\nHalf 2:\n\n')
# # Perform Dunn's test  
# dunn_results2 <- dunn.test( as.numeric(r2$rank), r2$strategy,  method = "bonferroni")  
# # Print the results  
# print(dunn_results2) 



dtlong    <- ldply(dunn.tests, as.data.frame) ## list of lists; apply for each test list in replication list 
dtlong.h1 <- ldply(dunn.tests.h1, as.data.frame) ## list of lists; apply for each test list in replication list 
dtlong.h2 <- ldply(dunn.tests.h2, as.data.frame) ## list of lists; apply for each test list in replication list 




## Zscores from each replications strategy comparison (across utility vectors)
zsigtab <- dtlong %>% group_by(comparisons, M, N, BI_PROB) %>% 
  summarize(sig_abs_z2 = sum(abs(Z) > 2)/ n())
#
grp.ids <- which(names(zsigtab) %in% c('comparisons','M','N'))
zcnts <- apply(zsigtab[, -grp.ids], 1, function(x)paste(c('Prop(|Z|>2) = '), round(x,3), sep=''))
zcnts <- rbind(zsigtab$comparisons, zcnts)
#
# cat(
#   paste(apply(zcnts, 2, function(x)paste(x[1],':', paste(x[-1], collapse = ';  '))), collapse = '\n' , '\n\n')
# )
#






##====================== DUNN TEST (PAIRWISE Wilcoxon RANK TEST) of RANK ORDERING COUNTS  =============================




for (facet_scale in c('free','free_x','free_y'))
{
  

  dzplots <- list()
  
  for (jj in 1:length(BI_PROBs)) {
    jj_BI_PROB <- BI_PROBs[ jj ]
    cat(sprintf('\n   BI_PROB = %s \n', jj_BI_PROB))
  
  
  
      dzplot0 <- dtlong %>%  filter(BI_PROB == jj_BI_PROB) %>%
        ggplot(aes(x=abs(Z), color=comparisons, fill=comparisons, linetype=comparisons)) +
        geom_histogram(aes(y=..density..), position='identity', alpha=.2, linewidth=1) +
        geom_density(alpha=.2, linewidth=1) +
        geom_vline(aes(xintercept = mean, color=comparisons, linetype = comparisons),
                   data=dtlong%>%group_by(comparisons, N, M)%>%summarize(mean=mean(abs(Z)))) +
        geom_vline(xintercept = c(2)) +
        # color_palette(palette = 'Dark2') +
        scale_color_manual(values = RColorBrewer::brewer.pal( length(unique(utdf$strategy)), 'Dark2' ) ) +
        facet_grid(N ~ M, scales = facet_scale) +
        labs(x='Z-Score Abs. Val. (Bonferroni Adj.)') +
        theme_bw()  +  ggtitle(sprintf('Pairwise Dunn Test of Avg. Utility by Strategy\nBI_PROB = %s (facet scales %s )',
                                       jj_BI_PROB,  facet_scale))
        # ggtitle(sprintf('Pairwise Dunn Test of Avg. Utility by Strategy from %s replications\n%s',  nreps,
        #                 paste(apply(zcnts, 2, function(x)paste(x[1],':', paste(x[-1], collapse = ';  '))), collapse = '\n' )))
      print(dzplot0)
      ggsave(sprintf('%s_sensitivity-gridsearch_start0_dunn_test_abs_Z_scores_by_strategy_BI_PROB%s__%s.png',
                     .ts, jj_BI_PROB,  facet_scale ),
             dzplot0,
             height = 6, width = 9, units = 'in', dpi=600)
  
  
      # dzplot1 <- dtlong %>%  filter(BI_PROB < .1) %>%
      #   ggplot(aes(x=abs(Z), color=comparisons, fill=comparisons, linetype=comparisons)) +
      #   geom_histogram(aes(y=..density..), position='identity', alpha=.3, linewidth=1) +
      #   geom_density(alpha=.2, linewidth=1) +
      #   geom_vline(aes(xintercept = mean, color=comparisons, linetype = comparisons),
      #              data=dtlong%>%group_by(comparisons, M, N)%>%summarize(mean=mean(abs(Z)))) +
      #   geom_vline(xintercept = c(2)) +
      #   # color_palette(palette = 'Dark2') +
      #   scale_color_manual(values = RColorBrewer::brewer.pal( length(unique(utdf$strategy)), 'Dark2' ) ) +
      #   facet_grid(M ~ N, facet_scale) +
      #   labs(x='Z-Score Abs. Val. (Bonferroni Adj.)') +
      #   theme_bw()  +  ggtitle(sprintf('Pairwise Dunn Test of Avg. Utility by Strategy BI_PROB0.01_%s', facet_scale))
      # # ggtitle(sprintf('Pairwise Dunn Test of Avg. Utility by Strategy from %s replications\n%s',  nreps,
      # #                 paste(apply(zcnts, 2, function(x)paste(x[1],':', paste(x[-1], collapse = ';  '))), collapse = '\n' )))
      # print(dzplot1)
      # ggsave(sprintf('%s_sensitivity-gridsearch_start0_dunn_test_abs_Z_scores_by_strategy_BI_PROB0.99_%s.png', .ts, facet_scale ),
      #        dzplot1,
      #        height = 6, width = 9, units = 'in', dpi=600)
  
  
      dzplots[[ length(dzplots) + 1 ]] <- dzplot0
  
  }
  
  ##
  combined_dzplots <- ggarrange(
    plotlist = dzplots, # lplot[[1]], lplot[[2]],  lplot[[3]],  lplot[[4]],
    ncol = 1, nrow = length(BI_PROBs),  
    # widths = c(4.1,0.9), # Adjust column widths
    common.legend = TRUE, # Share a common legend if needed
    legend = "bottom"#,     # Place legend at the bottom
    # align = 'h'
  ) 
  # print( combined_tran_profile_plot )
  ggsave(sprintf('%s_sensitivity-gridsearch_start0_dunn_test_abs_Z_scores_by_strategy_COMBINEDFACETGRID_%s.png', 
                 .ts,  facet_scale ), 
         combined_dzplots,  device = 'png',
         height = 7*length(BI_PROBs), width = 8, units = 'in', dpi=800)



}

# ## P-vals from each replications strategy comparison (across utility vectors)
# sigtab <- dtlong %>% group_by(comparisons) %>% 
#   summarize(sig_05 = sum(P.adjusted < .05)/ n(),
#             sig_01 = sum(P.adjusted < .01) / n(),
#             sig_001 = sum(P.adjusted < .001)/ n())
# 
# testcnts <- apply(sigtab[,-1], 1, function(x)paste(c('Prop(<.05) = ','Prop(<.01) = ','Prop(<.001) = '), round(x,3), sep=''))
# testcnts <- rbind(sigtab$comparisons,testcnts)
# 
# cat(
#   paste(apply(testcnts, 2, function(x)paste(x[1],':', paste(x[-1], collapse = ';  '))), collapse = '\n' )
# )
# 
# nbins <- 31
# hist_data <- hist(dtlong$P.adjusted, breaks = nbins, plot = FALSE)  
# hist_data$density <- hist_data$counts / sum(hist_data$counts)  # Normalize counts  
# 
# 
# dtplot <- dtlong %>% ggplot(aes(x=P.adjusted, color=comparisons, fill=comparisons, linetype=comparisons))+ 
#   geom_histogram(aes(y=..density..), position='identity', bins = nbins , alpha=.1, linewidth=1) +
#   geom_density(alpha=.2, linewidth=2) +
#   geom_vline(xintercept = c(.05, .01, .001), color=c('black'), linetype=c(3,2,1), linewidth=1) +
#   scale_color_manual(values = RColorBrewer::brewer.pal( length(unique(utdf$strategy)), 'Dark2' ) ) +
#   scale_x_log10() +
#   theme_bw() + 
#   ggtitle(sprintf('Dunn Test of Utility by Strategy (Bonferroni Adj. P-Values) from %s replications:\n%s', 
#                    nreps, 
#                   paste(apply(testcnts, 2, function(x)paste(x[1],':', paste(x[-1], collapse = ';  '))), collapse = '\n' )
#                    ))
# print(dtplot)
# ggsave(sprintf('%s_dunn_test_pvals_by_strategy.png', .ts ),
#        dtplot, 
#        height = 6, width = 9, units = 'in', dpi=600)





cat('\n\nEND.')
####
## END SINK
####
sink()


















