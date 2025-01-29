rm(list=ls())  ## uncomment to clear environ before run
########################
##
##   SAOM-NK Runscript:
##
##   1. Test 1
##
##
#######################

## Directories
dir_proj <- 'C://Users//sdr8y//OneDrive - University of Missouri//Research//Search_networks//SaoMNK//R'
dir_proj_base <- 'C://Users//sdr8y//OneDrive - University of Missouri//Research//Search_networks'
dir_data <- 'D://Search_networks'

## default settings: Users do not change; TODO: implment within restricted class attributes
DV_NAME <- 'self$bipartite_rsienaDV'




# readxl::read_excel(file.path(dir_proj_base, 'SAOMNK_gridsearch_params_base_v1.xlsx'))
# 
# 
# read.csv(file.path(dir_proj, 'SAOM_NK_gridsearch_params_base_v1.csv'), stringsAsFactors = F)



self <- readRDS('D:\\Search_networks\\__SAOMNK_gridsearch_input_basic_M3__173663250316762\\grid_results\\6_SAOMNK_gridsearch_input_basic_M3_173663250316762_dc6c23c9-6d2d-4d7c-9506-c2b7ec924378.rds')




library(ggridges) 
library(dplyr)
library(forcats)

# ## https://wilkelab.org/ggridges/articles/gallery.html
# ##  EXAMPLE OF density ridges like a movie of evolving density
# # text labels
# ggplot(lincoln_weather, aes(x = `Mean Wind Speed[MPH]`, y = `Month`, fill = ..x..)) +
#   geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
#   scale_y_discrete(expand = c(0, 0, 0.2, 0)) + # add top margin
#   scale_x_continuous(expand = c(0, 0)) +
#   theme(legend.position = "none") +
#   scale_fill_viridis(option = "C", direction = -1) +
#   labs(title = "Daily Mean Wind Speeds in Lincoln, NE in 2016") +
#   xlab("Mean Wind Speed (MPH)")
# 
# ##  EXAMPLE OF multi-group (multicolor) density ridges like a movie of evolving group density




##****************************************************##
## INPUT FILE
grid.input.file <- 'SAOMNK_gridsearch_input_basic_M3.xlsx'
## ITERATIONS
ITERATIONS <- 10000
##****************************************************##

inputs <- read_excel(file.path(dir_proj_base, grid.input.file))

##
environ_seed_params <- unique(inputs$name[inputs$type %in% c('environment','seed')])
eff_cov_params      <- unique(inputs$name[inputs$type %in% c('effects','coCovars','varCovars','dyCovars')])
eff_struct_params   <- unique(inputs$name[inputs$type %in% c('effects')])
eff_component_params<- unique(inputs$name[inputs$type %in% c('coCovars','varCovars','dyCovars') & grepl('self\\$compon.+',inputs$interaction1) ])
eff_strategy_params <- unique(inputs$name[inputs$type %in% c('coCovars','varCovars','dyCovars') & grepl('self\\$strat.+',inputs$interaction1) ])




# get_actor_strategies = function() {
#   coveffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
#   covDvTypes <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$interaction1)
#   stratDV_ids     <- grep('self\\$strat_\\d{1,2}_coCovar', covDvTypes ) ## ex: "self$strat_1_coCovar" 
#   strat_effs <- coveffs[ stratDV_ids ]
#   nstrat <- length(strat_effs)
#   if (nstrat == 1) return(self$strat_1_coCovar)
#   if (nstrat == 2) return(paste(self$strat_1_coCovar, 
#                                 self$strat_2_coCovar, sep='_'))
#   if (nstrat == 3) return(paste(self$strat_1_coCovar, 
#                                 self$strat_2_coCovar, 
#                                 self$strat_3_coCovar, sep='_'))
#   if (nstrat == 4) return(paste(self$strat_1_coCovar, 
#                                 self$strat_2_coCovar, 
#                                 self$strat_3_coCovar, 
#                                 self$strat_4_coCovar, sep='_'))
#   if (nstrat == 5) return(paste(self$strat_1_coCovar, 
#                                 self$strat_2_coCovar, 
#                                 self$strat_3_coCovar, 
#                                 self$strat_4_coCovar, 
#                                 self$strat_5_coCovar, sep='_'))
#   # if (nstrat == 1) return(paste(strat_effs, self$strat_1_coCovar, sep='= '))
#   # if (nstrat == 2) return(paste(strat_effs, paste(self$strat_1_coCovar, 
#   #                                                 self$strat_2_coCovar, sep='_'), sep='= '))
#   # if (nstrat == 3) return(paste(strat_effs, paste(self$strat_1_coCovar, 
#   #                                                 self$strat_2_coCovar, 
#   #                                                 self$strat_3_coCovar, sep='_'), sep='= '))
#   # if (nstrat == 4) return(paste(strat_effs, paste(self$strat_1_coCovar, 
#   #                                                 self$strat_2_coCovar, 
#   #                                                 self$strat_3_coCovar, 
#   #                                                 self$strat_4_coCovar, sep='_'), sep='= '))
#   # if (nstrat == 5) return(paste(strat_effs, paste(self$strat_1_coCovar, 
#   #                                                 self$strat_2_coCovar, 
#   #                                                 self$strat_3_coCovar, 
#   #                                                 self$strat_4_coCovar, 
#   #                                                 self$strat_5_coCovar, sep='_'), sep='= '))
# }




actor_ids=c() 
wave_ids=c()
thin_factor=2
thin_wave_factor=1
smooth_method='loess'  ##"lm", "glm", "gam", "loess","auto"
show_utility_points=T
scale_utility=TRUE
return_plot=FALSE
plot_file=NA
plot_dir=NA


actor_strat <- as.factor( self$get_actor_strategies() )
#
nstep <- sum(!self$chain_stats$stability)
#
coveffs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$effect)
covparams <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)x$parameter)
covfixs   <- sapply(self$config_structure_model$dv_bipartite$coCovars, function(x)ifelse(x$fix,'','(var)'))
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
compoeffs   <- coveffs[ componentDV_ids ]
compoparams <- covparams[ componentDV_ids ]
compofixs   <- covfixs[ componentDV_ids ]
#
strateffs   <- coveffs[ stratDV_ids ]
stratparams <- covparams[ stratDV_ids ]
stratfixs   <- covfixs[ stratDV_ids ]
#
actor_component_period <- self$M * self$N
#
transition_pds <- 5
density_ridges_rel_min_height = 1e-07  ## prevents density ridges colored lines from covering full x-axis (clarifies group separation)
## Compare 2 actors utilty
dat <- self$actor_wave_util %>% 
  filter(chain_step_id %% thin_factor == 0) %>% 
  filter(wave_id %% thin_wave_factor == 0 ) %>% 
  mutate(
    strategy = actor_strat[ actor_id ],
    chain_below_med =  chain_step_id < median(chain_step_id),
    actor_component_period = 1 + floor( chain_step_id / actor_component_period )
  ) %>% 
  mutate(
    stabilization_summary_period = ifelse(actor_component_period <= transition_pds, 
                                          actor_component_period, 
                                          sprintf('%s+\n(%s-%s)',transition_pds+1,transition_pds+1,max(actor_component_period)))
  )
# dat$chain_half <- factor(ifelse(dat$chain_below_med, '1st Half', '2nd Half'))

#
util_lab <- 'Actor Utility'
if(scale_utility) {
  util_sc <- scale(dat$utility)
  util_lab <- sprintf('Actor Utility\n(Standardized Center = %.2f; Scale = %.2f)',
                      attr(util_sc, 'scaled:center'), 
                      attr(util_sc, 'scaled:scale'))
  if (!all(dat$utility == 0))
    dat <- dat %>% mutate(utility = c(scale(utility)))
}
#
density_rng <- range(dat$utility, na.rm=T)
density_absdiff_scale <- abs(diff(density_rng)) * 0.15
util_lim <- c(density_rng[1] - density_absdiff_scale,  density_rng[2] + density_absdiff_scale)
#
point_size <- 10 / log( nstep )
point_alpha <- min( 1,  1/log10( nstep ) )
#
if(length(actor_ids))
  dat <- dat %>% filter(actor_id %in% actor_ids)
if(length(wave_ids))
  dat <- dat %>% filter(wave_id %in% wave_ids)
dat_acp_stabil_means <- dat %>% group_by(stabilization_summary_period, strategy) %>% 
  dplyr::summarize(mean=mean(utility, na.rm=T)) %>%
  mutate(PeriodFct = fct_rev(as.factor(stabilization_summary_period)))
#
# plt <- ggplot(dat, aes(x=chain_step_id, y=utility)) + 
#   geom_hline(data=dat_acp_means, aes(yintercept=mean), linetype=3, col='black' ) +
#   facet_grid(wave_id ~ .) 


##==============================================
strat_legend_title <- sprintf("Strategy (%s) :  ", paste(strateffs, collapse = '_'))
strat_break <- levels(actor_strat) 
strat_labs <- sapply(1:length(levels(actor_strat)), function(i) {
  a <- levels(actor_strat)[i]
  names(a) <- a ## # names(a) <- sprintf('%s: %s', i, a)
  return(a)
}) 
dat_dens_rigde <- dat %>%
  mutate(PeriodFct = fct_rev(as.factor(stabilization_summary_period))) 
  # mutate(PeriodFct = as.factor(stabilization_summary_period))

group_dens_means <- dat_dens_rigde %>% ungroup() %>% group_by(strategy) %>% 
  dplyr::summarize(n=n(),mean=mean(utility,na.rm=T))

# dat_ridge_lines_density <-  dat_dens_rigde %>% ungroup() %>%
#   group_by(stabilization_summary_period, strategy) %>%
#   summarise(density = density(utility), 
#             strategy=paste(unique(strategy),collapse = "|"),
#             stabilization_summary_period=paste(unique(stabilization_summary_period),collapse = "|") ) #%>%
#   # unnest_wider(density) %>%
#   # rename(density_value = y, group = y) %>%
#   # unnest(cols = c(x, y)) %>%
#   # mutate(y = as.numeric(factor(group))) # Assign numeric y-values for alignment


plt.dr <- ggplot(dat_dens_rigde, aes(y = PeriodFct, x = utility, color=strategy, fill=strategy)) +
  # geom_density_ridges(stat = "binline", alpha=.15,
  #                     bins = 90, draw_baseline = FALSE) +
  ######
  stat_density_ridges(aes(point_color = strategy, point_fill = strategy, point_shape = strategy),
                      quantile_lines = TRUE, alpha = .3, rel_min_height = density_ridges_rel_min_height,
                      point_size=.4,
                      jittered_points = T, 
                      position = position_raincloud(adjust_vlines = FALSE, ygap = -.1, height = .15),# "raincloud",
                      # position = position_points_sina(rel_min = 0.1, rel_max = 0.9, seed = NULL),
                      quantiles = c(0.5), linewidth=.75 ) +
  # stat_density_ridges(aes(point_color = strategy, point_fill = strategy, point_shape = strategy),
  #                     quantile_lines = TRUE, alpha = 0.25,
  #                     point_size=.2,
  #                     quantiles = c(0.5), linewidth=.7 ) +
  ######
  ## geom_density_ridges(
  ##   aes(x = utility , fill = strategy),
  ##   # jittered_points = TRUE, scale = .95, rel_min_height = 0.2,
  ##   # point_shape = "|", point_size=1,
  ##   # position =  position_points_jitter(height = 0),
  ##   alpha = .4, #color = "darkgray",
  ##   from = min(dat$utility)-0.5, to = max(dat$utility)+0.5
  ## ) +
  # geom_line(data = dat_ridge_lines_density, aes(x = x, y = y + ..dat_ridge_lines_density..), 
  #           inherit.aes = FALSE, size = 0.8) +
  ###
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_cyclical(
    breaks = strat_break,
    labels = strat_labs,
    values = hue_pal()(length(levels(actor_strat))), # c("#ff0000", "#0000ff", "#ff8080", "#8080ff"),
    # name = sprintf("strategy\n(%s)", paste(strateffs, collapse = '_')), 
    guide = "legend"
  ) +
  labs(
    x = util_lab,
    y = sprintf(" Time Period \n(Actor-Component-Period = %s decision steps)", actor_component_period ),
    title = "Actor Utility Stabilization Paths",
    subtitle = sprintf("(Decision Chain Iterations: %s )", nstep), 
    # caption = "(plot design adapted from source: Marc Belzunces (@marcbeldata))",
    color = strat_legend_title,
    fill = strat_legend_title,
    point_color = strat_legend_title,
    point_fill = strat_legend_title,
    point_shape = strat_legend_title
  ) +
  geom_vline(xintercept = 0, linetype=1) +
  geom_vline(data = group_dens_means, 
             aes(xintercept = mean, color=strategy),
             linetype=3, linewidth=1.3) +
  geom_text(data = group_dens_means, 
            aes(x = mean, y = 1+length(unique(dat_dens_rigde$stabilization_summary_period)), label = round(mean, 2), color=strategy), 
            inherit.aes = FALSE, size = 5, nudge_x=-.1, nudge_y=.35  ) +
  coord_cartesian(clip = "off") +
  theme_ridges(grid = T, center=T) + theme(legend.position = 'bottom')
print(plt.dr) 

# #
# plt.dr.dat <- ggplot_build(plt.dr)$data[[1]]
# plt.dr.b$data
##==============================================

# plyr::count(plt.dr.dat$group)
# nstrats <- length(levels(actor_strat))
# npds <- length(levels(as.factor(dat$stabilization_summary_period)))
# 
# getGrpId <- function(strat, pd, nstrats, npds) (strat-1)*npds + pd #(npds-pd+1)
# getGrpId(strat=2, pd=1, nstrats=nstrats, npds=npds)
# 
# getPdFromPlotDf <- function(strat, pd, nstrats, npds)  
# getStratFromPlotDf <- function(strat, pd, nstrats, npds)  
# 
# plt.dr.dat %>% 
#   mutate(group_fct = as.factor(group)) %>%
#   filter(!is.na(density)) %>%   
#   group_by(group = group_fct) %>% 
#   dplyr::summarize(
#     n=n(), x_min= min(x, na.rm=T), x_max=max(x, na.rm=T), y_pos=first(y)
#   )
# plt.dr.dat$strategy <- NA
# plt.dr.dat$pd <- NA
# for (i in 1:nrow(plt.dr.dat)) {
#   plt.dr.dat$strategy[i] <- plt.dr.dat$group[i]
#   plt.dr.dat$pd[i]       <- 
# }


head(plt.dr.dat)


# 
# line_data <- dat_dens_rigde %>%
#   group_by(group, as.factor(y) ) %>%
#   summarise(
#     x_min = min(x[!is.na(density)]),  # Minimum x with density > 0
#     x_max = max(x[!is.na(density)]),  # Maximum x with density > 0
#     y_pos = first(y)                 # y-position for the group
#   )
# # Add connecting lines to the plot
# (plt.dr.lined <- plt.dr +
#   geom_segment(data = line_data,
#                aes(x = x_min, xend = x_max, y = y_pos, yend = y_pos),
#                inherit.aes = FALSE, size = 0.8)
# )




###






















##---------------------------------------------------------------------------------------------
##
search_rsiena_multiwave_plot_actor_utility_strategy_summary = function(actor_ids=c(), 
                                                                       wave_ids=c(),
                                                                       thin_factor=1, 
                                                                       thin_wave_factor=1,
                                                                       smooth_method='loess',  ##"lm", "glm", "gam", "loess","auto"
                                                                       show_utility_points=T,
                                                                       scale_utility=TRUE,
                                                                       return_plot=FALSE,
                                                                       plot_file=NA, plot_dir=NA
) {
  ## actor strategy
  if ( attr(self$strat_1_coCovar, 'nodeSet') != 'ACTORS' )
    stop("Actor Strategy self$strat_1_coCovar not set.")
  actor_strat <- as.factor( self$get_actor_strategies() )
  #
  nstep <- sum(!self$chain_stats$stability)
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
  actor_component_period <- self$M * self$N
  ## Compare 2 actors utilty
  dat <- self$actor_wave_util %>% 
    filter(chain_step_id %% thin_factor == 0) %>% 
    filter(wave_id %% thin_wave_factor == 0 ) %>% 
    mutate(
      strategy = actor_strat[ actor_id ],
      chain_below_med =  chain_step_id < median(chain_step_id),
      actor_component_period = 1 + floor( chain_step_id / actor_component_period )
    )
  dat$chain_half <- factor(ifelse(dat$chain_below_med, '1st Half', '2nd Half'))

  #
  util_lab <- 'Actor Utility'
  if(scale_utility) {
    util_sc <- scale(dat$utility)
    util_lab <- sprintf('Actor Utility\n(Standardized Center = %.2f; Scale = %.2f)',
                        attr(util_sc, 'scaled:center'), 
                        attr(util_sc, 'scaled:scale'))
    if (!all(dat$utility == 0))
      dat <- dat %>% mutate(utility = c(scale(utility)))
  }
  #
  density_rng <- range(dat$utility, na.rm=T)
  density_absdiff_scale <- abs(diff(density_rng)) * 0.15
  util_lim <- c(density_rng[1] - density_absdiff_scale,  density_rng[2] + density_absdiff_scale)
  #
  point_size <- 10 / log( nstep )
  point_alpha <- min( 1,  1/log10( nstep ) )
  #
  if(length(actor_ids))
    dat <- dat %>% filter(actor_id %in% actor_ids)
  if(length(wave_ids))
    dat <- dat %>% filter(wave_id %in% wave_ids)
  dat_wave_means <- dat %>% group_by(wave_id) %>% 
    dplyr::summarize(mean=mean(utility, na.rm=T))
  #
  plt <- ggplot(dat, aes(x=chain_step_id, y=utility)) + 
    geom_hline(data=dat_wave_means, aes(yintercept=mean), linetype=3, col='black' ) +
    facet_grid(wave_id ~ .) 
  if(show_utility_points)
    plt <- plt + geom_point(aes(color=strategy), alpha=point_alpha, shape=1, size=point_size)  # geom_line(alpha=.2) +#geom_smooth(method='loess', alpha=.1) + 
  if(self$exists(smooth_method))
    plt <- plt + geom_smooth(aes(linetype=actor_id, color=strategy, fill=strategy), method = smooth_method, linewidth=1, alpha=.09)
  #
  plt <- plt + geom_vline(xintercept = which(1:max(dat$chain_step_id) %% actor_component_period==0) , linetype = 3, color='gray')
  # 
  plt <- plt + theme_bw() + 
    ylim(util_lim) + 
    ylab(util_lab) +
    xlab('Actor Decision Chain Ministep') +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.key.width = unit(0.8, "cm"),  # Adjust legend key width
      legend.spacing.x = unit(0.3, "cm")   # Adjust spacing between keys
    )
  
  ACTORS     <- sienaNodeSet(self$M, nodeSetName="ACTORS")
  COMPONENTS <- sienaNodeSet(self$N, nodeSetName="COMPONENTS")
  # strat_mat_varCovar <- varCovar(strat_eff$x, nodeSet = 'ACTORS')
  if (! 'coCovars' %in% names(self$config_structure_model$dv_bipartite) ) {
    rsiena_data <- sienaDataCreate(list(self$bipartite_rsienaDV), nodeSets = list(ACTORS, COMPONENTS))
    return(rsiena_data)
  }
  
  plt <- plt + ggtitle(sprintf('Environment: Actors (M) = %s, Components (N) = %s, Init.Prob. = %.2f\nActor Strategy:  %s\nComponent Payoff:  %s\nStructure:  %s', 
                               self$M, self$N, self$BI_PROB,
                               paste( paste(paste(strateffs[stratDV_ids], stratparams[stratDV_ids], sep='= '), stratfixs[stratDV_ids], sep='' ), collapse = ';  '),
                               paste( paste(paste(strateffs[componentDV_ids], stratparams[componentDV_ids], sep='= '), stratfixs[componentDV_ids], sep=''), collapse = ';  '),
                               paste( paste(paste(structeffs, structparams, sep='= '), structfixs, sep=''), collapse = ';  ')
  ))
  plt <- plt +  guides(color = guide_legend(nrow = 1))
  
  #### Density
  stratmeans <- dat %>% group_by(strategy, wave_id) %>% 
    dplyr::summarize(n=n(), mean=mean(utility, na.rm=T), sd=sd(utility, na.rm=T))
  ## Actor density fact plots comparing H1 to H2 utility distribution
  plt2 <- ggplot(dat, aes(x=utility, color=strategy, fill=strategy)) + ##linetype=chain_half
    geom_density(alpha=.1, linewidth=1)  +
    # geom_histogram(alpha=.1, position = 'dodge') +
    geom_vline(data = stratmeans, aes(xintercept = mean, color=strategy), linetype=2, linewidth=.9) +
    geom_vline(data = dat_wave_means,  aes(xintercept=mean), linetype=3, col='black' ) +
    labs(y='', x='') +
    # xlim(c(ggplot_build(plt)$layout$panel_params[[1]]$y.range)) + 
    xlim(util_lim) +
    coord_flip() +
    facet_grid(wave_id ~ .) +
    ylab('Actor Utility Density') +
    # geom_vline(xintercept = 0, linetype=2, linewidth=.9, color='gray')+
    theme_bw() + theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = 'none', 
      plot.margin=unit(c(5.5, 5.5, 5.5, -23), 'pt'),
      axis.text.y = element_blank(),
      axis.ticks.y=element_blank()#,
      # axis.text.x = element_blank(),
      # axis.ticks.x=element_blank()#,
    ) + ggtitle('\n\n\n')
  
  # combined_plot <- plot_grid(plt ,#+ theme(panel.spacing = unit(0, "lines")), 
  #                            plt2 ,# + theme(panel.spacing = unit(0, "lines")), 
  #                            ncol = 2, rel_widths = c(4, 1)) 
  # print(combined_plot)
  
  combined_plot <- ggarrange(
    plt, plt2, 
    ncol = 2, 
    widths = c(4.1,0.9), # Adjust column widths
    common.legend = TRUE, # Share a common legend if needed
    legend = "bottom"#,     # Place legend at the bottom
    # align = 'h'
  ) 
  # print(combined_plot)
  #
  # print(plt)
  #
  if(!is.na(plot_file))
    ggsave(file = file.path(ifelse(is.na(plot_dir),getwd(),plot_dir), sprintf("%s_%s.png", self$config_environ_params$name, plot_file)), 
           combined_plot, 
           width = 10, height = 8, units = 'in', dpi = 600)
  
  if(return_plot)
    return(combined_plot)
}
##---------------------------------------------------------------------------------------------











































# outlist <- list()
# Ns <- 2 * 3^c(1:5)  ## 6 12 24 48 96
# for (i in 1:length(Ns)) 
# {
  
  
  ##****************************************##
  ## I. SIM SETUP 
  ##****************************************##
  ##----------------------
  ## 1. Environment:  determines what types of entities and strategies are permissible
  ##----------------------
  environ_params <- list(
    M = 3,        ## Actors
    N = 12,       ## Components
    BI_PROB = 0, ## Environmental Density (DGP hyperparameter)
    component_matrix_start = 'rand', ##**TODO** Implement: 'rand','modular','semi-modular',...
    rand_seed = 123,
    visualize_init = F,
    name = sprintf('_testSAOMNK_4_loopN%s', 12 )#,
    ## use starting matrix param that is character ('modular',etc) or random (random seed 12345)
  )
  
  
  ##----------------------
  ## 2. ACTOR STRATEGY:  constant (for now)
  ##----------------------
  
  # 2.a. Actor Strategies List (covariates X of SAOM structural covariate effecs like outActX,inPopX,simX,... )
  strategies <- list(
    egoX   =  c(-1,0, 1), #c(0),
    inPopX =  c(1,0, -1)  #c(0),
  )
  # strategies <- list(
  #   egoX   =  c(0),
  #   inPopX =  c(0)
  # )
  
  ## 2.b. Component Payoffs vector
  component_payoffs <-  runif(environ_params$N, min = 0, max = 1)
  
  
  ## 2. Strategies sets the objective function as a linear combination of network stats across DVs
  #
  actor_strats <- lapply(strategies, function(strat) rep(strat,  environ_params$M/length(strat)) )
  #
  structure_model <- list(
    dv_bipartite = list(
      name = 'self$bipartite_rsienaDV',
      effects = list( ##**STRUCTURAL EFFECTS -- dyadic/network endogeneity sources**
        list(effect='density', parameter= -10, fix=T, dv_name=DV_NAME), ##interaction1 = NULL
        list(effect='inPop',   parameter=  -1,  fix=T, dv_name=DV_NAME), #interaction1 = NUL
        list(effect='outAct',  parameter= -1, fix=T, dv_name=DV_NAME)#, #interaction1 = NULL
        # list(effect='outInAss', parameter=0, fix=F, dv_name=DV_NAME), #interaction1 = NULL
        # list(effect='cycle4', parameter=.5, fix=T, dv_name=DV_NAME)#, #interaction1 = NULL
      ),
      coCovars = list( ##**STRATEGY -- MONADIC CONSTANT COVARIATE EFFECTS **
        list(effect='altX',   parameter= 15, fix=T,dv_name=DV_NAME,interaction1='self$component_1_coCovar', x= component_payoffs ),
        # list(effect='outActX',parameter= .5,fix=T,dv_name=DV_NAME,interaction1='self$component_1_coCovar', x= component_payoffs ), #interaction1 = NULL
        list(effect='egoX',   parameter= 1,fix=T,dv_name=DV_NAME,interaction1='self$strat_1_coCovar', x= actor_strats[[1]] ), #interaction1 = NULL
        list(effect='inPopX', parameter= 1, fix=T,dv_name=DV_NAME,interaction1='self$strat_2_coCovar', x= actor_strats[[2]] )#, #interaction1 = NULL
        # list(effect=c('egoX','inPopX'), parameter= -2, fix=T,dv_name=DV_NAME,interaction1=c('self$strat_1_coCovar','self$strat_2_coCovar'), x= actor_strats[[2]] )#, #interaction1 = NULL
        # list(effect=c('outActX','inPopX'),parameter= -.3, fix=T, dv_name=DV_NAME,interaction1=c('self$strat_1_coCovar','self$strat_2_coCovar'),
        #      x = c(1,0,-1,1,0,-1,1,0,-1) )#, #interaction1 = NULL
      ),
      varCovars = list() ##**MONADIC TIME-VARYING COVARIATE EFFECTS -- DYNAMIC STRATEGY PROGRAMS**
    )
  )
  
  
  
  
  
  
  
  ##****************************************##
  ## II. SIM ANALYSIS 
  ##****************************************##
  ###############  Load R6 Class DEPENDENCIES ############################
  ## Biparite Environment Search Simulation Class
  SaomNkRSienaBiEnv <- source(file.path(dir_proj, 'SAOM_NK_R6_model.R'))$value
  # ## RSiena search Class
  # SaomNkRSienaBiEnv_search_rsiena <- source(file.path(dir_proj, 'SAOM_NK_R6_search_rsiena_model.R'))$value
  ###########
  ## Working director
  setwd(dir_proj)
  
  
  
  
  ## INIT SIM ENVIRONMENT: 
  env1 <- SaomNkRSienaBiEnv$new(environ_params)
  ## 1.1. Search 1. SHORT Run
  env1$search_rsiena_multiwave_run(
    structure_model, 
    waves=1, ##id='_short_run',
    iterations = env1$M * env1$N, 
    rand_seed = 12345
  )
  
  
  
  ## 1.2. Search 2. LONG Run
  env2 <- SaomNkRSienaBiEnv$new(environ_params)
  env2$search_rsiena_multiwave_run(
    structure_model, 
    waves=1,
    iterations = 8000,  ## 50 * env2$M * round(sqrt(env2$N)), 
    rand_seed = 12345
  )
  
  # Process results
  env1$search_rsiena_multiwave_process_results()
  env2$search_rsiena_multiwave_process_results()
  
  # Plot Utility
  env1$search_rsiena_multiwave_plot('utility_strategy_summary', 
                                    plot_file = as.character(as.numeric(Sys.time())) )  ## thin_factor = 1
  env2$search_rsiena_multiwave_plot('utility_strategy_summary',  thin_factor = 5, 
                                    plot_file = as.character(as.numeric(Sys.time())) ) ## thin_factor = 1
  
  # Plot K 4panel (degree evolution)
  env1$search_rsiena_multiwave_plot('K_4panel', 
                                    plot_file = as.character(as.numeric(Sys.time())) )  ## thin_factor = 1
  env2$search_rsiena_multiwave_plot('K_4panel', thin_factor = 5, 
                                    plot_file = as.character(as.numeric(Sys.time())) ) ## thin_factor = 1
  
#   
#   
#   outlist[[ sprintf('N%s',Ns[i]) ]] <- list(env1=env1, env2=env2, N=Ns[i])
# }
# 
# 
# 
# saveRDS(outlist, 
#         file = file.path(dir_data, sprintf('%s_%s.rds',environ_params$name,as.character(as.integer(Sys.time())))) )
# ## end




