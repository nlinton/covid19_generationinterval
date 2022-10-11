
# Result functions ----

run_custom_stan <- function(stan_file=NULL, iter_nums=c(1250, 2000), data_list=NULL, init_list=NULL, fixed_param=F) {
  
  options(mc.cores = parallel::detectCores())
  cmdstan_path = "~/miniconda3/bin/cmdstan"
  
  data_file = tempfile(fileext=".json")
  cmdstanr::write_stan_json(data_list, data_file)
  
  recalc = TRUE
  output_dir = "stan_output"
  if(!dir.exists(output_dir))
    dir.create(output_dir)
  if(recalc==TRUE) {
    start_time <- Sys.time()
    set_mod <- cmdstan_model(stan_file)
    mod <- set_mod$sample(
      data=data_file,
      # init=init_file,
      iter_sampling=iter_nums[1],
      iter_warmup=iter_nums[2],
      chains=4,
      parallel_chains=4,
      adapt_delta=0.98,
      seed=123,
      fixed_param=fixed_param)
    end_time <- Sys.time()
    
    out_dir <- "stan_output"
    dir.create(out_dir)
    mod$save_output_files(dir=out_dir, 
                          basename="variantRt",
                          timestamp=FALSE,
                          random=FALSE)
    mod$save_data_file(dir=out_dir,
                       basename="data",
                       timestamp=FALSE,
                       random=FALSE)
    # print(end_time-start_time)
    # mod$cmdstan_diagnose() # check everything is ok
  }
  else {
    output_files = c()
    for (fl_ in list.files(output_dir, "*.csv", full.names=TRUE))
      if(!grepl("diagnostic", fl_, fixed=TRUE))
        output_files=c(output_files, fl_)
    mod <- as_cmdstan_fit(output_files)
  }
  return(mod)
}


sel_vals <- function(mod) {
  selected_pars = c("index_copula", "index_GI", "index_IP")
  list(rbindlist(lapply(seq_along(selected_pars), function(x) { mod$draws(selected_pars[x]) %>% unlist %>% table %>% {./sum(.)} %>% as.data.frame %>% rename(value := ".") %>%
      arrange(desc(Freq)) %>% slice(1L) %>% dplyr::select(value) })) %>% pull() %>% as.character() %>% as.numeric())
}


# Plot functions ----

## Figure 1 ----

cont_plot <- function(mod) {
  
  selected_vals <- sel_vals(mod)
  plot_pars <- c(paste0("param1[1,",selected_vals[[1]][2],"]"), # [GI, selected dist]
                 paste0("param2[1,",selected_vals[[1]][2],"]"), # [GI, selected dist]
                 paste0("param1[2,",selected_vals[[1]][3],"]"), # [IP, selected dist]
                 paste0("param2[2,",selected_vals[[1]][3],"]"), # [IP, selected dist]
                 "tau")
  fitted_pars <- mod$summary(plot_pars, ~quantile(.x, probs = c(0.5)))
  fitted_pars <- data.frame(fitted_pars[,2]) %>% pull()
  
  library(copula)
  copcol = scales::hue_pal()(7)[5]
  xlim = c(0,10)
  ylim = c(0,10)
  zguide = seq(0,0.18,by=0.01)
  ncuts = 10
  ngrid = 200
  
  tau = as.numeric(fitted_pars[5])
  clayton_theta = (2*tau)/(1-tau)
  fitted_cop <- mvdc(copula=claytonCopula(eval(parse(text=paste0(cops[selected_vals[[1]][1]],"_theta"))), dim=2), 
                     margins=c(dists_short[selected_vals[[1]][2]],dists_short[selected_vals[[1]][3]]),
                     paramMargins=list(list(shape=fitted_pars[1], scale=fitted_pars[2]),
                                       list(meanlog=fitted_pars[3], sdlog=fitted_pars[4])))
  ggplotify::as.grob(contourplot2(fitted_cop, dMvdc,
                                  at=zguide,
                                  xlim = xlim, ylim = ylim, cuts = ncuts, n.grid = ngrid,
                                  col.regions = colorRampPalette(c("white", copcol), space="Lab"),
                                  xlab = list("Generation interval (days)",cex=1.6),
                                  ylab = list("Incubation period (days)",cex=1.6),
                                  scales=list(x=list(cex=1.5),y=list(cex=1.5))))
}

## Figure 2 ----

slabs <- c("All cases", "", "Infector under 30", "Infector aged 30-59", "Infector ages 60+", "Infectee under 30", "Infectee aged 30-59", "Infectee ages 60+", "", 
           "Female infectors", "Male infectors", "Female infectees", "Male infectees", "", "Household", "Social contact", "Community", "", 
           "Wave 1", "Wave 2", "Wave 3", "", "Cluster-related", "Contact chain-related", "Domestic travel-related", "Import-related", "", "Asymptomatic infectors")

strat_l <- function(dat,lab) {
  
  options(repr.plot.width=6,repr.plot.height=6)
  cols = scales::hue_pal()(7)[c(6,3,2,4,5,1,7)] 
  vline_med = dat$median[1]
  vline_lci = dat$lci[1]
  vline_uci = dat$uci[1]
  ggplot(dat) + 
    geom_rect(aes(x=median, y=desc(y)), xmin=vline_lci, xmax=vline_uci, ymin=21, ymax=-Inf, fill="#f2f3f4") +
    geom_vline(xintercept=vline_med, color="#dedede", size=2) +
    geom_errorbarh(aes(xmin=lci, xmax=uci, y=desc(y), color=data_cat), height=0, size=2) +
    geom_point(aes(x=median, y=desc(y)), color="black") +
    scale_y_continuous("", breaks=-(1:length(slabs)), labels=slabs[1:length(slabs)]) +
    scale_x_continuous(lab, breaks=seq(0,7,1), expand=c(0,0), limits=c(0,7)) +
    scale_color_manual(values=cols, guide="none") +
    theme_classic() +
    theme(text = element_text(size=16),
          axis.text = element_text(size=16),
          plot.margin = margin(1, 0.2, 1, 0.2, "cm"))
}
strat_r <- function(dat,lab) {
  options(repr.plot.width=4,repr.plot.height=6)
  cols = scales::hue_pal()(7)[c(6,3,2,4,5,1,7)] 
  vline_med = dat$median[1]
  vline_lci = dat$lci[1]
  vline_uci = dat$uci[1]
  ggplot(dat) + 
    geom_rect(aes(x=median, y=desc(y)), xmin=vline_lci, xmax=vline_uci, ymin=21, ymax=-Inf, fill="#f2f3f4") +
    geom_vline(xintercept=vline_med, color="#dedede", size=2) +
    geom_errorbarh(aes(xmin=lci, xmax=uci, y=desc(y), color=data_cat), height=0, size=2) +
    geom_point(aes(x=median, y=desc(y)), color="black") +
    scale_y_continuous("", breaks=-(1:length(slabs)), labels=NULL) +
    scale_x_continuous(lab, breaks=seq(0,10,1), expand=c(0,0), limits=c(0,8)) +
    scale_color_manual(values=cols, guide="none") +
    theme_classic() +
    theme(text = element_text(size=16),
          axis.text.x = element_text(size=16),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(1, 0.2, 1, 0.2, "cm")) #tlbr
  
}



post_prior <- function(mod) {
  
  options(repr.plot.width=10, repr.plot.height=5)
  days_lab <- function(x) paste0(x, " (days)")
  set.seed(123)
  par = c('mu[1]', 'mu[2]')
  var = c("Generation interval", "Incubation period")
  posterior = rbindlist(lapply(seq_along(par), function(x) { mod$draws(par[x], format = "df") %>% select(-contains('.')) %>% rename(mu=1) %>% mutate(variable = var[x], type = "posterior") } ))
  prior = sapply(seq_along(par), function(x) { rnorm(nrow(posterior)/length(par), mu_mean_priors[x], mu_sigma_priors[x]) } ) %>% data.frame() %>% rename_all(~ var) %>% 
    pivot_longer(names(.), names_to = "variable", values_to = "mu") %>% relocate(mu) %>% mutate(type = "prior")
  yrep = bind_rows(posterior, prior)
  
  ggplot(data = yrep) + 
    geom_density(aes(mu, fill=type), alpha=.5) + 
    facet_wrap(~variable, strip.position = c("bottom"), labeller=as_labeller(days_lab)) +
    scale_y_continuous("Density") +
    scale_x_continuous("", limits=c(0,9), expand=c(0,0), breaks=seq(0,8,2)) +
    scale_fill_discrete("") +
    theme_classic() +
    theme(text = element_text(size=16),
          strip.background=element_rect(color="white", fill="white"),
          strip.text = element_text(size = 16, margin = margin(0.2, 0, 0.2, 0, "cm"), hjust=0.5), #
          strip.placement = "outside",
          panel.spacing.x = unit(2, "lines"))
}



