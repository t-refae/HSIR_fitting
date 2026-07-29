library(targets)
library(tarchetypes)      # tar_map()
library(stantargets)
library(dplyr)

tar_source()

cfg       <- load_params("config.yml")

w <- Sys.getenv("HSIR_WARMUP");   if (nzchar(w)) cfg$iter_warmup   <- as.integer(w)
s <- Sys.getenv("HSIR_SAMPLING"); if (nzchar(s)) cfg$iter_sampling <- as.integer(s)

vw <- Sys.getenv("VOI_WARMUP");   if (nzchar(vw)) cfg$voi_iter_warmup   <- as.integer(vw)
vs <- Sys.getenv("VOI_SAMPLING"); if (nzchar(vs)) cfg$voi_iter_sampling <- as.integer(vs)

pw <- Sys.getenv("PTS_WARMUP");   if (nzchar(pw)) cfg$pts_iter_warmup   <- as.integer(pw)
ps <- Sys.getenv("PTS_SAMPLING"); if (nzchar(ps)) cfg$pts_iter_sampling <- as.integer(ps)

stan_file <- select_stan_file(cfg)
scenarios <- read_scenarios(cfg)          # one row per param combo

## subset testing
# scenarios <- scenarios %>% filter(cv != 0)
# print(scenarios)

tar_option_set(
  packages = c("cmdstanr", "posterior", "deSolve", "dplyr", "yaml", "crew", "gtools"),
  format   = "rds",
  controller = crew::crew_controller_local(workers = cfg$workers, seconds_idle = 30),
  storage    = "worker",   # workers write draws straight to the store
  retrieval  = "worker"    # avoids shipping large draw objects back to main
)

#### FTS ####

fits <- tar_map(
  values = scenarios,
  names  = id,
  
  tar_target(scen_params, scenario_params(cfg, beta, gamma, cv)),
  tar_target(sim_data,    simulate_epidemic(scen_params)),
  
  tar_stan_mcmc(
    fit,
    stan_files      = stan_file,
    data            = build_stan_data(sim_data, scen_params),
    seed            = cfg$seed,
    chains          = cfg$chains,
    parallel_chains = cfg$parallel_chains,
    iter_warmup     = cfg$iter_warmup,
    iter_sampling   = cfg$iter_sampling,
    init            = make_inits(cfg$model_type, cfg$chains, cfg$seed),
    refresh         = cfg$refresh
  )
)

model_name <- tools::file_path_sans_ext(basename(stan_file))

export_targets <- lapply(scenarios$id, function(sid) {
  draws_sym <- as.symbol(sprintf("fit_draws_%s_%s",       model_name, sid))
  summ_sym  <- as.symbol(sprintf("fit_summary_%s_%s",     model_name, sid))
  diag_sym  <- as.symbol(sprintf("fit_diagnostics_%s_%s", model_name, sid))
  out_path  <- file.path("outputs/FTS", sprintf("fit_%s_%s.rds", model_name, sid))
  tar_target_raw(
    sprintf("fit_export_%s", sid),
    substitute(
      write_fit_bundle(D, S, G, P),
      list(D = draws_sym, S = summ_sym, G = diag_sym, P = out_path)
    ),
    format = "file"
  )
})

#### VOI ####
voi_scenarios  <- voi_read_scenarios(cfg)         # id, beta, gamma, cv, regime
voi_mcmc       <- voi_mcmc_settings(cfg)
voi_stan_sir   <- voi_stan_file(cfg, "SIR")
voi_stan_hsir  <- voi_stan_file(cfg, "HSIR")
voi_sir_model  <- tools::file_path_sans_ext(basename(voi_stan_sir))
voi_hsir_model <- tools::file_path_sans_ext(basename(voi_stan_hsir))

# Homogeneous (SIR) fits to the HSIR-generated data
voi_fits_sir <- tar_map(
  values = voi_scenarios,
  names  = id,
  
  tar_target(voi_scen_params_sir, voi_scenario_params(cfg, beta, gamma, cv, regime)),
  tar_target(voi_sim_data_sir,    voi_simulate_truth(voi_scen_params_sir)),
  
  tar_stan_mcmc(
    voi_sir,
    stan_files      = voi_stan_sir,
    data            = voi_build_stan_data(voi_sim_data_sir, voi_scen_params_sir),
    seed            = voi_mcmc$seed,
    chains          = voi_mcmc$chains,
    parallel_chains = voi_mcmc$parallel_chains,
    iter_warmup     = voi_mcmc$iter_warmup,
    iter_sampling   = voi_mcmc$iter_sampling,
    refresh         = voi_mcmc$refresh
  )
)

# Heterogeneous (HSIR) fits to the same HSIR-generated data
voi_fits_hsir <- tar_map(
  values = voi_scenarios,
  names  = id,
  
  tar_target(voi_scen_params_hsir, voi_scenario_params(cfg, beta, gamma, cv, regime)),
  tar_target(voi_sim_data_hsir,    voi_simulate_truth(voi_scen_params_hsir)),
  
  tar_stan_mcmc(
    voi_hsir,
    stan_files      = voi_stan_hsir,
    data            = voi_build_stan_data(voi_sim_data_hsir, voi_scen_params_hsir),
    seed            = voi_mcmc$seed,
    chains          = voi_mcmc$chains,
    parallel_chains = voi_mcmc$parallel_chains,
    iter_warmup     = voi_mcmc$iter_warmup,
    iter_sampling   = voi_mcmc$iter_sampling,
    refresh         = voi_mcmc$refresh
  )
)

voi_export_sir <- lapply(voi_scenarios$id, function(sid) {
  draws_sym <- as.symbol(sprintf("voi_sir_draws_%s_%s",       voi_sir_model, sid))
  summ_sym  <- as.symbol(sprintf("voi_sir_summary_%s_%s",     voi_sir_model, sid))
  diag_sym  <- as.symbol(sprintf("voi_sir_diagnostics_%s_%s", voi_sir_model, sid))
  out_path  <- file.path("outputs/VOI", sprintf("voi_fit_%s_%s.rds", voi_sir_model, sid))
  tar_target_raw(
    sprintf("voi_export_sir_%s", sid),
    substitute(
      write_fit_bundle(D, S, G, P),
      list(D = draws_sym, S = summ_sym, G = diag_sym, P = out_path)
    ),
    format = "file"
  )
})

voi_export_hsir <- lapply(voi_scenarios$id, function(sid) {
  draws_sym <- as.symbol(sprintf("voi_hsir_draws_%s_%s",       voi_hsir_model, sid))
  summ_sym  <- as.symbol(sprintf("voi_hsir_summary_%s_%s",     voi_hsir_model, sid))
  diag_sym  <- as.symbol(sprintf("voi_hsir_diagnostics_%s_%s", voi_hsir_model, sid))
  out_path  <- file.path("outputs/VOI", sprintf("voi_fit_%s_%s.rds", voi_hsir_model, sid))
  tar_target_raw(
    sprintf("voi_export_hsir_%s", sid),
    substitute(
      write_fit_bundle(D, S, G, P),
      list(D = draws_sym, S = summ_sym, G = diag_sym, P = out_path)
    ),
    format = "file"
  )
})

#### PTS ####
# Partial time series: FTS combos (cv != 0) x {-2GI, -1GI, at peak}, fitted with
# the ORIGINAL HSIR Stan file (same as FTS, NOT the VOI Stan files). Reuses the
# FTS machinery so the simulation and Stan-data indexing match FTS; only the
# built Stan-data list is truncated to the horizon. Own MCMC settings (pts_mcmc).
pts_scenarios <- pts_make_scenarios(scenarios)   # 27 rows for the default grid
pts_mcmc      <- pts_mcmc_settings(cfg)
pts_stan_file <- stan_file                       # original HSIR Stan, as FTS
pts_model     <- model_name

pts_fits <- tar_map(
  values = pts_scenarios,
  names  = id,
  
  tar_target(pts_scen_params, scenario_params(cfg, beta, gamma, cv)),
  tar_target(pts_sim_data,    simulate_epidemic(pts_scen_params)),
  
  tar_stan_mcmc(
    pts_fit,
    stan_files      = pts_stan_file,
    data            = pts_build_stan_data(pts_sim_data, pts_scen_params, gamma, n_before),
    seed            = pts_mcmc$seed,
    chains          = pts_mcmc$chains,
    parallel_chains = pts_mcmc$parallel_chains,
    iter_warmup     = pts_mcmc$iter_warmup,
    iter_sampling   = pts_mcmc$iter_sampling,
    init            = make_inits(cfg$model_type, pts_mcmc$chains, pts_mcmc$seed),
    refresh         = pts_mcmc$refresh,
    
    adapt_delta   = pts_mcmc$adapt_delta,
    max_treedepth = pts_mcmc$max_treedepth,
    metric        = pts_mcmc$metric,
  )
)

pts_export <- lapply(pts_scenarios$id, function(sid) {
  draws_sym <- as.symbol(sprintf("pts_fit_draws_%s_%s",       pts_model, sid))
  summ_sym  <- as.symbol(sprintf("pts_fit_summary_%s_%s",     pts_model, sid))
  diag_sym  <- as.symbol(sprintf("pts_fit_diagnostics_%s_%s", pts_model, sid))
  out_path  <- file.path("outputs/PTS", sprintf("pts_fit_%s_%s.rds", pts_model, sid))
  tar_target_raw(
    sprintf("pts_export_%s", sid),
    substitute(
      write_fit_bundle(D, S, G, P),
      list(D = draws_sym, S = summ_sym, G = diag_sym, P = out_path)
    ),
    format = "file"
  )
})


list(
  #### FTS ####
  tar_target(scenarios_manifest, scenarios),
  
  fits,
  
  export_targets,
  
  tarchetypes::tar_combine(
    conv_summary,
    tarchetypes::tar_select_targets(fits, starts_with("fit_summary")),
    command = dplyr::bind_rows(!!!.x, .id = "target")
  ),
  
  tarchetypes::tar_combine(
    conv_diag,
    tarchetypes::tar_select_targets(fits, starts_with("fit_diagnostics")),
    command = dplyr::bind_rows(!!!.x, .id = "target")
  ),
  
  #### VOI ####
  tar_target(voi_scenarios_manifest, voi_scenarios),
  
  voi_fits_sir,
  voi_fits_hsir,
  
  voi_export_sir,
  voi_export_hsir,
  
  tarchetypes::tar_combine(
    voi_conv_summary,
    tarchetypes::tar_select_targets(voi_fits_sir,  starts_with("voi_sir_summary")),
    tarchetypes::tar_select_targets(voi_fits_hsir, starts_with("voi_hsir_summary")),
    command = dplyr::bind_rows(!!!.x, .id = "target")
  ),
  
  tarchetypes::tar_combine(
    voi_conv_diag,
    tarchetypes::tar_select_targets(voi_fits_sir,  starts_with("voi_sir_diagnostics")),
    tarchetypes::tar_select_targets(voi_fits_hsir, starts_with("voi_hsir_diagnostics")),
    command = dplyr::bind_rows(!!!.x, .id = "target")
  ),
  
  #### PTS ####
  tar_target(pts_scenarios_manifest, pts_scenarios),
  
  pts_fits,
  
  pts_export,
  
  tarchetypes::tar_combine(
    pts_conv_summary,
    tarchetypes::tar_select_targets(pts_fits, starts_with("pts_fit_summary")),
    command = dplyr::bind_rows(!!!.x, .id = "target")
  ),
  
  tarchetypes::tar_combine(
    pts_conv_diag,
    tarchetypes::tar_select_targets(pts_fits, starts_with("pts_fit_diagnostics")),
    command = dplyr::bind_rows(!!!.x, .id = "target")
  )
)