library(targets)
library(tarchetypes)      # tar_map()
library(stantargets)

tar_source()

cfg       <- load_params("config.yml")
stan_file <- select_stan_file(cfg)
scenarios <- read_scenarios(cfg)          # one row per param combo
print(scenarios)

# testing
scenarios <- scenarios[1,]

tar_option_set(
  packages = c("cmdstanr", "posterior", "deSolve", "dplyr", "yaml"),
  format   = "rds"
)

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

list(
  tar_target(scenarios_manifest, scenarios),
  fits
)