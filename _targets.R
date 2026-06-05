# _targets.R

library(targets)
library(stantargets)

tar_source()

# Read + validate config ONCE at construction time. tar_stan_mcmc() needs the
# .stan path here (it names targets from it), so the model switch is resolved
# now. `cfg` is a tracked global: editing config.yml invalidates the affected
# targets on the next tar_make().
cfg       <- load_params("config.yml")
stan_file <- select_stan_file(cfg)

tar_option_set(
  packages = c("cmdstanr", "posterior", "deSolve", "dplyr", "yaml"),
  format   = "rds"
)

list(
  tar_target(sim_data, simulate_epidemic(cfg)),
  
  tar_stan_mcmc(
    fit,
    stan_files      = stan_file,
    data            = build_stan_data(sim_data, cfg),
    seed            = cfg$seed,
    chains          = cfg$chains,
    parallel_chains = cfg$parallel_chains,
    iter_warmup     = cfg$iter_warmup,
    iter_sampling   = cfg$iter_sampling,
    init            = make_inits(cfg$model_type, cfg$chains, cfg$seed),
    refresh         = cfg$refresh
  )
)