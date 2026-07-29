# ============================================================================
# PTS fits — partial time series, fitted with the ORIGINAL HSIR Stan file
# ----------------------------------------------------------------------------
# Same (beta, gamma, cv) combinations as the full-time-series (FTS) fits,
# excluding cv = 0, each fitted at three observation horizons relative to the
# incidence peak (in generation intervals, GI = 1/gamma):
#
#       -2GI     two generation intervals before the peak
#       -1GI     one generation interval before the peak
#       at peak  up to the peak
#
# PTS reuses the FTS machinery (scenario_params, simulate_epidemic,
# build_stan_data, select_stan_file, make_inits). The HSIR Stan model fits a
# partial window via its `n_fit` field (the likelihood uses cases[1:n_fit]),
# so PTS leaves the full Stan-data list untouched and simply sets n_fit to the
# horizon length. The full epidemic is still simulated and integrated; only the
# likelihood window changes. PTS has its own MCMC settings (cfg$pts_*), falling
# back to the main values when unset.
# ============================================================================

.pts_or <- function(x, default) if (is.null(x)) default else x

# --- PTS's own MCMC settings (independent of FTS) ----------------------------
pts_mcmc_settings <- function(cfg) {
  list(
    chains          = .pts_or(cfg$pts_chains,          cfg$chains),
    parallel_chains = .pts_or(cfg$pts_parallel_chains, cfg$parallel_chains),
    iter_warmup     = .pts_or(cfg$pts_iter_warmup,     cfg$iter_warmup),
    iter_sampling   = .pts_or(cfg$pts_iter_sampling,   cfg$iter_sampling),
    refresh         = .pts_or(cfg$pts_refresh,         cfg$refresh),
    seed            = .pts_or(cfg$pts_seed,            cfg$seed),
    
    adapt_delta   = .pts_or(cfg$pts_adapt_delta,   0.8),
    max_treedepth = .pts_or(cfg$pts_max_treedepth, 10L),
    metric        = .pts_or(cfg$pts_metric,        "diag_e")
  )
}

# --- Scenario grid: FTS combos (cv != 0) x three horizons --------------------
pts_make_scenarios <- function(scenarios) {
  horizons <- data.frame(
    code     = c("m2GI", "m1GI", "peak"),
    n_before = c(2, 1, 0),
    horizon  = c("-2GI", "-1GI", "at peak"),
    stringsAsFactors = FALSE
  )
  base <- scenarios[scenarios$cv != 0, c("id", "beta", "gamma", "cv"), drop = FALSE]
  grid <- merge(base, horizons, by = character(0))          # cross join
  grid$id   <- paste(grid$id, grid$code, sep = "_")
  grid$code <- NULL
  grid <- grid[order(grid$id), c("id", "beta", "gamma", "cv", "n_before", "horizon")]
  rownames(grid) <- NULL
  if (anyDuplicated(grid$id)) {
    stop("PTS scenarios: `id` values are not unique.", call. = FALSE)
  }
  grid
}

# --- Set the partial-fit window via n_fit ------------------------------------
# Leaves cases/ts/n_days/y0/t0/N at full length (so all Stan-data indexing is
# identical to FTS) and sets n_fit to the number of days up to the horizon.
# n_before = generation intervals before the peak (0 = at the peak); the peak is
# taken from the observed case series.
pts_set_fit_window <- function(stan_data, gamma, n_before) {
  peak <- which.max(stan_data$cases)
  keep <- round(peak - n_before * (1 / gamma))
  keep <- min(max(keep, 3L), length(stan_data$cases))
  stan_data$n_fit <- as.integer(keep)
  stan_data
}

# --- Build PTS Stan data: full FTS build, then restrict the fit window -------
pts_build_stan_data <- function(sim_data, scen_params, gamma, n_before) {
  full <- build_stan_data(sim_data, scen_params)
  pts_set_fit_window(full, gamma, n_before)
}

.pts_num <- function(d, nm) {
  if (is.null(d) || is.null(d[[nm]])) return(NA_real_)
  v <- suppressWarnings(as.numeric(d[[nm]]))
  if (!length(v) || all(is.na(v))) NA_real_ else v
}
.pts_min <- function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
.pts_sum <- function(v) if (all(is.na(v))) NA_real_ else sum(v, na.rm = TRUE)

.pts_ebfmi <- function(dr) {
  if (!"energy__" %in% colnames(dr)) return(NA_real_)
  v <- tapply(dr$energy__, dr$.chain, function(E) {
    if (length(E) < 2L || stats::var(E) == 0) NA_real_
    else sum(diff(E)^2) / length(E) / stats::var(E)
  })
  .pts_min(as.numeric(v))
}

pts_summarise_convergence <- function(cfg = load_params("config.yml"),
                                      pts_dir = "outputs/PTS") {
  grid  <- pts_make_scenarios(read_scenarios(cfg))
  model <- tools::file_path_sans_ext(basename(select_stan_file(cfg)))
  
  rows <- list()
  for (i in seq_len(nrow(grid))) {
    sid  <- grid$id[i]
    path <- file.path(pts_dir, sprintf("pts_fit_%s_%s.rds", model, sid))
    if (!file.exists(path)) next
    
    b  <- readRDS(path)
    dr <- posterior::as_draws_df(as.data.frame(b$draws))
    s  <- posterior::summarise_draws(dr, posterior::default_convergence_measures())
    d  <- b$diagnostics
    
    eb <- .pts_min(.pts_num(d, "ebfmi"))
    if (is.na(eb)) eb <- .pts_ebfmi(dr)
    
    rows[[length(rows) + 1L]] <- data.frame(
      target             = sprintf("pts_fit_mcmc_%s_%s", model, sid),
      id                 = sid,
      horizon            = grid$horizon[i],
      beta               = grid$beta[i],
      gamma              = grid$gamma[i],
      cv                 = grid$cv[i],
      max_rhat           = max(s$rhat,     na.rm = TRUE),
      min_ess_bulk       = min(s$ess_bulk, na.rm = TRUE),
      min_ess_tail       = min(s$ess_tail, na.rm = TRUE),
      divergences        = .pts_sum(.pts_num(d, "num_divergent")),
      max_treedepth_hits = .pts_sum(.pts_num(d, "num_max_treedepth")),
      min_ebfmi          = eb,
      respliced          = !is.null(b$rerun),
      stringsAsFactors   = FALSE
    )
  }
  if (!length(rows)) stop("No PTS bundles found in ", pts_dir, call. = FALSE)
  res <- dplyr::bind_rows(rows)
  
  meta <- targets::tar_meta(fields = "seconds")
  res$seconds <- meta$seconds[match(res$target, meta$name)]
  
  res$ok <- round(res$max_rhat, digits = 3) <= 1.02 &
    (is.na(res$divergences) | res$divergences <= 100) &
    (is.na(res$min_ebfmi)   | res$min_ebfmi >= 0.3)
  
  res[order(res$id), c("id", "horizon", "beta", "gamma", "cv",
                       "max_rhat", "min_ess_bulk", "min_ess_tail",
                       "divergences", "max_treedepth_hits", "min_ebfmi",
                       "seconds", "respliced", "ok")]
}


#### REMOVE FROM HERE LATER ####
pts_posterior_draws <- function(cfg = load_params("config.yml"),
                                params = c("beta", "cv", "R0"),
                                pts_dir = "outputs/PTS") {
  grid  <- pts_make_scenarios(read_scenarios(cfg))
  model <- tools::file_path_sans_ext(basename(select_stan_file(cfg)))
  hlev  <- c("-2GI", "-1GI", "at peak")
  
  true_value <- function(p, beta, gamma, cv) switch(
    p, beta = beta, gamma = gamma, cv = cv,
    R0 = beta / gamma, D = 1 / gamma, NA_real_
  )
  
  out <- list()
  for (i in seq_len(nrow(grid))) {
    sid  <- grid$id[i]
    path <- file.path(pts_dir, sprintf("pts_fit_%s_%s.rds", model, sid))
    if (!file.exists(path)) next
    
    dr <- as.data.frame(readRDS(path)$draws)
    if (!"R0"    %in% names(dr) && all(c("beta", "D") %in% names(dr))) dr$R0    <- dr$beta * dr$D
    if (!"gamma" %in% names(dr) && "D" %in% names(dr))                 dr$gamma <- 1 / dr$D
    
    for (p in intersect(params, names(dr))) {
      out[[length(out) + 1L]] <- data.frame(
        combo     = sub("_(m2GI|m1GI|peak)$", "", sid),
        horizon   = factor(grid$horizon[i], levels = hlev),
        parameter = p,
        value     = dr[[p]],
        true      = true_value(p, grid$beta[i], grid$gamma[i], grid$cv[i]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out)) stop("No PTS bundles found in ", pts_dir, call. = FALSE)
  dplyr::bind_rows(out)
}

pts_ridge_plot <- function(params = c("beta", "cv", "R0"),
                           draws = NULL, show_truth = TRUE,
                           cfg = load_params("config.yml")) {
  if (is.null(draws)) draws <- pts_posterior_draws(cfg, params)
  
  panel <- function(p) {
    d <- draws[draws$parameter == p, , drop = FALSE]
    g <- ggplot2::ggplot(d, ggplot2::aes(x = value, y = combo, fill = combo)) +
      ggridges::geom_density_ridges(alpha = 0.7, scale = 1.05,
                                    rel_min_height = 0.01, colour = "white") +
      ggplot2::facet_wrap(~ horizon, nrow = 1) +
      ggplot2::labs(x = p, y = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "none")
    if (show_truth) {
      tv <- unique(d[, c("combo", "horizon", "true")])
      g <- g + ggplot2::geom_point(
        data = tv, inherit.aes = FALSE,
        ggplot2::aes(x = true, y = combo), shape = 124, size = 3
      )
    }
    g
  }
  
  panels <- lapply(params, panel)
  if (length(panels) == 1) panels[[1]] else patchwork::wrap_plots(panels, ncol = 1)
}
