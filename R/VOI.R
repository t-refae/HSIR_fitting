# R/VOI.R

.voi_or <- function(x, default) if (is.null(x)) default else x

# --- Stan file per fitted model (SIR = homogeneous, HSIR = heterogeneous) ---
voi_stan_file <- function(cfg, model) {
  switch(
    model,
    SIR  = .voi_or(cfg$voi_stan_sir,  "Stan/VOI_SIR.stan"),
    HSIR = .voi_or(cfg$voi_stan_hsir, "Stan/VOI_HSIR.stan"),
    stop("VOI: unknown model '", model, "' (use 'SIR' or 'HSIR').", call. = FALSE)
  )
}

# --- VOI's own MCMC settings -------------------------------------------------
# Independent of the main fits. Each key falls back to the shared cfg value only
# if its voi_* counterpart is unset, so setting the voi_* keys (or VOI_WARMUP /
# VOI_SAMPLING env vars, handled in _targets.R) makes VOI fully independent.
voi_mcmc_settings <- function(cfg) {
  list(
    chains          = .voi_or(cfg$voi_chains,          cfg$chains),
    parallel_chains = .voi_or(cfg$voi_parallel_chains, cfg$parallel_chains),
    iter_warmup     = .voi_or(cfg$voi_iter_warmup,     cfg$iter_warmup),
    iter_sampling   = .voi_or(cfg$voi_iter_sampling,   cfg$iter_sampling),
    refresh         = .voi_or(cfg$voi_refresh,         cfg$refresh),
    seed            = .voi_or(cfg$voi_seed,            cfg$seed)
  )
}

# --- Per-scenario parameters, resolved from the grid row + shared cfg --------
voi_scenario_params <- function(cfg, beta, gamma, cv, regime) {
  list(
    beta          = beta,
    gamma         = gamma,
    cv            = cv,
    regime        = as.character(regime),
    P             = .voi_or(cfg$voi_P,  1e4),
    i0            = .voi_or(cfg$voi_i0, 1 / .voi_or(cfg$voi_P, 1e4)),
    seed          = cfg$seed,
    max_days      = cfg$max_days,
    end_threshold = cfg$end_threshold
  )
}

# --- Data-generating process: heterogeneous-susceptibility SIR ---------------
voi_hsir_ode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    foi <- beta * I * S^(1 + cv^2)
    list(c(-foi, foi - gamma * I, gamma * I, foi))
  })
}

voi_simulate_truth <- function(params) {
  y0 <- c(S = 1 - params$i0, I = params$i0, R = 0, C = 0)
  solution <- as.data.frame(deSolve::ode(
    y      = y0,
    times  = seq_len(params$max_days),
    func   = voi_hsir_ode,
    parms  = c(beta = params$beta, gamma = params$gamma, cv = params$cv),
    method = "lsoda"
  ))
  mean_incidence <- diff(solution$C) * params$P
  active   <- which(mean_incidence >= params$end_threshold)   # trim sub-threshold tail
  last_day <- if (length(active)) max(active) else length(mean_incidence)
  mean_incidence <- mean_incidence[seq_len(max(last_day, 3L))]
  set.seed(params$seed)
  observed_cases <- rpois(length(mean_incidence), mean_incidence)
  list(
    params         = params,
    P              = params$P,
    y0             = y0,
    mean_incidence = mean_incidence,
    observed_cases = observed_cases
  )
}

# --- Crop to the observation horizon -----------------------------------------
voi_crop_by_regime <- function(observed_cases, mean_incidence, gamma, regime) {
  n_before <- switch(
    as.character(regime),
    full = NA_real_, "0GI" = 0, "1GI" = 1, "2GI" = 2,
    stop("VOI: unknown regime '", regime, "' (use full/0GI/1GI/2GI).",
         call. = FALSE)
  )
  if (is.na(n_before)) return(observed_cases)        # whole epidemic
  peak_day <- which.max(mean_incidence)
  horizon  <- round(peak_day - n_before * (1 / gamma))
  horizon  <- min(max(horizon, 3L), length(observed_cases))
  observed_cases[seq_len(horizon)]
}

# --- Assemble the Stan data list ---------------------------------------------
# Project contract (consistent with the main fits): t0 = 1, y0 is the day-1
# state, ts = 2:n_days, length(ts) == length(cases) == n_days - 1. The two VOI
# Stan files are written to this same contract.
voi_build_stan_data <- function(sim_data, params) {
  cases  <- voi_crop_by_regime(sim_data$observed_cases, sim_data$mean_incidence,
                               params$gamma, params$regime)
  n_days <- length(cases) + 1L
  ts     <- 2:n_days
  stopifnot(n_days >= 2L, length(ts) == n_days - 1L, length(cases) == n_days - 1L)
  list(
    n_days = n_days,
    y0     = unname(sim_data$y0),
    t0     = 1,
    ts     = ts,
    N      = sim_data$P,
    cases  = cases
  )
}

# --- Scenario grid: one row per (parameter set x regime) ---------------------
voi_default_scenarios <- function() {
  data.frame(
    id     = c("pts_2GI", "pts_1GI", "pts_0GI", "pts_full"),
    beta   = 1.2,
    gamma  = 0.4,
    cv     = 1,
    regime = c("2GI", "1GI", "0GI", "full"),
    stringsAsFactors = FALSE
  )
}

.voi_safe_id <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  ifelse(grepl("^[0-9]", x), paste0("s_", x), x)
}

voi_read_scenarios <- function(cfg) {
  csv  <- .voi_or(cfg$voi_scenarios_csv, "Data/voi_scenarios.csv")
  grid <- if (!is.null(csv) && file.exists(csv)) {
    utils::read.csv(csv, stringsAsFactors = FALSE)
  } else {
    voi_default_scenarios()
  }
  
  required <- c("id", "beta", "gamma", "cv", "regime")
  missing  <- setdiff(required, names(grid))
  if (length(missing)) {
    stop("VOI scenarios are missing required columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  for (col in c("beta", "gamma", "cv")) grid[[col]] <- as.numeric(grid[[col]])
  
  grid$regime <- as.character(grid$regime)
  bad <- !grid$regime %in% c("full", "0GI", "1GI", "2GI")
  if (any(bad)) {
    stop("VOI scenarios: `regime` must be full/0GI/1GI/2GI, but found: ",
         paste(unique(grid$regime[bad]), collapse = ", "), call. = FALSE)
  }
  
  grid$id <- .voi_safe_id(grid$id)
  if (anyDuplicated(grid$id)) {
    stop("VOI scenarios: `id` values must be unique after sanitisation.",
         call. = FALSE)
  }
  grid[, required]
}

voi_summarise_convergence <- function(cfg = load_params("config.yml"),
                                      voi_dir = "outputs/VOI") {
  grid  <- voi_read_scenarios(cfg)
  stans <- c(SIR  = tools::file_path_sans_ext(basename(voi_stan_file(cfg, "SIR"))),
             HSIR = tools::file_path_sans_ext(basename(voi_stan_file(cfg, "HSIR"))))
  
  rows <- list()
  for (m in names(stans)) {
    for (i in seq_len(nrow(grid))) {
      sid  <- grid$id[i]
      path <- file.path(voi_dir, sprintf("voi_fit_%s_%s.rds", stans[[m]], sid))
      if (!file.exists(path)) next
      
      b  <- readRDS(path)
      dr <- posterior::as_draws_df(as.data.frame(b$draws))
      s  <- posterior::summarise_draws(dr, posterior::default_convergence_measures())
      d  <- b$diagnostics
      
      d  <- b$diagnostics
      dv <- if (!is.null(d) && !is.null(d$divergent__)) as.numeric(d$divergent__) else NA_real_
      td <- if (!is.null(d) && !is.null(d$treedepth__)) as.numeric(d$treedepth__) else NA_real_
      eb <- if (!is.null(d) && !is.null(d$energy__) && !is.null(d$.chain)) {
        .pts_min(as.numeric(tapply(as.numeric(d$energy__), d$.chain, function(E) {
          if (length(E) < 2L || stats::var(E) == 0) NA_real_
          else sum(diff(E)^2) / length(E) / stats::var(E)
        })))
      } else NA_real_
      
      rows[[length(rows) + 1L]] <- data.frame(
        target             = sprintf("%s_mcmc_%s_%s",
                                     if (m == "SIR") "voi_sir" else "voi_hsir",
                                     stans[[m]], sid),
        model              = m,
        id                 = sid,
        regime             = grid$regime[i],
        max_rhat           = max(s$rhat,     na.rm = TRUE),
        min_ess_bulk       = min(s$ess_bulk, na.rm = TRUE),
        min_ess_tail       = min(s$ess_tail, na.rm = TRUE),
        divergences        = .pts_sum(dv), 
        max_treedepth_hits = .pts_sum(td >= 10),
        min_ebfmi          = eb,
        stringsAsFactors   = FALSE
      )
    }
  }
  if (!length(rows)) stop("No VOI bundles found in ", voi_dir, call. = FALSE)
  res <- dplyr::bind_rows(rows)
  
  meta <- targets::tar_meta(fields = "seconds")
  res$seconds <- meta$seconds[match(res$target, meta$name)]
  
  res$ok <- round(res$max_rhat, digits = 3) <= 1.02 &
    (is.na(res$divergences) | res$divergences <= 100) &
    (is.na(res$min_ebfmi)   | res$min_ebfmi >= 0.3)
  
  res[order(res$model, res$id), c("model", "id", "regime",
                                  "max_rhat", "min_ess_bulk", "min_ess_tail",
                                  "divergences", "max_treedepth_hits",
                                  "min_ebfmi", "seconds", "ok")]
}

#### REMOVE FROM HERE LATER ####
voi_posterior_draws <- function(cfg = load_params("config.yml"),
                                params = c("beta", "R0", "cv"),
                                voi_dir = "outputs/VOI") {
  grid <- voi_read_scenarios(cfg)
  rlev <- c("2GI", "1GI", "0GI", "full")
  stans <- c(SIR  = tools::file_path_sans_ext(basename(voi_stan_file(cfg, "SIR"))),
             HSIR = tools::file_path_sans_ext(basename(voi_stan_file(cfg, "HSIR"))))
  
  true_value <- function(p, beta, gamma, cv) switch(
    p, beta = beta, gamma = gamma, cv = cv,
    R0 = beta / gamma, D = 1 / gamma, NA_real_
  )
  
  out <- list()
  for (m in names(stans)) {
    for (i in seq_len(nrow(grid))) {
      sid  <- grid$id[i]
      path <- file.path(voi_dir, sprintf("voi_fit_%s_%s.rds", stans[[m]], sid))
      if (!file.exists(path)) next
      
      dr <- as.data.frame(readRDS(path)$draws)
      if (!"gamma" %in% names(dr) && "D" %in% names(dr)) dr$gamma <- 1 / dr$D
      if (!"R0" %in% names(dr) && all(c("beta", "D") %in% names(dr))) dr$R0 <- dr$beta * dr$D
      
      for (p in intersect(params, names(dr))) {
        out[[length(out) + 1L]] <- data.frame(
          model     = m,
          combo     = sub("_(2GI|1GI|0GI|full)$", "", sid),
          regime    = factor(grid$regime[i], levels = rlev),
          parameter = p,
          value     = dr[[p]],
          true      = true_value(p, grid$beta[i], grid$gamma[i], grid$cv[i]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(out)) stop("No VOI bundles found in ", voi_dir, call. = FALSE)
  res <- dplyr::bind_rows(out)
  res$model <- factor(res$model, levels = c("HSIR", "SIR"))
  res
}

voi_ridge_plot <- function(params = c("beta", "R0", "cv"),
                           draws = NULL, show_truth = TRUE,
                           cfg = load_params("config.yml")) {
  if (is.null(draws)) draws <- voi_posterior_draws(cfg, params)
  multi_combo <- length(unique(draws$combo)) > 1
  
  panel <- function(p) {
    d <- draws[draws$parameter == p, , drop = FALSE]
    g <- ggplot2::ggplot(d, ggplot2::aes(x = value, y = model, fill = model)) +
      ggridges::geom_density_ridges(alpha = 0.7, scale = 1.0,
                                    rel_min_height = 0.01, colour = "white") +
      (if (multi_combo) ggplot2::facet_grid(combo ~ regime)
       else             ggplot2::facet_wrap(~ regime, nrow = 1)) +
      ggplot2::labs(x = p, y = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "none")
    if (show_truth) {
      tv <- unique(d[, c("model", "combo", "regime", "true")])
      g <- g + ggplot2::geom_point(
        data = tv, inherit.aes = FALSE,
        ggplot2::aes(x = true, y = model), shape = 124, size = 3
      )
    }
    g
  }
  
  panels <- lapply(params, panel)
  if (length(panels) == 1) panels[[1]] else patchwork::wrap_plots(panels, ncol = 1)
}
