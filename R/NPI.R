# R/NPI.R

.npi_or <- function(x, default) if (is.null(x)) default else x

npi_stan_file <- function(cfg) .npi_or(cfg$npi_stan, "Stan/NPI_HSIR.stan")

npi_mcmc_settings <- function(cfg) {
  list(
    chains          = .npi_or(cfg$npi_chains,          cfg$chains),
    parallel_chains = .npi_or(cfg$npi_parallel_chains, cfg$parallel_chains),
    iter_warmup     = .npi_or(cfg$npi_iter_warmup,     cfg$iter_warmup),
    iter_sampling   = .npi_or(cfg$npi_iter_sampling,   cfg$iter_sampling),
    refresh         = .npi_or(cfg$npi_refresh,         cfg$refresh),
    seed            = .npi_or(cfg$npi_seed,            cfg$seed),
    adapt_delta     = .npi_or(cfg$npi_adapt_delta,     0.9),
    max_treedepth   = .npi_or(cfg$npi_max_treedepth,   12L),
    metric          = .npi_or(cfg$npi_metric,          "diag_e")
  )
}

npi_read_scenarios <- function(cfg) {
  path <- .npi_or(cfg$npi_scenarios_csv, "Data/npi_scenarios.csv")
  df <- utils::read.csv(path, stringsAsFactors = FALSE, strip.white = TRUE)
  names(df) <- tolower(trimws(names(df)))
  need <- c("theta", "beta", "gamma", "cv", "eff", "n_intervals")
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("npi_scenarios_csv missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  df$id <- if (!is.null(df$id)) .id_safe(df$id) else sprintf("npi_theta%d", df$theta)
  if (anyDuplicated(df$id)) stop("NPI scenarios: `id` values are not unique.", call. = FALSE)
  if (is.null(df$npi_start)) df$npi_start <- NA_real_
  df[order(df$theta), , drop = FALSE]
}

.npi_rhs <- function(t, state, parms) {
  S <- state[["S"]]; I <- state[["I"]]
  b <- if (t > parms$npi_start) parms$beta * (1 - parms$eff) else parms$beta
  foi <- b * I * S^(1 + (parms$cv)^2)
  list(c(S = -foi,
         I =  foi - parms$gamma * I,
         R =  parms$gamma * I,
         C =  foi))
}

.npi_solve <- function(beta, gamma, cv, eff, npi_start, i0, n_days) {
  y0 <- c(S = 1 - i0, I = i0, R = 0, C = 0)
  parms <- list(beta = beta, gamma = gamma, cv = cv, eff = eff, npi_start = npi_start)
  sol <- deSolve::ode(y = y0, times = 0:(n_days - 1), func = .npi_rhs,
                      parms = parms, method = "lsoda")
  list(y0 = as.numeric(y0), sol = as.data.frame(sol))
}

npi_peak_day <- function(cfg, beta, gamma, cv) {
  n_days <- .npi_or(cfg$npi_n_days, 60)
  i0     <- .npi_or(cfg$npi_i0, cfg$i0)
  s <- .npi_solve(beta, gamma, cv, eff = 0, npi_start = n_days + 1,
                  i0 = i0, n_days = n_days)
  which.max(diff(s$sol$C))
}

npi_start_day <- function(cfg, beta, gamma, cv, override = NA_real_) {
  if (!is.na(override)) return(as.numeric(override))
  gi <- 1 / gamma
  round(npi_peak_day(cfg, beta, gamma, cv) - .npi_or(cfg$npi_gi_before_peak, 1) * gi)
}

npi_build_stan_data <- function(cfg, row) {
  n_days <- .npi_or(cfg$npi_n_days, 60)
  P      <- .npi_or(cfg$npi_P,  cfg$P)
  i0     <- .npi_or(cfg$npi_i0, cfg$i0)
  seed   <- .npi_or(cfg$npi_seed, cfg$seed)
  
  start <- npi_start_day(cfg, row$beta, row$gamma, row$cv, row$npi_start)
  i_start <- as.integer(start)
  if (i_start < 1L || i_start > n_days - 2L) {
    stop("NPI start day ", start, " outside 1..", n_days - 2,
         " for scenario ", row$id, call. = FALSE)
  }
  
  s <- .npi_solve(row$beta, row$gamma, row$cv, row$eff, start, i0, n_days)
  incidence <- pmax(diff(s$sol$C), 1e-12)
  
  set.seed(seed)
  cases <- stats::rpois(length(incidence), incidence * P)
  
  gi     <- 1 / row$gamma
  n_fit  <- as.integer(floor(start + row$n_intervals * gi))
  n_fit  <- max(2L, min(n_fit, n_days - 1L))
  
  list(
    n_days  = n_days,
    y0      = s$y0,
    t0      = 0,
    ts      = as.array(as.numeric(seq_len(n_days - 1))),
    N       = P,
    cases   = as.array(as.integer(cases)),
    n_fit   = n_fit,
    i_start = i_start,
    eff     = row$eff
  )
}

npi_summarise_convergence <- function(cfg = load_params("config.yml"),
                                      npi_dir = "outputs/NPI") {
  grid  <- npi_read_scenarios(cfg)
  model <- tools::file_path_sans_ext(basename(npi_stan_file(cfg)))
  
  rows <- list()
  for (i in seq_len(nrow(grid))) {
    sid  <- grid$id[i]
    path <- file.path(npi_dir, sprintf("npi_fit_%s_%s.rds", model, sid))
    if (!file.exists(path)) next
    
    b  <- readRDS(path)
    dr <- posterior::as_draws_df(as.data.frame(b$draws))
    keep <- intersect(c("beta", "D", "cv", "R0", "gamma"), posterior::variables(dr))
    s <- posterior::summarise_draws(
      if (length(keep)) posterior::subset_draws(dr, variable = keep) else dr,
      posterior::default_convergence_measures()
    )
    
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
      target             = sprintf("npi_fit_mcmc_%s_%s", model, sid),
      theta              = grid$theta[i],
      id                 = sid,
      eff                = grid$eff[i],
      n_intervals        = grid$n_intervals[i],
      max_rhat           = max(s$rhat,     na.rm = TRUE),
      min_ess_bulk       = min(s$ess_bulk, na.rm = TRUE),
      min_ess_tail       = min(s$ess_tail, na.rm = TRUE),
      divergences        = .pts_sum(dv),
      max_treedepth_hits = .pts_sum(td >= 10),
      min_ebfmi          = eb,
      stringsAsFactors   = FALSE
    )
  }
  if (!length(rows)) stop("No NPI bundles found in ", npi_dir, call. = FALSE)
  res <- dplyr::bind_rows(rows)
  
  meta <- targets::tar_meta(fields = "seconds")
  res$seconds <- meta$seconds[match(res$target, meta$name)]
  
  res$ok <- round(res$max_rhat, digits = 3) <= 1.02 &
    (is.na(res$divergences) | res$divergences <= 100) &
    (is.na(res$min_ebfmi)   | res$min_ebfmi >= 0.3)
  
  res[order(res$theta), c("theta", "id", "eff", "n_intervals", "max_rhat",
                          "min_ess_bulk", "min_ess_tail", "divergences",
                          "max_treedepth_hits", "min_ebfmi", "seconds", "ok")]
}

npi_posterior_draws <- function(cfg = load_params("config.yml"),
                                params = c("beta", "cv", "R0"),
                                npi_dir = "outputs/NPI") {
  grid  <- npi_read_scenarios(cfg)
  model <- tools::file_path_sans_ext(basename(npi_stan_file(cfg)))
  
  true_value <- function(p, beta, gamma, cv) switch(
    p, beta = beta, gamma = gamma, cv = cv,
    R0 = beta / gamma, D = 1 / gamma, NA_real_
  )
  
  out <- list()
  for (i in seq_len(nrow(grid))) {
    sid  <- grid$id[i]
    path <- file.path(npi_dir, sprintf("npi_fit_%s_%s.rds", model, sid))
    if (!file.exists(path)) next
    
    dr <- as.data.frame(readRDS(path)$draws)
    if (!"R0"    %in% names(dr) && all(c("beta", "D") %in% names(dr))) dr$R0    <- dr$beta * dr$D
    if (!"gamma" %in% names(dr) && "D" %in% names(dr))                 dr$gamma <- 1 / dr$D
    
    for (p in intersect(params, names(dr))) {
      out[[length(out) + 1L]] <- data.frame(
        id        = sid,
        theta     = grid$theta[i],
        eff       = grid$eff[i],
        GI        = grid$n_intervals[i],
        parameter = p,
        value     = dr[[p]],
        true      = true_value(p, grid$beta[i], grid$gamma[i], grid$cv[i]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out)) stop("No NPI bundles found in ", npi_dir, call. = FALSE)
  dplyr::bind_rows(out)
}

npi_ridge_plot <- function(params = c("beta", "cv", "R0"),
                           draws = NULL, show_truth = TRUE,
                           cfg = load_params("config.yml")) {
  if (is.null(draws)) draws <- npi_posterior_draws(cfg, params)
  draws$GI_f  <- factor(draws$GI, levels = rev(sort(unique(draws$GI))),
                        labels = paste0("+", rev(sort(unique(draws$GI))), " GI"))
  draws$eff_f <- factor(draws$eff)
  
  panel <- function(p) {
    d <- draws[draws$parameter == p, , drop = FALSE]
    g <- ggplot2::ggplot(d, ggplot2::aes(x = value, y = GI_f, fill = eff_f)) +
      ggridges::geom_density_ridges(alpha = 0.6, scale = 1.1,
                                    rel_min_height = 0.01, panel_scaling = FALSE) +
      ggplot2::facet_wrap(~ eff_f, ncol = 1, labeller = ggplot2::label_both) +
      ggplot2::labs(x = p, y = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "none")
    if (show_truth) {
      tv <- unique(d[, c("eff_f", "GI_f", "true")])
      g <- g + ggplot2::geom_point(
        data = tv, inherit.aes = FALSE,
        ggplot2::aes(x = true, y = GI_f), shape = 124, size = 3
      )
    }
    g
  }
  
  panels <- lapply(params, panel)
  if (length(panels) == 1) panels[[1]] else patchwork::wrap_plots(panels, ncol = length(params))
}
