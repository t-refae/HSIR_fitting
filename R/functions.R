# R/functions.R
# helpers for H/SIR fitting targets pipeline

parse_fit_regime <- function(regime) {
  r <- tolower(trimws(as.character(regime)))
  if (r == "full") return(list(type = "full", k = NA_integer_))
  if (r == "peak") return(list(type = "gi",   k = 0L))
  m <- regmatches(r, regexec("^([0-9]+)\\s*gi$", r))[[1]]
  if (length(m) == 2L) return(list(type = "gi", k = as.integer(m[[2]])))
  stop("fit_regime must be 'full', 'peak', or '<k>GI' (e.g. '0GI','1GI','2GI'); got '",
       regime, "'.")
}


load_params <- function(config_file) {
  p <- yaml::read_yaml(config_file)
  p$model_type <- toupper(trimws(as.character(p$model_type)))
  if (!p$model_type %in% c("HSIR", "SIR"))
    stop("model_type must be 'HSIR' or 'SIR' (got '", p$model_type, "').")
  parse_fit_regime(p$fit_regime)
  p
}

select_stan_file <- function(params) {
  switch(params$model_type,
         HSIR = "Stan/HSIR.stan",
         SIR  = "Stan/SIR.stan"
  )
}

.epi_rhs <- function(t, state, parms) {
  S <- state[["S"]]; I <- state[["I"]]
  foi <- parms$beta * I * S^(1 + (parms$cv)^2)
  list(c(S = -foi,
         I =  foi - parms$gamma * I,
         R =  parms$gamma * I,
         C =  foi))
}


# simulated incidence with added Poisson noise
simulate_epidemic <- function(params) {
  cv_sim <- if (params$model_type == "SIR") 0 else params$cv
  parms  <- list(beta = params$beta, gamma = params$gamma, cv = cv_sim)
  state0 <- c(S = 1 - params$i0, I = params$i0, R = 0, C = 0)  # initial state
  

  full_out <- as.data.frame(deSolve::ode(
    y = state0, times = seq_len(params$max_days), func = .epi_rhs, parms = parms
  ))
  inc_full    <- diff(full_out$C)          # incidence days 2 .. max_days
  counts_full <- inc_full * params$P       # expected new cases / day
  
  # calculate peak deterministic incidence (index i -> calendar day i + 1)
  peak_index <- which.max(counts_full)
  
  # dynamic horizon: first POST-peak day with counts < end_threshold
  tail_below <- which(counts_full[peak_index:length(counts_full)] < params$end_threshold)
  if (length(tail_below) == 0L) {
    warning("Epidemic did not fall below end_threshold within max_days; using max_days.")
    end_index <- length(counts_full)
  } else {
    end_index <- peak_index - 1L + tail_below[1]
  }
  n_days <- as.integer(end_index + 1L)     # incidence index end_index -> day n_days
  
  # truncate and observe
  out <- full_out[seq_len(n_days), , drop = FALSE]
  inc <- inc_full[seq_len(n_days - 1L)]
  set.seed(params$seed)
  cases <- rpois(length(inc), inc * params$P)
  
  stopifnot(n_days >= 2L, length(cases) == n_days - 1L)
  
  message(sprintf("[simulate] dynamic n_days = %d (peak ~ day %d, GI = %.2f d, cv = %g)",
                  n_days, peak_index + 1L, 1 / params$gamma, cv_sim))
  
  list(
    states     = out,
    cases      = as.integer(cases),        # length (n_days - 1)
    incidence  = inc,
    n_days     = n_days,
    peak_index = as.integer(peak_index),   # index into incidence (1 .. n_days-1)
    GI_days    = 1 / params$gamma,         # generation interval = 1/gamma
    cv_used    = cv_sim
  )
}


build_stan_data <- function(sim, params) {
  n_days <- sim$n_days
  M      <- n_days - 1L                    # number of incidence observations
  fr     <- parse_fit_regime(params$fit_regime)
  
  if (fr$type == "full") {
    n_fit <- M
    message(sprintf("[fit window] full | n_fit = %d / %d days | n_days = %d",
                    n_fit, M, n_days))
  } else {
    gi    <- sim$GI_days
    i_cut <- sim$peak_index - fr$k * gi    # cutoff in incidence-index units
    n_fit <- as.integer(min(M, max(1L, floor(i_cut))))
    if (i_cut < 1)
      warning(sprintf("Peak (incidence day %d) is < %d GI (%.1f d) from the start; clamping n_fit to 1.",
                      sim$peak_index, fr$k, fr$k * gi))
    message(sprintf("[fit window] %dGI pre-peak | peak ~ day %d | GI = %.2f d | n_fit = %d / %d days | n_days = %d",
                    fr$k, sim$peak_index + 1L, gi, n_fit, M, n_days))
  }
  
  days <- seq_len(n_days)                  # 1, 2, ..., n_days
  ts   <- as.numeric(days[-1])             # 2, ..., n_days   (length n_days - 1)
  
  stopifnot(
    n_days >= 2L,
    length(ts)        == n_days - 1L,
    length(sim$cases) == n_days - 1L,
    n_fit >= 1L, n_fit <= n_days - 1L
  )
  
  list(
    n_days = n_days,
    y0     = as.numeric(c(1 - params$i0, params$i0, 0, 0)),  # [S,I,R,C] at day 1
    t0     = as.numeric(days[1]),                            # = 1 (day 1)
    ts     = ts,                                             # days 2 .. n_days
    N      = as.integer(params$P),
    cases  = sim$cases,                                      # length n_days - 1
    n_fit  = n_fit
  )
}


# per-chain initial values, jittered around prior medians

make_inits <- function(model_type, chains, seed = 1L) {
  set.seed(seed)
  lapply(seq_len(chains), function(ch) {
    init <- list(
      beta = runif(1, 0.6, 1.4),
      D    = runif(1, 3.5, 4.5)
    )
    if (toupper(trimws(model_type)) == "HSIR")
      init$cv <- runif(1, 1e-3, 0.25)
    init
  })
}
