#### full parameter grid to explore ####

rm(list=ls())

library(dplyr)
library(deSolve)
library(ggplot2)

## parameter ranges
# beta ~ U(0, 2)
# D = 1/gamma ~ U(2, 6)
# cv ~ U(0, 4)

beta_range <- data.frame(beta = c(seq(0.3, 0.6, by=0.15), 
                                  seq(0.6, 1.2, by=0.3), 
                                  seq(0.9, 1.8, by=0.45))
              )

gamma_range <- data.frame(gamma = c(0.2, 0.3, 0.4))

R0_range <- c(1.5, 3, 4.5)

beta_gamma_join <- cross_join(beta_range, gamma_range) %>%
  mutate(beta = round(beta, digits=2),
         gamma= round(gamma, digits=2),
         R0 = round(beta / gamma, digits=2)) %>%
  filter(R0 %in% R0_range) %>%
  arrange(R0, beta, gamma) %>%
  distinct()

cv_range <- data.frame(cv = c(0, 0.5, 1, 2))

# full range of parameter combinations explored
bgcv_join <- cross_join(beta_gamma_join, cv_range)

bgcv_join$param_combo_label <- NA

bgcv_join <- bgcv_join %>%
  mutate(
    param_combo_label = case_when(
      (R0 == 1.5 & beta == 0.30) | (R0 == 3 & beta == 0.6) | (R0 == 4.5 & beta == 0.90) ~ "slow",
      (R0 == 1.5 & beta == 0.45) | (R0 == 3 & beta == 0.9) | (R0 == 4.5 & beta == 1.35) ~ "reference",
      TRUE ~ "fast"
    ),
    id = seq_len(nrow(bgcv_join))
  ) %>%
  relocate(param_combo_label, .before = beta) %>%
  relocate(id, .before = param_combo_label)

write.csv(bgcv_join, file = "Data/parameter_grid.csv", row.names = FALSE)

# different types of fitting to perform

# to quantify wrt n_runs: VOI? NPI (NOT NPI_release) ?
fitting_types <- c("full", "partial", "NPI_release")
n_fitting_types <- c(1, 3, 4)

n_to_runs <- nrow(bgcv_join)*sum(n_fitting_types)

avg_length_in_hours <- 5

n_cores_avail <- 1

n_to_runs * avg_length_in_hours / (24 * n_cores_avail) # expected number of days to run

#### visualiation ####

# comparing I dynamics across a grid of param combos


HSIR_out <- function(beta, gamma, cv, P = 1e6, tmax = 800, dt = 1) {
  I0 <- 1/P  
  derivs <- function(t, y, p) {
    S <- y[["S"]]
    I <- y[["I"]]
    R <- y[["R"]] 
    C <- y[["C"]]
    
    dS <- -beta*I*S^(1+cv^2)
    dI <-  beta*I*S^(1+cv^2) - gamma*I
    dR <-  gamma*I
    dC <-  beta*I*S^(1+cv^2)
    
    list(c(dS, dI, dR, dC))
  }
  out <- ode(y     = c(S = 1 - I0, I = I0, R = 0, C = 0),
             times = seq(0, tmax, by = dt),
             func  = derivs)
  as.data.frame(out)
}

## run one combination, then trim the dead tail
## Keep time up to when 99.5% of the final epidemic size is reached, so
## each panel's x-axis frames the active epidemic instead of a long flat
## line.  (Works with the per-column free x-scale below.)
run_one <- function(p) {
  d   <- HSIR_out(p$beta, p$gamma, p$cv)
  C   <- d$S[1] - d$S                         # cumulative incidence
  idx <- which(C >= 0.995 * C[length(C)])[1] + 14
  tcut <- if (is.na(idx)) max(d$time) else d$time[idx]
  d   <- d[d$time <= tcut, c("time", "I")]
  data.frame(d, beta = p$beta, gamma = p$gamma, R0 = p$R0, cv = p$cv)
}

traj <- do.call(rbind, lapply(seq_len(nrow(bgcv_join)),
                              function(i) run_one(bgcv_join[i, ])))
traj$cv <- factor(traj$cv, levels = sort(unique(traj$cv)))

## Rows = R0  (free y  -> each R0 gets a y-range matched to its peak, so
##             the cv effect is legible in every panel)
## Cols = gamma (free x -> each gamma gets its own timescale)
## Colour = cv (the comparison of interest, easy within each panel)

p <- ggplot(traj, aes(time, I, colour = cv)) +
  geom_line(linewidth = 0.8) +
  facet_grid(R0 ~ gamma, scales = "free",
             labeller = label_bquote(rows = R[0] == .(R0),
                                     cols = gamma == .(gamma))) +
  scale_colour_viridis_d(name = expression(nu), end = 0.9) +
  labs(
    title    = "Infectious-compartment dynamics across parameter combinations",
    subtitle = expression("rows: "*R[0]*"      columns: "*gamma*
                            "      ("*beta == R[0]*gamma*")      colour: "* nu),
    x = "Time (Days)",
    y = expression(I(t))
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background  = element_rect(fill = "grey92"),
        legend.position   = "right",
        axis.title.y = element_text(angle=0, vjust=0.5))

print(p)

ggsave("I_dynamics.png", p, width = 11, height = 8, dpi = 300)

