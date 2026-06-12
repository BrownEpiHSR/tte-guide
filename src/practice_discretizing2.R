library(here)
library(tidyverse)
library(survival)
library(survminer)
library(broom)

# Function ----
make_pp <- function(dat_cont, rule = c("ceiling", "floor", "round"), t_cuts = 30) {
  
  rule <- match.arg(rule)
  
  # Pick the discretized time variable
  rule_var <- paste0(rule, "_t")
  
  # Build interval table
  max_t <- max(dat_cont[[rule_var]])
  cuts <- seq(0, max_t, by = t_cuts)
  
  interval_df <- tibble(
    interval = 1:(length(cuts) - 1),
    t_start = cuts[-length(cuts)],
    t_end   = cuts[-1]
  )
  
  # Merge rule-based event time
  dat <- dat_cont %>%
    mutate(event_time_disc = .data[[rule_var]])
  
  # Expand to person-period
  dat_pp <- dat %>%
    crossing(interval_df) %>%
    mutate(
      # at risk in [t_start, t_end) if follow-up time extends beyond t_start
      at_risk = time > t_start,
      # event occurs in the (single) discrete interval containing event_time_disc
      y = as.integer(outcome == 1 &
                       event_time_disc > t_start &
                       event_time_disc <= t_end)
    ) %>%
    # only keep intervals where subject is actually at risk
    filter(at_risk) %>%
    arrange(id, interval)
  
  return(dat_pp)
}

estimate_disc_hazard <- function(dat_pp) {
  dat_pp %>%
    group_by(interval, t_start, t_end) %>%
    summarize(
      events = sum(y),
      risk   = sum(at_risk),
      h      = events / risk,
      .groups = "drop"
    ) %>%
    mutate(H = cumsum(h))
}

estimate_cont_hazard <- function(dat_cont) {
  fit <- survfit(Surv(time, outcome) ~ 1, type = "fh")
  
  tibble(
    time = fit$time,
    dN   = fit$n.event,
    Y    = fit$n.risk,
    h    = dN / Y,
    H    = cumsum(h)
  )
}

compute_mae <- function(haz_disc, haz_cont) {
  haz_disc %>%
    select(t_end, H_disc = H) %>%
    left_join(
      haz_cont %>% select(time, H_cont = H),
      by = c("t_end" = "time")
    ) %>%
    drop_na() %>%
    summarize(mae = mean(abs(H_disc - H_cont))) %>%
    pull(mae)
}

compare_rules <- function(dat_cont, t_cuts = 30) {
  
  haz_cont <- estimate_cont_hazard(dat_cont)
  
  map_df(c("ceiling", "floor", "round"), function(rule) {
    
    dat_pp <- make_pp(dat_cont, rule = rule, t_cuts = t_cuts)
    haz_disc <- estimate_disc_hazard(dat_pp)
    mae <- compute_mae(haz_disc, haz_cont)
    
    tibble(rule = rule, mae = mae)
  })
}

# Synthetic data ----
  set.seed(2026)
  
  n <- 1000
  
  # Weibull parameters
  shape <- 1.5      # hazard increases over time
  scale <- 80       # roughly median ~ 80 * (log 2)^(1/shape)
  
  # Generate true event times
  t_event <- ceiling(rweibull(n, shape = shape, scale = scale))
  
  # Random censoring distribution
  # Choose rate to get ~30–40% censoring
  t_cens <- ceiling(rexp(n, rate = 1/120))
  
  # Administrative censoring at 180 days
  t_admin <- 180
  t_cuts <- 30
  
  # Observed time and event indicator
  time <- pmin(t_event, t_cens, t_admin)
  outcome <- as.integer(t_event <= t_cens & t_event <= t_admin)
  censor <- as.integer(t_cens <= t_event & t_cens <= t_admin)
    
  dat_cont <- tibble(
    id = 1:n,
    t_event,
    t_cens,
    time,
    outcome,
    censor
  ) %>%
  mutate(ceiling_t = ceiling(time/t_cuts)*t_cuts,
         floor_t = if_else(outcome==1, ceiling(time/t_cuts)*t_cuts, floor(time/t_cuts)*t_cuts),
         round_t = if_else(outcome==1, ceiling(time/t_cuts)*t_cuts, round(time/t_cuts)*t_cuts)
         )
  
  write.csv(dat_cont, here('dta', 'syn_survdta.csv'))

# Convert to person-period ----
  compare_rules(dat_cont, t_cuts = 30)
  