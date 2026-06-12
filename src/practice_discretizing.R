library(here)
library(tidyverse)
library(survival)
library(survminer)
library(broom)

# Function ----
assign_floor <- function(time, interval_df) {
  if (time == max(interval_df$t_end)) {
    return(max(interval_df$interval))
  }
  interval_df %>%
    filter(time >= t_start, time < t_end) %>%
    pull(interval)
}

assign_ceiling <- function(time, interval_df) {
  if (time == min(interval_df$t_start)) {
    return(min(interval_df$interval))
  }
  interval_df %>%
    filter(time > t_start, time <= t_end) %>%
    pull(interval)
}

assign_nearest <- function(time, interval_df) {
  interval_df %>%
    mutate(dist = abs(time - midpoint)) %>%
    pull(interval) %>%
    sample(., size=1L, replace=F) 
}

make_discrete <- function(dat_cont, cuts, rule = c("floor", "ceiling", "nearest")) {
  
  rule <- match.arg(rule)
  
  # Build interval lookup table
  interval_df <- tibble(
    interval = 1:(length(cuts) - 1),
    t_start = cuts[-length(cuts)],
    t_end   = cuts[-1],
    midpoint = (t_start + t_end) / 2
  )
  
  # Choose the assignment function
  assign_fun <- switch(
    rule,
    floor   = function(t) assign_floor(t, interval_df),
    ceiling = function(t) assign_ceiling(t, interval_df),
    nearest = function(t) assign_nearest(t, interval_df)
  )
  
  # Apply rule to get event interval for each subject
  dat_with_rule <- dat_cont %>%
    mutate(event_interval = map_int(time, assign_fun))
  
  # Build person-period dataset
  dat_pp <- dat_with_rule %>%
    crossing(interval = interval_df$interval) %>%
    left_join(interval_df, by = "interval") %>%
    mutate(
      y = as.integer(event == 1 & interval == event_interval)
    ) %>%
    group_by(id) %>%
    filter(interval<=event_interval) %>%
    ungroup() %>%
    arrange(id, interval)
  
  return(dat_pp)
}

make_disc_haz =function(x) {
  x %>%
  group_by(interval, t_start, t_end) %>%
  summarize(
    events = sum(y),
    risk   = n(),
    h      = events / risk,
    .groups = "drop"
  ) %>%
  mutate(
    H = cumsum(h)
  )
}

est_logitfun = function(x) {
  fit_logit <- glm(
    y ~ factor(interval), 
    data = x,
    family = binomial()
  )
  
  tidy(fit_logit)
  
  tibble(interval = unique(x$interval)) %>%
    mutate(
      est_h     = plogis(predict(fit_logit, newdata = ., type = "link")),
      est_H     = cumsum(est_h)
    )
}

summ_disc_haz = function(x) {
    inner_join(make_disc_haz(x), 
               est_logitfun(x), 
               by = join_by(interval)
               )
}

# Synthetic data ----
  set.seed(2026)
  
  n <- 1000
  
  # Weibull parameters
  shape <- 1.5      # hazard increases over time
  scale <- 80       # roughly median ~ 80 * (log 2)^(1/shape)
  
  # Generate true event times
  t_event <- round(rweibull(n, shape = shape, scale = scale),0)
  
  # Random censoring distribution
  # Choose rate to get ~30–40% censoring
  t_cens <- round(rexp(n, rate = 1/120), 0)
  
  # Administrative censoring at 180 days
  t_admin <- 180
  
  # Observed time and event indicator
  time <- pmin(t_event, t_cens, t_admin)
  event <- as.integer(t_event <= t_cens & t_event <= t_admin)
  
  dat_cont <- tibble(
    id = 1:n,
    t_event,
    t_cens,
    time,
    event
  )
  
  write.csv(dat_cont, here('dta', 'syn_survdta.csv'))

# Estimate hazard ----
  fit = survfit(Surv(time, event) ~ 1, type='fh')
  
  dat_haz = tibble(
    time = fit$time,
    dN = fit$n.event,
    Y = fit$n.risk,
    h_t = dN / Y,
    H_t = cumsum(h_t)
  )

## Discretized data ----
  # Define interval cutpoints
  cuts <- seq(0, 180, by = 30)
  
  # 3 different rules for counting
  dat_floor <- make_discrete(dat_cont, cuts, rule = "floor") 
  dat_ceiling  <- make_discrete(dat_cont, cuts, rule = "ceiling")  
  dat_nearest  <- make_discrete(dat_cont, cuts, rule = "nearest") 
  
  disc_haz_floor = summ_disc_haz(dat_floor)
  disc_haz_ceiling = summ_disc_haz(dat_ceiling)
  disc_haz_near = summ_disc_haz(dat_nearest)
  
# Parametric function ----
  
  pred_df
  