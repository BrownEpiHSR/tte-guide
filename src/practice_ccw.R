
# Setup ----
  library(tidyverse)
  library(survival)
  library(broom)
  library(future)
  library(furrr)
  library(data.table)

## CCW simulation; 
# 1. generates data
# 2. Estimates weights
# 3. Models outcome

run_sim = function(grace=90, disc_cuts=30, n=1000, t_cens = F, t_admin=360, stabilized=F) {

# Synthetic data
  # L affects BOTH event time and treatment time
    L <- rbinom(n, 1, 0.5)   # 50/50 binary confounder

  # People with L=1 have higher hazard (shorter survival)
    t_event <- if_else(
      L == 1,
      rweibull(n, shape = 1.1, scale = 300), # short survival
      rweibull(n, shape = 1.1, scale = 600)) |> # long survival
      ceiling()

  # People with L=1 get treated earlier (confounding)
     t_treat <- if_else(
      L == 0,
      rweibull(n, shape = 1.0, scale = 180), 
      rweibull(n, shape = 1.0, scale = 60)) |> # early initiator if L=1
      ceiling()
     
     if (t_cens) {  t_cens <- rweibull(n, shape = 1.0, scale = 900) |>
       ceiling() } else t_cens = Inf

  # observed time
    time  <- pmin(t_event, t_cens, t_admin)
    # event if time prior to censoring or admin end
    # outcome > t_admin > natural censoring (death as last)
    event <- as.integer(t_event == time)
  
    dat <- tibble(
      id = 1:n,
      L,
      t_event,
      t_treat,
      t_cens,
      time,
      event
    )

# cloning
  dat_clone <- dat %>%
    #duplicate data
    crossing(strategy = c("treat_by_grace", "never_treat")) %>%
    # assign treatment rules
    mutate(
      t_art = case_when(
        strategy == "treat_by_grace" & t_treat > grace  ~ grace,
        strategy == "treat_by_grace" & t_treat <= grace ~ Inf,
        strategy == "never_treat"  & is.finite(t_treat) ~ t_treat,
        strategy == "never_treat"  & !is.finite(t_treat) ~ Inf
      ),
      # adjust time/event for artificial censoring
      t_obs = pmin(time, t_art),
      c_art = t_art==t_obs,
      event_obs = as.integer(event == 1 & t_event==t_obs)
    )

  # --- KM-based estimators (unexpanded data) -----------------------------
  
  dat_clone_km <- dat_clone
  
  # Fit censoring model for IPCW (natural + artificial censoring)
  mod_km_ipcw <- glm(
    c_art ~ L * strategy,
    data = dat_clone_km,
    family = binomial()
  )
  
  # IPCW weights: 1 / P(C=0 | L, strategy)
  dat_clone_km$ipcw <- 1 / (1 - mod_km_ipcw$fitted.values)
  
  # Unweighted KM
  km_u <- survfit(
    Surv(t_obs, event_obs) ~ strategy,
    data = dat_clone_km
  )
  
  # Weighted KM
  km_w <- survfit(
    Surv(t_obs, event_obs) ~ strategy,
    data = dat_clone_km,
    weights = ipcw
  )
  
  # Extract cumulative incidence at each unique event time
  km_u_df <- broom::tidy(km_u) %>%
    filter(strata %in% c("strategy=treat_by_grace", "strategy=never_treat")) %>%
    mutate(risk = 1 - estimate) %>%
    select(time, strata, risk) 
  
  km_w_df <- broom::tidy(km_w) %>%
    filter(strata %in% c("strategy=treat_by_grace", "strategy=never_treat")) %>%
    mutate(risk = 1 - estimate) %>%
    select(time, strata, risk) 
    
  # Reshape wide
  km_u_wide <- km_u_df %>%
    mutate(strata = recode(strata,
                           "strategy=treat_by_grace" = "treat",
                           "strategy=never_treat" = "never")) %>%
    pivot_wider(names_from = strata, values_from = risk) %>%
    mutate(RD_km_u = treat - never)
  
  km_w_wide <- km_w_df %>%
    mutate(strata = recode(strata,
                           "strategy=treat_by_grace" = "treat",
                           "strategy=never_treat" = "never")) %>%
    pivot_wider(names_from = strata, values_from = risk) %>%
    mutate(RD_km_w = treat - never)
  
# Discretize data
  # Define intervals `disc_cuts`
    cuts <- seq(0, t_admin, by = disc_cuts)
    
    interval_df <- tibble(
      interval = 1:(length(cuts)-1),
      t_start = cuts[-length(cuts)],
      t_end   = cuts[-1]
    )

  # construct person-period dataset
    dat_pp <- dat_clone %>%
      crossing(interval_df) %>%
      mutate(
        # Clone-specific risk
        at_risk = t_obs > t_start,
        
        # Event in this interval
        y = as.integer(
          event_obs == 1 &
            t_obs > t_start &
            t_obs <= t_end
        ),
        
        # Censoring in this interval (natural or artificial)
        cens = as.integer(
          event_obs==0 & c_art==1 &
            t_obs > t_start &
            t_obs <= t_end
        )
      ) %>%
      # KEEP the censoring interval — do NOT drop it
      filter(at_risk | cens == 1) %>%
      arrange(id, strategy, interval) %>%
      # note if observation already treated in time < interval start
      mutate(treated = as.integer(t_treat <= t_start))

# Estimate censoring weights 
    p_cens <- rep(0, nrow(dat_pp))
    
    # Flag for relevant person-periods
    flg_1 = which(
      dat_pp$strategy == 'treat_by_grace' & 
        dat_pp$t_end == grace & 
        dat_pp$treated==0)
    
    flg_0 = which(dat_pp$strategy == 'never_treat')
    
    fit_cens1 <- glm(
      cens ~ L,
      family = binomial(),
      data = dat_pp[flg_1,]
    )
    
    fit_cens0 <- glm(
      cens ~ factor(interval)*L,
      family = binomial(),
      data = dat_pp[flg_0,])
    
    p_cens[flg_1] <- predict(fit_cens1, type = "response")
    p_cens[flg_0] <- predict(fit_cens0, type = "response")
    
    # Switch to data.table for speed
    setDT(dat_pp)
    
    dat_pp[, `:=`(
      p_cens = p_cens,
      p_nocens = 1 - p_cens # Pr(C=0)
    )]
    
    # Pr(C=0 | A=a, Time=t) - numerator for stabilized weights
    dat_pp[, p_num := (1 - mean(cens)), by = .(strategy, interval)]
    
    # weights
    dat_pp[, `:=`(w = (1 / cumprod(p_nocens)),
                 sw  = cumprod(p_num) / cumprod(p_nocens)), 
          by=.(strategy, id)]
    
    # apply weights for following period
    # in the grace period we are counting up until the end
    # so the weight applies to folks in interval +1 (after censored after dropped)
    dat_pp[order(id, strategy, interval),
      `:=`(
        w = shift(w, type='lag', fill = 1L),
        sw = shift(sw, type = 'lag', fill = 1L)
      ), by =.(id, strategy)]
    
    setDF(dat_pp)
  
# Estimate outcome models 
    # glm is slow, switch to parglm() with parallel computing for speed
    # using time + time^2; but splines or cut-points are better
    # Unweighted for comparison
    fit_u <- glm(
      y ~ strategy*factor(interval),
      family = binomial(),
      data = dat_pp
    )
    
    # unstabilized weights
    fit_w <- glm(
      y ~ strategy*factor(interval),
      family = quasibinomial(), #prevents warning about non-integer
      data = dat_pp,
      weights = w
    )
    
    #stabilized weights
    fit_sw <- glm(
      y ~ strategy*factor(interval),
      family = quasibinomial(), #prevents warning about non-integer
      data = dat_pp,
      weights = sw
    )

  # Conditional Y | Interval, strategy, IPCW(L) `weighted`
  # Conditional Y | Interval, strategy `unweighted`
  dat_pp$h_u = predict(fit_u, newdata = dat_pp, type = "response")
  dat_pp$h_w = predict(fit_w, newdata = dat_pp, type = "response")
  dat_pp$h_sw = predict(fit_sw, newdata = dat_pp, type = "response")
  
  dat_summ <- dat_pp %>%
    group_by(id, strategy) %>%
    arrange(interval) %>%
    mutate(
      # survival from your hazards
      S_u = cumprod(1 - h_u),
      S_w = cumprod(1 - h_w),
      S_sw = cumprod(1 - h_sw),
      
      # cumulative incidence = 1 - survival
      risk_u = 1 - S_u,
      risk_w = 1 - S_w,
      risk_sw = 1 - S_sw
    ) %>%
    ungroup() %>%
    group_by(strategy, interval, t_end) %>%
    summarize(
      risk_u = mean(risk_u),
      risk_w = mean(risk_w),
      risk_sw = mean(risk_sw),
      .groups = "drop"
    ) %>%
    # reshape wide → long to compute contrasts
    pivot_wider(
      names_from = strategy,
      values_from = c(risk_u, risk_w, risk_sw)
    ) %>%
    mutate(
      RD_u = risk_u_treat_by_grace - risk_u_never_treat,
      RD_w = risk_w_treat_by_grace - risk_w_never_treat,
      RD_sw = risk_sw_treat_by_grace - risk_sw_never_treat
    )

  # Merge KM estimates by nearest time point
  dat_km <- full_join(
    km_u_wide %>% select(time, RD_km_u),
    km_w_wide %>% select(time, RD_km_w),
    by = "time"
  )
  dat_km <- dat_km %>%
    mutate(interval = findInterval(time, cuts),
           t_end = cuts[pmin(interval + 1, length(cuts))]) %>%
    group_by(t_end) %>%
    summarize(
      RD_km_u = mean(RD_km_u, na.rm = TRUE),
      RD_km_w = mean(RD_km_w, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Final merge with parametric results
  out <- dat_summ %>%
    left_join(dat_km, by = "t_end")
  
return(out)
}


  d_try = run_sim(n=10000)
  
  # quick look
  d_try %>%
    select(RD_km_u, RD_km_w, RD_u, RD_w, RD_sw) %>%
    slice_tail()
  
# Run simulations for illustration ----

  d_simtable = crossing(disc_cuts = c(10, 15, 30), n=c(500, 1000, 2000))
  
  future::plan(multisession)
  set.seed(2026)
  
  d_simtable <- d_simtable %>%
    mutate(
      sims = map2(
        disc_cuts,
        n,
        \(dc, nn)
        future_map(
          1:50,
          \(i) run_sim(
            grace = 90,
            disc_cuts = dc,
            n = nn,
            t_cens = FALSE,
            t_admin = 360,
          ),
          .options = furrr_options(seed = TRUE)
        )
      )
    )
      
  d_long <- d_simtable %>%
    mutate(
      sims = map(sims, ~ bind_rows(.x, .id = "rep"))
    ) %>%
    unnest(sims)
  
  d_rd_mc <- d_long %>%
    group_by(disc_cuts, n, t_end) %>%
    summarize(
      # na.rm=T is because some Km estimators fail
      RD_sw_mean = mean(RD_sw, na.rm=T),
      RD_sw_lo   = quantile(RD_sw, 0.025, na.rm=T),
      RD_sw_hi   = quantile(RD_sw, 0.975, na.rm=T),
      RD_w_mean = mean(RD_w, na.rm=T),
      RD_w_lo   = quantile(RD_w, 0.025, na.rm=T),
      RD_w_hi   = quantile(RD_w, 0.975, na.rm=T),
      RD_u_mean = mean(RD_u, na.rm=T),
      RD_u_lo   = quantile(RD_u, 0.025, na.rm=T),
      RD_u_hi   = quantile(RD_u, 0.975), na.rm=T,
      RD_km_u_mean = mean(RD_km_u, na.rm=T),
      RD_km_u_lo   = quantile(RD_km_u, 0.025, na.rm=T),
      RD_km_u_hi   = quantile(RD_km_u, 0.975, na.rm=T),
      RD_km_w_mean = mean(RD_km_w, na.rm=T),
      RD_km_w_lo   = quantile(RD_km_w, 0.025, na.rm=T),
      RD_km_w_hi   = quantile(RD_km_w, 0.975, na.rm=T),
      .groups = "drop"
    ) 
  
  d_plotmeans = d_rd_mc %>%
    pivot_longer(cols = c(contains('mean')),
                 names_to = 'model',
                 names_pattern = '_(sw|u|w|km_u|km_w)_', 
                 values_to = 'RD') %>%
    select(t_end, disc_cuts, n, model, RD) 
  
  ggplot(d_plotmeans, aes(x = t_end)) +
    geom_line(aes(y = RD, linetype = model, color = model), linewidth=1.1) +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_grid(disc_cuts ~ n, labeller = label_both) +
    labs(
      title = "Risk Difference with Monte Carlo 95% Intervals",
      x = "Days",
      y = "Risk Difference",
    ) +
    theme_minimal()
  
  ggplot(d_rd_mc, aes(x = t_end)) +
    # Stabilized ribbon
    geom_ribbon(aes(ymin = RD_sw_lo, ymax = RD_sw_hi),
                fill = "green", alpha = 0.15) +
    # Weighted ribbon
    geom_ribbon(aes(ymin = RD_w_lo, ymax = RD_w_hi),
                fill = "red", alpha = 0.15) +
    # Unweighted ribbon
    geom_ribbon(aes(ymin = RD_u_lo, ymax = RD_u_hi),
                fill = "blue", alpha = 0.10) +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_grid(disc_cuts ~ n, labeller = label_both) +
    labs(
      title = "Risk Difference with Monte Carlo 95% Intervals",
      x = "Days",
      y = "Risk Difference",
      caption = "Green = Stabilized; Red = weighted; Blue = unweighted"
    ) +
    theme_minimal()
    
    d_ci_long <- d_long %>%
      select(disc_cuts, n, rep, interval, t_end,
             risk_w_treat_by_grace, risk_w_never_treat,
             risk_u_treat_by_grace, risk_u_never_treat,
             risk_sw_treat_by_grace, risk_sw_never_treat) %>%
      pivot_longer(
        cols = starts_with("risk_"),
        names_to = c("type", "strategy"),
        names_pattern = "risk_(sw|w|u)_(.*)",
        values_to = "risk"
      ) %>%
      mutate(
        type = case_when(type == "w" ~"Weighted", 
                         type == 'sw' ~ "Stabilized",
                         type == 'u' ~ 'Unweighted'),
        strategy = recode(strategy,
                          "treat_by_grace" = "Treat by Grace",
                          "never_treat" = "Never Treat")
      )
    
    d_ci_mc <- d_ci_long %>%
      group_by(disc_cuts, n, t_end, type, strategy) %>%
      summarize(
        risk_mean = mean(risk),
        risk_lo   = quantile(risk, 0.025),
        risk_hi   = quantile(risk, 0.975),
        .groups = "drop"
      )
    
    ggplot(d_ci_mc, aes(x = t_end)) +
      # Weighted ribbon
      geom_ribbon(
        data = subset(d_ci_mc, type == "Stabilized"),
        aes(ymin = risk_lo, ymax = risk_hi),
        fill = "green", alpha = 0.15
      ) +
    # Weighted ribbon
      geom_ribbon(
        data = subset(d_ci_mc, type == "Weighted"),
        aes(ymin = risk_lo, ymax = risk_hi),
        fill = "red", alpha = 0.15
      ) +
      # Unweighted ribbon
      geom_ribbon(
        data = subset(d_ci_mc, type == "Unweighted"),
        aes(ymin = risk_lo, ymax = risk_hi),
        fill = "blue", alpha = 0.10
      ) +
      # Mean curves
      geom_line(aes(y = risk_mean, color = type), size = 1.2) +
      
      facet_grid(strategy ~ n + disc_cuts, labeller = label_both) +
      
      labs(
        title = "Cumulative Incidence with Monte Carlo 95% Intervals",
        x = "Days",
        y = "Cumulative Incidence",
        caption = "Red = weighted; Blue = unweighted"
      ) +
      
      scale_color_manual(values = c("Stabilized" = "green", "Weighted" = "red", "Unweighted" = "blue")) +
      theme_minimal()

