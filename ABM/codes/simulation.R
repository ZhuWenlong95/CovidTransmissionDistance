require(tidyverse)
require(data.table)
# require(parallel)
require(ggplot2)
require(ggpubr)
require(ggsci)
require(gridExtra)
require(scales)
require(furrr)
require(sf)
require(gghalves)

source("codes/model_function_S.R")
source("codes/model_setup_S.R")
source("codes/model_step_S.R")
source("codes/model_outbreak_S.R")

theme_set(
  theme_bw(base_size = 22) +
    theme(
      axis.text = element_text(colour = "black"),
      # panel.grid = element_line(linewidth = 0.5),
      panel.grid = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "lines")
    )
)

# generation time sampling function
GT.mean <- 2.99
GT.sd <- 1.31
GTfn <- gamma_setup(dist_shape = GT.mean^2 / GT.sd^2, dist_scale = GT.sd^2/GT.mean)

# incubation period sampling function
inc.mean <- 4.06
inc.sd <- 1.91
incfn <- weibull_setup(dist_shape = (inc.sd / inc.mean)^(-1.086), dist_scale = inc.mean / gamma(1 + 1 / (inc.sd / inc.mean)^(-1.086)))
# summary(incfn(10000))

# onset to isolation delay sampling function
delay.mean <- 3
delay.sd <- 2
# delay.shape
(delay.sd / delay.mean)^(-1.086)
# delay.scale
delay.mean / gamma(1 + 1 / (delay.sd / delay.mean)^(-1.086))
delayfn_short <- weibull_setup(dist_shape = 0.4710656, dist_scale = 0.4451066) # mean=1,sd=2
delayfn_long <- weibull_setup(dist_shape = 1.553228, dist_scale = 3.336312) # mean=3, sd=2

# exposure to positive delay sampling function
posfn <- norm_setup(dist_mean = 1, dist_sd = 2)

# distributiong of transmission distance
mean_distance <- 20.17
distancefn <- exp_setup(dist_mean = mean_distance)
# summary(distancefn(1000))

# single simulation
# temp2 <- outbreak_model(num.initial.cases = 50,
#                cap_max_days = 100,
#                cap_cases = 26500,
#                r0isolated = 0,
#                r0community = 7.9,
#                disp.iso = 1,
#                disp.com = 0.1,
#                p.vaccinated = 0,
#                p.ascertain = 0.596,
#                incfn = incfn, delayfn_long = delayfn_long, delayfn_short = delayfn_short, posfn = posfn, GTfn = GTfn, distancefn = distancefn, inf_fn = inf_fn, k = 0.7,
#                NAAT=T,NAAT_distance=mean_distance/2,prob_NAAT = 0.1)
# temp2[[2]] %>%
#   ggplot() +
#   geom_line(aes(exposure, cumul_total))


# no testing -------------------------------------------------------------------------

plan(multisession, workers = 10)

df.nonNAAT <- future_map(rep(20, 10),
  .options = furrr_options(seed = NULL),
  .f = function(x) {
    out <- purrr::map(1:x, ~ outbreak_model(
      num.initial.cases = 50,
      cap_max_days = 100,
      cap_cases = 26500,
      r0isolated = 0,
      r0community = 7.9,
      disp.iso = 1,
      disp.com = 0.1,
      p.vaccinated = 0,
      p.ascertain = 0.596,
      incfn = incfn, delayfn_long = delayfn_long, delayfn_short = delayfn_short, posfn = posfn, GTfn = GTfn, distancefn = distancefn, inf_fn = inf_fn, k = 0.7,
      NAAT = F, NAAT_distance = mean_distance, prob_NAAT = 1
    ))
    return(out)
  }
)

plan("default")

df.nonNAAT <- do.call("rbind", df.nonNAAT)
df.nonNAAT.daily <- df_combine(data = df.nonNAAT) %>% group_by(exposure) %>% df_curve(type = "P")

df.nonNAAT.daily %>%
  ggplot(aes(exposure, y = cum_mean, ymin = cum_lower, ymax = cum_upper)) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  scale_x_continuous(limits = c(0, 50))

df.nonNAAT.daily %>%
  ggplot(aes(exposure, y = daily_mean, ymin = daily_lower, ymax = daily_upper)) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  scale_x_continuous(limits = c(0, 50))

# testing ----------------------------------------------------------------

## test proportion ####

df.NAAT.daily_prob <- df.NAAT.combine_prob <- NULL
plan(multisession, workers = 10)

for (v_prob in seq(0.1, 0.9, 0.1)) {
  df.NAAT <- future_map(rep(20, 10),
    .options = furrr_options(seed = NULL),
    .f = function(x) {
      out <- purrr::map(1:x, ~ outbreak_model(
        num.initial.cases = 50,
        cap_max_days = 100,
        cap_cases = 26500,
        r0isolated = 0,
        r0community = 7.9,
        disp.iso = 1,
        disp.com = 0.1,
        p.vaccinated = 0,
        p.ascertain = 0.596,
        incfn = incfn, delayfn_long = delayfn_long, delayfn_short = delayfn_short, posfn = posfn, GTfn = GTfn, distancefn = distancefn, inf_fn = inf_fn, k = 0.7,
        NAAT = T, NAAT_distance = mean_distance, prob_NAAT = v_prob
      ))
      return(out)
    }
  )
  print(v_prob)

  df.NAAT <- do.call("rbind", df.NAAT)
  temp_combine <- df_combine(data = df.NAAT) %>% mutate(prob_NAAT = v_prob)
  df.NAAT.combine_prob <- rbind(df.NAAT.combine_prob, temp_combine)

  temp_daily <- df_combine(data = df.NAAT) %>%
    group_by(exposure) %>%
    df_curve(type = "P") %>%
    mutate(prob_NAAT = v_prob)
  df.NAAT.daily_prob <- rbind(df.NAAT.daily_prob, temp_daily)
}

plan("default")

df.NAAT.daily <- df.NAAT.daily_prob %>%
  ungroup() %>%
  filter(prob_NAAT == 0.9) %>%
  select(-prob_NAAT)

### A curve ####

df.all.daily <- df.nonNAAT.daily %>%
  mutate(group = "No testing") %>%
  rbind(df.NAAT.daily %>%
    mutate(group = "With testing")) %>%
  ungroup()

df.all.daily %>%
  group_by(group) %>%
  filter(exposure == max(exposure)) %>%
  select(cum_lower, cum_mean, cum_upper) %>%
  view()

df.all.daily %>%
  group_by(group) %>%
  filter(daily_mean == max(daily_mean)) %>%
  select(daily_mean) %>%
  ungroup() %>%
  mutate(prop = 100 * (max(daily_mean) - min(daily_mean)) / max(daily_mean)) %>%
  view()

fig.all.cum <- df.all.daily %>%
  ggplot(aes(exposure, y = cum_mean, ymin = cum_lower, ymax = cum_upper, colour = group, fill = group)) +
  geom_line(linewidth = 2) +
  geom_ribbon(colour = NA, alpha = 0.1) +
  scale_x_continuous(name = "Infection time", expand = c(0, 0), limits = c(0, 50)) +
  scale_y_continuous(name = expression("Cumulative No. of cases （" %*% ~ 10^4 * ")"), expand = expansion(mult = c(0.01, 0.1)), breaks = seq(0, 150000, 30000), labels = seq(0, 15, 3)) +
  scale_color_manual(values = c("#811400", "#006d81")) +
  scale_fill_manual(values = c("#811400", "#006d81")) +
  theme(
    legend.position = c(0.1, 0.9),
    legend.justification = c(0, 1),
    legend.title = element_blank()
  )
fig.all.cum

fig.all.daily <- df.all.daily %>%
  ggplot(aes(exposure, y = daily_mean, ymin = daily_lower, ymax = daily_upper, colour = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(colour = NA, alpha = 0.1) +
  scale_x_continuous(name = "Infection time", expand = c(0, 0), limits = c(0, 50)) +
  scale_y_continuous(name = expression(atop("Daily No. of", "cases (" * 10^3 * ")")), expand = expansion(mult = c(0.01, 0.1)), breaks = seq(0, 20000, 4000), labels = seq(0, 20, 4)) +
  scale_color_manual(values = c("#811400", "#006d81")) +
  scale_fill_manual(values = c("#811400", "#006d81")) +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA)
  )
fig.all.daily

# fig.all <- fig.all.cum + annotation_custom(ggplotGrob(fig.all.daily), xmin = 16, xmax = Inf, ymin = -Inf, ymax = 52500)
# fig.all

df.all.daily %>%
  group_by(group) %>%
  filter(cum_mean == max(cum_mean)) %>%
  filter(exposure == min(exposure)) %>%
  mutate(cum_cases = paste0(prettyNum(ceiling(cum_mean), big.mark = ","), " (95% CI: ", prettyNum(ceiling(cum_lower), big.mark = ","), ", ", prettyNum(ceiling(cum_upper), big.mark = ","), ")")) %>%
  select(group, cum_cases) %>%
  view()

### B proportion ####

df.NAAT.daily_prob <- df.NAAT.combine_prob %>% 
  group_by(exposure,prob_NAAT) %>% 
  df_curve(type = "P")

fig.NAAT_prob <- df.NAAT.daily_prob %>%
  rbind(df.nonNAAT.daily %>% mutate(prob_NAAT=0)) %>% 
  ungroup() %>%
  filter(exposure == max(exposure)) %>%
  mutate(prob_NAAT = factor(prob_NAAT)) %>%
  # view()
  ggplot(aes(prob_NAAT, cum_mean, ymin = cum_lower, ymax = cum_upper)) +
  geom_col(fill = NA, linewidth = 1, colour="black",width = 0.5) +
  geom_errorbar(linewidth = 1, width = 0.25) +
  # scale_alpha_discrete(range = c(0.9, 0.1)) +
  scale_y_continuous(name = expression("Cumulative No. of cases (" %*% ~ 10^4 * ")"), breaks = seq(0, 150000, 30000), labels = seq(0, 15, 3), limits = c(0, 150000), expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(name = "Proportion of population tested (%)",labels=seq(0,90,10)) +
  theme(
    legend.position = "none"
  )
fig.NAAT_prob

## distance ####

# early testing
df.NAAT.daily_dis <- df.NAAT.combine_dis <- df.NAAT.propTest_dis <- NULL
plan(multisession, workers = 10)
for (i in qexp(rate = 1 / mean_distance, p = seq(0.1, 0.9, 0.1))) {
  NAAT_distance_temp <- i
  df.NAAT <- future_map(rep(20, 10),
    .options = furrr_options(seed = NULL),
    .f = function(x) {
      out <- purrr::map(1:x, ~ outbreak_model(
        num.initial.cases = 50,
        cap_max_days = 100,
        cap_cases = 26500,
        r0isolated = 0,
        r0community = 7.9,
        disp.iso = 1,
        disp.com = 0.1,
        p.vaccinated = 0,
        p.ascertain = 0.596,
        incfn = incfn, delayfn_long = delayfn_long, delayfn_short = delayfn_short, posfn = posfn, GTfn = GTfn, distancefn = distancefn, inf_fn = inf_fn, k = 0.7,
        NAAT = T, NAAT_distance = NAAT_distance_temp, prob_NAAT = 0.9
      ))
      return(out)
    }
  )
  df.NAAT <- do.call("rbind", df.NAAT)
  
  df.NAAT.propTest_dis <- rbind(df.NAAT.propTest_dis,prop_PopTest(data = df.NAAT,v.NAAT_distance = NAAT_distance_temp))

  temp <- df_combine(data = df.NAAT) %>% mutate(group = round(i, 2))
  df.NAAT.combine_dis <- rbind(df.NAAT.combine_dis, temp)

  temp <- df_combine(data = df.NAAT) %>%
    group_by(exposure) %>%
    df_curve(type = "P") %>%
    mutate(group = round(i, 2))
  df.NAAT.daily_dis <- rbind(df.NAAT.daily_dis, temp)

  rm(list = c("df.NAAT", "temp"))
  gc()
  print(i)
}
plan("default")

df.NAAT.daily_dis <- df.NAAT.combine_dis %>% 
  group_by(exposure,group) %>% 
  df_curve(type = "P")

fig.NAAT_dis <-
  df.NAAT.daily_dis %>%
  ungroup() %>%
  filter(exposure == max(exposure)) %>%
  mutate(
    group = as.numeric(str_remove(group, "MTD\\*")),
    # group = factor(group, labels = sprintf("%.1f", group))
    group = factor(group, labels = paste0("P",seq(10,90,10)))
  ) %>%
  ggplot(aes(group, cum_mean, ymin = cum_lower, ymax = cum_upper)) +
  geom_col(fill = NA, linewidth = 1, colour="black",width = 0.5) +
  geom_errorbar(linewidth = 1, width = 0.25) +
  # geom_hline(yintercept = 54641)+
  scale_y_continuous(name = expression("Cumulative No. of cases (" %*% ~ 10^4 * ")"), breaks = seq(0, 150000, 30000), labels = seq(0, 15, 3), limits = c(0, 150000), expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(name = "Distance used to construct\npotential transmission areas") +
  theme(
    legend.position = "none"
  )
fig.NAAT_dis

# df.NAAT.propTest_dis %>% 
#   mutate(
#     prop.in=100*N.in/N.total,
#     NAAT_distance=factor(NAAT_distance)) %>% 
#   group_by(NAAT_distance) %>% 
#   # reframe(
#   #   prop.in.upper=Rmisc::CI(prop.in)[3],
#   #   prop.in.mean=Rmisc::CI(prop.in)[2],
#   #   prop.in.lower=Rmisc::CI(prop.in)[1]
#   #   ) %>% 
#   ggplot(aes(NAAT_distance,prop.in))+
#   geom_boxplot()+
#   scale_y_continuous(breaks = seq(0,100,10))+
#   theme(panel.grid.major.y = element_line("darkgrey"))
#   # view()

#### output ####

# ggsave(
#   ggarrange(fig.all.cum,
#     ggarrange(fig.NAAT_prob, fig.NAAT_dis, ncol = 2, labels = c("B", "C"), font.label = list(size = 20), align = "hv"),
#     nrow = 2, labels = c("A", NA), font.label = list(size = 20)
#   ),
#   width = 16 * 0.8, height = 16 * 0.8, dpi = 300, filename = paste0("results/TD Figure 4 model NAAT ", Sys.Date(), ".png")
# )

# save image --------------------------------------------------------------

# save.image(file = paste0("results/TD model analysis ", Sys.Date(), ".RData"))
