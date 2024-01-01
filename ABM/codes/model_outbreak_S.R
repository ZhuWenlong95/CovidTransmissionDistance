outbreak_model <- function(num.initial.cases = NULL,
                           cap_max_days = NULL, cap_cases = NULL,
                           r0isolated = NULL, r0community = NULL,
                           disp.iso = NULL, disp.com = NULL,
                           # p.asym = NULL, p.asym.from.asym = NULL,
                           p.vaccinated = NULL, p.ascertain = NULL,
                           # p.test = NULL,
                           incfn = NULL, delayfn_long = NULL, delayfn_short = NULL,
                           posfn = NULL, GTfn = NULL, distancefn = NULL,
                           NAAT = NULL, NAAT_distance = NULL,
                           inf_fn = NULL, k = NULL, prob_NAAT = NULL
                           # , t = NULL
) {
  # Set initial values for loop indices
  total.cases <- num.initial.cases
  latest.onset <- 0
  extinct <- FALSE

  # Initial setup
  case_data_raw <- outbreak_setup(
    num.initial.cases = num.initial.cases,
    incfn = incfn,
    # p.asym = p.asym,
    delayfn_long = delayfn_long,
    posfn = posfn
    # p.test = p.test
  )


  # t <- t
  # test.start <- min(case_data_raw$detected_time)
  # start.diff = start.diff

  case_data <- copy(case_data_raw)[, `:=`(isolated_time = detected_time, new_cases = NA)]


  # Preallocate
  effective_r0_vect <- c()
  cases_in_gen_vect <- c()
  positive_in_gen_vect <- c()
  cases_avoid_vect <- c()

  # Model loop
  while (latest.onset < cap_max_days & total.cases < cap_cases & !extinct) {
    out <- outbreak_step(
      case_data = case_data,
      disp.iso = disp.iso, disp.com = disp.com,
      r0isolated = r0isolated, r0community = r0community,
      # p.asym = p.asym, p.asym.from.asym = p.asym.from.asym,
      p.vaccinated = p.vaccinated, p.ascertain = p.ascertain,
      # p.test = p.test,
      incfn = incfn, delayfn_long = delayfn_long, delayfn_short = delayfn_long,
      posfn = posfn, GTfu = GTfn, distancefn = distancefn,
      NAAT = NAAT, NAAT_distance = NAAT_distance,
      inf_fn = inf_fn, k = k, prob_NAAT = prob_NAAT
      # , tests = tests, t = t
    )


    case_data <- out[[1]]
    effective_r0_vect <- c(effective_r0_vect, out[[2]])
    cases_in_gen_vect <- c(cases_in_gen_vect, out[[3]])
    positive_in_gen_vect <- c(positive_in_gen_vect, out[[4]])
    cases_avoid_vect <- c(cases_avoid_vect, out[[5]])
    total.cases <- nrow(case_data) + sum(cases_avoid_vect, na.rm = T)
    latest.onset <- max(case_data$onset)
    extinct <- all(case_data$isolated)
  }

  # results clean -----------------------------------------------------------


  # df1 <- case_data %>%
  #   mutate(
  #     exposure = round(exposure)
  #   ) %>%
  #   group_by(exposure) %>%
  #   summarise(daily_total = n())

  df1 <- case_data %>%
    mutate(
      exposure = round(onset)
    ) %>%
    group_by(exposure) %>%
    summarise(daily_total = n())

  df.out <- data.frame(exposure = 1:cap_max_days) %>%
    left_join(df1, by = "exposure") %>%
    mutate(
      daily_total = ifelse(is.na(daily_total), 0, daily_total),
      cumul_total = cumsum(daily_total)
    )

  return(list(case_data, df.out))
}
