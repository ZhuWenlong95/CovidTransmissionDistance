# exposure: time of exposure
# onset: time of symptom onset, onset = exposure + incubfn_sample
# pos time: agents are undetectable before the time, equal to min(exposure + pos_samples, onset)
# detected_time: the time agents are detected
# test_time: the time of agents are tested during massive NAATs
# isolated_time: time of agents are isolated


outbreak_step <- function(case_data = NULL, disp.iso = NULL, disp.com = NULL, r0isolated = NULL, r0community = NULL,
                          # p.asym = NULL, p.asym.from.asym = NULL,
                          p.vaccinated = NULL, p.ascertain = NULL,
                          # p.test = NULL,
                          incfn = NULL, delayfn_long = NULL, delayfn_short = NULL, posfn = NULL, GTfu = NULL, distancefn = NULL,
                          NAAT = NAAT, NAAT_distance = NAAT_distance,
                          inf_fn = NULL, k = NULL, prob_NAAT = NULL
                          # tests = NULL, t = NULL
) {
  # For each case in case_data, draw new_cases from a negative binomial distribution
  # with an R0 and dispersion dependent on if isolated=TRUE
  # new_cases <- c()
  cases_avoid <- 0

  # browser()

  case_data_old <- case_data[!(!generated_NC & is.na(new_cases)), ]
  case_data_new <- case_data[!generated_NC & is.na(new_cases), ]

  case_data_new <- case_data_new[isolated == T, ][, new_cases := rnbinom(.N, size = disp.iso, mu = r0isolated)] %>%
    rbind(
      case_data_new[isolated == F, ][, new_cases := rnbinom(.N, size = disp.com, mu = r0community)]
    )

  case_data <- rbind(case_data_old, case_data_new)

  # Select cases that have generated any new cases
  new_case_data <- case_data[new_cases > 0 & !generated_NC]
  # The total new cases generated
  total_new_cases <- sum(new_case_data$new_cases, na.rm = T)

  # If no new cases drawn, outbreak is over so return case_data
  if (total_new_cases == 0) {
    # If everyone is isolated it means that either control has worked or everyone has had a chance to infect but didn't
    case_data$isolated <- TRUE
    case_data$generated_NC <- TRUE

    effective_r0 <- 0
    cases_in_gen <- 0
    positive_in_gen <- 0
    cases_avoid <- 0
    out <- list(case_data, effective_r0, cases_in_gen, positive_in_gen, cases_avoid)
    names(out) <- c("cases", "effective_r0", "cases_in_gen", "positive_in_gen", "cases_avoid")

    return(out)
  }

  # Compile a data.table for all new cases, new_cases is the amount of people that each infector has infected
  inc_samples <- incfn(total_new_cases)

  prob_samples_raw <- data.table(
    # time when new cases were exposed, a draw from serial interval based on infector's onset
    # exposure = unlist(purrr::map2(
    #   new_case_data$new_cases, new_case_data$onset,
    #   function(x, y) {
    #     inf_fn(rep(y, x), k)
    #   }
    # )),
    # time when new cases were exposed, a draw from generation time based on infector's exposure
    exposure = rep(new_case_data$exposure, new_case_data$new_cases) + GTfn(sum(new_case_data$new_cases)),

    # records the infector of each new person
    infector = unlist(purrr::map2(
      new_case_data$caseid, new_case_data$new_cases,
      function(x, y) {
        rep(as.integer(x), as.integer(y))
      }
    )),
    # records when infector exposure
    infector_exposure = unlist(purrr::map2(
      new_case_data$exposure, new_case_data$new_cases,
      function(x, y) {
        rep(x, as.integer(y))
      }
    )),
    # records when infector onset
    infector_onset_time = unlist(purrr::map2(
      new_case_data$onset, new_case_data$new_cases,
      function(x, y) {
        rep(x, as.integer(y))
      }
    )),
    # records when infector isolated
    infector_isolated_time = unlist(purrr::map2(
      new_case_data$isolated_time, new_case_data$new_cases,
      function(x, y) {
        rep(x, as.integer(y))
      }
    )),
    infector_NAAT_time = unlist(purrr::map2(
      new_case_data$NAAT_time, new_case_data$new_cases,
      function(x, y) {
        rep(x, as.integer(y))
      }
    )),
    infector_detected_time = unlist(purrr::map2(
      new_case_data$detected_time, new_case_data$new_cases,
      function(x, y) {
        rep(x, as.integer(y))
      }
    )),
    # records if infector asymptomatic
    # infector_asym = unlist(purrr::map2(new_case_data$asym, new_case_data$new_cases,
    #                                    function(x, y) {
    #                                      rep(x, y)
    #                                    })),
    # draws a sample to see if this person is vaccinated
    vaccinated = purrr::rbernoulli(n = total_new_cases, p = p.vaccinated),
    distance = distancefn(total_new_cases),
    # draws a sample to see if this person is missed
    # missed = purrr::rbernoulli(n = total_new_cases, p = 1 - p.ascertain),
    # draws a sample to see if this person is tested positive
    # tested = purrr::rbernoulli(n = total_new_cases, p = p.test),
    # sample from the incubation period for each new person
    incubfn_sample = inc_samples,
    NAAT_selected = 0,
    NAAT_time = Inf,
    isolated = FALSE,
    generated_NC = FALSE,
    new_cases = NA
  )

  # filter out new cases prevented by vaccination and isolation
  prob_samples <- prob_samples_raw[exposure < infector_exposure + 14, ][, `:=`(onset = exposure + incubfn_sample, incubfn_sample = NULL)]

  cases_avoid <- cases_avoid + prob_samples[vaccinated == T, .N]
  prob_samples <- prob_samples[vaccinated == FALSE, ]

  # If no new cases remained, outbreak is over so return case_data
  if (nrow(prob_samples) == 0) {
    # If everyone is isolated it means that either control has worked or everyone has had a chance to infect but didn't
    case_data$isolated <- TRUE
    case_data$generated_NC <- TRUE

    effective_r0 <- 0
    cases_in_gen <- 0
    positive_in_gen <- 0
    # cases_avoid <- 0
    out <- list(case_data, effective_r0, cases_in_gen, positive_in_gen, cases_avoid)
    names(out) <- c("cases", "effective_r0", "cases_in_gen", "positive_in_gen", "cases_avoid")

    return(out)
  }


  # cases whose parents are asymptomatic are more likely to be asymptomatic
  # prob_samples <- prob_samples %>%
  #   filter(infector_asym) %>%
  #   mutate(asym=rbernoulli(nrow(.),p=p.asym.from.asym)) %>%
  #   rbind(
  #     prob_samples %>%
  #       filter(!infector_asym) %>%
  #       mutate(asym=rbernoulli(nrow(.),p=p.asym))
  #   )

  # sample from the positive period for each new person
  pos_samples <- posfn(nrow(prob_samples)) # interval between exposrure and detectable
  prob_samples[, pos_time := vect_min(prob_samples$exposure + pos_samples, prob_samples$onset)]

  # the time you should be detected
  prob_samples[, delay := delayfn_long(.N)][, detected_time := onset + delay][, delay := NULL]

  # the time you should be tested

  # cases_avoid <- 0
  if (NAAT) {
    prob_samples[distance <= NAAT_distance, NAAT_selected := rbinom(.N, size = 1, prob = prob_NAAT)][NAAT_selected == 1, NAAT_delay := delayfn_long(.N)][, NAAT_time := fifelse(is.na(NAAT_delay), pos_time + Inf, pos_time + NAAT_delay)] # [, NAAT_selected := NULL]
  }
  # browser()
  cases_avoid <- cases_avoid + prob_samples[infector_detected_time > infector_NAAT_time & exposure > infector_NAAT_time, .N]
  prob_samples <- prob_samples[exposure <= infector_isolated_time, ]


  if (nrow(prob_samples) == 0) {
    # If everyone is isolated it means that either control has worked or everyone has had a chance to infect but didn't
    case_data$isolated <- TRUE
    case_data$generated_NC <- TRUE

    effective_r0 <- 0
    cases_in_gen <- 0
    positive_in_gen <- 0
    out <- list(case_data, effective_r0, cases_in_gen, positive_in_gen, cases_avoid)
    names(out) <- c("cases", "effective_r0", "cases_in_gen", "positive_in_gen", "cases_avoid")

    return(out)
  }

  # Set new case ids for new people
  prob_samples$caseid <- (nrow(case_data) + 1):(nrow(case_data) + nrow(prob_samples))

  samples <- prob_samples %>%
    # the time you are isolated
    mutate(isolated_time = vect_min(detected_time, NAAT_time)) %>%
    # select the columns needed
    select(exposure, caseid, infector, onset, pos_time, vaccinated, distance, detected_time, NAAT_selected, NAAT_time, isolated, generated_NC, isolated_time, new_cases)

  ## Number of new cases
  positive_in_gen <- nrow(samples)
  # cases_in_gen <- sum(!vect_isTRUE(samples$asym))
  cases_in_gen <- nrow(samples)

  ## Estimate the effective r0
  effective_r0 <- nrow(samples) / nrow(case_data[!vect_isTRUE(case_data$generated_NC)])
  # effective_r0 <- nrow(samples) / nrow(case_data[!vect_isTRUE(case_data$isolated)])

  # Everyone in case_data so far has had their chance to infect and are therefore considered isolated
  case_data$isolated <- TRUE
  case_data$generated_NC <- TRUE

  # bind original cases + new secondary cases
  case_data <- data.table::rbindlist(list(case_data, samples),
    use.names = TRUE
  )

  # Return
  out <- list(case_data, effective_r0, cases_in_gen, positive_in_gen, cases_avoid)
  names(out) <- c("cases", "effective_r0", "cases_in_gen", "positive_in_gen", "cases_avoid")

  return(out)
}
