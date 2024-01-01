outbreak_setup <- function(num.initial.cases, incfn, delayfn_long, posfn) {
  # Set up table of initial cases
  exposure_samples <- sample(1, num.initial.cases, replace = TRUE)
  inc_samples <- incfn(num.initial.cases) + exposure_samples
  pos_samples <- posfn(num.initial.cases) + exposure_samples

  case_data <- data.table(
    exposure = exposure_samples, # Exposure time of 0 for all initial cases
    # asym = F,#purrr::rbernoulli(num.initial.cases, p.asym),
    caseid = 1:(num.initial.cases), # set case id
    infector = 0,
    onset = inc_samples,
    pos_time = vect_min(inc_samples, pos_samples),
    vaccinated = FALSE,
    distance = 0
    # missed = TRUE
  )

  # set isolation time for cluster to minimum time of onset of symptoms + draw from delay distribution
  case_data <- case_data[, detected_time := onset + delayfn_long(1)][, `:=`(NAAT_selected = 0, NAAT_time = detected_time, isolated = FALSE, generated_NC = FALSE)]
  # [, NAAT_time := detected_time][, isolated := FALSE][, generated_NC := FALSE]

  # case_data$detected_time[case_data$asym] <- Inf
  # browser()
  # case_data$tested <- purrr::rbernoulli(num.initial.cases, p.test)
  # return
  return(case_data)
}
