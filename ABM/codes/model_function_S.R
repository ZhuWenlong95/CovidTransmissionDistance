weibull_setup <- function(dist_shape = NULL, dist_scale = NULL) {
  out <- purrr::partial(rweibull,
    shape = dist_shape,
    scale = dist_scale
  )
  return(out)
}

lognorm_setup <- function(meanlog = NULL, sdlog = NULL) {
  out <- purrr::partial(rlnorm,
    meanlog = meanlog,
    sdlog = sdlog
  )
  return(out)
}

gamma_setup <- function(dist_shape = NULL, dist_scale = NULL) {
  out <- purrr::partial(rgamma,
    shape = dist_shape,
    scale = dist_scale
  )
  return(out)
}

norm_setup <- function(dist_mean = NULL, dist_sd = NULL) {
  out <- purrr::partial(rnorm,
    mean = dist_mean,
    sd = dist_sd
  )
  return(out)
}

exp_setup <- function(dist_mean = NULL) {
  out <- purrr::partial(rexp,
    rate = 1 / dist_mean
  )
  return(out)
}

# # incubation period sampling function
# incfn <- weibull_setup(dist_shape = 3.34,
#                     dist_scale = 9.28)
#
# # onset to isolation delay sampling function
# delayfn <- weibull_setup(dist_shape = 2.5,
#                       dist_scale = 5)
#
# # exposure to positive delay sampling function
# posfn <- norm_setup(dist_mean = 3,
#                     dist_sd = 1)


vect_isTRUE <- function(x) {
  purrr::map_lgl(x, isTRUE)
}

vect_max <- function(x, y) {
  purrr::map2_dbl(x, y, max)
}

vect_min <- function(x, y) {
  purrr::map2_dbl(x, y, min)
}


inf_fn <- function(inc_samp = NULL, k = NULL) {
  out <- sn::rsn(
    n = length(inc_samp),
    xi = inc_samp,
    omega = 2,
    alpha = k
  )

  out <- ifelse(out < 1, 1, out)

  return(out)
}


# my function -------------------------------------------------------------------------

RI <- function(x, p = 0.975) {
  Mean <- mean(x)
  # upper <- HDInterval::hdi(x)[2]
  # lower <- HDInterval::hdi(x)[1]
  upper <- mean(x) + qnorm(p) * sd(x)
  lower <- mean(x) - qnorm(p) * sd(x)
  return(c(upper = upper, mean = Mean, lower = lower))
}

df_combine <- function(data = NULL) {
  df.daily <- NULL
  for (i in 1:length(data)) {
    df.daily <- df.daily %>%
      rbind(data[[i]][[2]] %>% mutate(id = i))
  }
  return(df.daily)
}

df_curve <- function(data = NULL, type = "CI") {
  df.daily <- data
  if (type == "CI") {
    df.daily2 <- df.daily %>%
      reframe(
        daily_temp = list(Rmisc::CI(daily_total)),
        cum_temp = list(Rmisc::CI(cumul_total)),
      )
  } else {
    df.daily2 <- df.daily %>%
      reframe(
        daily_temp = list(c(quantile(daily_total, 0.75, na.rm = T), mean(daily_total, na.rm = T), quantile(daily_total, 0.25, na.rm = T))),
        cum_temp = list(c(quantile(cumul_total, 0.75, na.rm = T), mean(cumul_total, na.rm = T), quantile(cumul_total, 0.25, na.rm = T)))
      )
  }
  df.daily2 %>%
    rowwise() %>%
    mutate(
      daily_lower = daily_temp[3],
      daily_mean = daily_temp[2],
      daily_upper = daily_temp[1],
      cum_lower = cum_temp[3],
      cum_mean = cum_temp[2],
      cum_upper = cum_temp[1],
    ) %>%
    return()
}

df_location <- function(data, df_grid, .) {
  df.location <- NULL
  for (i in unique(data$infector)) {
    temp <- data %>% filter(infector == i)
    # browser()

    if (i == 0) {
      n_temp <- rpois(1, 16)
      if (n_temp > 20) {
        n_temp <- 20
      }

      temp <- temp %>%
        # arrange(onset) %>%
        mutate(
          infector_temp = c(0, 0, rep(1, n_temp), rep(2, nrow(.) - 2 - n_temp)),
          distance_temp = distancefn(nrow(.))
        )

      i_temp <- 0
      while (i_temp <= 2) {
        if (i_temp == 0) {
          temp1 <- temp %>%
            filter(infector_temp == i_temp) %>%
            mutate(X = 0, Y = 0) %>%
            select(caseid, X, Y)
        } else {
          infector_XY <- df.location %>%
            filter(caseid == i_temp) %>%
            select(X, Y)

          center.point <- st_point(x = c(as.numeric(infector_XY$X), as.numeric(infector_XY$Y)))

          temp_temp <- temp %>% filter(infector_temp == i_temp)

          temp1 <- sapply(temp_temp$distance_temp,
            FUN = function(v.distance) {
              temp <- st_buffer(center.point, v.distance) %>%
                # 取圆的边界
                st_boundary() %>%
                # 取相交部分
                st_intersection(df_grid)

              # 判断是否有相交的部分
              if (length(temp) == 0) {
                # next
                temp <- st_bbox(df_grid)
                out <- list(
                  X = runif(n = 1, min = temp["xmin"], max = temp["xmax"]),
                  Y = runif(n = 1, min = temp["ymin"], max = temp["ymax"])
                )
              } else {
                temp <- temp %>%
                  # 合并线段为多段线对象
                  st_combine() %>%
                  # 将MULTILINESTRING转换为LINESTRING
                  st_cast(., "LINESTRING") %>%
                  # 取点
                  st_line_sample(n = 100, type = "random") %>%
                  st_coordinates()

                out <- temp[sample(1:nrow(temp), 1), 1:2] %>% as.list()
              }
              return(out)
            }
          ) %>%
            t() %>%
            as.data.table() %>%
            mutate(caseid = temp_temp$caseid) %>%
            select(caseid, everything())
        }
        i_temp <- i_temp + 1
        df.location <- rbind(df.location, temp1)
      }

      # temp1 <- temp %>%
      #   mutate(X = runif(nrow(.),-20,20), Y = runif(nrow(.),-20,20)) %>%
      #   select(caseid, X, Y)
    } else {
      infector_XY <- df.location %>%
        filter(caseid == i) %>%
        select(X, Y)

      center.point <- st_point(x = c(as.numeric(infector_XY$X), as.numeric(infector_XY$Y)))


      temp1 <- sapply(temp$distance,
        FUN = function(v.distance) {
          temp <- st_buffer(center.point, v.distance) %>%
            # 取圆的边界
            st_boundary() %>%
            # 取相交部分
            st_intersection(df_grid)

          # 判断是否有相交的部分
          if (length(temp) == 0) {
            # next
            temp <- st_bbox(df_grid)
            out <- list(
              X = runif(n = 1, min = temp["xmin"], max = temp["xmax"]),
              Y = runif(n = 1, min = temp["ymin"], max = temp["ymax"])
            )
          } else {
            temp <- temp %>%
              # 合并线段为多段线对象
              st_combine() %>%
              # 将MULTILINESTRING转换为LINESTRING
              st_cast(., "LINESTRING") %>%
              # 取点
              st_line_sample(n = 100, type = "random") %>%
              st_coordinates()

            out <- temp[sample(1:nrow(temp), 1), 1:2] %>% as.list()
          }
          return(out)
        }
      ) %>%
        t() %>%
        as.data.table() %>%
        mutate(caseid = temp$caseid) %>%
        select(caseid, everything())
      df.location <- rbind(df.location, temp1)
    }
  }
  return(df.location)
}

prop_PopTest <- function(data, v.NAAT_distance) {
  df_case <- NULL

  for (i in 1:length(data)) {
    df_case <- df_case %>% rbind(data[[i]][[1]] %>% mutate(id = i))
  }

  N.total <- df_case[, .N, id][, `:=`(N.total = N, N = NULL)]
  N.in <- df_case[distance <= v.NAAT_distance, .N, id][, `:=`(N.in = N, N = NULL)]
  N.selected <- df_case[NAAT_selected == 1, .N, id][, `:=`(N.selected = N, N = NULL)]

  N.total %>%
    left_join(N.in, by = "id") %>%
    left_join(N.selected, by = "id") %>%
    mutate(NAAT_distance = v.NAAT_distance) %>%
    return()
}
