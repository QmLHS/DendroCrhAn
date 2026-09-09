#' Calculate Annual Accumulated Climate Anomalies from Monthly Series
#'
#' @description
#' Computes monthly climate anomalies (deviations from long-term monthly means)
#' for numeric variables in a time series and aggregates them into annual total sums.
#' Optionally applies an absolute value transformation to assess overall anomaly magnitude.
#'
#' @param df A \code{data.frame} where the first column contains character dates
#' formatted strictly as \code{"mmm-YYYY"} (e.g., \code{"gen-1901"}, exactly 8 characters long),
#' followed by one or more numeric climate variables.
#' @param mod Logical. If \code{TRUE}, takes the absolute value of monthly deviations
#' before summing them annually, measuring total variability regardless of sign.
#' If \code{FALSE} (default), preserves signed deviations (positive/negative balance).
#'
#' @return A \code{data.frame} with an \code{anno} column (numeric year) and aggregated
#' columns prefixed with \code{sum_dif_} containing annual anomaly sums for each numeric input variable.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(stringr)
#'
#' # Create sample monthly dataset
#' dates <- paste0(rep(c("gen", "feb", "mar", "apr", "mag", "giu",
#'                       "lug", "ago", "set", "ott", "nov", "dic"), 2),
#'                 "-", rep(c(2020, 2021), each = 12))
#'
#' sample_df <- data.frame(
#'   data = dates,
#'   prec = runif(24, 10, 100),
#'   temp = runif(24, 5, 30)
#' )
#'
#' # Calculate annual climate anomalies
#' res_anomalies <- scostamento(df = sample_df, mod = FALSE)
#' head(res_anomalies)
#' }
#'
#' @export
scostamento <- function(df, mod = FALSE){
  if (nchar(df[1,1]) != 8) {stop("The first column of dataframe must be date on format: gen-1901!")}
  df_due <- df
  colnames(df_due[1]) <- "data"
  df_due <- df_due %>%
    mutate (mese = as.factor(stringr::str_sub(data, start = 1, end = 3))) %>%
    group_by (mese) %>%
    mutate(across(where(is.numeric), ~ .x - mean(.x, na.rm = TRUE),
                  .names = "dif_{.col}")) %>%
    ungroup()
  if (mod) {df_tre <- df_due %>%
    mutate(across(starts_with("dif_"), abs))
  } else {df_tre <- df_due}
  df_qtr <- df_tre %>%
    mutate(anno = as.numeric(stringr::str_sub(data, -4))) %>%
    group_by(anno) %>%
    summarise(across(starts_with("dif_"), ~ sum(.x, na.rm = TRUE),
                     .names = "sum_{.col}")) %>%
    ungroup()
  return(df_qtr)}

#' Prepare Data for Tree-Ring and Climate Visualization
#'
#' @description
#' Processes Ring-Width Index (RWI) data and SPEI climate series to prepare all necessary
#' data structures for plotting. Performs detrending, chronology calculation, pointer year detection
#' (via BSGC and Interval Trend methods), SPEI normalization, and cross-correlation analysis (CCF).
#'
#' @param RWI A \code{data.frame} or \code{rwl} object containing raw or indexed tree-ring width series,
#' with years as row names and tree/core IDs as column names.
#' @param spei_df A \code{data.frame} with 2 columns: the first containing years (\code{Anno})
#' and the second containing SPEI (Standardised Precipitation-Evapotranspiration Index) values.
#'
#' @return A named \code{list} containing five elements:
#' \describe{
#'   \item{\code{rwi}}{A long-format \code{data.frame} with columns \code{Anno}, \code{serie},
#'   and \code{valore} containing the spline-detrended individual tree-ring series.}
#'   \item{\code{master}}{A \code{data.frame} with columns \code{Anno} and \code{TRW}
#'   containing the biweight robust mean master chronology.}
#'   \item{\code{spei}}{A \code{data.frame} with columns \code{Anno} and \code{SPEI}
#'   containing scaled SPEI values matched to the common RWI time period.}
#'   \item{\code{df}}{A \code{data.frame} of identified negative pointer years with columns
#'   \code{year}, \code{Methods} (\code{"BSGC"} or \code{"IT"}), and \code{y_plot} for graphic positioning.}
#'   \item{\code{corr_df}}{A \code{data.frame} summarizing maximum cross-correlation values (\code{Corr})
#'   and corresponding lags (\code{Lag}) between the master chronology, SPEI, and pointer year indicators.}
#' }
#'
#' @examples
#' \dontrun{
#' library(dplR)
#' library(pointRes)
#' library(dplyr)
#' library(tidyr)
#'
#' # Example setup using dplR sample data
#' data(ca533)
#'
#' # Mock SPEI data frame for the same period
#' years <- as.numeric(rownames(ca533))
#' spei_sample <- data.frame(
#'   Anno = years,
#'   SPEI = runif(length(years), -2.5, 2.5)
#' )
#'
#' # Process data for plotting
#' plot_data <- plot_preparation(RWI = ca533, spei_df = spei_sample)
#'
#' # Inspect prepared master chronology
#' head(plot_data$master)
#' }
#'
#' @export
plot_preparation <- function(RWI, spei_df) {

  ccf_multiple <- function(s1, altre_serie, nomi) {
    correlazioni_max <- c()
    lag_max <- c()
    for (s_corrente in altre_serie) {
      risultato_ccf <- ccf(s1, s_corrente, plot = FALSE, na.action = na.omit)
      valori_acf <- as.numeric(risultato_ccf$acf)
      valori_lag <- as.numeric(risultato_ccf$lag)
      indice_max <- which.max(abs(valori_acf))
      correlazioni_max <- c(correlazioni_max, round(valori_acf[indice_max], 3))
      lag_max <- c(lag_max, valori_lag[indice_max])
    }
    matrice_risultati <- matrix(c(correlazioni_max, lag_max), nrow = 2, byrow = TRUE)
    df_finale <- as.data.frame(matrice_risultati)
    rownames(df_finale) <- c("Corr", "Lag")
    colnames(df_finale) <- nomi
    return(df_finale)
  }

  rwl <- common_interval(RWI)$rwl
  rwl_bsgc <- common.interval(RWI, make.plot = FALSE)
  det <- detrend(rwl, method = "Spline", nyrs = 30, f = 0.5)
  it <- interval.trend(rwl, trend.thresh = 0, IT.thresh = 95, make.plot = FALSE)
  BSGC <- bsgc(rwl_bsgc, make.plot = FALSE, maxlag = 5)
  n_rep <- length(names(BSGC$neg))
  df_bsgc <- data.frame(year = names(BSGC$neg), Methods = rep("BSGC", n_rep))
  pointer_years_df <- it$out[it$out$nature < 0, 1, drop = FALSE]
  pointer_years_df$Methods <- rep("IT", nrow(pointer_years_df))
  df_py <- rbind(df_bsgc, pointer_years_df)
  df_py$year <- as.numeric(df_py$year)
  cron <- chron(det, biweight = TRUE, prewhiten = TRUE)
  det <- tibble::rownames_to_column(det, var = "Anno")
  det_long <- det %>%
    tidyr::pivot_longer(
      cols = -c(Anno),
      names_to = "serie",
      values_to = "valore"
    ) %>% as.data.frame()
  det_long$Anno <- as.numeric(det_long$Anno)
  det_long <- arrange(det_long, serie)

  cron_df <- data.frame(
    Anno = as.numeric(rownames(cron)),
    TRW = as.numeric(cron[, 1])
  )
  min_old <- min(spei_df[, 2], na.rm = TRUE)
  max_old <- max(spei_df[, 2], na.rm = TRUE)
  spei_df[, 2] <- ((spei_df[, 2] - min_old) / (max_old - min_old)) + 1.5
  colnames(spei_df) <- c("Anno", "SPEI")
  spei_df <- spei_df[spei_df$Anno %in% time(rwl), ]
  max_y <- max(c(spei_df$SPEI, cron_df$TRW), na.rm = TRUE) * 1.05
  min_y <- min(c(spei_df$SPEI, cron_df$TRW), na.rm = TRUE) * 0.8
  df_py <- df_py %>% mutate(
    y_plot = case_when(Methods == "BSGC" ~ 2.60, Methods == "IT" ~ 2.75)
  )
  corr <- left_join(cron_df, spei_df, by = "Anno")

  if (nrow(df_py) != 0) {
    for (py_meth in unique(df_py$Methods)) {
      num <- ncol(corr) + 1
      corr[, num] <- rep(0, nrow(corr))
      corr[corr$Anno %in% df_py[df_py$Methods == py_meth, "year"], num] <- 1
      colnames(corr)[num] <- py_meth
    }
    s1 <- corr[, 2]
    if (ncol(corr) < 5) {
      altre_serie <- list(s2 = corr[, 3], s3 = corr[, 4])
    } else {
      altre_serie <- list(s2 = corr[, 3], s3 = corr[, 4], s4 = corr[, 5])
    }
    nomi <- colnames(corr)[c(3:ncol(corr))]
  } else {
    altre_serie <- list(s2 = corr[, 3])
    nomi <- colnames(corr)[c(3)]
  }
  df_corr <- ccf_multiple(s1, altre_serie, nomi)
  df_corr$Misure <- c(max(det_long$Anno) - 5, min_y)
  return(list(rwi = det_long, master = cron_df, spei = spei_df, df = df_py, corr_df = df_corr))
}
