#' Find the Longest Common Interval in an RWL Object
#'
#' Identifies and extracts the longest continuous time period within an RWL
#' dataset that meets constraints on time range (1901-2025), minimum duration (12 years),
#' and minimum sample depth (tree-based evaluation prioritized, falling back to simple series count >= 6).
#'
#' @param rwl A data frame or `rwl` object where rows represent years (with year numbers
#'   as row names) and columns represent individual measurement series.
#' @param min.trees Numeric. Minimum number of active series/trees required per year (default is 6).
#' @param start.yr Numeric. The earliest year to consider for the common interval (default is 1901).
#' @param end.yr Numeric. The latest year to consider for the common interval (default is 2025).
#' @param min.period.len Numeric. Minimum required length (in years) for a continuous valid period
#'   to be extracted (default is 12).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{found}{Logical. `TRUE` if a continuous interval matching all criteria was found; `FALSE` otherwise.}
#'   \item{rwl}{A subsetted data frame of the input `rwl` containing only the years within
#'     the common interval and removing empty columns. Returns `NULL` if no valid interval was found.}
#' }
#'
#' @export
common_interval <- function(rwl, min.trees = 6, start.yr = 1901,
                            end.yr = 2025, min.period.len = 12) {

  if (!is.data.frame(rwl)) stop("'rwl' must be a data.frame")

  yrs <- as.numeric(row.names(rwl))
  rwlNotNA <- !is.na(rwl)

  use_tree_count <- FALSE
  samp.depth <- numeric(nrow(rwl))

  ids <- tryCatch({
    dplR::autoread.ids(rwl)
  }, error = function(e) NULL)

  if (!is.null(ids) && "tree" %in% colnames(ids)) {
    tree_vec <- ids$tree
    samp.depth.trees <- apply(rwlNotNA, 1, function(row_mask) {
      length(unique(tree_vec[row_mask]))
    })

    potential_trees <- samp.depth.trees >= min.trees & yrs >= start.yr & yrs <= end.yr
    runs_trees <- rle(potential_trees)

    # Se tramite alberi unici troviamo almeno un blocco valido della lunghezza minima
    if (any(runs_trees$values & runs_trees$lengths >= min.period.len)) {
      use_tree_count <- TRUE
      samp.depth <- samp.depth.trees
    }
  }

  if (!use_tree_count) {
    samp.depth <- rowSums(rwlNotNA)
  }

  potential.years <- samp.depth >= min.trees & yrs >= start.yr & yrs <= end.yr

  runs <- rle(potential.years)
  if (!any(runs$values)) {
    return(list(found = FALSE, rwl = NULL))
  }

  max.run.len <- max(runs$lengths[runs$values == TRUE], 0)
  output <- list(found = FALSE, rwl = NULL)

  if (max.run.len >= min.period.len) {
    end.indices <- cumsum(runs$lengths)
    start.indices <- c(1, end.indices[-length(end.indices)] + 1)

    valid.runs <- which(runs$values == TRUE & runs$lengths >= min.period.len)

    # Seleziona l'intervallo continuo più lungo
    target.run <- valid.runs[which.max(runs$lengths[valid.runs])]
    row.indices <- start.indices[target.run]:end.indices[target.run]

    res <- rwl[row.indices, , drop = FALSE]

    # Mantiene solo le colonne che hanno dati nell'intervallo ritagliato
    keep.col.mask <- colSums(!is.na(res)) > 0

    output$found <- TRUE
    output$rwl <- res[, keep.col.mask, drop = FALSE]
  }

  return(output)
}

#' Process RWL Files to Build Chronologies and Availability Summary
#'
#' Batch-processes all Tucson-format `.rwl` files within a target directory.
#' For each file, it performs detrending, common interval extraction, and chronology
#' building, while also computing site-level sample depth coverage (threshold-based
#' vs. standard common interval). Returns individual chronology data frames, a summary
#' data frame of site availability over time, and execution statistics.
#'
#' @param dir_file A character string specifying the local directory path containing
#'   the `.rwl` files to be processed.
#' @param start_year Numeric. The starting year for the global time window analysis (default is 1901).
#' @param end_year Numeric. The ending year for the global time window analysis (default is 2025).
#'
#' @return A list containing three elements:
#' \describe{
#'   \item{cronologie}{A named list of data frames (one per successfully processed file).
#'     Each data frame contains two columns: \code{anno} (Year) and \code{rwl} (Chronology value).}
#'   \item{df_disponibilita}{A long-format data frame suitable for plotting, containing
#'     \code{Year}, \code{threshold_Type}, and \code{Relative_Frequency} (percentage of available sites).}
#'   \item{stand_read}{Integer. The total number of `.rwl` files successfully read and parsed.}
#' }
#'
#' @importFrom dplR read.rwl detrend chron common.interval
#' @importFrom dplyr %>% filter
#' @importFrom tidyr pivot_longer
#' @export
elabora_cronologie <- function(dir_file, start_year = 1901, end_year = 2025) {

  nomi_file <- list.files(dir_file, pattern = "\\.rwl$", full.names = FALSE)
  totale_stands <- length(nomi_file)
  results_base <- list()
  anni_globali <- start_year:end_year
  storage_soglia <- setNames(numeric(length(anni_globali)), anni_globali)
  storage_common <- setNames(numeric(length(anni_globali)), anni_globali)
  siti_letti_con_successo <- 0

  # --- CICLO SUI FILE ---
  for (i in seq_along(nomi_file)) {
    stand <- nomi_file[i]
    file_path <- file.path(dir_file, stand)

    # Soppressione dei Warning di lettura per pulizia e velocità
    RWL <- tryCatch({
      suppressWarnings(
        capture.output(
          res <- dplR::read.rwl(file_path, format = "tucson"),
          type = "output"
        )
      )
      res
    }, error = function(e) NULL)

    if (!is.null(RWL) && ncol(RWL) > 0) {

      # A) DENDROCHRONOLOGICAL ANALYSIS
      cron_df <- tryCatch({
        DET <- dplR::detrend(RWL, method = "Spline", nyrs = 30, f = 0.5)

        # Uso esplicito della nuova funzione con parametri allineati (min 6 serie, min 12 anni)
        det <- common_interval(DET, min.trees = 6, start.yr = start_year, end.yr = end_year, min.period.len = 12)
        if (!det$found) stop("Common Interval not found")

        cron <- dplR::chron(det$rwl, biweight = TRUE, prewhiten = TRUE)
        data.frame(
          anno = as.numeric(rownames(cron)),
          rwl  = as.numeric(cron[, 1])
        )
      }, error = function(e) NULL)

      if (!is.null(cron_df)) {
        results_base[[stand]] <- cron_df
      }

      # B) SAMPLE DEPTH (Logica Semplificata e Allineata alla Nuova common_interval)
      siti_letti_con_successo <- siti_letti_con_successo + 1
      anni <- as.numeric(rownames(RWL))

      # CONTEGGIO SEMPLICE COLONNE (senza autoread.ids)
      count <- rowSums(!is.na(RWL))

      # Standard dplR common interval per il secondo parametro del grafico
      ci_rwl <- tryCatch(dplR::common.interval(RWL, make.plot = FALSE), error = function(e) NULL)
      anni_ci <- if (!is.null(ci_rwl)) as.numeric(rownames(ci_rwl)) else NULL

      dati_sito <- data.frame(
        Anno = anni,
        Supera_Soglia = count >= 6,  # SOGLIA FISSA A 6 SERIE
        In_Common_Interval = anni %in% anni_ci
      )

      # C) AGGIORNAMENTO STORAGE
      anni_ok_soglia <- as.character(dati_sito$Anno[dati_sito$Supera_Soglia &
                                                      dati_sito$Anno >= start_year &
                                                      dati_sito$Anno <= end_year])
      anni_ok_common <- as.character(dati_sito$Anno[dati_sito$In_Common_Interval &
                                                      dati_sito$Anno >= start_year &
                                                      dati_sito$Anno <= end_year])

      storage_soglia[anni_ok_soglia] <- storage_soglia[anni_ok_soglia] + 1
      storage_common[anni_ok_common] <- storage_common[anni_ok_common] + 1
    }

    if (i %% 100 == 0) cat("Processati:", i, "/", totale_stands, "\n")
  }

  # --- COSTRUZIONE DATAFRAME PER GRAFICO ---
  df_disponibilita <- data.frame(
    Year = anni_globali,
    Threshold_6_Tree = (storage_soglia / siti_letti_con_successo) * 100,
    Common_Interval_100 = (storage_common / siti_letti_con_successo) * 100
  )

  # Utilizzo namespaces espliciti di dplyr e tidyr per evitare errori
  df_disponibilita <- df_disponibilita %>%
    dplyr::filter(Threshold_6_Tree > 0 | Common_Interval_100 > 0) %>%
    tidyr::pivot_longer(cols = -Year, names_to = "threshold_Type", values_to = "Relative_Frequency")

  message("Process completed.\nTotal Stand: ", totale_stands,
          "\nSuccessfully read: ", siti_letti_con_successo, "\n")

  return(list(
    cronologie = results_base,
    df_disponibilita = df_disponibilita,
    stand_read = siti_letti_con_successo
  ))
}
