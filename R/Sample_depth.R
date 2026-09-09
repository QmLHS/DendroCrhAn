#' Calculate Sample Depth and Determine Thresholds for RWL Objects
#'
#' Evaluates the sample depth over time for a Ring-Width Length (RWL) object.
#' It attempts to extract tree identities using `autoread.ids`. If tree IDs are
#' successfully retrieved, sample depth is calculated in terms of unique trees per year
#' with a recommended threshold of 6 trees. If ID extraction fails, it falls back to
#' counting active series per year with a threshold of 12 series.
#'
#' @param RWL A data frame or `rwl` object (typically from the `dplR` package)
#'   where rows represent years (with year numbers as row names) and columns represent
#'   individual measurement series.
#'
#' @return A list containing three elements:
#' \describe{
#'   \item{df}{A data frame with columns \code{Anno} (Year), \code{n_presenti} (number of active series),
#'     and \code{n_valore} (number of trees or series depending on ID availability).}
#'   \item{soglia}{A numeric value indicating the suggested minimum sample depth threshold (6 for trees, 12 for series).}
#'   \item{label}{A character string specifying the y-axis label suitable for plotting
#'     (\code{"Number of trees"} or \code{"Number of series"}).}
#' }
#'
#' @importFrom dplR autoread.ids
#' @export
get_sample_depth <- function(RWL) {
  # 1. Prova a leggere gli ID degli alberi
  ids <- tryCatch(dplR::autoread.ids(RWL), error = function(e) NULL)

  # 2. Prepariamo il conteggio anno per anno
  coverage_df <- data.frame(
    Anno = as.numeric(rownames(RWL)),
    n_presenti = rowSums(!is.na(RWL))
  )

  # 3. Logica della soglia
  if (!is.null(ids)) {
    # Se abbiamo gli ID, la soglia è 6 alberi.
    n_alberi_per_anno <- sapply(1:nrow(RWL), function(i) {
      length(unique(ids$tree[which(!is.na(RWL[i, ]))]))
    })
    coverage_df$n_valore <- n_alberi_per_anno
    soglia_val <- 6
    label_y <- "Number of trees"
  } else {
    # Fallback: se autoread.ids fallisce, usiamo 12 serie
    soglia_val <- 12
    coverage_df$n_valore <- coverage_df$n_presenti
    label_y <- "Number of series"
  }

  return(list(df = coverage_df, soglia = soglia_val, label = label_y))
}
