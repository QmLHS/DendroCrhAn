#' Construct a Distributed Lag Matrix from Monthly Climate Data
#'
#' @description
#' Generates a multi-year distributed lag matrix for a specified target variable
#' from monthly climate time series. For each year in the series, it extracts a backwards
#' sequence of \code{lag} months ending at a chosen reference month (\code{mese_n}), capturing
#' seasonal and multi-year legacy effects on ecological processes (e.g., tree growth).
#'
#' @param df A \code{data.frame} containing monthly climate data. Must include a \code{data} column
#' with date strings formatted as \code{"mmm-YYYY"} (e.g., \code{"jan-1901"}) and a numeric variable column.
#' @param val Character string specifying the name of the climate variable column to extract (e.g., \code{"temp"}, \code{"pre"}, \code{"bal"}).
#' @param mese_n Integer (1 to 12). The target reference month that terminates each lag sequence
#' (e.g., \code{9} for September, \code{2} for February). Default is \code{9}.
#' @param lag Integer. Total number of monthly lags to include in the sequence going backwards from \code{mese_n}. Default is \code{12}.
#'
#' @return A numeric \code{matrix} where each row represents a target year (derived from \code{date_finali})
#' and each column represents a monthly lag step going backwards from \code{mese_n}, named in the format \code{"month_M_year_Y"}.
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(stringr)
#' library(lubridate)
#'
#' # Create sample monthly time series (2 years)
#' dates <- paste0(rep(c("jan", "feb", "mar", "apr", "may", "jun",
#'                       "jul", "aug", "sep", "oct", "nov", "dec"), 2),
#'                 "-", rep(c(2019, 2020), each = 12))
#'
#' sample_meteo <- data.frame(
#'   data = dates,
#'   temp = runif(24, 10, 25)
#' )
#'
#' # Build a 12-month lag matrix ending in September (month 9)
#' lag_mat <- lag_mese(df = sample_meteo, val = "temp", mese_n = 9, lag = 12)
#' head(lag_mat)
#' }
#'
#' @export
lag_mese <- function(df, val, mese_n=9, lag = 12) {
  mesi_ita <- c("gen", "feb", "mar", "apr", "mag", "giu",
                "lug", "ago", "set", "ott", "nov", "dic")
  # 1. Preparazione (Ordine cronologico per estrarre correttamente)
  df_clean <- df %>%
    mutate(
      mese_str = str_split_fixed(data, "-", 2)[,1],
      anno_val = as.numeric(str_split_fixed(data, "-", 2)[,2]),
      mese_num = match(mese_str, mesi_ita),
      data_obj = make_date(year = anno_val, month = mese_num, day = 1)
    ) %>%
    arrange(data_obj)
  # 2. Identifichiamo i Settembre (punti di partenza di ogni riga)
  punti_settembre <- df_clean %>%
    filter(month(data_obj) == mese_n) %>%
    pull(data_obj)
  # 3. Estrazione a ritroso
  estrai_sequenza <- function(data_rif) {
    res <- df_clean %>%
      filter(data_obj <= data_rif) %>% # Prendo tutto ciò che viene prima del Settembre di riferimento
      slice_tail(n = lag) %>%     # Prendo gli ultimi n mesi
      arrange(desc(data_obj)) %>%      # Inverto l'ordine: Settembre, Agosto, Luglio...
      pull(val)
    if(length(res) == lag) return(res) else return(NULL)}
  # 4. Creazione lista e pulizia anni incompleti
  lista_risultati <- lapply(punti_settembre, estrai_sequenza)
  indici_validi <- sapply(lista_risultati, function(x) !is.null(x))
  lista_finale <- lista_risultati[indici_validi]
  date_finali <- punti_settembre[indici_validi]
  # 5. Costruzione Matrice
  matrice_finale <- do.call(rbind, lista_finale)
  # Nomi righe: l'anno del settembre scelto
  rownames(matrice_finale) <- year(date_finali)
  # Nomi colonne: mostriamo il ritardo (Lag)
  mesi_ripetuti <- rep(c(mese_n:1, 12:(mese_n+1)), length.out = lag)
  offset_anni <- floor((0:(lag - 1) + 3) / 12)
  colnames(matrice_finale) <- paste0("month_", mesi_ripetuti, "_year_", offset_anni)
  return(matrice_finale)}

#' Prepare Data Structures for Mixed-Effects Distributed Lag Non-linear Models (DLNM)
#'
#' @description
#' Filters tree-ring stands based on climate or species groups and prepares synchronized datasets
#' for mixed-effects DLNM analysis. For each selected stand, it extracts tree-ring chronologies,
#' constructs a multi-year monthly climate lag matrix, computes annual climate control anomalies,
#' and binds all valid sites (with >10 common years) into a unified data frame and lag matrix.
#'
#' @param df_tax A \code{data.frame} containing taxonomic and climate group metadata per stand.
#' Must include \code{stand_id} and \code{pixel_id} in the first two columns, followed by climate group
#' and species group in columns 3 and 4.
#' @param var_meteo Character string specifying the target climate variable to construct the lag matrix
#' (e.g., \code{"temp"}, \code{"pre"}, or \code{"bal"} for water balance).
#' @param var_ctrl Character string specifying the climate control variable to be processed via \code{scostamento()}.
#' @param month_n Integer (1 to 12). Target end-month for the backward lag matrix calculation in \code{lag_mese()}.
#' @param g_clima Optional character string to filter stands by climate group (\code{gruppo_climatico}).
#' If \code{NULL}, no climate group filter is applied.
#' @param g_specie Optional character string to filter stands by species group (\code{gruppo_specie}).
#' If \code{NULL}, no species group filter is applied.
#' @param lag_n Integer. Number of backwards monthly lags to include in the climate matrix. Default is \code{24}.
#'
#' @return A named \code{list} containing three elements:
#' \describe{
#'   \item{\code{df}}{A unqiue \code{data.frame} merging tree-ring chronologies, annual control variables,
#'   and stand identifiers (\code{ID}) across all selected sites.}
#'   \item{\code{mat}}{A combined numeric \code{matrix} binding the monthly climate lag matrices for all selected sites.}
#'   \item{\code{var_ctrl_name}}{Character string echoing the input \code{var_ctrl} variable name.}
#' }
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#'
#' # Mock taxonomic metadata
#' tax_metadata <- data.frame(
#'   stand_id = c("STD_01", "STD_02"),
#'   pixel_id = c("PXL_A", "PXL_B"),
#'   g_clim = c("Med", "Med"),
#'   g_spec = c("Fagus", "Pinus")
#' )
#'
#' # Assuming external objects 'crono' and 'clima' are defined in the environment:
#' # dlnm_input <- prepara_dlnm_misto(
#' #   df_tax = tax_metadata,
#' #   var_meteo = "temp",
#' #   var_ctrl = "pre",
#' #   month_n = 9,
#' #   g_clima = "Med",
#' #   g_specie = NULL,
#' #   lag_n = 24
#' # )
#' }
#'
#' @export
prepara_dlnm_misto <- function(df_tax, var_meteo, var_ctrl, month_n,
                               g_clima, g_specie, lag_n = 24, crono, clima) {
  raga <- df_tax
  colnames(raga)[c(2,3)] <- c("gruppo_climatico","gruppo_specie")
  if (!is.null(g_clima)) {raga <- raga %>% filter(gruppo_climatico == g_clima)}
  if (!is.null(g_specie)) {raga <- raga %>% filter(gruppo_specie == g_specie)}
  stand_selezionati <- raga %>% select(stand_id, pixel_id)
  big_df <- data.frame()
  big_X_meteo <- matrix(nrow = 0, ncol = lag_n)
  for(i in 1:nrow(stand_selezionati)) {
    s_id <- stand_selezionati$stand_id[i]
    p_id <- stand_selezionati$pixel_id[i]
    s_id_comp <- paste0(s_id,".rwl")
    y_stand <- crono[[s_id_comp]]
    df_meteo <- clima[[p_id]]
    if(var_meteo == "bal") {df_meteo$bal <- df_meteo$pre - df_meteo$pet}
    X_mat <- lag_mese(df_meteo, val = var_meteo, mese_n = month_n,lag = lag_n)
    df_per_scostamento <- df_meteo[, c("data", var_ctrl)]
    ctrl_annuale <- as.data.frame(scostamento(df_per_scostamento))
    colnames(ctrl_annuale) <- c("anno", var_ctrl)
    anni_comuni <- intersect(rownames(X_mat), as.character(y_stand$anno))
    anni_comuni <- intersect(anni_comuni, as.character(ctrl_annuale$anno))
    if(length(anni_comuni) > 10) {
      df_tmp <- y_stand %>%
        filter(anno %in% anni_comuni) %>%
        left_join(ctrl_annuale %>% mutate(anno = as.numeric(anno)), by = "anno") %>%
        mutate(ID = s_id)
      X_tmp <- X_mat[anni_comuni, ]
      big_df <- rbind(big_df, df_tmp)
      big_X_meteo <- rbind(big_X_meteo, X_tmp)}}
  return(list(df = big_df, mat = big_X_meteo,
              var_ctrl_name = var_ctrl))}

#' Add an Autoregressive Growth Term and Synchronize Climate Lag Matrix
#'
#' @description
#' Computes a 1-year lagged tree-ring growth variable (\code{rwl_prev}) within each stand/tree
#' sequence to account for temporal autocorrelation. Automatically prevents boundary-crossing
#' leakage between distinct tree/stand IDs, filters out incomplete initial observations (\code{NA}s),
#' aligns the climate lag matrix accordingly, and performs an internal verification check to confirm temporal accuracy.
#'
#' @param df A \code{data.frame} containing tree-ring data. Must include columns \code{ID} (stand or tree identifier),
#' \code{anno} (year), and \code{rwl} (ring-width index or growth measurement).
#' @param mat A numeric \code{matrix} representing the climate lag data, whose rows directly correspond
#' to the rows of \code{df}.
#'
#' @return A named \code{list} containing two aligned elements:
#' \describe{
#'   \item{\code{df}}{A \code{data.frame} filtered for complete cases with the added \code{rwl_prev} autoregressive column.}
#'   \item{\code{mat}}{A numeric \code{matrix} synchronized with the filtered data frame rows.}
#' }
#'
#' @examples
#' \dontrun{
#' # Sample tree-ring data for two stands
#' sample_df <- data.frame(
#'   ID = rep(c("S01", "S02"), each = 5),
#'   anno = rep(2010:2014, 2),
#'   rwl = c(1.2, 1.0, 0.9, 1.1, 1.3, 0.8, 0.7, 0.9, 1.0, 1.1)
#' )
#'
#' # Mock matching climate matrix (10 rows)
#' sample_mat <- matrix(rnorm(10 * 12), nrow = 10, ncol = 12)
#'
#' # Add autoregressive term and synchronize
#' res <- termine_auto(df = sample_df, mat = sample_mat)
#' head(res$df)
#' }
#'
#' @export
termine_auto <- function(df, mat){
  indici_teste <- c(1, which(df$ID[-1] != df$ID[-nrow(df)]) + 1)
  rwl_slittato <- c(NA, df$rwl[-nrow(df)])
  rwl_slittato[indici_teste] <- NA
  df_aggiornato <- df
  df_aggiornato$rwl_prev <- rwl_slittato
  righe_da_tenere <- which(!is.na(df_aggiornato$rwl_prev))
  df_final <- df_aggiornato[righe_da_tenere,]
  mat_final <- mat[righe_da_tenere,]
  test_idx <- sample(2:nrow(df_final), 1)
  cor_id <- df_final$ID[test_idx]
  cor_anno <- df_final$anno[test_idx]
  valore_reale_anno_prec <- df$rwl[df$ID == cor_id & df$anno == (cor_anno - 1)]
  if(valore_reale_anno_prec == df_final$rwl_prev[test_idx]) {
    cat("✅ Allineamento Verificato!")
  } else {
    stop("❌ ERRORE DI ALLINEAMENTO!")}
  return(list(df=df_final,mat=mat_final))}

#' Verify Synchronized Temporal Alignment across Chronologies, Data Frames, and Lag Matrices
#'
#' @description
#' Performs a multi-level audit to verify exact temporal synchronization among original site
#' chronologies (\code{crono}), the extracted DLNM data frame (\code{df}), and the corresponding climate lag matrix (\code{mat}).
#' Checks boundary alignment for start and end years across all sites, accounts for initial
#' lag truncation offsets (e.g., historical climate start thresholds or autoregressive term adjustments),
#' and prints a detailed alignment diagnostic report.
#'
#' @param output_dlnm A named \code{list} produced by DLNM data preparation functions (e.g., \code{prepara_dlnm_misto}
#' or \code{termine_auto}). Must contain at least \code{df} (data frame with stand \code{ID} and \code{anno})
#' and \code{mat} (climate lag matrix with year row names).
#' @param report Logical. If \code{TRUE} (default), detailed site-by-site diagnostic reports and step-transition boundaries
#' are printed to the console. If \code{FALSE}, only final summary counts are displayed.
#' @param rwl_prev Logical. Set to \code{TRUE} if an autoregressive prior-year growth term was added
#' (e.g., via \code{termine_auto}), shifting the expected start year forward by +1 year. Default is \code{FALSE}.
#'
#' @return Invisibly returns \code{NULL}. Outputs verification messages, warnings for temporal discrepancies,
#' and summary counts (total stands evaluated and discrepancy counts) directly to the console.
#'
#' @examples
#' \dontrun{
#' # Mock output list matching function expectations
#' mock_df <- data.frame(
#'   ID = rep("STD_01", 5),
#'   anno = 2005:2009
#' )
#' mock_mat <- matrix(1:20, nrow = 5, ncol = 4)
#' rownames(mock_mat) <- 2005:2009
#'
#' mock_output <- list(df = mock_df, mat = mock_mat)
#'
#' # Assuming external chronology object 'crono' exists in the environment:
#' # crono <- list(STD_01 = data.frame(anno = 2005:2009, rwl = runif(5)))
#' # controllo(output_dlnm = mock_output, report = TRUE, rwl_prev = FALSE)
#' }
#'
#' @export
controllo <- function(output_dlnm, report = TRUE, rwl_prev = FALSE) {
  df_out <- output_dlnm$df
  mat_out <- output_dlnm$mat
  stands_out <- unique(as.character(df_out$ID))
  n_cumulato <- 0
  cont_errori <- 0
  cont_serie <- length(stands_out)
  cat("============================================================\n")
  cat("   VERIFICA COERENZA TRIPLA: CRONO <-> DF <-> MATRICE\n")
  cat("============================================================\n\n")
  for (i in seq_along(stands_out)) {
    s_id <- stands_out[i]
    s_id_c <- paste0(s_id,".rwl")
    df_sub <- df_out[df_out$ID == s_id, ]
    n_righe <- nrow(df_sub)
    idx_i <- n_cumulato + 1
    idx_f <- n_cumulato + n_righe
    anno_i_out <- as.numeric(df_sub$anno[1])
    anno_f_out <- as.numeric(df_sub$anno[n_righe])
    anno_i_mat <- as.numeric(rownames(mat_out)[idx_i])
    anno_f_mat <- as.numeric(rownames(mat_out)[idx_f])
    anno_i_crono <- as.numeric(crono[[s_id_c]]$anno[1])
    anno_f_crono <- as.numeric(crono[[s_id_c]]$anno[nrow(crono[[s_id_c]])])
    if(report){cat(sprintf("SITO: [%s] (Righe: %d)\n", s_id, n_righe))
      cat(sprintf("  • Posizione Matrice: righe %d a %d\n", idx_i, idx_f))
      cat(sprintf("  • Range Anni Matrice: %d --> %d\n", anno_i_mat, anno_f_mat))
      cat(sprintf("  • Range Anni Dataframe: %d --> %d\n", anno_i_out, anno_f_out))
      cat(sprintf("  • Range Anni Crono:   %d --> %d\n", anno_i_crono, anno_f_crono))}
    anno_x <- ifelse(ncol(mat_out)>9 & ncol(mat_out)<22, 1902, 1903)
    cont_i_crono <- ifelse(anno_i_crono>1902, anno_i_crono, anno_x)
    if(rwl_prev) {cont_i_crono <- cont_i_crono+1}
    cont_f_crono <- anno_f_crono
    inizio_ok <- (anno_i_out == cont_i_crono && anno_i_mat == cont_i_crono)
    fine_ok <- (anno_f_out == anno_f_crono && anno_f_mat == anno_f_crono)
    if (inizio_ok && fine_ok) {
      if(report){cat("  ✅ Allineamento perfetto (o deroga 1901-03 rispettata)\n")}
    } else { cont_errori <- cont_errori+1
    if(report){cat("  ⚠️ ATTENZIONE: Discrepanza rilevata!\n")
      if(!inizio_ok) cat(sprintf("     - Errore INIZIO: Crono(%d) vs Output(%d)\n", anno_i_crono, anno_i_out))
      if(!fine_ok)   cat(sprintf("     - Errore FINE:   Crono(%d) vs Output(%d)\n", anno_f_crono, anno_f_out))}}
    if (i < length(stands_out)) {
      anno_f_attuale_mat <- rownames(mat_out)[idx_f]
      anno_i_prossimo_mat <- rownames(mat_out)[idx_f + 1]
      prossimo_id <- stands_out[i+1]
      if(report){cat(sprintf("  • GRADINO: Fine [%s]: %s | Inizio [%s]: %s\n",
                             s_id, anno_f_attuale_mat, prossimo_id, anno_i_prossimo_mat))
        cat("    🔗 Salto temporale rilevato (corretto per cambio sito).\n")
        cat("------------------------------------------------------------\n")}}
    n_cumulato <- idx_f}
  cat(sprintf("NUMERO SITI: %s\n", cont_serie))
  cat(sprintf("NUMERO DISCREPANSE: %s\n", cont_errori))}
