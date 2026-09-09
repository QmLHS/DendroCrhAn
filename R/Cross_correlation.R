#' Compute and Plot Significant Species-Specific Cross-Correlations across Lags
#'
#' @description
#' Calculates cross-correlation functions (CCF) between annual climate anomalies (SPEI)
#' and tree-ring width index (RWI) series for selected species groups across specified time lags.
#' Filters for statistically significant positive correlations ($r > 1.96/sqrt(N)$),
#' computes species-level response frequencies, and builds a dual-axis \code{ggplot2} visualization.
#'
#' @param df_tax A \code{data.frame} containing stand metadata. Must include columns \code{stand_id},
#' \code{pixel_id}, and \code{gruppo_specie}.
#' @param sel_orders Optional character vector specifying target species groups (\code{gruppo_specie})
#' to filter from \code{df_tax}. If \code{NULL}, all species groups are processed.
#' @param max_lag Integer. Maximum time lag (in years) to evaluate in the cross-correlation function.
#' Default is \code{4}.
#'
#' @return A named \code{list} containing three elements:
#' \describe{
#'   \item{\code{dati}}{A \code{data.frame} of significant positive correlations per stand and lag
#'   with columns \code{gruppo_specie}, \code{stand_id}, \code{lag}, and \code{cor}.}
#'   \item{\code{freq}}{A \code{data.frame} summarizing site-level response frequencies per species
#'   and lag, including absolute (\code{freq_ass}) and relative percentage (\code{freq_rel}) frequencies.}
#'   \item{\code{grafico}}{A \code{ggplot} object displaying jittered point-linerange correlation
#'   values on the primary y-axis and relative response frequency bars on the secondary y-axis, faceted by species.}
#' }
#'
#' @import ggplot2
#' @importFrom dplyr %>% filter group_by group_modify ungroup summarise n_distinct n left_join mutate bind_rows
#' @importFrom paletteer paletteer_d scale_color_paletteer_d
#' @importFrom stats ccf complete.cases
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(ggplot2)
#' library(paletteer)
#'
#' # Mock taxonomy data
#' tax_data <- data.frame(
#'   stand_id = c("S01", "S02", "S03"),
#'   pixel_id = c("P01", "P02", "P03"),
#'   gruppo_specie = c("Fagus sylvatica", "Fagus sylvatica", "Picea abies")
#' )
#'
#' # Execute function for specified species
#' res <- correlazione_significativa_lag(
#'   df_tax = tax_data,
#'   sel_orders = c("Fagus sylvatica"),
#'   max_lag = 4
#' )
#'
#' # Print plot
#' print(res$grafico)
#' }
#'
#' @export
correlazione_significativa_lag <- function(df_tax, sel_orders, crono, clima, max_lag = 4) {
  # 1. Filtro per Genere
  df_lavoro <- df_tax
  if (!is.null(sel_orders)) {
    df_lavoro <- df_lavoro %>% filter(gruppo_specie %in% sel_orders)}
  # 2. Funzione interna: Estrae solo i dati POSITIVI e SIGNIFICATIVI
  get_sig_cor_data <- function(id_stand, id_pixel) {
    id_s <- paste0(id_stand, ".rwl"); id_p <- as.character(id_pixel)
    if (!(id_s %in% names(crono)) | !(id_p %in% names(clima))) return(NULL)
    tryCatch({
      df_a <- crono[[id_s]]; df_b <- clima[[id_p]]
      df_spei_t <- as.data.frame(scostamento(df_b[, c("data", "spei")]))
      df_m <- merge(df_a, df_spei_t, by = "anno")
      if (nrow(df_m) < 10) return(NULL)
      ccf_res <- ccf(df_m$sum_dif_spei, df_m$rwl, lag.max = max_lag, plot = FALSE, na.action = na.pass)
      n_obs <- sum(complete.cases(df_m))
      soglia <- 1.96 / sqrt(n_obs)
      # MODIFICA CRUCIALE: Teniamo solo cor > soglia (esclude negative e non sig)
      return(data.frame(
        stand_id = id_s,
        lag = abs(as.numeric(ccf_res$lag)),
        cor = as.numeric(ccf_res$acf)
      ) %>% filter(as.numeric(ccf_res$lag) <= 0, cor > soglia))
    }, error = function(e) return(NULL))}
  # 3. Elaborazione Dati
  res_sig <- df_lavoro %>%
    group_by(gruppo_specie) %>%
    group_modify(~ {
      bind_rows(mapply(get_sig_cor_data, .x$stand_id, .x$pixel_id, SIMPLIFY = FALSE))
    }) %>% ungroup()
  # 4. Calcolo Frequenze (per le colonne sienna)
  n_siti_tot <- df_lavoro %>% group_by(gruppo_specie) %>% summarise(n_tot = n_distinct(stand_id), .groups = 'drop')
  df_frequenze <- res_sig %>%
    group_by(gruppo_specie, lag) %>%
    summarise(freq_ass = n(), .groups = 'drop') %>%
    left_join(n_siti_tot, by = "gruppo_specie") %>%
    mutate(freq_rel = (freq_ass / n_tot) * 100)
  # 5. Parametri Grafici
  my_scale <- 50
  pos_jitter <- position_jitter(width = 0.25, seed = 123)
  pal <- paletteer::paletteer_d("fishualize::Anisotremus_virginicus")[c(1,3)]
  # 6. Costruzione Grafico
  p <- ggplot(res_sig, aes(x = as.numeric(as.character(lag)), y = cor, color = gruppo_specie)) +
    # Barre (senza jitter)
    geom_col(data = df_frequenze,
             aes(x = as.numeric(as.character(lag)) + 0.4, y = freq_rel / my_scale),
             inherit.aes = FALSE, fill = pal[2], alpha = 0.5, width = 0.15) +
    # Linee e Punti (stesso jitter, stessi dati)
    geom_linerange(aes(ymin = 0, ymax = cor, y = cor),
                   position = pos_jitter, alpha = 0.4, linewidth = 0.4) +
    geom_point(position = pos_jitter, size = 1.5, alpha = 0.7) +
    scale_y_continuous(
      name = "Correlation Coefficient (r)",
      expand = expansion(mult = c(0, 0.05)),
      sec.axis = sec_axis(trans = ~ . * my_scale, name = "Significant Responses (%)")
    ) +
    # Definiamo i limiti qui: questo NON elimina i dati
    coord_cartesian(ylim = c(0, 0.8)) +
    scale_x_continuous(breaks = seq(0, max_lag, 1), name = "Lag (Years)") +
    paletteer::scale_color_paletteer_d("fishualize::Anisotremus_virginicus") +
    facet_wrap(~gruppo_specie, ncol = 1) +
    theme_minimal() +
    labs(
      title = "Species-Specific Patterns of Significant Cross-Correlations") +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.text.x = element_text(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.line.y.left = element_line(color = pal[1]),
      axis.text.y.left = element_text(color = pal[1]), # Numeri etichette
      axis.ticks.y.left = element_line(color = pal[1]),
      axis.title.y.left = element_text(color = pal[1], face = "bold"),
      axis.line.y.right = element_line(color = pal[2]),
      axis.text.y.right = element_text(color = pal[2]),
      axis.ticks.y.right = element_line(color = pal[2]),
      axis.title.y.right = element_text(color = pal[2],
                                        face = "bold", margin = margin(l = 10)),
      strip.text = element_text(face = "bold.italic", size = 12))
  return(list(dati = res_sig, freq = df_frequenze, grafico = p))
}
