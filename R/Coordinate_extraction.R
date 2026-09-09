#' Compare Tree-Ring File Headers with Metadata
#'
#' @description
#' Extracts and compares geographic coordinates and species information between a
#' tree-ring data file (e.g., `.rwl`) and a reference metadata file. It extracts
#' header information using `ecfh` and matches it against an ITRDB
#' metadata database.
#'
#' @param file_path A character string. The path to the tree-ring data file (e.g., a `.rwl` file).
#' @param metadata_path A character string. The path to the tab-delimited metadata file containing ITRDB codes.
#'
#' @return A named \code{list} containing three elements:
#' \itemize{
#'   \item{metadata}{A list containing the extracted metadata coordinates (Longitude, Latitude) and site info.}
#'   \item{header}{A list containing the extracted header coordinates (Y, X), species, and raw header content.}
#'   \item{comparison_report}{A \code{data.frame} summarizing the comparison results (see Details).}
#' }
#'
#' @details
#' The \code{comparison_report} data frame contains the following columns:
#' \itemize{
#'   \item{\code{stand_id}:}{ Character. The standardized lowercase identifier of the stand extracted from the filename.}
#'   \item{\code{header_text}:}{ Character. The squished, single-string representation of the original file header lines preserved from the extraction phase.}
#'   \item{\code{meta_lat}:}{ Numeric. Latitude extracted from the metadata file (\code{NA_real_} if not found).}
#'   \item{\code{meta_lon}:}{ Numeric. Longitude extracted from the metadata file (\code{NA_real_} if not found).}
#'   \item{\code{meta_species}:}{ Character. Species code extracted from the metadata file (\code{NA_character_} if not found).}
#'   \item{\code{diff_lat}:}{ Numeric. Absolute difference between header latitude (\code{coord_X}) and metadata latitude. \code{NA_real_} if data is missing.}
#'   \item{\code{diff_lon}:}{ Numeric. Absolute difference between header longitude (\code{coord_Y}) and metadata longitude. \code{NA_real_} if data is missing.}
#'   \item{\code{correspond_coord}:}{ Logical. \code{TRUE} if both \code{diff_lat} and \code{diff_lon} are strictly less than 1; \code{FALSE} otherwise; \code{NA} if inputs are missing.}
#'   \item{\emph{Other columns}:}{ Any additional columns originally returned by the internal \code{ecfh$report} data frame.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage within a package context:
#' report <- ec(
#'   file_path = "path/to/site.rwl",
#'   metadata_path = "path/to/metadata.txt"
#' )
#' print(report$comparison_report)
#' }
ec <- function(file_path, metadata_path) {
  result <- list(metadata = list(), header = list(), comparison_report = NULL)
  metadata <- read.delim(metadata_path)
  nome_stand_up <- toupper(gsub("\\.rwl$", "", basename(file_path), ignore.case = TRUE))

  # 1. Estrazione Metadati
  meta_lon <- NA_real_
  meta_lat <- NA_real_
  meta_species <- NA_character_
  if (nome_stand_up %in% metadata$ITRDB_Code) {
    coordinate <- metadata[metadata$ITRDB_Code == nome_stand_up, c("Latitude", "Longitude")]
    info_stand <- metadata[metadata$ITRDB_Code == nome_stand_up, c("Elevation", "Species")]
    meta_lon <- coordinate$Longitude[1]
    meta_lat <- coordinate$Latitude[1]
    meta_species <- info_stand$Species[1]
    result$metadata$coordinate <- c(meta_lon, meta_lat)
    result$metadata$info <- info_stand}

  # 2. Estrazione Header da ecfh
  coordinate_estratte <- ecfh(file_path, lista_specie = unique(metadata$Species))
  if (!is.null(coordinate_estratte$coord_Y) && !is.null(coordinate_estratte$coord_X)) {
    result$header$coordinate <- c(coordinate_estratte$coord_Y, coordinate_estratte$coord_X)}
  if (!is.null(coordinate_estratte$Specie)) {
    result$header$specie <- coordinate_estratte$Specie}
  result$header$head <- coordinate_estratte$header

  # 3. CONFRONTO E CREAZIONE REPORT UNIFICATO
  rep <- coordinate_estratte$report
  rep$stand_id <- tolower(nome_stand_up)
  rep <- rep[, c("stand_id", setdiff(names(rep), "stand_id"))]
  rep$meta_lat <- meta_lat
  rep$meta_lon <- meta_lon
  rep$meta_species <- meta_species

  # Estraiamo i valori assicurandoci che siano scalari singoli
  head_lat <- coordinate_estratte$coord_X[1]
  head_lon <- coordinate_estratte$coord_Y[1]

  # Usiamo isTRUE() per evitare che NA interrompa l'if
  ha_head <- !is.null(head_lat) && !is.na(head_lat)
  ha_meta <- !is.null(meta_lat) && !is.na(meta_lat)
  if (ha_head && ha_meta) {
    rep$diff_lat <- abs(head_lat - meta_lat)
    rep$diff_lon <- abs(head_lon - meta_lon)
    rep$correspond_coord <- (rep$diff_lat < 1 && rep$diff_lon < 1)
  } else {
    rep$diff_lat <- NA_real_
    rep$diff_lon <- NA_real_
    rep$correspond_coord <- NA}
  result$comparison_report <- rep
  return(result)}

#' Parse Tree-Ring File Headers for Species and Coordinates
#'
#' @description
#' Reads the first few lines of a tree-ring data file to extract species codes
#' (4-letter acronyms) and geographic coordinates. The function standardizes the
#' header text and sequentially tests multiple regex patterns to identify
#' coordinate formats (e.g., M-compact, explicit N/W characters, decimals).
#'
#' @param file_path A character string. The path to the file to be processed.
#' @param lista_specie A character vector. A reference list of valid species codes
#'   used to cross-reference and validate the extracted 4-letter species code.
#'
#' @return A named \code{list} containing five elements:
#' \describe{
#'   \item{header}{A character string containing the squished, combined text of the first 3 lines of the file (\code{NULL} if the file cannot be read).}
#'   \item{Specie}{A character string. The validated 4-letter species code if a unique match is found; \code{NULL} otherwise.}
#'   \item{coord_X}{Numeric. The parsed X coordinate (typically Longitude, or processed Westing) scaled by 100 or 1; \code{NULL} if not found.}
#'   \item{coord_Y}{Numeric. The parsed Y coordinate (typically Latitude, or processed Northing) scaled by 100 or 1; \code{NULL} if not found.}
#'   \item{report}{A \code{data.frame} containing processing logs and parsed metrics (see Details).}
#' }
#'
#' @details
#' The \code{report} data frame contains a single row with the following columns:
#' \itemize{
#'   \item{\code{header_text}:}{ Character. The squished, single-string representation of the file header lines (\code{NA_character_} if file reading failed).}
#'   \item{\code{pattern_matched}:}{ Character. The name of the specific regex pattern that matched the coordinate string (e.g., \code{"PATTERN M-COMPATTO"}, \code{"PATTERN 3 (N...W Compatto)"}). Defaults to \code{"NOT MATCH"}.}
#'   \item{\code{specie_estratta}:}{ Character. The 4-letter species string extracted from the text (\code{NA} if no match could be identified).}
#'   \item{\code{coord_X_raw}:}{ Character. The raw, unparsed substring isolated for the X coordinate component.}
#'   \item{\code{coord_Y_raw}:}{ Character. The raw, unparsed substring isolated for the Y coordinate component.}
#'   \item{\code{coord_X}:}{ Numeric. The final calculated and scaled numeric X value (\code{NA_real_} if missing or parsing failed).}
#'   \item{\code{coord_Y}:}{ Numeric. The final calculated and scaled numeric Y value (\code{NA_real_} if missing or parsing failed).}
#'   \item{\code{na_coercion}:}{ Logical. \code{TRUE} if a numeric conversion error occurred during string processing, forcing a fallback to \code{NA_real_}; \code{FALSE} otherwise.}
#'   \item{\code{note}:}{ Character. A text field appending error messages or detailed logs regarding NA conversions (e.g., \code{"[NA coercion su: ...]"} if applicable).}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' valid_species <- c("PCAB", "FASY", "LADE")
#' header_data <- ecfh(
#'   file_path = "path/to/sample.rwl",
#'   lista_specie = valid_species
#' )
#' if (!is.null(header_data$Specie)) {
#'   message("Found species: ", header_data$Specie)
#' }
#' }

ecfh <- function(file_path, lista_specie) {
  report_line <- data.frame(
    header_text = NA_character_,
    pattern_matched = "NOT MATCH",
    specie_estratta = NA_character_,
    coord_X_raw = NA_character_,
    coord_Y_raw = NA_character_,
    coord_X = NA_real_,
    coord_Y = NA_real_,
    na_coercion = FALSE,
    note = "",
    stringsAsFactors = FALSE)
  result <- list(header = NULL, Specie = NULL, coord_X = NULL, coord_Y = NULL, report = NULL)
  header_lines <- tryCatch({
    readLines(file_path, n = 3, warn = FALSE)
  }, error = function(e) {
    warning("Cannot read file: ", file_path)
    return(NULL)})

  if (is.null(header_lines)) {
    report_line$note <- "Error read file"
    result$report <- report_line
    return(result)}
  header_text <- stringr::str_squish(paste(header_lines, collapse = " "))
  result$header <- header_text
  report_line$header_text <- header_text

  # Estrazione Specie
  tutte_quattro <- unlist(stringr::str_extract_all(header_text, "\\b[A-Z]{4}\\b"))
  if (length(tutte_quattro) == 1) result$Specie <- tutte_quattro
  match_lista <- intersect(tutte_quattro, lista_specie)
  if (length(match_lista) == 1) result$Specie <- match_lista
  report_line$specie_estratta <- ifelse(is.null(result$Specie), NA, result$Specie)
  parse_num <- function(val_str, scale = 100) {
    if (is.na(val_str) || val_str == "") return(NA_real_)
    cleaned_str <- gsub("[^0-9-]", "", val_str)
    val_num <- suppressWarnings(as.numeric(cleaned_str))
    if (is.na(val_num)) {
      report_line$na_coercion <<- TRUE
      report_line$note <<- paste(report_line$note, "[NA coercion su:", val_str, "]")
      return(NA_real_)}
    return(val_num / scale)}

  # --- NUOVO PATTERN: M COMPATTO (es: "M 522701326 __") ---
  mM <- stringr::str_match(header_text, "M\\s*([[:digit:]]{8,9})\\s+_{1,}")
  if (!is.na(mM[1, 1])) {
    report_line$pattern_matched <- "PATTERN M-COMPATTO"
    raw_str <- mM[1, 2]
    len <- nchar(raw_str)
    raw_lon <- substr(raw_str, 1, 4)
    raw_lat <- substr(raw_str, 5, len)
    report_line$coord_X_raw <- raw_lon
    report_line$coord_Y_raw <- raw_lat
    result$coord_Y <- parse_num(raw_lat, scale = 100)
    result$coord_X <- parse_num(raw_lon, scale = 100)
    report_line$coord_X <- result$coord_X
    report_line$coord_Y <- result$coord_Y
    result$report <- report_line
    return(result)}

  # --- PATTERN 3 CORRETTO (es: "5256N00355W") ---
  m3 <- stringr::str_match(header_text, "([[:digit:]]{4})\\s*N\\s*([[:digit:]]{4,5})\\s*W")
  if (!is.na(m3[1, 1])) {
    report_line$pattern_matched <- "PATTERN 3 (N...W Compatto)"
    report_line$coord_X_raw <- m3[1, 2]
    report_line$coord_Y_raw <- m3[1, 3]
    result$coord_Y <- parse_num(m3[1, 3], scale = 100)
    lon_val <- parse_num(m3[1, 2], scale = 100)
    result$coord_X <- if(!is.na(lon_val)) -abs(lon_val) else NA_real_
    report_line$coord_X <- result$coord_X
    report_line$coord_Y <- result$coord_Y
    result$report <- report_line
    return(result)}

  # --- PATTERN 1 (M N E standard con spazi) ---
  m1 <- stringr::str_match(header_text, "(-?[[:digit:]]+)\\s*M\\s*(-?[[:digit:]]+)\\s*N?(-?[[:digit:]]+)\\s*E?")
  if (!is.na(m1[1, 1])) {
    report_line$pattern_matched <- "PATTERN 1 (M N E)"
    report_line$coord_X_raw <- m1[1, 3]
    report_line$coord_Y_raw <- m1[1, 4]
    result$coord_X <- parse_num(m1[1, 3], scale = 100)
    result$coord_Y <- parse_num(m1[1, 4], scale = 100)
    report_line$coord_X <- result$coord_X
    report_line$coord_Y <- result$coord_Y
    result$report <- report_line
    return(result)}

  # --- PATTERN 4 (signed DDMM-DDDMM) ---
  m4 <- stringr::str_match(header_text, "([+-][[:digit:]]{3,5})([+-][[:digit:]]{3,6})")
  if (!is.na(m4[1, 1])) {
    report_line$pattern_matched <- "PATTERN 4 (signed DDMM)"
    report_line$coord_X_raw <- m4[1, 2]
    report_line$coord_Y_raw <- m4[1, 3]
    result$coord_X <- parse_num(m4[1, 2], scale = 100)
    result$coord_Y <- parse_num(m4[1, 3], scale = 100)
    report_line$coord_X <- result$coord_X
    report_line$coord_Y <- result$coord_Y
    result$report <- report_line
    return(result)}

  # --- PATTERN 2 (decimali seguiti da _) ---
  m2 <- stringr::str_match(header_text, "(-?[[:digit:]]+\\.?[[:digit:]]*)\\s+(-?[[:digit:]]+\\.?[[:digit:]]*)\\s+_{1,}")
  if (!is.na(m2[1, 1])) {
    report_line$pattern_matched <- "PATTERN 2 (decimals with _)"
    report_line$coord_X_raw <- m2[1, 2]
    report_line$coord_Y_raw <- m2[1, 3]
    val2 <- parse_num(m2[1, 2], scale = 1)
    val3 <- parse_num(m2[1, 3], scale = 1)
    if (!is.na(val2) && val2 < 1000) {
      result$coord_X <- val3
      result$coord_Y <- val2
    } else {
      result$coord_X <- val2 / 100
      result$coord_Y <- val3 / 100}
    report_line$coord_X <- result$coord_X
    report_line$coord_Y <- result$coord_Y
    result$report <- report_line
    return(result)}

  # --- PATTERN 10 CORRETTO (Rimozione trattini, inversione X/Y e segno meno) ---
  m10 <- stringr::str_match(header_text, "(.{10})\\s+_")
  if (!is.na(m10[1, 1])) {
    report_line$pattern_matched <- "PATTERN 10 (fixed width)"
    raw_val <- m10[1, 2]
    raw_x <- substr(raw_val, 1, 5)
    raw_y <- substr(raw_val, 6, 10)
    report_line$coord_X_raw <- raw_x
    report_line$coord_Y_raw <- raw_y
    ha_trattino <- grepl("-", raw_x) || grepl("-", raw_y)
    clean_x <- gsub("-", "", raw_x)
    clean_y <- gsub("-", "", raw_y)
    val_x <- parse_num(clean_x, scale = 100)
    val_y <- parse_num(clean_y, scale = 100)
    if (ha_trattino && !is.na(val_y)) {
      val_y <- -abs(val_y)}
    result$coord_X <- val_x
    result$coord_Y <- val_y
    report_line$coord_X <- result$coord_X
    report_line$coord_Y <- result$coord_Y
    result$report <- report_line
    return(result)}

  result$report <- report_line
  return(result)}

#' Extract Metadata and Spatialize Master Tree-Ring Chronologies
#'
#' @description
#' Extracts site metadata and coordinates from Tucson header files and/or a metadata text file,
#' resolves coordinate choices based on user preferences, corrects scaling/formatting issues,
#' and returns both a combined data frame and a spatial vector object (`terra::SpatVector`).
#'
#' @param stand_validi A character vector containing the names/IDs of the valid stand files to process.
#' @param path_metadata A character string specifying the file path to the metadata text file (`.txt`).
#' @param dir_file A character string specifying the directory path where the stand RWL files are stored.
#' @param header Logical. If \code{TRUE}, primary priority for latitude and longitude is given to
#'   the Tucson header values (\code{coord_X} and \code{coord_Y}). If \code{FALSE} (default), priority
#'   is given to the text metadata file (\code{meta_lat} and \code{meta_lon}). In both cases, if the
#'   primary source contains \code{NA}, the secondary source is used as a fallback.
#'
#' @return A named \code{list} with two elements:
#' \describe{
#'   \item{\code{df_coor}}{A \code{data.frame} containing all extracted and processed metadata,
#'     including final coordinates (\code{lat_definitiva}, \code{lon_definitiva}) and species information.}
#'   \item{\code{punti_vect}}{A \code{terra::SpatVector} object containing the spatial points with
#'     WGS84 projection (+proj=longlat +datum=WGS84) for valid geometries, or \code{NULL} if no valid
#'     coordinates are found.}
#' }
#'
#' @details
#' The function compares coordinates between headers and metadata files using an internal
#' helper function \code{ec()}. It handles common digitization anomalies such as un-scaled coordinates
#' (values > 100) and sign/hemisphere inversions.
#'
#' @importFrom terra vect
#' @importFrom dplyr bind_rows
#'
#' @export

estrai_metadati_e_master <- function(stand_validi, path_metadata, dir_file, header = FALSE) {
  metadata_txt <- read.delim(path_metadata, stringsAsFactors = FALSE)
  lista_report <- list()

  for (i in seq_along(stand_validi)) {
    stand <- stand_validi[i]
    file_path <- file.path(dir_file, stand)

    ec_res <- tryCatch({
      ec(file_path, path_metadata)
    }, error = function(e) NULL)

    if (is.null(ec_res) || is.null(ec_res$comparison_report)) next

    rep <- ec_res$comparison_report
    st_id_clean <- rep$stand_id[1]

    # Recupero Altitudine dai metadati
    idx_meta <- which(toupper(metadata_txt$ITRDB_Code) == toupper(st_id_clean))
    rep$altitude <- if (length(idx_meta) > 0) metadata_txt$Elevation[idx_meta[1]] else NA_real_

    # --- ASSEGNAZIONE PULITA E DIRETTA DI LATITUDINE E LONGITUDINE ---
    if (header) {
      # Caso HEADER = TRUE: priorità a Tucson (coord_X = Lat, coord_Y = Lon)
      lat_base <- ifelse(!is.na(rep$coord_X), rep$coord_X, rep$meta_lat)
      lon_base <- ifelse(!is.na(rep$coord_Y), rep$coord_Y, rep$meta_lon)
    } else {
      # Caso HEADER = FALSE: priorità ai Metadati .txt
      lat_base <- ifelse(!is.na(rep$meta_lat), rep$meta_lat, rep$coord_X)
      lon_base <- ifelse(!is.na(rep$meta_lon), rep$meta_lon, rep$coord_Y)
    }

    # --- Controllo numeri esageratamente grandi (es. gradi senza punto decimale) ---
    if (!is.na(lon_base) && !is.na(lat_base) && lon_base > 100 && lat_base > 100) {
      lon_base <- lon_base / 100
      lat_base <- lat_base / 100
    }

    # --- VERIFICA CORRISPONDENZA SEGNI ---
    ha_X <- !is.na(rep$coord_X) && !is.na(rep$meta_lat) && !is.na(rep$diff_lat)
    ha_Y <- !is.na(rep$coord_Y) && !is.na(rep$meta_lon) && !is.na(rep$diff_lon)

    if (ha_X && rep$diff_lat > 2) {
      lat_prov <- abs((-rep$coord_X) - rep$meta_lat)
      if (lat_prov <= 2) {
        rep$diff_lat <- lat_prov
      }
    }
    if (ha_Y && rep$diff_lon > 2) {
      lon_prov <- abs((-rep$coord_Y) - rep$meta_lon)
      if (lon_prov <= 2) {
        rep$diff_lon <- lon_prov
      }
    }
    if (ha_X && ha_Y && rep$diff_lat <= 2 && rep$diff_lon <= 2) {
      rep$correspond_coord <- TRUE
    }

    # Assegnazione definitiva delle coordinate selezionate
    rep$lat_definitiva <- lat_base
    rep$lon_definitiva <- lon_base

    # Specie definitiva (Priorità metadati)
    rep$specie_definitiva <- ifelse(!is.na(rep$meta_species) & rep$meta_species != "",
                                    rep$meta_species, rep$specie_estratta)
    lista_report[[stand]] <- rep
    if (i %% 100 == 0) message("Extracted Metadata : ", i, " / ", length(stand_validi), "\n")
  }

  message("End Extraction. Total Metadata : ", length(lista_report), " / ", length(stand_validi), "\n")

  if (length(lista_report) == 0) return(list(df_coor = NULL, punti_vect = NULL))

  # Unione dei report
  df_master_totale <- if (requireNamespace("dplyr", quietly = TRUE)) {
    dplyr::bind_rows(lista_report)
  } else {
    do.call(rbind, lista_report)
  }

  # Filtraggio righe con coordinate valide per la spazializzazione
  df_valid_geo <- df_master_totale[!is.na(df_master_totale$lon_definitiva) &
                                     !is.na(df_master_totale$lat_definitiva), ]

  punti_vect <- NULL
  if (nrow(df_valid_geo) > 0) {
    punti_vect <- terra::vect(
      df_valid_geo,
      geom = c("lon_definitiva", "lat_definitiva"),
      crs = "+proj=longlat +datum=WGS84"
    )
  }

  return(list(df_coor = df_master_totale, punti_vect = punti_vect))
}
