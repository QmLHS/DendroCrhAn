#' Estrazione rapida e unione di serie temporali climatiche da file NetCDF
#'
#' @description
#' Legge in modo ottimizzato (in RAM) serie temporali da uno o piu file NetCDF per
#' coordinate geografiche specificate. La funzione filtra automaticamente i pixel unici,
#' calcola gli indici spaziali dedicati per ciascuna griglia e unisce le variabili
#' climatiche in una lista di \code{data.frame} indicizzati per \code{pixel_id}.
#'
#' @param df_unici A \code{data.frame} contenente i punti geografici di interesse.
#'   Deve obbligatoriamente includere le colonne:
#'   \itemize{
#'     \item \code{pixel_id}: Identificatore unico del pixel (es. "399_265").
#'     \item \code{longitude}: Longitudine in gradi decimali.
#'     \item \code{latitude}: Latitudine in gradi decimali.
#'   }
#' @param v_file_clima Vettore di caratteri (\code{character}) con i nomi dei file NetCDF
#'   da cui estrarre i dati (es. \code{c("spei01.nc", "CRU_TS_pre.nc")}).
#' @param v_nomi Vettore facoltativo di caratteri (\code{character}) con i nomi da assegnare
#'   alle variabili nel dataset finale (es. \code{c("spei", "pre")}). Se \code{NULL} (default),
#'   i nomi vengono ricavati automaticamente rimuovendo l'estensione dai nomi dei file.
#' @param cartella_dati Stringa (\code{character}) indicante il percorso della cartella
#'   in cui sono salvati i file NetCDF. Il valore predefinito e \code{"dati"}.
#'
#' @return Una \code{list} con un elemento per ciascun \code{pixel_id} unico.
#'   Ogni elemento e un \code{data.frame} contenente:
#'   \itemize{
#'     \item \code{data}: Vettore temporale (estratto dalle unità di tempo del NetCDF).
#'     \item Una o piu colonne corrispondenti ai nomi definiti in \code{v_nomi}
#'           (es. \code{spei}, \code{pre}) contenenti le serie temporali estratte.
#'   }
#'
#' @import ncdf4
#' @export
#'
#' @examples
#' \dontrun{
#' clima <- estrazione_clima(
#'   df_unici = df_pixel,
#'   v_file_clima = c("spei01.nc", "CRU_TS_pre.nc"),
#'   v_nomi = c("spei", "pre"),
#'   cartella_dati = "dati"
#' )
#' }
estrazione_clima <- function(df_unici, v_file_clima, v_nomi = NULL, cartella_dati = "dati") {

  if (is.null(v_nomi)) {
    v_nomi <- tools::file_path_sans_ext(basename(v_file_clima))}

  # 1. Filtriamo per avere solo 1 riga per ciascun pixel_id unico (es. i 467 pixel)
  df_pixel_unici <- df_unici[!duplicated(df_unici$pixel_id), ]
  dati_file <- list()

  # 2. Caricamento in RAM veloce file per file usando le coordinate
  for (j in seq_along(v_file_clima)) {
    nome_var <- v_nomi[j]
    percorso <- file.path(cartella_dati, v_file_clima[j])
    nc <- nc_open(percorso)
    map_lon <- ncvar_get(nc, "lon")
    map_lat <- ncvar_get(nc, "lat")
    # Trova gli indici lon/lat specifici per QUESTO file ma per i pixel UNICI
    find_idx <- function(lon, lat) {
      c(idx_lon = which.min(abs(map_lon - lon)),
        idx_lat = which.min(abs(map_lat - lat)))}
    indices <- t(mapply(find_idx, df_pixel_unici$longitude, df_pixel_unici$latitude))
    # Estrazione date
    time_raw <- ncvar_get(nc, "time")
    time_units <- ncatt_get(nc, "time", "units")$value
    if (!is.null(time_units) && grepl("since", time_units)) {
      origin_str <- sub(".*since ", "", time_units)
      if (grepl("days", time_units, ignore.case = TRUE)) {
        vett_date <- as.Date(time_raw, origin = origin_str)
      } else if (grepl("hours", time_units, ignore.case = TRUE)) {
        vett_date <- as.POSIXct(time_raw * 3600, origin = origin_str, tz = "UTC")
      } else {
        vett_date <- time_raw}
    } else {
      vett_date <- time_raw}
    var_nc_name <- names(nc$var)[names(nc$var) %in% c(nome_var, "spei", "pre", "tmp", "pet")][1]
    if (is.na(var_nc_name)) var_nc_name <- names(nc$var)[1]
    array_3d <- ncvar_get(nc, var_nc_name)
    nc_close(nc)
    dati_file[[nome_var]] <- list(
      indices = indices,
      date = vett_date,
      array = array_3d)}
  # 3. Assemblaggio della lista denominata esattamente per pixel_id (es. "399_265")
  lista_meteo_raw <- list()
  for (i in 1:nrow(df_pixel_unici)) {
    p_id <- as.character(df_pixel_unici$pixel_id[i]) # "399_265"
    df_pixel_accumulato <- NULL
    for (j in seq_along(v_file_clima)) {
      nome_var <- v_nomi[j]
      idx_lon <- dati_file[[nome_var]]$indices[i, "idx_lon"]
      idx_lat <- dati_file[[nome_var]]$indices[i, "idx_lat"]
      vett_date <- dati_file[[nome_var]]$date
      valori <- dati_file[[nome_var]]$array[idx_lon, idx_lat, ]
      df_temp <- data.frame(
        data_raw = as.Date(vett_date),
        valore = valori)
      colnames(df_temp)[2] <- nome_var
      if (is.null(df_pixel_accumulato)) {
        df_pixel_accumulato <- df_temp
      } else {
        df_pixel_accumulato <- merge(df_pixel_accumulato, df_temp, by = "data_raw", all = TRUE)}}
    df_pixel_accumulato <- df_pixel_accumulato[order(df_pixel_accumulato$data_raw), ]
    df_pixel_accumulato$data <- tolower(format(df_pixel_accumulato$data_raw, "%b-%Y"))
    df_pixel_accumulato$data_raw <- NULL
    df_pixel_accumulato <- df_pixel_accumulato[, c("data", v_nomi)]
    lista_meteo_raw[[p_id]] <- df_pixel_accumulato}
  return(lista_meteo_raw)}

#' Reassign Ocean Climate Class to Nearest Land Cell
#'
#' @description
#' Identifies spatial points within a metadata data frame that fall into the "Ocean"
#' climate zone (based on a Köppen-Geiger climate raster map) and reassigns them to
#' the climate zone of the nearest terrestrial (non-Ocean) cell using a local spatial crop.
#'
#' @param df_meta A data frame containing spatial points with metadata. Must include
#'   columns \code{longitude}, \code{latitude}, and \code{climate_zona}.
#' @param kg_map A \code{SpatRaster} object representing the Köppen-Geiger climate classification map.
#' @param rat A data frame acting as the Raster Attribute Table (RAT), containing at least an
#'   \code{ID} column and a climate name column (\code{climate} or \code{climate_zona}).
#' @param buffer_deg Numeric. The geographic bounding box extension (in degrees) used
#'   to crop the raster locally around each Ocean point to search for land cells. Default is \code{0.5}.
#'
#' @return A modified version of \code{df_meta} where entries with \code{climate_zona == "Ocean"}
#'   have been updated with the nearest land-based climate zone classification.
#'
#' @details
#' To optimize RAM usage and execution speed, the function avoids converting the entire global
#' raster into vector points. Instead, it creates a local bounding box around each Ocean point
#' to extract adjacent terrestrial cells and determines spatial proximity via \code{terra::nearest()}.
#'
#' @importFrom terra vect crs ext crop hasValues values as.points nearest
#' @importFrom dplyr left_join
#'
#' @export
Ocean_mask <- function(df_meta, kg_map, rat, buffer_deg = 0.5) {

  # Check required columns in df_meta
  req_cols <- c("longitude", "latitude", "climate_zona")
  if (!all(req_cols %in% colnames(df_meta))) {
    stop("df_meta must contain columns: 'longitude', 'latitude', and 'climate_zona'")
  }

  # Identify mask for points falling in 'Ocean'
  ocean_mask <- !is.na(df_meta$climate_zona) & df_meta$climate_zona == "Ocean"

  if (any(ocean_mask)) {

    # Identify the column name used for climate label in RAT
    col_rat_name <- if ("climate" %in% colnames(rat)) "climate" else "climate_zona"
    ocean_id <- rat$ID[rat[[col_rat_name]] == "Ocean"]

    # Mask out Ocean cells from the raster
    kg_land <- kg_map
    kg_land[kg_land == ocean_id] <- NA

    # Convert Ocean points to SpatVector
    punti_ocean <- terra::vect(
      df_meta[ocean_mask, c("longitude", "latitude")],
      geom = c("longitude", "latitude"),
      crs = terra::crs(kg_map)
    )

    nuovi_id <- numeric(nrow(punti_ocean))

    # Process each Ocean point with a local bounding box crop
    for (i in seq_len(nrow(punti_ocean))) {
      pt <- punti_ocean[i]

      e <- terra::ext(pt) + buffer_deg
      sub_raster <- terra::crop(kg_land, e)

      if (terra::hasValues(sub_raster) && any(!is.na(terra::values(sub_raster)))) {
        pts_sub <- terra::as.points(sub_raster)
        nearest_idx <- terra::nearest(pt, pts_sub)

        nuovi_id[i] <- pts_sub[[1]][nearest_idx$to_id, 1]
      } else {
        nuovi_id[i] <- NA_real_
      }
    }

    # Map new IDs to RAT climate names
    nuovi_climi <- data.frame(ID = nuovi_id) %>%
      dplyr::left_join(rat, by = "ID")

    # Update df_meta
    df_meta$climate_zona[ocean_mask] <- nuovi_climi[[col_rat_name]]
  }

  return(df_meta)
}
