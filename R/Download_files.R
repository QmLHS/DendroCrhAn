#' Download .rwl Files from a Web Directory
#'
#' This function Connects to the NOAA ITRDB server for European tree-ring measurements, parses the HTML directory listing, filters for valid '.rwl' dataset files, and optionally downloads them to a specified local directory.
#'
#' @param url A character string specifying the base URL of the web directory
#'   containing the `.rwl` files. Must include a trailing slash (e.g., `"https://example.com/data/"`).
#' @param name_folder A character string specifying the path to the local directory
#'   where the downloaded files should be saved.
#' @param download Logical. If `TRUE` (default), downloads the matching files
#'   to `name_folder`. If `FALSE`, only parses and returns the file list without downloading.
#'
#' @return A character vector containing the names of the identified `.rwl` files.
#'
#' @importFrom httr GET status_code content
#' @importFrom rvest read_html html_nodes html_attr
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#'
#' @export
download_files <- function(url,name_folder,download=TRUE) {
  # Fetch the directory listing HTML
  response = GET(url)
  # Check if the request was successful
  if (status_code(response) != 200) {
    stop("Failed to retrieve the directory listing.")}
  # Parse the HTML content
  page = content(response, "text")
  parsed_html = read_html(page)
  # Extract all anchor (<a>) elements, this includes the rwl and the txt
  links = parsed_html %>% html_nodes("a") %>% html_attr("href")
  # Filter out parent directory links and non-rwl file links
  file_links = links[grepl("\\.rwl", links)]
  # Delete the element ".rwl" from the string
  file_links_spill = str_split(file_links, ".rwl") #Split the string
  int=c()
  for (i in 1:length(file_links_spill)){
    # delete the tail and pick up only the name
    int=c(int,file_links_spill[[i]][1])}
  final.list=c()
  for (i in int){
    # If the string end with number, have to record it
    if (grepl("\\s*[[:digit:]]$",i)==TRUE){
      final.list=c(final.list,paste0(i, ".rwl"))}}
  int=final.list
  # Download the file list
  if (download){
    if (!dir.exists(name_folder)) {
      dir.create(name_folder, recursive = TRUE)}
    for (file_link in int) {
      # Construct the full URL for each file
      file_url = paste0(url, file_link)
      # Define the local destination file path
      destfile = file.path(name_folder, basename(file_link))
      # Download the file
      download.file(file_url, destfile, mode = "wb")
      # Print status
      cat("Downloaded:", file_link, "\n")}}
  return(final.list)}
