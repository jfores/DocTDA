#' parse_max_energy_vina
#'
#' Retrieves max binding energy from a vina log output file.
#'
#' @param log_file_path
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' parse_max_energy_vina(log_file_path)
#' }
parse_max_energy_vina <- function(log_file_path){
  con <- file(log_file_path)
  read_temp <- readLines(con)
  max_bin_en <- as.numeric(strsplit(read_temp[which(grepl("kcal/mol",read_temp)) + 2]," ")[[1]][13])
  close(con)
  return(max_bin_en)
}
