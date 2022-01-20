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
  read_temp <- readLines(log_file_path)
  max_bin_en <- as.numeric(strsplit(read_temp[which(grepl("kcal/mol",read_temp)) + 2]," ")[[1]][13])
  return(max_bin_en)
}
