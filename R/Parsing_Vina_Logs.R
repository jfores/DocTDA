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


#' read_doc_res_from_list_dirs
#'
#' Parses a list of docking results from a vector of directories.
#'
#' @param x vector of directories to vina log files.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' read_doc_res_from_list_dirs(x)
#' }
read_doc_res_from_list_dirs <- function(x){
  directories <- x
  x <- gsub(".*/","",x)
  x <- data.frame(do.call("rbind",strsplit(x,"_AND_")))
  x[,2] <- gsub(".log","",x[,2])
  x$Dir_file <- directories
  x$Max_Energy <- rep(NA,nrow(x))
  for(i in 1:nrow(x)){
    print(i)
    try({
      x[i,4] <- parse_max_energy_vina(x[i,3])
    },silent=TRUE)
  }
  return(x)
}




#' aux_same_drug_diff_struc
#'
#' Auxiliary function to compute diferences in docking energies between different structures docked to the same drug.
#'
#' @param x dataframe of different structures docked to the same ligand.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' aux_same_drug_diff_struc(x)
#' }
aux_same_drug_diff_struc <- function(x){
  counter_test <- 1
  list_out <- list()
  for(i in 1:nrow(x)){
    for(j in (i+1):nrow(x)){
      list_out[[counter_test]] <- c(x[i,1],x[j,1],abs(x[i,4]-x[j,4]))
      counter_test <- counter_test + 1
    }
  }
  df_diffs  <- data.frame(do.call("rbind",list_out))
  df_diffs$X3 <- as.numeric(df_diffs$X3)
  colnames(df_diffs) <- c("Str_1","Str_2","Abs_Diff")
  return(df_diffs)
}




#' same_drug_diff_struc
#'
#' @param x Complete dataset of docking results.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' same_drug_diff_struc(x)
#' }
same_drug_diff_struc <- function(x){
  x <- split(x,as.factor(x[,2]))
  x <- lapply(x,aux_same_drug_diff_struc)
  x <- do.call("rbind",x)
  return(x)
}
