#' Read multiple pdb files
#'
#' This function reads multiple pdb files containing multiple structures.
#'
#' @param paths_to_structures directory were the structures are stored.
#'
#' @return a list of the loaded structures.
#' @export
#'
#' @examples
read_multi_pdb <- function(paths_to_structures){
  structures <- lapply(paths_to_structures,read.pdb)
  names(structures) <- gsub(".pdb","",gsub(".*/","",paths_to_structures))
  return(structures)
}
