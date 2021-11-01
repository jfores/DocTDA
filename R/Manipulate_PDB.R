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
#' \dontrun{
#' read_multi_pdb(paths_to_structures)
#' }
read_multi_pdb <- function(paths_to_structures){
  structures <- lapply(paths_to_structures,bio3d::read.pdb)
  names(structures) <- gsub(".pdb","",gsub(".*/","",paths_to_structures))
  return(structures)
}


#' aux_for_extract_multiple_alpha
#'
#' Auxiliary function for the extract_multiple_alpha function. Given a pdb structure it extracts the coordinats from its alpha carbons.
#'
#' @param struc A particular pdb structure loaded using the read.pdb function from the bio3d package.
#'
#' @return Returns a matrix containing the coordinates of the alpha carbons.
#' @export
#'
#' @examples
#' \dontrun{
#' aux_for_extract_multiple_alpha(struc)
#' }
aux_for_extract_multiple_alpha <- function(struc){
  ca.inds <- bio3d::atom.select(struc,string = "calpha")
  struc_filt <- bio3d::trim.pdb(struc,ca.inds)$atom[,c("x","y","z")]
  return(struc_filt)
}


#' extract_multiple_alpha
#'
#' This function extracts the alpha carbon coordintes from a list containing multiple structures.
#'
#' @param structures A list of pdb structures loaded employing the bio3d read.pdb function.
#'
#' @return Returns a list with multiple matrices of alfa carbon coordinates. One for each structure in the input lits.
#' @export
#'
#' @examples
#' \dontrun{
#' aux_for_extract_multiple_alpha(struc)
#' }
extract_multiple_alpha <- function(structures){
  alpha_multi <- lapply(structures,aux_for_extract_multiple_alpha)
  names(alpha_multi) <- names(structures)
  return(alpha_multi)
}



