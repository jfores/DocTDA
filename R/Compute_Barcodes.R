
#' Computes max distance in a point cloud.
#'
#' @param alpha_coord Matrix including the alpha carbon coordinates of a particular structure.
#' @param dimension_betti Maximum Betti dimension to compute.
#'
#' @return returns the diameter or maximum distance of the point cloud.
#' @export
#'
#' @examples
#' \dontrun{
#' compute_max_d(alpha_coord)
#' }
compute_max_d <- function(alpha_coord,dimension_betti = 2){
  diameter_to_hom <- base::max(as.matrix(stats::dist(as.matrix(alpha_coord))))
  return(diameter_to_hom)
}


#' Compute maximum distance multi
#'
#' Compute the maximum distance for a list of alpha carbon coordinate matrices.
#'
#' @param alpha_multi A list including a single or multiple matrices of alpha carbon coordinates generated using the extract_multiple_alpha function.
#'
#' @return Returns a list of max distances.
#' @export
#'
#' @examples
#' \dontrun{
#' compute_max_d_multi(alpha_multi)
#' }
compute_max_d_multi <- function(alpha_multi){
  max_d_multi <- lapply(alpha_multi,compute_max_d)
  names(max_d_multi) <- names(alpha_multi)
  return(max_d_multi)
}



#' Compute horology for one point cloud.
#'
#' @param alpha_coord a matrix containing the alpha carbon coordinates.
#' @param dimension_betti the maximum bettin dimension to compute.
#'
#' @return returns the barcode of the input point cloud data.
#' @export
#'
#' @examples
#' \dontrun{
#' compute_homology(alpha_coord,2)
#' }
compute_homology <- function(alpha_coord,dimension_betti = 2){
  diameter_to_hom <- base::max(base::as.matrix(stats::dist(base::as.matrix(alpha_coord))))
  barcode_out <- TDAstats::calculate_homology(alpha_coord,dim = dimension_betti,threshold = diameter_to_hom,standardize = FALSE)
  return(barcode_out)
}



#' Computes homology for multiple alfa carbon coordinate matrices stored into a list.
#'
#' @param alpha_multi A list containing the matrices storing the alpha carbon possitions for different structures.
#' @param dimension_betti Max Betti dimension to compute.
#'
#' @return returns a lists with the barcodes of the multiple point cloud data
#' @export
#'
#' @examples
#' \dontrun{
#' compute_homology_multi(alpha_multi,dimension_betti)
#' }
compute_homology_multi <- function(alpha_multi,dimension_betti = 2){
  barcodes_multi <- lapply(alpha_multi,compute_homology,dimension_betti)
  names(barcodes_multi) <- names(alpha_multi)
  return(barcodes_multi)
}
