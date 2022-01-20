

#' compute_sim_barcode_pair_aux
#'
#' @param pbf_a_one_dim persistent betty function for a particular dimenson derived from barcode A
#' @param pbf_b_one_dim persistent betty function for a particular dimenson derived from barcode B
#' @param range_x_axis Comparison range.
#'
#' @return returns the persistent betty similarity for the dimensions of
#' @export
#'
#' @examples
#' \dontrun{
#' compute_sim_barcode_pair_aux(pbf_a_one_dim,pbf_b_one_dim,range_x_axis)
#' }
compute_sim_barcode_pair_aux <- function(pbf_a_one_dim,pbf_b_one_dim,range_x_axis){
  f_1 <- approxfun(range_x_axis, pmin(pbf_a_one_dim,pbf_b_one_dim))
  f_1_int <- integrate(f_1, min(range_x_axis), max(range_x_axis),subdivisions = 1000L)
  f_2 <- approxfun(range_x_axis, pmax(pbf_a_one_dim,pbf_b_one_dim))
  f_2_int <- integrate(f_2, min(range_x_axis), max(range_x_axis),subdivisions = 1000L)
  sim_dim <- f_1_int$value/f_2_int$value
  return(sim_dim)
}



#' compute_sim_barcode_pair
#'
#' Compute the barcode betti similarity function for all the studied dimension in two given structures.
#'
#' @param barcode_a barcode structure a
#' @param barcode_b barcode structure b
#' @param resolution resolution parameter that indicates de width of the intervals used to approximate the function.
#'
#' @return returns a list with the similarities between the barcodes of
#' @export
#'
#' @examples
compute_sim_barcode_pair <- function(barcode_a,barcode_b,resolution = 0.01){
  max_a <- max(barcode_a[,2:3])
  max_b <- max(barcode_b[,2:3])
  range_x_axis <- seq(from = 0, to = max(max_a,max_b)+5, by = resolution )
  pbf_a <- compute_pbf(range_x_axis,barcode_a,k=2)
  pbf_b <- compute_pbf(range_x_axis,barcode_b,k=2)
  return(mapply(compute_sim_barcode_pair_aux, pbf_a,pbf_b,MoreArgs = list(range_x_axis)))
}
