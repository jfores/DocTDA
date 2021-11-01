
compute_max_d <- function(alpha_coord,dimension_betti = 2){
  diameter_to_hom <- max(as.matrix(dist(as.matrix(alpha_coord))))
  return(diameter_to_hom)
}
