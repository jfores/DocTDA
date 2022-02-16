#' get_level_one
#'
#' Carries out PAM clustering of  a matrix of topological-based distances.
#'
#' @param dist_mat Matrix of topological distances among biological assemblies.
#' @param n_clust  Number of clusters to retrieve.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' get_level_one(dist_mat,n_clust)
#' }
get_level_one <- function(dist_mat,n_clust){
  pam_clust <- cluster::pam(dist_mat,diss = TRUE,k = n_clust)
  medoids_data <- pam_clust$medoids
  df_data <- data.frame(pam_clust$clustering)
  df_data$Str <- rownames(df_data)
  splitted_clust <- split(df_data,as.factor(df_data[,1]))
  names(splitted_clust) <- medoids_data
  return(splitted_clust)
}

#' get_level_two
#'
#' Carries out the second level of PAM clustering.
#'
#' @param df_data Data.frame output from get_level_one
#' @param dist_mat Topological distance matrix
#' @param n_clust Number of clusters to produce in level two clustering.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' get_level_two(dist_mat,n_clust)
#' }
get_level_two <- function(df_data,dist_mat,n_clust){
  dist_mat_filt <- dist_mat[df_data[,2],df_data[,2]]
  pam_clust_filt <- cluster::pam(dist_mat_filt,diss = TRUE,k = n_clust)
  medoids_data_filt <- pam_clust_filt$medoids
  df_data_filt <- data.frame(pam_clust_filt$clustering)
  df_data_filt$Str <- rownames(df_data_filt)
  splitted_clust_filt <- split(df_data_filt,as.factor(df_data_filt[,1]))
  names(splitted_clust_filt) <- medoids_data_filt
  return(splitted_clust_filt)
}
