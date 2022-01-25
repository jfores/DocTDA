#' format_matrices
#'
#'
#'
#' @param dist_mat list of similarity matrices computed by compute_sim_two_dir_bar function.
#' @param db_data  Data frame containing the infromation from the reference database.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' format_matrices(dist_mat,db_data)
#' }
format_matrices <- function(dist_mat,db_data){
  df_temp <- cbind(reshape2::melt(dist_mat[[1]]),reshape2::melt(dist_mat[[2]]),reshape2::melt(dist_mat[[3]]),reshape2::melt(dist_mat[[4]]))
  df_temp <- df_temp[,c(1,2,3,6,9,12)]
  colnames(df_temp) <- c("DB_ID","Q_ID","B_0","B_1","B_2","B_Av")
  df_temp$DB_PDB <- gsub("_.*","",df_temp[,1])
  merged_df <- merge(test,db_data,by.x = "DB_PDB",by.y = "PDB.ID")
  return(merged_df)
}
