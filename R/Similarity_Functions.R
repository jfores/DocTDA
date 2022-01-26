

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
  f_1_int <- integrate(f_1, min(range_x_axis), max(range_x_axis),subdivisions = 1000L,stop.on.error = FALSE)
  f_2 <- approxfun(range_x_axis, pmax(pbf_a_one_dim,pbf_b_one_dim))
  f_2_int <- integrate(f_2, min(range_x_axis), max(range_x_axis),subdivisions = 1000L,stop.on.error = FALSE)
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
#' \dontrun{
#' compute_sim_barcode_pair(barcode_a,barcode_b,resolution)
#' }
compute_sim_barcode_pair <- function(barcode_a,barcode_b,resolution = 0.01){
  try({
  max_a <- max(barcode_a[,2:3])
  max_b <- max(barcode_b[,2:3])
  range_x_axis <- seq(from = 0, to = max(max_a,max_b)+5, by = resolution )
  pbf_a <- compute_pbf(range_x_axis,barcode_a,k=2)
  pbf_b <- compute_pbf(range_x_axis,barcode_b,k=2)
  return(mapply(compute_sim_barcode_pair_aux, pbf_a,pbf_b,MoreArgs = list(range_x_axis)))
  })
  return(c(0,0,0))
}






#' compute_sim_two_dir_bar
#'
#' @param dir_a Paths to the first set of barcodes.
#' @param dir_b Paths to the second set of barcodes.
#' @param dir_out Output directory.
#' @param n_cores Number of processors for parallel computing.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' compute_sim_two_dir_bar(dir_a,dir_b,dir_out,n_cores)
#' }
compute_sim_two_dir_bar <- function(dir_a,dir_b,dir_out,n_cores = 30){
  doParallel::registerDoParallel(n_cores)
  if(file.exists(paste(dir_out,"/list_out.Rda",sep=""))){
    list_out <- get(load(paste(dir_out,"/list_out.Rda",sep="")))
    start_iter <- get(load(paste(dir_out,"/start_iter.Rda",sep="")))
  }else{
    list_out <- list()
    start_iter <- 1
  }
  for(i in start_iter:length(dir_a)){
    print(i)
    temp_1 <- get(load(file = dir_a[i]))
    list_out[[i]] <- foreach::foreach(j = 1:length(dir_b)) %dopar%
      {
        temp_2 <- get(load(file = dir_b[j]))
        temp_3 <- compute_sim_barcode_pair(temp_1,temp_2)
      }
    save(file = paste(dir_out,"/list_out.Rda",sep=""),list_out)
    save(file = paste(dir_out,"/start_iter.Rda",sep=""),i)
  }
  foreach::registerDoSEQ()
  row_names_mat <- gsub("_BC.Rda","",gsub(".*/","",dir_a))
  col_names_mat <- gsub("_BC.Rda","",gsub(".*/","",dir_b))
  b_cero <- do.call("rbind",lapply(list_out,function(x) unlist(lapply(x,function(x)x[1]))))
  rownames(b_cero) <- row_names_mat
  colnames(b_cero) <- col_names_mat
  b_one <- do.call("rbind",lapply(list_out,function(x) unlist(lapply(x,function(x)x[2]))))
  rownames(b_one) <- row_names_mat
  colnames(b_one) <- col_names_mat
  b_two <- do.call("rbind",lapply(list_out,function(x) unlist(lapply(x,function(x)x[3]))))
  rownames(b_two) <- row_names_mat
  colnames(b_two) <- col_names_mat
  b_mean <- (b_cero + b_one + b_two)/3
  rownames(b_mean) <- row_names_mat
  colnames(b_mean) <- col_names_mat
  list_mat <- list(b_cero,b_one,b_two,b_mean)
  return(list_mat)
}



#' compute_sim_one_dir_bar
#'
#' @param path_to_bcs vector including the complete paths to all barcode files.
#' @param dir_out directory where the output data must be written.
#' @param n_cores Number of processores used in for each
#' @param tail_to_remove Tail tagg to remove from barcode name.
#'
#' @importFrom foreach %dopar%
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' compute_sim_one_dir_bar(path_to_bcs,dir_out,n_cores,tail_to_remove)
#' }
compute_sim_one_dir_bar <- function(path_to_bcs,dir_out,n_cores = 30,tail_to_remove = ".Rda"){
  #`%dopar%` <- foreach::`%dopar%`
  if(file.exists(paste(dir_out,"/list_out.Rda",sep=""))){
    list_out <- get(load(paste(dir_out,"/list_out.Rda",sep="")))
    start_iter <- get(load(paste(dir_out,"/start_iter.Rda",sep="")))
    start_iter <- start_iter + 1
  }
  else{
    list_out <- list()
    start_iter <- 1
  }
  for(i in start_iter:length(path_to_bcs)){
    doParallel::registerDoParallel(n_cores)
    print(i)
    temp_1 <- get(load(file = path_to_bcs[i]))
    list_out[[i]] <- foreach::foreach(j = i:length(path_to_bcs)) %dopar%
      {
        temp_2 <- get(load(file = path_to_bcs[j]))
        temp_3 <- compute_sim_barcode_pair(temp_1,temp_2)
      }
    save(file = paste(dir_out,"/list_out.Rda",sep=""),list_out)
    save(file = paste(dir_out,"/start_iter.Rda",sep=""),i)
    foreach::registerDoSEQ()
    gc()
  }
  for(i in 1:length(list_out)){
    #print(i -1)
    temp <- rep(list(c(0,0,0)),i-1)
    list_out[[i]] <- c(temp,list_out[[i]])
  }
  row_col_names <- gsub(tail_to_remove,"",gsub(".*/","",path_to_bcs))
  b_cero <- do.call("rbind",lapply(list_out,function(x) unlist(lapply(x,function(x)x[1]))))
  rownames(b_cero) <- row_col_names
  colnames(b_cero) <- row_col_names
  b_one <- do.call("rbind",lapply(list_out,function(x) unlist(lapply(x,function(x)x[2]))))
  rownames(b_one) <- row_col_names
  colnames(b_one) <- row_col_names
  b_two <- do.call("rbind",lapply(list_out,function(x) unlist(lapply(x,function(x)x[3]))))
  rownames(b_two) <- row_col_names
  colnames(b_two) <- row_col_names
  b_mean <- (b_cero + b_one + b_two)/3
  rownames(b_mean) <- row_col_names
  colnames(b_mean) <- row_col_names
  list_mat <- list(b_cero,b_one,b_two,b_mean)
  return(list_mat)
}



#' mount_matrix
#'
#' Mounts distance matrix from partially computed distances.
#'
#' @param dir_files directory where the output intermediate files from compute_sim_one_dir_bar or compute_sim_two_dir_bar are placed
#' @param structures_ini original list of barcodes to compute used to compute similarities.
#' @param tail_to_remove  character string to remove from the tail of the barcode file paths.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' mount_matrix(dir_files,structures_ini,tail_to_remove)
#' }
mount_matrix <- function(dir_files,structures_ini,tail_to_remove = ".Rda"){
  list_out <- get(load(paste(dir_files,"/list_out.Rda",sep="")))
  for(i in 1:length(list_out)){
    #print(i -1)
    temp <- rep(list(c(0,0,0)),i-1)
    list_out[[i]] <- c(temp,list_out[[i]])
  }
  row_col_names <- gsub(tail_to_remove,"",gsub(".*/","",structures_ini))
  b_cero <- do.call("rbind",lapply(list_out,function(x) unlist(lapply(x,function(x)x[1]))))
  rownames(b_cero) <- row_col_names[1:nrow(b_cero)]
  colnames(b_cero) <- row_col_names
  b_one <- do.call("rbind",lapply(list_out,function(x) unlist(lapply(x,function(x)x[2]))))
  rownames(b_one) <- row_col_names[1:nrow(b_one)]
  colnames(b_one) <- row_col_names
  b_two <- do.call("rbind",lapply(list_out,function(x) unlist(lapply(x,function(x)x[3]))))
  rownames(b_two) <- row_col_names[1:nrow(b_two)]
  colnames(b_two) <- row_col_names
  b_mean <- (b_cero + b_one + b_two)/3
  rownames(b_mean) <- row_col_names[1:nrow(b_mean)]
  colnames(b_mean) <- row_col_names
  list_mat <- list(b_cero,b_one,b_two,b_mean)
  return(list_mat)
}



#' aux_fun_complete_matrices
#'
#' @param x upper diagonal matrix to complete.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' aux_fun_complete_matrices(x)
#' }
aux_fun_complete_matrices <- function(x){
  x_t <- t(x)
  x[lower.tri(x)] <- x_t[lower.tri(x_t)]
  return(x)

}


#' complete_matrices
#'
#' @param sim_mat_list list of similarity matrices derived from compute_sim_one_dir_bar or compute_sim_two_dir_bar
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' complete_matrices(sim_mat_list)
#' }
complete_matrices <- function(sim_mat_list){
  return(lapply(sim_mat_list,aux_fun_complete_matrices))
}


