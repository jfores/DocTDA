

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
  max_a <- max(barcode_a[,2:3])
  max_b <- max(barcode_b[,2:3])
  range_x_axis <- seq(from = 0, to = max(max_a,max_b)+5, by = resolution )
  pbf_a <- compute_pbf(range_x_axis,barcode_a,k=2)
  pbf_b <- compute_pbf(range_x_axis,barcode_b,k=2)
  return(mapply(compute_sim_barcode_pair_aux, pbf_a,pbf_b,MoreArgs = list(range_x_axis)))
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
    list_out[[i]] <- foreach(j = 1:length(dir_b)) %dopar%
      {
        temp_2 <- get(load(file = dir_b[j]))
        temp_3 <- compute_sim_barcode_pair(temp_1,temp_2)
      }
    save(file = paste(dir_out,"/list_out.Rda",sep=""),list_out)
    save(file = paste(dir_out,"/start_iter.Rda",sep=""),i)
  }
  doParallel::registerDoSEQ()
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
