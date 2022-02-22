

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
  if(nrow(x) != ncol(x)){
    x <- rbind(x,matrix(0,nrow = ncol(x) - nrow(x),ncol = ncol(x)))
  }
  x_t <- t(x)
  x[lower.tri(x)] <- x_t[lower.tri(x_t)]
  rownames(x) <- colnames(x)
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


#' transform_sim_mat
#'
#' Gets completed similarity matrices and generates a clean data.frame
#'
#'
#' @param sim_matrices_comp list of similarity matrices
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' transform_sim_mat(sim_matrices_comp)
#' }
transform_sim_mat <- function(sim_matrices_comp){

  #Getting similarity matrices removing similarities between the same structures and setting the lower diagonal to cero...

  sim_cero <-  sim_matrices_comp[[1]]
  diag(sim_cero) <- 0
  sim_cero[lower.tri(sim_cero)] <- 0
  sim_cero_melt <- reshape2::melt(sim_cero)

  sim_one <-  sim_matrices_comp[[2]]
  diag(sim_one) <- 0
  sim_one[lower.tri(sim_one)] <- 0
  sim_one_melt <- reshape2::melt(sim_one)

  sim_two <-  sim_matrices_comp[[3]]
  diag(sim_two) <- 0
  sim_two[lower.tri(sim_two)] <- 0
  sim_two_melt <- reshape2::melt(sim_two)

  sim_av <- sim_matrices_comp[[4]]
  diag(sim_av) <- 0
  sim_av[lower.tri(sim_av)] <- 0
  sim_av_melt <- reshape2::melt(sim_av)

  #merging simmilarities for all dimensions...

  all_sim_melt <- cbind(sim_cero_melt,sim_one_melt,sim_two_melt,sim_av_melt)
  all_sim_melt <- all_sim_melt[,-c(4,5,7,8,10,11)]
  colnames(all_sim_melt) <- c("Str_1","Str_2","D_Zero","D_One","D_Two","D_Av")

  #Removing similarities between identical structures...

  all_sim_melt <- all_sim_melt[!all_sim_melt$Str_1 == all_sim_melt$Str_2,]

  #Removing similarities with AV_sim = 0

  all_sim_melt <- all_sim_melt[!all_sim_melt$D_Av == 0,]

  #Ordering matrix based on the average similarities with decreasing values...

  all_sim_melt <- all_sim_melt[order(all_sim_melt$D_Av,decreasing = TRUE),]

  return(all_sim_melt)
}


#' sample_structures_a
#'
#' Add group labels based on av similarity values.
#'
#' @param x output dataframe from transform_sim_mat.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' sample_structures_a(x)
#' }
sample_structures_a <- function(x){
  x$group <- rep(NA, nrow(x))
  x[x[,3] > 0.95,7] <- 1
  x[x[,3] > 0.9 & x[,6] <= 0.95,7] <- 2
  x[x[,3] > 0.85 & x[,6] <= 0.9,7] <- 3
  x[x[,3] > 0.8 & x[,6] <= 0.85,7] <- 4
  x[x[,3] > 0.75 & x[,6] <= 0.8,7] <- 5
  x[x[,3] > 0.70 & x[,6] <= 0.75,7] <- 6
  x[x[,3] > 0.65 & x[,6] <= 0.70,7] <- 7
  x[x[,3] > 0.60 & x[,6] <= 0.65,7] <- 8
  x[x[,3] > 0.55 & x[,6] <= 0.60,7] <- 9
  x[x[,3] > 0.50 & x[,6] <= 0.55,7] <- 10
  x[x[,3] > 0.45 & x[,6] <= 0.50,7] <- 11
  x[x[,3] > 0.40 & x[,6] <= 0.45,7] <- 12
  x[x[,3] > 0.35 & x[,6] <= 0.40,7] <- 13
  x[x[,3] > 0.30 & x[,6] <= 0.35,7] <- 14
  x[x[,3] > 0.25 & x[,6] <= 0.30,7] <- 15
  x[x[,3] > 0.20 & x[,6] <= 0.25,7] <- 16
  x[x[,3] > 0.15 & x[,6] <= 0.20,7] <- 17
  x[x[,3] > 0.10 & x[,6] <= 0.15,7] <- 18
  x[x[,3] > 0.05 & x[,6] <= 0.10,7] <- 19
  x[x[,3] <= 0.05,7] <- 20
  return(x)
}



#' aux_function_samp_b
#'
#' Auxiliar function for sample_structures_b
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' aux_function_samp_b(x)
#' }
aux_function_samp_b <- function(x){
  try({
    x[min(which(x[,8] == "No")),8] <- "Yes"
    return(x)
  },silent = TRUE)
  return(x)
}


#' sample_structures_b
#'
#' Sample n structures from each group. Uses sample_structures_a output to select n similarities from each group.
#'
#'
#' @param x output form sample_structures_a
#' @param rounds Number of elements to be selected from each group.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' sample_structures_b(x,rounds = 2)
#' }
sample_structures_b <- function(x,rounds = 2){
  x$sampled <- rep("No",nrow(x))
  x_splitted <- split(x,as.factor(x[,7]))
  for(i in 1:rounds){
    x_splitted <- lapply(x_splitted,aux_function_samp_b)
  }
  x <- do.call("rbind",x_splitted)
  return(x)
}


#' remove_ceros
#'
#' Function that removes rows and columns that only contain zeros from the similiarity matrices.
#'
#' @param x List of similarity matrices derived from the complete_matrices functions.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' remove_ceros(x)
#' }
remove_ceros <- function(x){
  y <- x[[4]]
  y <- which(rowSums(y == 0) == ncol(y))
  list_out <- lapply(x,function(x) x[-y,-y])
  return(list_out)
}





#' search_bcs
#'
#' Given a specific barcode, the function carries out searches in the database and returns the barcodes displaying the higher simmilarity with the query barcode.
#'
#' @param bc_to_test Problem barcode.
#' @param bc_to_test_name Probelm barcode name.
#' @param level_one Level one clustering data of the database derived from get_level_one function
#' @param level_two Level two clustering data from the database derived from the get_level_two funciton.
#' @param dir_bcs Directory where the barcodes for the datbase structures are stored.
#' @param n_cores Number of cores available for parallel colmuting.
#' @param lev_one_thr Similarity threshold for level one clustering medioids.
#' @param lev_two_thr Similarity threshold for level two clustering medioids.
#' @param lev_three_thr Similarity thresuld for level three results.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' search_bcs(bc_to_test,bc_to_test_name,level_one,level_two,dir_bcs = "/media/data/jaume/DocTDA/Data/Round_1/Barcodes",n_cores,lev_one_thr = 0.8,lev_two_thr = 0.85,lev_three_thr = 0.90)
#' }
search_bcs <- function(bc_to_test,bc_to_test_name,level_one,level_two,dir_bcs = "/media/data/jaume/DocTDA/Data/Round_1/Barcodes",n_cores,lev_one_thr = 0.8,lev_two_thr = 0.85,lev_three_thr = 0.90){

  # Prepare medioids to compute similarities.

  dir_bcs <- dir(dir_bcs,full.names = TRUE)
  bio_ass <- substr(gsub(".*/","",dir_bcs),1,10)
  df_data <- data.frame(bio_ass,dir_bcs)
  df_data_filt <- df_data[df_data$bio_ass %in% names(level_one),]
  df_data_filt_ord <- df_data_filt[match(names(level_one),df_data_filt$bio_ass),]

  # Compute distances with medioids of level 1 clusters.

  doParallel::registerDoParallel(n_cores)
  list_out <- foreach::foreach(j = 1:nrow(df_data_filt_ord)) %dopar%
    {
      temp_bc <- get(load(file = df_data_filt_ord[j,2]))
      temp_sim <- compute_sim_barcode_pair(bc_to_test,temp_bc)
    }
  df_lev_one <- cbind(df_data_filt_ord,do.call("rbind",list_out))
  df_lev_one$Av_Sim <- rowSums(df_lev_one[,c(3,4,5)])/3
  colnames(df_lev_one)[3] <- "Sim_0"
  colnames(df_lev_one)[4] <- "Sim_1"
  colnames(df_lev_one)[5] <- "Sim_2"

  foreach::registerDoSEQ()

  # Prepare level two mediods to compute similarities.

  df_lev_one_filt <- df_lev_one[df_lev_one$Av_Sim > lev_one_thr,,drop = F]
  lev_two_medioids <- c()
  list_out <- list()
  for(i in 1:nrow(df_lev_one_filt)){
    list_out[[i]] <- data.frame(rep(df_lev_one_filt[i,1],length(unlist(lapply(level_two,names)[[df_lev_one_filt[i,1]]]))),unlist(lapply(level_two,names)[[df_lev_one_filt[i,1]]]))
  }
  df_lev_two <- do.call("rbind",list_out)
  colnames(df_lev_two) <- c("lev_one_meds","lev_two_meds")
  df_data_filt_2 <- df_data[df_data[,1] %in% df_lev_two[,2],]
  df_data_filt_2_ord <- df_data_filt_2[match(df_lev_two[,2],df_data_filt_2[,1]),]
  df_lev_two <- cbind(df_lev_two,df_data_filt_2_ord)

  # Compute distances to level two mediods.

  doParallel::registerDoParallel(n_cores)
  list_out <- foreach::foreach(j = 1:nrow(df_lev_two)) %dopar%
    {
      temp_bc <- get(load(file = df_lev_two[j,4]))
      temp_sim <- compute_sim_barcode_pair(bc_to_test,temp_bc)
    }

  df_lev_two <- cbind(df_lev_two,do.call("rbind",list_out))
  df_lev_two$Av_Sim <- rowSums(df_lev_two[,c(5,6,7)])/3
  colnames(df_lev_two)[5] <- "Sim_0"
  colnames(df_lev_two)[6] <- "Sim_1"
  colnames(df_lev_two)[7] <- "Sim_2"

  foreach::registerDoSEQ()

  # Prepare level three data for computing similarities.

  lev_three_temp <- df_lev_two
  lev_three_temp_filt <- lev_three_temp[lev_three_temp$Av_Sim > lev_two_thr,]
  list_out <- list()
  for(i in 1:nrow(lev_three_temp_filt)){
    lev_one_id <- lev_three_temp_filt[i,1]
    lev_two_id <- lev_three_temp_filt[i,2]
    df_to_save <- level_two[[lev_one_id]][[lev_two_id]]
    df_to_save <- cbind(rep(lev_one_id,nrow(df_to_save)),rep(lev_two_id,nrow(df_to_save)),df_to_save)
    df_to_save <- df_to_save[,-3]
    colnames(df_to_save) <- c("lev_one_meds","lev_two_meds","lev_three")
    list_out[[i]] <- df_to_save
  }
  lev_three_temp <- do.call("rbind",list_out)
  df_data_filt_3 <- df_data[df_data$bio_ass %in% lev_three_temp$lev_three,]
  df_data_filt_3_ord <- df_data_filt_3[match(lev_three_temp$lev_three,df_data_filt_3$bio_ass),]
  df_lev_three <- cbind(lev_three_temp,df_data_filt_3_ord)

  # Compute distances for level three data.

  doParallel::registerDoParallel(n_cores)
  list_out <- foreach::foreach(j = 1:nrow(df_lev_three)) %dopar%
    {
      temp_bc <- get(load(file = df_lev_three[j,5]))
      temp_sim <- compute_sim_barcode_pair(bc_to_test,temp_bc)
    }

  foreach::registerDoSEQ()


  df_lev_three <- cbind(df_lev_three,do.call("rbind",list_out))
  df_lev_three$Av_Sim <- rowSums(df_lev_three[,c(6,7,8)])/3
  colnames(df_lev_three)[6] <- "Sim_0"
  colnames(df_lev_three)[7] <- "Sim_1"
  colnames(df_lev_three)[8] <- "Sim_2"
  df_lev_three <- df_lev_three[df_lev_three$Av_Sim > lev_three_thr,]
  df_lev_three <- df_lev_three[,-4]
  query_Structure <- rep(bc_to_test_name,nrow(df_lev_three))
  df_lev_three <- cbind(query_Structure,df_lev_three)

  return(df_lev_three)
}

