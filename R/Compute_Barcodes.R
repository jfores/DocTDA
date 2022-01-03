
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


#' add_pseudo_barcode
#'
#' This function adds a pseudobarcode to each dimension of the generated barcode.
#'
#' @param barcode a barcode obtained from the compute_barcode or compute_barcode_multi functions
#'
#' @return returns a matrix with a barcode for which a pseudobarcoded has been added to each demiension.
#' @export
#'
#' @examples
#' \dontrun{
#' add_pseudo_barcode(barcode)
#' }
add_pseudo_barcode <- function(barcode){
  splited_barcode <- base::split(data.frame(barcode),as.factor(as.character(barcode[,"dimension"])))
  lacking_dimensions <- c(0,1,2)[!c(0,1,2) %in% unique(barcode[,"dimension"])]
  list_out <- list()
  counter <- 1
  counter_added <- 1
  for(i in 1:3){
    print(i)
    if((i -1) %in% lacking_dimensions){
      line_add_lack <- c(i-1,0,0.5)
      list_out[[counter]] <- line_add_lack
      counter <- counter + 1
    }else{
      df_add <- base::rbind(c(i-1,0,0.5),splited_barcode[[counter_added]])
      list_out[[counter]] <- df_add
      counter <- counter + 1
      counter_added <- counter_added + 1
    }
  }
  barcode_pseudo <- as.matrix(base::do.call("rbind",list_out))
  return(barcode_pseudo)
}


#' extract_alpha_chains
#'
#' Extract alpha chains from a pdb structure. Auxiliar function for compute_bc_round.
#'
#' @param structure
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' extract_alpha_chains(structure)
#' }
extract_alpha_chains <- function(structure){
  list_out <- list()
  chains_names <- unique(structure$atom$chain)
  for(i in 1:length(chains_names)){
    ca.inds <- bio3d::atom.select(structure,string = "calpha",chain = chains_names[i])
    struc_filt <- bio3d::trim.pdb(structure,ca.inds)$atom[,c("x","y","z")]
    list_out[[i]] <- struc_filt
  }
  names(list_out) <- chains_names
  return(list_out)
}



#' compute_bc_round
#'
#' @param x a dataframe with the whole database information.
#' @param round_n the number of the round to be computed.
#' @param dir_bc the directory where barcode files must be stored.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' compute_bc_round(x,round_n,dir_bc)
#' }
compute_bc_round <- function(x,round_n,dir_bc){
  print(x)
  round_data <- x[x$Assigned_Round == round_n,]
  round_data_filt <- unique(data.frame(round_data[,c("n_alpha","bio_ass_pahts")]))
  for(i in 1:nrow(round_data_filt)){
    print(i)
    try({
      if(round_data_filt[i,1] <= 1000){
        name_file <- round_data_filt[i,2]
        root_file <- gsub("_clean.pdb","",strsplit(name_file,"\\/")[[1]][length(strsplit(name_file,"\\/")[[1]])])
        if(sum(grepl(root_file,dir(dir_bc))) < 1){
          pdb_temp <- bio3d::read.pdb(round_data_filt[i,2])
          alpha_carbs <- extract_multiple_alpha(list(pdb_temp))
          barcode <- compute_homology(alpha_carbs[[1]])
          file_to_save <- paste(dir_bc,"/",root_file,".Rda",sep = "")
          save(file = file_to_save,barcode)}
      }else{
        pdb_temp <- bio3d::read.pdb(round_data_filt[i,2])
        list_of_chains <- extract_alpha_chains(pdb_temp)
        name_file <- round_data_filt[i,2]
        root_file <- gsub("_clean.pdb","",strsplit(name_file,"\\/")[[1]][length(strsplit(name_file,"\\/")[[1]])])
        if(sum(grepl(root_file,dir(dir_bc))) < 1){
          for(j in 1:length(list_of_chains)){
            barcode <- compute_homology(list_of_chains[[j]])
            file_to_save <- paste(dir_bc,"/",root_file,"_",names(list_of_chains)[j],".Rda",sep = "")
            save(file = file_to_save,barcode)
          }
        }
      }
      gc()
    })
  }
}
