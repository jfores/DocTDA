


#' create_config_pdb
#'
#' Create pdb configuration file.
#'
#' @param x pdb file.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' create_config_pdb(x)
#' }
#'
create_config_pdb <- function(x){
  pdb_readed <- bio3d::read.pdb(x)
  #elements <- bio3d::elements
  #elements <- bio3d::elements
  elements <<- bio3d::elements
  center_of_mass <- bio3d::com(pdb_readed)
  x_size <- range(pdb_readed$atom$x)[2] -range(pdb_readed$atom$x)[1]
  y_size <- range(pdb_readed$atom$y)[2] -range(pdb_readed$atom$y)[1]
  z_size <- range(pdb_readed$atom$z)[2] -range(pdb_readed$atom$z)[1]
  print(list(center_of_mass[1],center_of_mass[2],center_of_mass[3],x_size,y_size,z_size))
  return(list(center_of_mass[1],center_of_mass[2],center_of_mass[3],x_size,y_size,z_size))
}

#' create_vina_command
#'
#' Function to create vina commands.
#'
#' @param x Path to protein pdbqt.
#' @param y Path to ligand pdbqt.
#' @param z Output pathway.
#' @param w Path to receptor pdb files.
#' @param exhaust Exhaustiveness parameter for vina.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' create_vina_command(x,y,z,w,exhaust)
#' }
#'
create_vina_command <- function(x,y,z,w,exhaust = 8){
  temp_rec_one <- gsub("_clean.pdbqt","",strsplit(x,"\\/")[[1]][length(strsplit(x,"\\/")[[1]])])
  temp_lig_one <- gsub(".pdbqt","",strsplit(y,"\\/")[[1]][length(strsplit(y,"\\/")[[1]])])
  out_file <- paste(z,"/",temp_rec_one,"_AND_",temp_lig_one,".pdbqt",sep="")
  temp_pdb_file <- w
  print(temp_pdb_file)
  temp_log <- gsub("pdbqt","log",out_file)
  #print(temp_log)
  #print("Here...")
  list_out <- create_config_pdb(temp_pdb_file)
  #print("Here 2...")
  #print(list_out)
  vina_command <- paste("vina --receptor ",x, " --ligand ",y," --out ",out_file, " --center_x ", round(list_out[[1]],digits = 2)," --center_y ",round(list_out[[2]],digits = 2)," --center_z ",round(list_out[[3]],digits = 2)," --size_x ", round(list_out[[4]],digits = 2)," --size_y ",round(list_out[[5]],digits = 2)," --size_z ",round(list_out[[6]],digits = 2)," --exhaustiveness ",exhaust," --log ",temp_log,sep="")
  if(!(out_file %in% dir(z,full.names = TRUE))){
    print("Running vina...")
    system(vina_command)

  }
  #print(vina_command)
}

#' compute_by_square
#'
#' Performs docking completing a rectangle.
#'
#' @param bio_assay_pdbqt vector with the absolute paths to the pdbqt files of the proteic structures we want to dock.
#' @param drug_pdbqt_paths vector with the absolute paths to the drug pdbqt files.
#' @param output_path string indicating the folder where docking data should be stored.
#' @param bio_assay_paths Vector indicating the absolute paths to the pdb files for the protein structures.
#' @param exhaust Exhaustiveness parameter for vina.
#' @param path_save_control_matrix_path Path where the control matrix will be sotored or has been saved previouly.
#' @param verbose Logical value. Print information or not.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' compute_by_square(bio_assay_pdbqt,drug_pdbqt_paths,output_path,bio_assay_paths,exhaust = 40,path_save_control_matrix_path,verbose = TRUE)
#' }
#'
compute_by_square <- function(bio_assay_pdbqt,drug_pdbqt_paths,output_path,bio_assay_paths,exhaust = 40,save_control_matrix_path,verbose = TRUE){
  if(file.exists( paste(save_control_matrix_path,"/control_matrix.Rda",sep=""))){
    print(paste(save_control_matrix_path,"/control_matrix.Rda",sep=""))
    control_matrix <- get(load(file = paste(save_control_matrix_path,"/control_matrix.Rda",sep="")))
    row_count <- get(load(file =  paste(save_control_matrix_path,"/row_count.Rda",sep="")))
    col_count <- get(load(file =  paste(save_control_matrix_path,"/col_count.Rda",sep="")))
    len_row <- length(bio_assay_pdbqt)
    len_col <- length(drug_pdbqt_paths)
  }else{
    len_row <- length(bio_assay_pdbqt)
    len_col <- length(drug_pdbqt_paths)
    row_count <- 1
    col_count <- 1
    control_matrix <- matrix(FALSE,nrow = len_row,ncol = len_col)
    print(dim(control_matrix))
  }
  continue_process <- TRUE
  while(continue_process){
    if(verbose){
      print(paste("Row: ", row_count,sep = ""))
      print(paste("Col: ", col_count,sep = ""))
    }
    if((row_count * col_count) %% 100 == 0){
      save(file = paste(save_control_matrix_path,"/control_matrix.Rda",sep=""),control_matrix)
      save(file = paste(save_control_matrix_path,"/row_count.Rda",sep=""),row_count)
      save(file = paste(save_control_matrix_path,"/col_count.Rda",sep=""),col_count)
    }
    if(row_count == 1 & col_count == 1){
      try({
        create_vina_command(bio_assay_pdbqt[row_count],drug_pdbqt_paths[col_count],output_path,bio_assay_paths[row_count],exhaust = exhaust)
      })
      control_matrix[row_count,col_count] <- TRUE
      row_count <- row_count + 1
      col_count <- col_count + 1
    }else{
      for(i in 1:col_count){
        try({
          create_vina_command(bio_assay_pdbqt[row_count],drug_pdbqt_paths[i],output_path,bio_assay_paths[row_count],exhaust = exhaust)
        })
        control_matrix[row_count,i] <- TRUE
      }
      for(i in 1:(row_count-1)){
        try({
          create_vina_command(bio_assay_pdbqt[i],drug_pdbqt_paths[col_count],output_path,bio_assay_paths[i],exhaust = exhaust)
        })
        control_matrix[i,col_count] <- TRUE
      }
      row_count <- row_count + 1
      col_count <- col_count + 1
    }
    if((col_count > len_col) | (row_count > len_row)){
      continue_process <- FALSE
    }
  }
  if(len_row > len_col){
    print("Len_row larger len_col")
    continue_process <- TRUE
    while(continue_process){
      if(verbose){
        print(paste("Row: ", row_count,sep = ""))
        print(paste("Col: ", col_count,sep = ""))
      }
      if((row_count * col_count) %% 100 == 0){
        save(file = paste(save_control_matrix_path,"/control_matrix.Rda",sep=""),control_matrix)
        save(file = paste(save_control_matrix_path,"/row_count.Rda",sep=""),row_count)
        save(file = paste(save_control_matrix_path,"/col_count.Rda",sep=""),col_count)
      }
      for(i in 1:len_col){
        print(i)
        try({
          create_vina_command(bio_assay_pdbqt[row_count],drug_pdbqt_paths[i],output_path,bio_assay_paths[row_count],exhaust = exhaust)
        })
        control_matrix[row_count,i] <- TRUE
      }
      row_count <- row_count + 1
      print(row_count)
      print(len_row)
      if(row_count > len_row){
        print("Here we are")
        continue_process <- FALSE
      }
    }
  }
  if(len_col > len_row){
    continue_process <- TRUE
    while(continue_process){
      if(verbose){
        print(paste("Row: ", row_count,sep = ""))
        print(paste("Col: ", col_count,sep = ""))
      }
      if((row_count * col_count) %% 100 == 0){
        save(file = paste(save_control_matrix_path,"/control_matrix.Rda",sep=""),control_matrix)
        save(file = paste(save_control_matrix_path,"/row_count.Rda",sep=""),row_count)
        save(file = paste(save_control_matrix_path,"/col_count.Rda",sep=""),col_count)
      }
      for(i in 1:len_row){
        print(i)
        try({
          create_vina_command(bio_assay_pdbqt[i],drug_pdbqt_paths[col_count],output_path,bio_assay_paths[i],exhaust = exhaust)
        })
        control_matrix[i,col_count] <- TRUE
      }
      col_count <- col_count + 1
      print(row_count)
      print(len_row)
      if(col_count > len_col){
        print("here")
        continue_process <- FALSE
      }
    }
  }
  return(control_matrix)
}



#' parallel_drug_prep
#'
#' Function to prepare drug structures for docking.
#'
#' @param x Path to sdf drug file.
#' @param path_to_mk_prepare Absolute path to mk_prepare_ligand.py scrii
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' parallel_drug_prep(x,path_to_mk_prepare)
#' }
#'
parallel_drug_prep <- function(x,path_to_mk_prepare = "/home/antoniojr/anaconda3/envs/py36/scripts/mk_prepare_ligand.py"){
  command_to_run <- paste(path_to_mk_prepare," -i ",x," -o",gsub("Drug_Name","Drug_PDBQT",gsub("sdf","pdbqt",x))," --pH 7.4")
  print(command_to_run)
  system(command_to_run)
}
