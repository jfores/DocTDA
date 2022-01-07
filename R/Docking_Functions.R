


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
  pdb_readed <- read.pdb(x)
  center_of_mass <- com(pdb_readed)
  x_size <- range(pdb_readed$atom$x)[2] -range(pdb_readed$atom$x)[1]
  y_size <- range(pdb_readed$atom$y)[2] -range(pdb_readed$atom$y)[1]
  z_size <- range(pdb_readed$atom$z)[2] -range(pdb_readed$atom$z)[1]
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
create_vina_command <- function(x,y,z,w,exhaust){
  temp_rec_one <- gsub("_clean.pdbqt","",strsplit(x,"\\/")[[1]][length(strsplit(x,"\\/")[[1]])])
  temp_lig_one <- gsub(".pdbqt","",strsplit(y,"\\/")[[1]][length(strsplit(y,"\\/")[[1]])])
  out_file <- paste(out_path,temp_rec_one,"_AND_",temp_lig_one,".pdbqt",sep="")
  temp_pdb_file <- gsub("pdbqt","pdb",strsplit(x,"\\/")[[1]][length(strsplit(x,"\\/")[[1]])])
  temp_pdb_file <- paste(w,temp_pdb_file,sep="")
  print(temp_pdb_file)
  temp_log <- gsub("pdbqt","log",out_file)
  print(temp_log)
  list_out <- create_config_pdb(temp_pdb_file)
  print(list_out)
  vina_command <- paste("vina --receptor ",x, " --ligand ",y," --out ",out_file, " --center_x ", round(list_out[[1]],digits = 2)," --center_y ",round(list_out[[2]],digits = 2)," --center_z ",round(list_out[[3]],digits = 2)," --size_x ", round(list_out[[4]],digits = 2)," --size_y ",round(list_out[[5]],digits = 2)," --size_z ",round(list_out[[6]],digits = 2)," --exhaustiveness ",exhaust," --log ",temp_log,sep="")
  print(vina_command)
  #system(vina_command)
}
