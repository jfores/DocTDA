
#' Extract Information from DrugBank
#'
#' @param path_to_drugbank_xml Path to the file containing the xml format drugbank database.
#' @param path_to_data Path to the folder were the output data should be stored.
#' @param path_to_uniprot_data Path to the file containing the mapping between uniprot IDs and the PDB structures.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' extract_DrugBank_data(alpha_coord)
#' }
#'
extract_DrugBank_data <- function(path_to_drugbank_xml,path_to_data,path_to_uniprot_data,verbose = FALSE){
  if(verbose){
    print("Starting process...")
  }
  bool_A <- paste(path_to_data,"/Drugs.Rda",sep="") %in% dir(path_to_data,full.names = T)
  bool_B <- paste(path_to_data,"/drug_groups.Rda",sep="") %in% dir(path_to_data,full.names = T)
  bool_C <- paste(path_to_data,"/db_targets.Rda",sep="") %in% dir(path_to_data,full.names = T)
  bool_D <- paste(path_to_data,"/targets_polypeptides_db.Rda",sep="") %in% dir(path_to_data,full.names = T)
  bool_E <- paste(path_to_data,"/drugs_general_information.Rda",sep="") %in% dir(path_to_data,full.names = T)
  if(!(bool_A & bool_B & bool_C & bool_D & bool_E)){
    if(verbose){
      print("Parsing drugbank...")
    }
    db_xml <- dbparser::read_drugbank_xml_db(path_to_drugbank_xml)
    drugs <- dbparser::drugs()
    drugs_general_information <- drugs$general_information
    names(drugs_general_information)[1] <- "DB_ID"
    save(file = paste(path_to_data,"/Drugs.Rda",sep=""),drugs)
    save(file = paste(path_to_data,"/drugs_general_information.Rda",sep=""),drugs_general_information)
    drug_groups <- dbparser::drug_groups()
    save(file = paste(path_to_data,"/drug_groups.Rda",sep=""),drug_groups)
    db_targets <- dbparser::targets()
    names(db_targets)[1] <- "Protein_ID_2"
    save(file = paste(path_to_data,"/db_targets.Rda",sep=""),db_targets)
    targets_polypeptides_db <- dbparser::targets_polypeptides()
    save(file = paste(path_to_data,"/targets_polypeptides_db.Rda",sep=""),targets_polypeptides_db)
  }else{
    if(verbose){
      print("Loading data...")
    }
    drugs <- get(load(file = paste(path_to_data,"/Drugs.Rda",sep="")))
    drug_groups <- get(load(file = paste(path_to_data,"/drug_groups.Rda",sep="")))
    db_targets <- get(load(file = paste(path_to_data,"/db_targets.Rda",sep="")))
    targets_polypeptides_db <- get(load(file = paste(path_to_data,"/targets_polypeptides_db.Rda",sep="")))
    drugs_general_information <- get(load(file = paste(path_to_data,"/drugs_general_information.Rda",sep="")))
  }
  print(names(targets_polypeptides_db)[1])
  print(names(targets_polypeptides_db)[20])
  names(db_targets)[1] <- "Protein_ID_2"
  names(targets_polypeptides_db)[1] <- "Protein_ID"
  names(targets_polypeptides_db)[20] <- "Protein_ID_2"
  print(names(db_targets))
  print(names(targets_polypeptides_db))
  joined_one <- dplyr::left_join(db_targets,targets_polypeptides_db,by = "Protein_ID_2")
  names(joined_one)[6] <- "DB_ID"
  if(verbose){
    print("Joinig data...")
  }
  names(drug_groups)[2] <- "DB_ID"
  joined_two <- dplyr::left_join(drugs_general_information,drug_groups,by = "DB_ID")
  joined_trhee <- dplyr::left_join(joined_two,joined_one,by = "DB_ID")
  joined_trhee <- joined_trhee[!is.na(joined_trhee$Protein_ID),]
  joined_trhee <- joined_trhee[,c("DB_ID","other_keys","type","name","description","group","Protein_ID","gene_name","Protein_ID_2","name.x","organism.x","known_action","general_function","specific_function","cellular_location","transmembrane_regions","signal_regions","theoretical_pi","molecular_weight","chromosome_location","amino_acid_sequence","amino_acid_format")]
  if(verbose){
    print("Getting pdb structures from uniprot mapping...")
  }
  uniprot_to_pdb <- read.csv(path_to_uniprot_data,sep = "\t")
  names(uniprot_to_pdb)[9] <- "PDB"
  temp_A <- dplyr::mutate(uniprot_to_pdb,PDB = strsplit(as.character(PDB), ";"))
  uniprot_to_pdb <- tidyr::unnest(temp_A,PDB)
  # %>%   mutate(PDB = strsplit(as.character(PDB), ";")) %>%  unnest(PDB) -> uniprot_to_pdb
  names(uniprot_to_pdb)[1] <- "Protein_ID"
  if(verbose){
    print("Joining data...")
  }
  joined_four <- dplyr::left_join(joined_trhee,uniprot_to_pdb,by = "Protein_ID")
  joined_four <- joined_four[!is.na(joined_four$PDB),]
  save(file = paste(path_to_data,"/joined_four.Rda",sep=""),joined_four)
}


#' Download a large number of PDB files.
#'
#' @param pdb_ids_vec PDB ids for the structures to download.
#' @param path_to_PDB_folder Path to the folder were the structures must be downloaded.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' download_pdb_structures(pdb_ids_vec,path_to_PDB_folder)
#' }
#'
download_pdb_structures <- function(pdb_ids_vec,path_to_PDB_folder,verbose = FALSE){
  pdb_ids_vec <- base::unique(pdb_ids_vec)
  seq_to_download <- base::seq(1,base::length(pdb_ids_vec),1000)
  for(i in 1:base::length(seq_to_download)){
    if(i < base::length(seq_to_download)){
      from_id <- seq_to_download[i]
      to_id <-  seq_to_download[i+1] -1
      if(verbose){
        print(paste("Downloading structures from: ",from_id," to: ",to_id,".",sep=""))
      }
      bio3d::get.pdb(pdb_ids_vec[base::seq(from_id,to_id)],path = path_to_PDB_folder)
    }else{
      from_id <- seq_to_download[i]
      to_id <-  base::length(pdb_ids_vec)
      if(verbose){
        print(paste("Downloading structures from: ",from_id," to: ",to_id,".",sep=""))
      }
      bio3d::get.pdb(pdb_ids_vec[base::seq(from_id,to_id)],path = path_to_PDB_folder)
    }
  }
}


#' Get and clean biological assemblies.
#'
#' @param pdb_object A pdb class object.
#' @param directory_to_write Path to the directory where clean biological units should be written.
#' @param structure_name Name of the original pdb strcuture.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' get_bio_clean(pdb_object,directory_to_write,structure_name)
#' }
#'
get_bio_clean <- function(pdb_object,directory_to_write,structure_name){
  if(is.null(pdb_object$remark)){
    temp_clean <- bio3d::atom.select(pdb_object,"protein",value=TRUE)
    temp_clean <- bio3d::clean.pdb(temp_clean,rm.lig = TRUE,rm.wat = TRUE,fix.aa = TRUE)
    bio3d::write.pdb(temp_clean,file = paste(directory_to_write,"/",structure_name,"_","bio_",1,"_clean.pdb",sep=""))
  }
  else{
  biounit_list <- bio3d::biounit(pdb_object)
  for(i in 1:length(biounit_list)){
    temp_clean <- bio3d::atom.select(biounit_list[[i]],"protein",value=TRUE)
    temp_clean <- bio3d::clean.pdb(temp_clean,rm.lig = TRUE,rm.wat = TRUE,fix.aa = TRUE)
    biounit_list[[i]] <- temp_clean
    bio3d::write.pdb(biounit_list[[i]],file = paste(directory_to_write,"/",structure_name,"_","bio_",i,"_clean.pdb",sep=""))
    }
  }
}

#' Extract Bio And clean
#'
#' Extracts biological assemblies from pdb structures and removes waters, ligans, and non-proteic atoms.
#'
#' @param pdb_structures_paths a vector containing the complete paths to the pdb structures.
#' @param directory_to_write a string including the destination folder where clean biological assemblies must be stored.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' extract_and_clean_multi(pdb_structures_paths,directory_to_write)
#' }
#'
extract_and_clean_multi <- function(pdb_structures_paths,directory_to_write){
  print("Extracting biological assemblies and removing ligands, waters, and non-proteic atoms...")
  pb <- utils::txtProgressBar(min = 0, max = length(pdb_structures_paths))
  for(i in 1:length(pdb_structures_paths)){
    utils::setTxtProgressBar(pb, i)
    try({
    pdb_to_read <- pdb_structures_paths[i]
    structure_name <- strsplit(pdb_to_read,"\\/")[[1]][length(strsplit(pdb_to_read,"\\/")[[1]])]
    structure_name <- strsplit(structure_name,"\\.")[[1]][1]
    temp_pdb <- bio3d::read.pdb(file = pdb_to_read)
    get_bio_clean(temp_pdb,directory_to_write,structure_name)
    })
  }
  close(pb)
}


#' get_bio_assembly_information
#'
#' @param path_to_bio_ass Pathway to biological assemblies folder.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' get_bio_assembly_information(path_to_bio_ass)
#' }
#'
get_bio_assembly_information <- function(path_to_bio_ass){
  bio_ass_pahts <- dir(path_to_bio_ass,full.names = TRUE)
  original_structure <- unlist(lapply(strsplit(unlist(lapply(strsplit(bio_ass_pahts,"\\/"),function(x)x[length(x)])),"_"),function(x)x[1]))
  data_frame_structure_path <- data.frame(original_structure,bio_ass_pahts)
  list_out <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(bio_ass_pahts), style = 3)
  for(i in 1:nrow(data_frame_structure_path)){
    utils::setTxtProgressBar(pb, i)
    temp_pdb <- bio3d::read.pdb(data_frame_structure_path[i,2])
    n_chains <- length(unique(temp_pdb$atom$chain))
    name_chains <- paste(unique(temp_pdb$atom$chain),collapse = ", ")
    n_alpha_carbons <- sum(temp_pdb$calpha)
    n_alpha_carbons_by_chain <- paste(lapply(split(temp_pdb$atom[temp_pdb$calpha,],as.factor(temp_pdb$atom[temp_pdb$calpha,]$chain)),nrow),collapse = ", ")
    list_out[[i]] <- c(n_alpha_carbons,n_chains,name_chains,n_alpha_carbons_by_chain)
  }
  close(pb)
  data_frame_lengths <- base::do.call("rbind",list_out)
  colnames(data_frame_lengths) <- c("n_alpha","n_chains","name_chains","leng_n_chains")
  data_frame_paths_and_structures <- cbind(data_frame_structure_path,data_frame_lengths)
  data_frame_paths_and_structures$PDB_BioAss <- gsub("_clean.pdb","",unlist(lapply(strsplit(data_frame_paths_and_structures$bio_ass_pahts,"\\/"),function(x) x[length(x)])))
  data_frame_paths_and_structures$n_alpha <- as.numeric(data_frame_paths_and_structures$n_alpha)
  colnames(data_frame_paths_and_structures)[1] <- "PDB"
  return(data_frame_paths_and_structures)
}



#' rename_drug_sdf
#'
#' Renames and SDF file generated using the split_sdf.py script using its DrugBank ID.
#'
#' @param path_to_sdfs Path to the location of sdf 3D drug structures.
#' @param path_to_out Path were the new files must be written.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' rename_drug_sdf(path_to_sdfs,path_to_out)
#' }
#'
rename_drug_sdf <- function(path_to_sdfs,path_to_out){
  dir_drug_3D <- dir(path_to_sdfs,full.names = T)
  print(paste(path_to_sdfs,"/3D structures.sdf",sep=""))
  drug_files <- dir_drug_3D[grepl("3D.*sdf",dir_drug_3D)]
  print(table(grepl(paste(path_to_sdfs,"/3D structures.sdf",sep=""),drug_files)))
  drug_files <- drug_files[!grepl(paste(path_to_sdfs,"/3D structures.sdf",sep=""),drug_files)]
  drug_names <- list()
  new_files <- c()
  for(i in 1:length(drug_files)){
    fileName = drug_files[i]
    con=file(fileName,open="r")
    line=readLines(con)
    drug_name <- line[which(line == "> <DATABASE_ID>") + 1]
    drug_name <- paste(drug_name,collapse = ", ")
    drug_names[[i]] <- drug_name
    close.connection(con)
    new_file_name <- paste(path_to_out,"/",drug_name,".sdf",sep="")
    new_files <- c(new_files,new_file_name)
  }
  df_return <- data.frame(drug_files,new_files)
  file.rename(from = df_return$drug_files,to = df_return$new_files)
  return(df_return)
}


#' aux_func
#'
#' auxiliar function to carry out the random sampling to generate different project rounds.
#'
#' @param x a df with the same column structure than joined_5.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' aux_function(x)
#' }
#'
aux_function <- function(x){
  sampled <- sample(1:nrow(x),1)
  bio_sampled <- as.character(x[sampled,"PDB_BioAss"])
  #bool <- control_file_2$PDB_BioAss == bio_sampled
  return(bio_sampled)
}


#' assign_rounds_for_bc
#'
#' Generate randmly sampled round to carry out the project in different steps.
#'
#' @param path_to_joined_five Paht to the joined_five file.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' assign_rounds_for_bc(path_to_joined_five)
#' }
#'
assign_rounds_for_bc <- function(path_to_joined_five){
  joined_five <- get(load(file = path_to_joined_five))
  control_file <- joined_five
  control_file$BarCode_Computed <- rep("No",nrow(control_file))
  control_file$Assigned_Round <- rep(NA,nrow(control_file))
  control_file <- control_file[!is.na(control_file$PDB_BioAss),]
  control_file_2 <- control_file
  sampling_cont <- TRUE
  counter_num <- 1
  while(sampling_cont){
    print(counter_num)
    to_split <- control_file_2[is.na(control_file_2$Assigned_Round),]
    splitted <- split(to_split,as.factor(to_split$Protein_ID))
    temp_strucutres <- unlist(lapply(splitted,aux_function))
    control_file_2[control_file_2$PDB_BioAss %in% temp_strucutres,"Assigned_Round"] <- counter_num
    if(all(!is.na(control_file_2$Assigned_Round))){
      sampling_cont <- FALSE
    }
    counter_num <- counter_num + 1
  }
  return(control_file_2)
}
