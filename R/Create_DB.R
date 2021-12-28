
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
extract_DrugBank_data <- function(path_to_drugbank_xml,path_to_data,path_to_uniprot_data){
  db_xml <- dbparser::read_drugbank_xml_db(path_to_drugbank_xml)
  drugs <- dbparser::drugs()
  drugs_general_information <- drugs$general_information
  names(drugs_general_information)[1] <- "DB_ID"
  save(file = paste(path_to_data,"/Drugs.Rda",sep=""),drugs)
  drug_groups <- dbparser::drug_groups()
  save(file = paste(path_to_data,"/drug_groups.Rda",sep=""),drug_groups)
  db_targets <- dbparser::targets()
  save(file = paste(path_to_data,"/db_targets.Rda",sep=""),db_targets)
  names(db_targets)[1] <- "Protein_ID_2"
  targets_polypeptides_db <- dbparser::targets_polypeptides()
  save(file = paste(path_to_data,"/targets_polypeptides_db.Rda",sep=""),targets_polypeptides_db)
  names(targets_polypeptides_db)[1] <- "Protein_ID"
  names(targets_polypeptides_db)[20] <- "Protein_ID_2"
  joined_one <- dplyr::left_join(db_targets,targets_polypeptides_db,by = "Protein_ID_2")
  names(joined_one)[6] <- "DB_ID"
  joined_two <- dplyr::left_join(drugs_general_information,drug_groups,by = "DB_ID")
  joined_trhee <- dplyr::left_join(joined_two,joined_one,by = "DB_ID")
  joined_trhee <- joined_trhee[!is.na(joined_trhee$Protein_ID),]
  joined_trhee <- joined_trhee[,c("DB_ID","other_keys","type","name","description","group","Protein_ID","gene_name","Protein_ID_2","name.x","organism.x","known_action","general_function","specific_function","cellular_location","transmembrane_regions","signal_regions","theoretical_pi","molecular_weight","chromosome_location","amino_acid_sequence","amino_acid_format")]
  uniprot_to_pdb <- read.csv(path_to_uniprot_data,sep = "\t")
  names(uniprot_to_pdb)[9] <- "PDB"
  uniprot_to_pdb %>%   mutate(PDB = strsplit(as.character(PDB), ";")) %>%  unnest(PDB) -> uniprot_to_pdb
  names(uniprot_to_pdb)[1] <- "Protein_ID"
  joined_four <- dplyr::left_join(joined_trhee,uniprot_to_pdb,by = "Protein_ID")
  joined_four <- joined_four[!is.na(joined_four$PDB),]
  save(file = paste(paht_to_folder,"/joined_four.Rda",sep=""),joined_four)

}

