# rm(list = ls(all = T)); 

library(tidyverse);

####################################################################################################
### user sets working dir, folder names, and desired output file name
####################################################################################################
setwd("/Users/howardchoi/Desktop/optm_human_2020/");

folder_names = 
  c("MetSulfoxidation_MSFragger");

output_csv_name = "optm_occupancy_human_2020_msfragger.csv";

database_file_name = "_ip2_ip2_data_nhlbi_database__SwissProt_Human_contaminant_05-20-2020_reversed.fasta";
####################################################################################################
####################################################################################################
psm_dta = read_tsv("MetSulfoxidation_MSFragger/psm_116_2.tsv");
peptide_dta = read_tsv("MetSulfoxidation_MSFragger/peptide.tsv");

psm_dw0 = 
  psm_dta %>% 
  select(peptide = Peptide, 
         modified_peptide = `Modified Peptide`, 
         assigned_modification = `Assigned Modifications`, 
         protein = Protein, 
         uniprot = `Protein ID`);

# remove cys static modififications from assigned_modification column;
cys_static_modification_regexp = "[[:digit:]]+C\\(57\\.0215\\)[, ]*";
endpoint_regexp = ", $"

# remove Reverse or contaminant protein match
reverse_regexp = "Reverse"
contaminant_regexp = "contaminant"

psm_dw =
  psm_dw0 %>% 
  mutate(assigned_diff_modification = assigned_modification %>% str_remove_all(cys_static_modification_regexp) %>% str_remove(endpoint_regexp)) %>% 
  mutate(assigned_diff_modification = ifelse(assigned_diff_modification == "", NA, assigned_diff_modification)) %>% 
  mutate(num_diff_modified = (assigned_diff_modification %>% str_count("\\("))) %>% 
  arrange(desc(num_diff_modified)) %>% 
  filter(!((protein %>% str_detect(reverse_regexp)) | (protein %>% str_detect(contaminant_regexp))));

# ####################################################################################################
# ### read db fasta file
# ####################################################################################################
# library(seqinr)
# 
# db_dta = read.fasta(file = "_ip2_ip2_data_nhlbi_database__SwissProt_Human_contaminant_05-20-2020_reversed.fasta",
#             seqtype = "AA", as.string = TRUE);
# 
# db_dta_var_list = rownames(summary(db_dta)); # row names
# valid_protein_idx = !((db_dta_var_list %>% str_detect(reverse_regexp)) | (db_dta_var_list %>% str_detect(contaminant_regexp)))
# 
# db_dta =
#   db_dta[valid_protein_idx];
# 
# all_peptides =
#   psm_dw %>% 
#   select(Peptide, protein) %>% 
#   unique();
# 
# all_peptides_with_position_info = tibble();
# for (i in c(1:dim(all_peptides)[1])){
#   each_peptide = all_peptides$Peptide[i];
#   each_protein = all_peptides$protein[i];
#   
#   pro_idx = (db_dta_var_list == each_protein)
#   
#   each_tbl =
#     (db_dta[pro_idx] %>% str_locate_all(each_peptide))[[1]] %>% 
#     tbl_df() %>% 
#     add_column(
#       protein = each_protein,
#       peptide = each_peptide
#     )
#   
#   if ((each_tbl %>% dim())[1] <1){
#     
#     each_peptide_IL =
#       each_peptide %>%
#       str_replace_all("I|L", "[IL]")
#     
#     each_tbl =
#       (db_dta[pro_idx] %>% str_locate_all(each_peptide_IL))[[1]] %>% 
#       tbl_df() %>% 
#       add_column(
#         protein = each_protein,
#         peptide = each_peptide
#       )
#   }
#   
#   # if ((each_tbl %>% dim())[1] <1){
#   #   break
#   #   }
#   
#   all_peptides_with_position_info =
#     all_peptides_with_position_info %>%
#     bind_rows(each_tbl)
#   
#   if((i %% 1000) == 0){
#     print(paste(i,"/", dim(all_peptides)[1]))
#   }
# }
# 
# all_peptides_with_position_info %>% 
#   write_csv("all_peptides_with_position_info_msfragger.csv");

####################################################################################################
### read db fasta file
####################################################################################################
database_dta = read_csv("_ip2_ip2_data_nhlbi_database__SwissProt_Human_contaminant_05-20-2020_reversed.csv");

sep_cols = paste("sep_col_", (1:max(psm_dw$num_diff_modified, na.rm = T)), sep = "");

psm_global_position_dw = 
  psm_dw %>%
  filter(num_diff_modified > 1) %>% 
  separate(assigned_diff_modification, into = sep_cols, sep = ",") %>% 
  gather(key = seps, value = assigned_diff_modification, sep_cols) %>% 
  select(-seps) %>% 
  filter(!is.na(assigned_diff_modification)) %>% 
  bind_rows(
    psm_dw %>% 
      filter(num_diff_modified == 1 | is.na(num_diff_modified))
  ) %>% 
  left_join(database_dta,
            by = c("peptide", "protein")) %>% 
  mutate(diff_mod_local_position =  as.numeric(assigned_diff_modification %>% str_replace("[[:alpha:]]\\([0-9\\.]+\\)", ""))) %>% 
  mutate(diff_mod_global_position = (start + (diff_mod_local_position - 1)));
  

####################################################################################################
### occupancy calculations 
####################################################################################################
optm_data = 
  psm_global_position_dw %>% 
  filter(!is.na(num_diff_modified)) %>% 
  group_by(uniprot, diff_mod_global_position) %>% 
  summarize(spec_count_modified = n()) %>% 
  ungroup();

optm_dw = tibble();
for (i in c(1:dim(optm_data)[1])){
  each_optm = optm_data[i,];
  
  each_non_modified = 
    psm_global_position_dw %>% 
    filter(uniprot == each_optm$uniprot & 
             is.na(assigned_diff_modification) & 
             ((start <= each_optm$diff_mod_global_position) & (end >= each_optm$diff_mod_global_position)));
  each_non_modified
  optm_dw = 
    optm_dw %>% 
    bind_rows(each_optm %>% 
                add_column(spec_count_unmodified = dim(each_non_modified)[1]));
}

optm_dw %>% 
  add_column(type = folder_names[1],
             sample = "sample_116_2") %>% 
  write_csv("MetO_116_2.csv")




####################################################################################################
####################################################################################################
optm_data = tibble(); # store all optms and their occupancies in this data frame

# go over each folder
for(each_folder_name in folder_names){
  ###########
  # each_folder_name = folder_names[1]
  ###########
  
  print(each_folder_name)
  
  dta_files = list.files(path = each_folder_name, pattern = ".txt"); # read all text files in a folder
  
  ### Unfiltered DTA: get protein, peptide, and peptide position for all peptides
  ufd_dta_files = dta_files[dta_files %>% str_starts("UFD")];
  ufd_raw_peptide = tibble();
  
  for(each_ufd_file_name in ufd_dta_files){
    ###########
    # each_ufd_file_name = ufd_dta_files[1]
    ###########
    
    print(Sys.time())
    print(each_ufd_file_name)
    
    raw_ufd_data = read_lines(paste(each_folder_name, each_ufd_file_name, sep = "/"));
    
    ufd_protein_peptide_regexp = "[LD]\t"
    ufd_peptide_regexp = "D\t"
    
    raw_ufd_protine_peptide = 
      raw_ufd_data[(raw_ufd_data %>% str_starts(ufd_protein_peptide_regexp))];
    
    var_num = ((raw_ufd_protine_peptide[(raw_ufd_protine_peptide[1:50] %>% str_starts(ufd_peptide_regexp))][1] %>% str_count("\t")) + 1);
    
    raw_ufd_protein_peptide_dw = 
      raw_ufd_protine_peptide %>% 
      tbl_df() %>% 
      separate(value, into = paste("var", c(1:var_num), sep="_"), sep = "\t");
    
    ufd_pro_pep_dw =
      raw_ufd_protein_peptide_dw %>% 
      select(var_1, var_2, var_13, var_14)
    
    rm(raw_ufd_protein_peptide_dw)# to manage memory
    
    proteins = "";
    num_pep_line = -1;
    
    for(i in c(1:dim(ufd_pro_pep_dw)[1])){
      
      each_line_indicator = ufd_pro_pep_dw$var_1[i];
      
      if(each_line_indicator == "L"){ #protein line
        each_protein = ufd_pro_pep_dw$var_2[i];
        
        if(num_pep_line != 0){ # first line of protein
          proteins = each_protein
        }
        else{
          proteins = paste(proteins, each_protein, sep = ";")
        }
        num_pep_line = 0;
      }
      else { #peptide line
        ufd_raw_peptide = 
          ufd_raw_peptide %>% 
          bind_rows(
            tibble(peptide = ufd_pro_pep_dw$var_13[i], 
                   start_position = ufd_pro_pep_dw$var_14[i], 
                   Uniprot = proteins))
        
        num_pep_line =  num_pep_line + 1
      }
    }
    ufd_raw_peptide = unique(ufd_raw_peptide);
    print(paste("ended at: ", Sys.time(), sep = ""))
  }
}

# write_csv(ufd_raw_peptide, "peptide_position_info.csv")

###########################################################################
### wrangle data
###########################################################################
peptide_position_info = read_csv("peptide_position_info.csv");

peptide_position_info = (peptide_position_info %>% unique());

rev_idx = (peptide_position_info$Uniprot %>% str_detect("[Rr]everse")) | (peptide_position_info$Uniprot %>% str_detect("[Cc]ontaminant"))

write_csv(peptide_position_info[!rev_idx,], "valid_peptide_position_info.csv")


