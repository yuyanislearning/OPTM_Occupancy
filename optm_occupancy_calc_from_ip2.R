# rm(list = ls(all = T)); 


####################################################################################################
### call libraries
####################################################################################################
library(tidyverse); # for data wrangling
library(seqinr); # to read fasta file


####################################################################################################
### set parameters: working dir, folder names,  desired output file name, and database fasta file name 
####################################################################################################
setwd("/Users/howardchoi/Desktop/optm_mouse_2019/");

folder_names = 
  c("CysCarbonylation",
  "CysSulfonylation");
  # c("MetSulfoxidation",
  #   "ArgHydroxylation",
  #   "AsnHydroxylation",
  #   "AspHydroxylation",
  #   "CysCarbonylation",
  #   "CysSulfinylation",
  #   "CysSulfonylation",
  #   "LysCarbonylation",
  #   "LysHydroxylation",
  #   "MetSulfonation",
  #   "MetSulfoxidation",
  #   "PheHydroxylation",
  #   "ProCarbonylation",
  #   "TryHydroxylation",
  #   "TyrHydroxylation",
  #   "ValHydroxylation")
    
output_csv_name = "optm_occupancy_mouse_2019.csv"

database_file_name = "_ip2_ip2_data_nhlbi_database__SwissProt_mouse_06-12-2016_reversed.fasta";
# database_file_name = "_ip2_ip2_data_nhlbi_database__SwissProt_Human_contaminant_05-20-2020_reversed.fasta";


####################################################################################################
####################################################################################################
is_dta_protein_line = function(lines){
  
  protein_regexp = "\\w*sp\\|";
  contaminant_regexp = "contaminant";
  
  bool_lines = ((lines %>% str_starts(protein_regexp)) | 
              (lines %>% str_starts(contaminant_regexp)))
  
  return (bool_lines)
}

is_dta_peptide_line = function(lines){
  
  peptide_regexp = "[\\*\t]";
  summary_header_regexp = "\tProteins\tPeptide IDs\tSpectra";
  
  bool_lines = ((lines %>% str_starts(peptide_regexp)) & 
              !(lines %>% str_starts(summary_header_regexp)))
  
  return (bool_lines)
}

valid_cys_line = function(sequences){
  
  wrong_cys_regexp = "C[[:upper:]\\.]"
  
  output = !(sequences %>% str_detect(wrong_cys_regexp) %in% TRUE);
  
  return(output)
}

wrangle_dta_select_filter = function(a_dta_select_filter){
  
  dw_start_time = Sys.time();
  # a_dta_select_filter = paste(each_folder_name,each_file_name, sep = "/")
  
  raw_data = read_lines(a_dta_select_filter);
  
  # 1 means protein line; 2 means peptide line; 0 means neither
  protein_peptide_line_indicator = rep(0, length(raw_data));
  protein_peptide_line_indicator[(raw_data %>% is_dta_protein_line())] = 1;
  protein_peptide_line_indicator[(raw_data %>% is_dta_peptide_line())] = 2;
  
  protein_header_regexp = "Locus\tSequence Count\t"
  raw_protein_header = # protein header
    (raw_data[(raw_data %>% str_starts(protein_header_regexp))] %>%  
       str_split("\t"))[[1]];
  
  peptide_header_regexp = "Unique\tFileName\t"
  raw_peptide_header = # peptide header
    (raw_data[(raw_data %>% str_starts(peptide_header_regexp))] %>%  
    str_split("\t"))[[1]];
  
  # convert the dataset to table
  table_data = 
    tbl_df(raw_data) %>% 
    separate(value, into = raw_peptide_header, sep = "\t");
  
  table_data_unique = table_data$Unique;
  table_data_file_name = table_data$FileName;
  table_data_redundancy = table_data$Redundancy;
  table_data_sequence = table_data$Sequence;
  rm(table_data);

  # indicate valid cys peptide line
  if (a_dta_select_filter %>% str_detect("[Cc]ys")){
    valid_cys_line_indicator = valid_cys_line(table_data_sequence);
  }
  else{
    valid_cys_line_indicator = rep(FALSE, length(table_data_sequence));
  }
  
  proteins = "";
  num_peptide_lines = -1;
  # dta_dw = tibble();
  
  spec_num_vec = c();
  redundancy_vec = c();
  modified_peptide_vec = c();
  protein_vec = c();
  
  for(i in c(1:length(protein_peptide_line_indicator))){
    
    if(protein_peptide_line_indicator[i] == 1){ # protein line: 1
      each_protein = table_data_unique[i];
      
      if(num_peptide_lines != 0){ # first line of protein
        proteins = each_protein
        }
      else{
        proteins = paste(proteins, each_protein, sep = ";")
        }
      num_peptide_lines = 0;
      }
    else if((protein_peptide_line_indicator[i] == 2)){ # peptide line: 2
      
      if(valid_cys_line_indicator[i]){
        spec_num_vec = c(spec_num_vec, table_data_file_name[i]);
        redundancy_vec = c(redundancy_vec, table_data_redundancy[i]);
        modified_peptide_vec = c(modified_peptide_vec, table_data_sequence[i]);
        protein_vec = c(protein_vec, proteins);
      }
      
      num_peptide_lines =  num_peptide_lines + 1
      }
    
    # if(i%%10000 == 0){
    #   print(paste(i, "/", length(protein_peptide_line_indicator), sep = " "))
    #   }
    }
  
  dta_dw = 
    tibble(spec_num = spec_num_vec,
           redundancy = redundancy_vec,
           modified_peptide = modified_peptide_vec,
           protein = protein_vec);
  
  max_protein_num = (dta_dw$protein %>% str_count(";") %>% max()) + 1
  sep_protein_header = paste("sep", c(1:max_protein_num), sep = "_");
  
  dta_reformat = 
    dta_dw %>% 
    separate(protein, into = sep_protein_header, sep = ";") %>% 
    gather(seps, protein, sep_protein_header) %>% 
    filter(!is.na(protein)) %>% 
    select(-seps);
  
  dw_end_time = Sys.time();
  print(paste(a_dta_select_filter, "reading time:", round(dw_end_time - dw_start_time, digits = 2), sep = " "))
  
  return(dta_reformat)
}


prepare_db_fasta = function(a_db_file){
  fasta_dta = read.fasta(file = a_db_file,
                      seqtype = "AA", as.string = TRUE);
  
  db_dta_var_list = rownames(summary(fasta_dta)); # row names
  # remove reverse and contaminant
  valid_protein_idx = !((db_dta_var_list %>% str_detect("Reverse")) | (db_dta_var_list %>% str_detect("contaminant"))) 
  
  fasta_dta =
    fasta_dta[valid_protein_idx];
  
  return(fasta_dta)
}

peptide_mapper = function(a_all_peptides, a_db_dta){
  # a_all_peptides = all_peptides
  # a_db_dta = db_dta
  
  a_db_dta_var_list = rownames(summary(a_db_dta)); # row names
  
  start_positions_vec = c();
  end_positions_vec = c();
  proteins_vec = c();
  peptides_vec = c();
  
  # a_all_peptides_with_position_info = tibble();
  
  for (i in c(1:dim(a_all_peptides)[1])){
    each_peptide = a_all_peptides$peptide[i];
    each_protein = a_all_peptides$protein[i];
    
    pro_idx = (a_db_dta_var_list == each_protein)
    
    start_end_positions = (a_db_dta[pro_idx] %>% str_locate_all(each_peptide))[[1]];
    # each_tbl =
    #   (a_db_dta[pro_idx] %>% str_locate_all(each_peptide))[[1]] %>%
    #   tbl_df() %>%
    #   add_column(
    #     protein = each_protein,
    #     peptide = each_peptide
    #   )
    
    if ((start_end_positions %>% dim())[1] < 1){ # I L confusion
      each_peptide_IL =
        each_peptide %>%
        str_replace_all("I|L", "[IL]")
      
      start_end_positions = (a_db_dta[pro_idx] %>% str_locate_all(each_peptide_IL))[[1]];
      # start_positions = start_end_positions[,1];
      # end_positions = start_end_positions[,2];
      # 
      # each_tbl =
      #   (a_db_dta[pro_idx] %>% str_locate_all(each_peptide_IL))[[1]] %>%
      #   tbl_df() %>%
      #   add_column(
      #     protein = each_protein,
      #     peptide = each_peptide
      #   )
    }
    
    # if ((start_end_positions %>% dim())[1] < 1){ 
    #  break 
    # }
    start_positions = start_end_positions[,1];
    end_positions = start_end_positions[,2];
    proteins = rep(each_protein, length(start_positions));
    peptides = rep(each_peptide, length(start_positions));
    
    start_positions_vec = c(start_positions_vec, start_positions);
    end_positions_vec = c(end_positions_vec, end_positions);
    proteins_vec = c(proteins_vec, proteins);
    peptides_vec = c(peptides_vec, peptides);

    # a_all_peptides_with_position_info =
    #   a_all_peptides_with_position_info %>%
    #   bind_rows(each_tbl)
  }
  
  a_all_peptides_with_position_info = 
    tibble(start = start_positions_vec,
         end = end_positions_vec,
         protein = proteins_vec,
         peptide = peptides_vec);
  
  # each_tbl =
  #   (a_db_dta[pro_idx] %>% str_locate_all(each_peptide_IL))[[1]] %>%
  #   tbl_df() %>%
  #   add_column(
  #     protein = each_protein,
  #     peptide = each_peptide
  #   )
  
  return(a_all_peptides_with_position_info)
}


####################################################################################################
### re-format dta-select filter files
####################################################################################################
### get fasta db data 
db_dta = prepare_db_fasta(database_file_name)


### go over each optm folder and dta select fileter files in each folder
# optm_15_data = tibble(); # store all optms and their occupancies in this data frame

# go over each folder
for(each_folder_name in folder_names){
  
  print(each_folder_name)
  #print(Sys.time()) # check 
  
  dta_files = list.files(path = each_folder_name, pattern = ".txt"); # read all text files in a folder
  dta_formatted_data = tibble(); # store all dta filters for each optm and their occupancies in this data frame
  
  for(each_file_name in dta_files){
    
    print(each_file_name)
    
    dta_data = wrangle_dta_select_filter(paste(each_folder_name,each_file_name, sep = "/"))
    
    if(each_file_name %>% str_detect("[Cc]ys")){
      static_cys_regexp = "C\\(57\\.02146\\)"
      dta_data = 
        dta_data %>% 
        mutate(modified_peptide = (modified_peptide %>% str_replace_all(static_cys_regexp, "C")));
    }
    
    
    # regular expression for clean up
    modification_regexp = "\\([0-9\\.]*\\)";
    dot_regexp = "\\.";
    dash_regexp = "-";
    modification_localization_regexp = "[A-Z]\\([0-9\\.]*\\)";
    
    dta_dw = 
      dta_data %>% 
      #select(spec_num:modified_peptide) %>% 
      #add_column(protein = dta_data$protein$value) %>% 
      mutate(num_diff_modified = (modified_peptide %>% str_count("\\("))) %>% 
      mutate(peptide = (modified_peptide %>% str_replace_all(modification_regexp, ""))) %>% 
      mutate(peptide = (peptide %>% str_replace_all(dot_regexp, ""))) %>% 
      mutate(peptide = (peptide %>% str_replace_all(dash_regexp, ""))) %>% 
      mutate(peptide_for_localization = (modified_peptide %>% str_replace_all(modification_localization_regexp, "\\*"))) %>% 
      mutate(peptide_for_localization = (peptide_for_localization %>% str_replace_all(dot_regexp, ""))) %>% 
      mutate(peptide_for_localization = (peptide_for_localization %>% str_replace_all(dash_regexp, ""))) %>% 
      mutate(diff_mod_local_position = (peptide_for_localization %>% str_locate_all("\\*")));
      
    # unnest modified ones
    unnest_modified = 
      dta_dw %>% 
      filter(num_diff_modified >= 1) %>% 
      unnest(diff_mod_local_position);
    
    # remove Reverse or contaminant protein match
    reverse_regexp = "Reverse"
    contaminant_regexp = "contaminant"
    
    dta_local_position_dw =
      unnest_modified %>% 
      select(spec_num:peptide) %>% 
      add_column(diff_mod_local_position = unnest_modified$diff_mod_local_position[,1]) %>% 
      bind_rows(
        dta_dw %>% 
          filter(num_diff_modified == 0) %>% 
          select(spec_num:peptide) %>% 
          add_column(diff_mod_local_position = NA)) %>% 
      filter(!(protein %>% str_starts(reverse_regexp))) %>%   # remove Reverse
      filter(!(protein %>% str_starts(contaminant_regexp))) %>%   # remove contaminants
      add_column(sample_rep = each_file_name,
                 type = each_folder_name);
      
    dta_formatted_data = 
      dta_formatted_data %>% 
      bind_rows(dta_local_position_dw)
    }
  
  
  ### peptide mapping using db_dta
  all_peptides =
    dta_formatted_data %>% 
    select(peptide, protein) %>% 
    unique();
  
  all_peptides_with_position_info = peptide_mapper(all_peptides, db_dta)
    
  
  ### calculate occupancy
  dta_global_position_dw = 
    dta_formatted_data %>% 
    left_join(all_peptides_with_position_info,
              by = c("protein", "peptide")) %>% 
    mutate(diff_mod_global_position = start + diff_mod_local_position - 1);
  
  optm_data = 
    dta_global_position_dw %>% 
    filter(num_diff_modified > 0) %>% 
    group_by(protein, diff_mod_global_position, sample_rep, type) %>% 
    summarize(spec_count_modified = sum(as.numeric(redundancy))) %>% 
    ungroup();
  
  # optm_dw = tibble();
  protein_vec = c();
  diff_mod_global_position_vec = c();
  sample_rep_vec = c();
  type_vec = c();
  spec_count_modified_vec = c();
  spec_count_unmodified_vec = c();
  
  for (i in c(1:dim(optm_data)[1])){
    each_optm = optm_data[i,];
    
    each_non_modified = 
      dta_global_position_dw %>% 
      filter(protein == each_optm$protein & 
               num_diff_modified == 0 &
               sample_rep == each_optm$sample_rep &
               type == each_optm$type &
               ((start <= each_optm$diff_mod_global_position) & (end >= each_optm$diff_mod_global_position)));
    
    protein_vec = c(protein_vec, each_optm$protein);
    diff_mod_global_position_vec = c(diff_mod_global_position_vec, each_optm$diff_mod_global_position);
    sample_rep_vec = c(sample_rep_vec, each_optm$sample_rep);
    type_vec = c(type_vec, each_optm$type);
    spec_count_modified_vec = c(spec_count_modified_vec, each_optm$spec_count_modified);
    spec_count_unmodified_vec = c(spec_count_unmodified_vec, sum(as.numeric(each_non_modified$redundancy)));

    # optm_dw = 
    #   optm_dw %>% 
    #   bind_rows(each_optm %>% 
    #               add_column(spec_count_unmodified = sum(as.numeric(each_non_modified$redundancy))));
  }
  
  optm_dw = 
    tibble(protein = protein_vec,
           diff_mod_global_position = diff_mod_global_position_vec,
           sample_rep = sample_rep_vec,
           type = type_vec,
           spec_count_modified = spec_count_modified_vec,
           spec_count_unmodified = spec_count_unmodified_vec) %>% 
    filter(!is.na(diff_mod_global_position));
    
  
  # # store all optms
  # optm_15_data =
  #   optm_15_data %>% 
  #   bind_rows(optm_dw)
  optm_dw %>% write_csv(paste(each_folder_name, ".csv", sep = ""));
  # rm(optm_dw);
  #print(Sys.time()) # print the end time
}
    
# dta_formatted_data %>% write_csv("met_sul_dta_prep.csv")




