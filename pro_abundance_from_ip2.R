##############################################################
# EXTRACT NSAF VALUE

wrangle_dta_filter = function(a_dta_select_filter){
  # note the time
  dw_start_time = Sys.time();
  
  raw_data = read_lines(a_dta_select_filter);
  
  # 1 means protein line; 2 means peptide line; 0 means neither
  protein_peptide_line_indicator = rep(0, length(raw_data));
  protein_peptide_line_indicator[(raw_data %>% is_dta_protein_line())] = 1;

  protein_header_regexp = "Locus\tSequence Count\t"
  raw_protein_header = # protein header
    (raw_data[(raw_data %>% str_starts(protein_header_regexp))] %>%  
       str_split("\t"))[[1]];
  table_data = 
    tbl_df(raw_data) %>% 
    separate(value, into = raw_protein_header, sep = "\t");

  table_data_locus = table_data$Locus;
  table_data_spec = table_data$`Spectrum Count`;
  table_data_nsaf = table_data$NSAF;
  rm(table_data);
  
  
  # create a peptide central table
  proteins = "";
  num_peptide_lines = -1;
  # dta_dw = tibble();
  
  nsaf_vec = rep(NA, length(raw_data));
  spec_vec = rep(NA, length(raw_data));
  protein_vec = rep(NA, length(raw_data));
  count = 1
  
  for(i in 1:length(protein_peptide_line_indicator)){
    
    if(protein_peptide_line_indicator[i] == 1){ # protein line: 1
      protein_vec[count] = table_data_locus[i];
      nsaf_vec[count] = table_data_nsaf[i];
      spec_vec[count] =  table_data_spec[i];
      count = count +1
      num_peptide_lines =  num_peptide_lines + 1
    }
  }
  
  # remove NA
  rm_index = is.na(nsaf_vec)
  
  dta_dw = 
    tibble(spec_num = spec_vec[!rm_index],
           nsaf = nsaf_vec[!rm_index],
           protein = protein_vec[!rm_index]);
  return(dta_dw)
}


setwd("~/Desktop/Projects/peipei/OPTM/DATA/2019 Ding new mouse data/")

all_dw = tibble()
for(day in c('01h','01d','03d','07d', '14d')){
  for(rep in c('1','2','3')){
    dta_dw = wrangle_dta_filter(paste0('./ProCarbonylation/Proline Carbonylation_C57-',
                                       day,'-S',rep, '.txt'))
    dta_dw$day = day
    dta_dw$rep = rep
    dta_dw$strain = 'C57'
    all_dw = all_dw %>%
      bind_rows(dta_dw)
    
  }
}

write_csv(all_dw, '../../New_calculation/Data/C57_abundance.csv')











