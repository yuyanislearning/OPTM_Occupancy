# rm(list = ls(all = TRUE));

library(tidyverse);  # packages

setwd("/Users/howardchoi/Desktop/OPTMGit/") # set working directory

######################################################################################################################################################
### read 15 optm data
######################################################################################################################################################
optm_data_first_batch = read_csv("optm_15_first_batch.csv");

datasets_second_batch = 
  c("ArgHydroxylation.csv",
    "AsnHydroxylation.csv",
    "AspHydroxylation.csv",
    "LysCarbonylation.csv",
    "LysHydroxylation.csv",
    "MetSulfonation.csv",
    "MetSulfoxidation.csv",
    "PheHydroxylation.csv",
    "ProCarbonylation.csv",
    "TryHydroxylation.csv",
    "TyrHydroxylation.csv",
    "ValHydroxylation.csv",
    "CysCarbonylation.csv",
    "CysSulfinylation.csv",
    "CysSulfonylation.csv");

optm_data_second_batch0 = tibble();
for (each_file in datasets_second_batch){
  
  each_data = read_csv(paste("optm_15_second_batch/", each_file, sep = ""));
  
  optm_data_second_batch0 = 
    optm_data_second_batch0 %>% 
    bind_rows(each_data);
}

optm_data_second_batch_dw = 
  optm_data_second_batch0 %>% 
  separate(protein, into = c("sp","uniprot", "protein"), sep = "\\|") %>%
  mutate(sample_rep = (sample_rep %>% str_replace("\\.txt", ""))) %>% 
  mutate(sample_rep = (sample_rep %>% str_replace("-R", ""))) %>% 
  mutate(sample_rep = (sample_rep %>% str_replace("cPTM_", ""))) %>% 
  separate(sample_rep, into = c("fn", "strain", "day", "replicate"), sep = "[_-]") %>% 
  select(-sp, -protein,-fn)

optm_data_second_batch = 
  optm_data_second_batch_dw %>% 
  mutate(strain = ifelse(strain == "BABLc", "BALB", strain)) %>% 
  mutate(day = paste(day %>% str_sub(3), day %>% str_sub(1,2), sep = "")) %>% 
  mutate(occupancy = spec_count_modified/(spec_count_modified + spec_count_unmodified)) %>% 
  mutate(type = ifelse(type == "MetSulfoxidation", "MetO", type)) %>% 
  mutate(type = ifelse(type == "MetSulfonation", "MetO2", type)) %>% 
  
  mutate(type = ifelse(type == "CysCarbonylation", "Cys4HNE", type)) %>% 
  mutate(type = ifelse(type == "CysSulfinylation", "CysSO2H", type)) %>% 
  mutate(type = ifelse(type == "CysSulfonylation", "CysSO3H", type)) %>% 
  
  mutate(type = ifelse(type == "LysCarbonylation", "Lys2AAA" , type)) %>% 
  mutate(type = ifelse(type == "ProCarbonylation", "ProCH" , type)) %>% 
  
  mutate(type = ifelse(type == "ArgHydroxylation", "Arg-OH" , type)) %>% 
  mutate(type = ifelse(type == "AsnHydroxylation", "Asn-OH" , type)) %>% 
  mutate(type = ifelse(type == "AspHydroxylation", "Asp-OH" , type)) %>% 
  mutate(type = ifelse(type == "LysHydroxylation", "Lys-OH" , type)) %>% 
  mutate(type = ifelse(type == "PheHydroxylation", "Phe-OH" , type)) %>% 
  mutate(type = ifelse(type == "TryHydroxylation", "Trp-OH" , type)) %>% 
  mutate(type = ifelse(type == "TyrHydroxylation", "Tyr-OH" , type)) %>% 
  mutate(type = ifelse(type == "ValHydroxylation", "Val-OH" , type)) %>% 
  add_column(group = "ISO") %>% 
  select(uniprot, ptm_location = diff_mod_global_position, strain, 
         group, day, occupancy, optm = type, replicate);
  
    
######################################################################################################################################################
######################################################################################################################################################
optm_dw_first_batch_ini = 
  optm_data_first_batch %>% 
  filter(optm != "Pro-OH") %>% 
  filter(day == "d00" | day == "d01" | day == "d03" | day == "d05" | day == "d07" | day == "d10" | day == "d14") %>% 
  add_column(replicate = "S4");

optm_dw = 
  optm_dw_first_batch_ini %>% 
  bind_rows(optm_dw_first_batch_ini %>%
              filter(day == "d00") %>%
              mutate(group = "ISO")) %>% 
  bind_rows(optm_data_second_batch) %>% 
  group_by(uniprot, ptm_location, strain, group, day, optm, replicate) %>% 
  summarise(occupancy = mean(occupancy)) %>% 
  ungroup();

optm_reform = 
  optm_dw %>%
  spread(group, occupancy) %>% 
  gather(group, occupancy, CTRL:ISO) %>% 
  spread(strain, occupancy) %>% 
  gather(strain, occupancy, AJ:FVB);

optm_ftr = 
  optm_reform %>% 
  group_by(uniprot, ptm_location, group, day, optm, replicate) %>% 
  summarize(non_na_count = sum(!is.na(occupancy))) %>% 
  ungroup() %>% 
  spread(group,non_na_count);

the_min_val = min(optm_dw$occupancy)/2; # 0.001587301


optm_ftred = 
  optm_ftr %>%
  #filter( (CTRL >= 3 & ISO >= 3) | (CTRL == 0 & ISO >= 3) | (CTRL >= 3 & ISO == 0)) %>%
  # filter( (CTRL >= 6 & ISO >= 6)) %>%  # | (CTRL == 0 & ISO >= 3) | (CTRL >= 3 & ISO == 0)) %>%
  gather(group, non_na_count, CTRL:ISO) %>% 
  left_join(optm_reform,
            by = c("uniprot", "ptm_location", "day",   "optm", "group")) %>% 
  left_join(optm_ftr %>% 
              filter((CTRL == 0 & ISO >= 3) | (CTRL >= 3 & ISO == 0)) %>% 
              select(uniprot:optm) %>% 
              add_column(abolished_initiated = TRUE),
            by = c("uniprot", "ptm_location", "day",   "optm")) %>% 
  mutate(occupancy_val = ifelse(((non_na_count == 0) & (abolished_initiated == TRUE)), the_min_val, occupancy)) %>% 
  select(-non_na_count, -occupancy, -abolished_initiated);

optm_ftred %>% write_csv("total_mouse_15optm_abolished_initiated_no_filter.csv")
