setwd("~/Desktop/Projects/peipei/OPTM/DATA/2019 Ding new mouse data/")

#output_csv_name = "optm_occupancy_mouse_2019.csv"
database_file_name = "../fasta/_ip2_ip2_data_nhlbi_database__SwissProt_mouse_06-12-2016_reversed.fasta"
db_dta = prepare_db_fasta(database_file_name)

folders = list.files()
file_start = c('Arginine Hydroxylation','Asparagine Hydroxylation',
               'Aspartate Hydroxylation','Cysteine Carbonylation','Cysteine Sulfinylation',
               'Cystein Sulfonylation','Lysine Carbonylation',
               'Lysine Hydroxylation','Methionine Sulfonation',
               'cPTM_Methionine Sulfoxidation','Phenylalanine Hydroxylation',
               'Proline Carbonylation','Tryptophan Hydroxylation',
               'Tyrosine Hydroxylation','Valine Hydroxylation'
               )
optm_type = c('Arg-OH','Asn-OH','Asp-OH','Cys4HNE','CysSO2H','CysSO3H','Lys2AAA','Lys-OH',
              'MetO2','MetO','Phe-OH','ProCH','Trp-OH','Tyr-OH','Val-OH')

optm_all = tibble()
for(i in 1:15){
  for(day in c('01h','01d','03d','07d','14d')){
    filenames = paste0(folders[i],'/',file_start[i], '_C57-',day, '-S',c('1','2','3'),'.txt')
    optm_dw = optm_occupancy_cal(database_file_name,  filenames=filenames )
    optm_dw$type = optm_type[i]
    optm_dw$occupancy = optm_dw$spec_count_modified/(optm_dw$spec_count_modified + optm_dw$spec_count_unmodified)
    optm_dw$day = day
    optm_dw$species = 'C57'
    optm_all = optm_all %>% 
      bind_rows(optm_dw)
  }
}

write_csv(optm_all,'../../New_calculation/Data/C57.txt'  )
















