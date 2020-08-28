
rm(list = ls(all = T));

library(tidyverse);

setwd("D:/15-OPTM_Mouse_2019/DownloadDTA/");

exp_files_csv <-list.files (pattern="*.csv")

print (exp_files_csv)

for(each_exp_name in exp_files_csv){
  
  print (each_exp_name)
  
  search_list <- read.csv(file=paste(each_exp_name),header=FALSE)
  
  t_row <- nrow(search_list)
  
  nr <- 1
  
  for (nr in  1:t_row) {
    
    if (search_list[nr,2] !="" & search_list[nr,3] != "Not Finished.") {
      
    
      if (str_detect(each_exp_name,"Run_")){
        
      samplename <- paste(sub("Run_","",each_exp_name))
      samplename <- paste(sub(".csv","",samplename),paste("-R"),paste(".csv"),sep="")
      
       } else {
         
         samplename <- each_exp_name};
        
      print(samplename);
      print(search_list[nr,2]);
      print ((paste (search_list[nr,2], paste ("_"), paste (sub(".csv","",samplename)), paste(".txt"),sep="")));
      
        download.file (url=paste("http://goldfish.scripps.edu/ip2/ip2_data/nhlbi/",paste(sub(".*//nhlbi/","",search_list[nr,15])),paste("DTASelect-filter.txt"),sep=""), 
           destfile = paste (search_list[nr,2], paste ("_"), paste (sub(".csv","",samplename)), paste(".txt"),sep=""));
      
   
     

# download.file (url=paste("http://goldfish.scripps.edu/ip2/ip2_data/nhlbi/",paste(sub(".*//nhlbi/","",search_list[nr,15])),paste("DTASelect.txt"),sep=""), 
     #                destfile = paste ("UFD_", paste (sub(".csv","",samplename)), paste ("_"), paste (sub("No_","",search_list[nr,2])), paste(".txt"),sep=""));
                     
                   
        }
      
                   
  }

  
}
