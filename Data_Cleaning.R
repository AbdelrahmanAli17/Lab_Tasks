library(XML)

#load XML file 

csf_data <- xmlToDataFrame(xmlParse("csf_metabolites.xml"))


smdb_ids_all <- c()

for(i in 1:nrow(csf_data)){
  
  string <- csf_data$biological_properties[i]
  
  smdb_ids <- gregexpr("SMP\\d+", string, perl = TRUE)[[1]]
  
  smdb_ids_all <- c(smdb_ids_all, substring(string, smdb_ids, smdb_ids + 7))
  

}

smdp_id <- as.data.frame(smdb_ids_all)

n_metabolites <- length(smdp_id[smdp_id$smdb_ids_all == "SMP00044" ,])

n_metabolites
