expand_geno = function(dta){
  
  for(i in 1:ncol(dta)){
    for(j in 1:nrow(dta)){
      
      if(dta[[i]][j] == 0){dta[[i]][j] = 11}
      else if(dta[[i]][j] == 2){dta[[i]][j] = 22}
      else if(dta[[i]][j] == 1){r = sample(c(12,21),1); dta[[i]][j] = r}
      
      
    }
  }
  return(dta)
}

#setwd("C:/Users/eadey/OneDrive/Desktop/UMN/Thesis PhD/F5_GS/2018/f5_gs/randFor")

geno <- read.csv("2018_Wheat_F5_GStrain_SNPs_20pct_missing_RF_output.txt", sep = "\t", header = TRUE)
tmp = geno[c(500:504, 541:544),]
row.names(tmp)=NULL
tmp=tmp %>% tibble::column_to_rownames('ID')
tmp4 = expand_geno(tmp3)
