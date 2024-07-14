

coancenstry = function(geno_dta, type = 'G05N'){
  
  #'@param geno_dta - genotype dataframe of n individuals (rows) by m markers (columns)
  #'@param type - how realized relationship is calculated. 'G05N' (Forni et al., 2011) is default. Other option is 'GOFN' (VanRaden, 2008).
  #'@return dataframe with pairwise coa
  
  # make sure marker data is -1, 0, 1
  # ge_id in rownames
  # marker names in columns
  
  # for G05, the minor allele freq is set as 0.5 across all markers
  # this is because calculation is prone to error based on the num of markers and inds
  # GOF calc minor allele freq as observed freq in the population
    
  
  dta_mat = as.matrix(geno_dta)
  
  # this is for observed allele freq in the pop
  
  if (type == 'GOFN'){
    
    minor_af = round((apply(dta_mat,2,sum)+nrow(dta_mat))/(nrow(dta_mat)*2),3)
    
    af = 2*(minor_af-.5)
    
    afMat = matrix(af,byrow=T,nrow=nrow(dta_mat),ncol=ncol(dta_mat))
    
  }
  
  
  # with allele freq set at 0.5
  
  else if (type == 'G05N'){
    
    minor_af = array(0.5,ncol(dta_mat))
    
    af = minor_af
    
    afMat = matrix(af,byrow=T,nrow=nrow(dta_mat),ncol=ncol(dta_mat))
    
  }
  
  
  
  Z = as.matrix(dta_mat - afMat)
  
  ZZt = Z %*% t(Z)

  # also, the ZZt matrix is divided by tr(ZZ)/n instead of the normal 2*sum((q*maf)(1-maf))
  # using tr(ZZt)/n limits diag elements to 1. Therefore, G05N and GOFN
  
  trZZ = sum(diag(ZZt))/nrow(ZZt)
  
  ModifiedG = ZZt/trZZ
  
  out_coa = data.frame(Parent_1 = rownames(ModifiedG)[row(ModifiedG)[upper.tri(ModifiedG)]],
                       Parent_2 = colnames(ModifiedG)[col(ModifiedG)[upper.tri(ModifiedG)]],
                       coa = ModifiedG[upper.tri(ModifiedG)])
  
  # for those coa > 1, change to max out at 1
  
  out_coa = out_coa %>% mutate(coa = ifelse(coa > 1, 1.0, coa))
  
  return(out_coa)
  
  
}


tmp_minus1 = tmp - 1

outG05 = coancenstry(tmp_minus1)
outGOF = coancenstry(tmp_minus1, type = 'GOFN')
