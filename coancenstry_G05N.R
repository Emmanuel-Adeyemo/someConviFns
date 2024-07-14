

coancenstry_G05N = function(geno_dta){
  
  #'@return dataframe with pairwise coa
  
  # make sure marker data is -1, 0, 1
  # ge_id in rownames
  # marker names in columns
  
  # for G05, the minor allele freq is set as 0.5 across all markers
  # instead of calc it as observed freq in the population
  # this is because it is prone to error based on the num of markers and inds
  
  # also, the ZZt matrix is divided by tr(ZZ)/n instead of the normal 2*sum((q*maf)(1-maf))
  # using tr(ZZt)/n limits diag elements to 1. Therefore, G05N
  
  dta_mat = as.matrix(geno_dta)
  
  # this is for observed allele freq in the pop
  
  #minor_af = round((apply(dta_mat,2,sum)+nrow(dta_mat))/(nrow(dta_mat)*2),3)
  
  # with allele freq set at 0.5
  
  minor_af = array(0.5,ncol(dta_mat))
  af = minor_af
  afMat = matrix(af,byrow=T,nrow=nrow(dta_mat),ncol=ncol(dta_mat))
  
  Z = as.matrix(dta_mat - afMat)
  
  ZZt = Z %*% t(Z)
  
  trZZ = sum(diag(ZZt))/nrow(ZZt)
  
  G05N = ZZt/trZZ
  
  out_coa = data.frame(Parent_1 = rownames(G05N)[row(G05N)[upper.tri(G05N)]],
                       Parent_2 = colnames(G05N)[col(G05N)[upper.tri(G05N)]],
                       coa = G05N[upper.tri(G05N)])
  
  # for those coa > 1, change to max out at 1
  
  out_coa = out_coa %>% mutate(coa = ifelse(coa > 1, 1.0, coa))
  
  return(out_coa)
  
  
}

# tmp data is coded 0, 1, 2

tmp_minus1 = tmp - 1

out = coancenstry_G05N(tmp_minus1)
