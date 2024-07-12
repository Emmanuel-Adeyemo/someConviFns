
get_coa = function(dta){
  
  k = list()
  
  out = list()
  
  for(i in 1:ncol(dta)){
    
    if(dta[[i]][1] == dta[[i]][2]){k[[i]] = 1}
    else if(dta[[i]][1] != dta[[i]][2]){k[[i]] = 0}
    
  }
  
  out['Parent_1'] = rownames(dta)[1]
  out['Parent_2'] = rownames(dta)[2]
  
  out['coa'] = round(Reduce('+', k) / (ncol(dta)),3)
  
  out_call = do.call('cbind', out) %>% as.data.frame()
  
  return(out_call)
}
