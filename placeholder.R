xy <- as.data.frame(combn(rownames(tmp), 2))

yx = sapply(xy, function(x) as.character(x), simplify = FALSE)



coa_list = list()
for(pair in 1:length(yx)){
  
  
  tmp_geno = tmp %>% filter(rownames(tmp) %in% yx[[pair]])
  
  coa_list[[pair]] = get_coa(tmp_geno)
}

a = do.call('rbind', coa_list)
