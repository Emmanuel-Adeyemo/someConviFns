xy <- as.data.frame(combn(rownames(tmp), 2))

yx = sapply(xy, function(x) as.character(x), simplify = FALSE)


# pairwise matching

matching_list = list()
            
for(pair in 1:length(yx)){
  
  
  tmp_geno = tmp %>% filter(rownames(tmp) %in% yx[[pair]])
  
  matching_list[[pair]] = get_matching(tmp_geno)
}

a = do.call('rbind', matching_list)
