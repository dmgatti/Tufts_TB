# QTL heatmap function.
# lod: data.frame containing marker, chr, pos (in Mb) and the 
#      LOD scores for each phenotype.
qtl_heatmap = function(lod) {
  
  # Split markers and LOD.
  mkr = lod[,1:3]
  lod = as.matrix(lod[,-(1:3)])
  
  # Get unique choromosomes.
  unique_chr = distinct(qtl, chr) %>% pull(chr)
  
  # Apportion 500 markers proportionally to each chromosome.
  nm = 500
  max_pos  = sapply(map, max)
  prop_chr = max_pos / sum(max_pos)
  num_mkr  = round(nm * prop_chr)
  
  # Interpolate a new marker map.
  new_map = mapply(function(m, num) { approx(x = m, n = num)$y }, m = map, num = num_mkr)
  
  # Interpolated QTL map and LODs.
  new_qtl = NULL
  
  for(chr in names(new_map)) {
    
    # Get breakpoints and midpoints.
    brks = cut(map[[chr]], new_map[[chr]])
    mids = sapply(split(map[[chr]], brks), mean)
    
    # Get LOD for this chromosome.
    curr_lod = lod[mkr$chr == chr,]
    
    # Make new data frame for this chromosome.
    tmp = data.frame(chr = rep(chr, length(mids)),
                     pos = mids,
                     matrix(0, length(mids), ncol(curr_lod), dimnames = list(NULL, colnames(curr_lod))))
    
    for(j in colnames(curr_lod)) {
      tmp[[j]] = sapply(split(curr_lod[,j], brks), max)
    } # for(j)
    
    new_qtl = rbind(new_qtl, tmp)
    
  } # for(chr)
  
  # Cluster phenotypes.
  tmp = as.matrix(new_qtl[,-(1:2)])
  cl = hclust(as.dist(1.0 - cor(tmp, use = "pairwise")), method = 'average')
  rm(tmp)
  new_qtl = new_qtl[,c(1:2, cl$order + 2)]
    
  colors = c('white', 'grey95', 'grey90', 'grey85', 'grey50', 'grey20', 'black')
  
  new_qtl %>%
    mutate(chr = factor(chr, levels = c(1:19, 'X'))) %>%
    gather(pheno, lod, -(chr:pos)) %>%
    mutate(pheno = factor(pheno, levels = cl$labels[cl$order])) %>%
    ggplot() +
    geom_tile(aes(x = pos, y = 1, color = lod, fill = lod), width = 5) +
    scale_color_gradientn(colors = colors, values = scales::rescale(c(0, 5.5, 6, 6.5, 7.0, max(lod) + 0.1), to = c(0,1))) +
    scale_fill_gradientn(colors  = colors, values = scales::rescale(c(0, 5.5, 6, 6.5, 7.0, max(lod) + 0.1), to = c(0,1))) +
    facet_grid(pheno ~ chr, scales = 'free_x') +
    theme(panel.spacing = unit(0, 'lines'),
          axis.text.x = element_text(angle = 90, vjust = 0),
          axis.text.y = element_blank(),
          axis.ticks  = element_blank(),
          axis.title  = element_blank(),
          strip.text.x = element_text(size = 24),
          strip.text.y = element_text(size = 14, angle = 0))
  
} # qtl_heatmap()
