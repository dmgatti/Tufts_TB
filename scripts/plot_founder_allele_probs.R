library(qtl2convert)
library(qtl2)
library(tidyverse)

# probs: qtl2-style genoprobs object containing a list of 3 dimensional
#        founder allele probs with dimensions: samples x 8 founders x markers.
#        Each list element represents one chromosome.
#        
# map:   qtl2-style map object containing a list of numeric vectors, each 
#        of which contains marker locations for one chromosome.
# chr:   name of chromosome to plot. Must match the names on probs and map.
plot_fp = function(probs, map, chr) {
  
  if(!chr %in% names(probs)) {
    stop(paste('Chromosome', chr, 'not found in probs'))
  }

  founder_names = c('A/J', 'C57BL/6J', '129S1/SvImJ', 'NOD/ShiLtJ', 'NZO/HlLtJ', 'CAST/EiJ', 'PWK/PhJ', 'WSB/EiJ')
  n_founders = length(founder_names)
  
  # Keep the current chromosome.  
  pr = probs[[chr]]
  mp = map[[chr]]
  
  # Get founder allele probs
  fp = apply(pr, 2:3, mean, na.rm = TRUE)
  
  # Subset to keep common markers.
  common_markers = intersect(names(mp), colnames(fp))
  mp = mp[common_markers]
  fp = fp[,common_markers]
  rownames(fp) = founder_names
  
  # Make long data for ggplot2.
  fp_long = tibble(founder = rep(rownames(fp), each = ncol(fp)),
                   marker  = rep(names(mp), n_founders),
                   pos     = rep(mp, n_founders),
                   ap = as.vector(t(fp)))
  
  # Rename CC colors.
  cc_colors = qtl2::CCcolors
  names(cc_colors) = founder_names
  
  fp_long %>%
    ggplot(aes(pos, ap)) +
      geom_line(aes(color = founder), size = 1) +
      scale_color_manual(values = cc_colors) +
      facet_wrap(~founder, ncol = 1)
  
} # plot_fp()

# Example code.

# Set filenames.
probs_file = '/media/dmgatti/hdb/projects/TB/haplo_reconstr/fda_do_alleleprobs_qc_filtered.rds'
map_file   = '/media/dmgatti/hdb/projects/TB/haplo_reconstr/qtl2/pmap.csv'

# Load in founder allele probs and marker map.
probs = readRDS(probs_file)
map   = readr::read_csv(map_file, col_types = 'ccn') %>%
          as.data.frame()
map   = qtl2convert::map_df_to_list(map = map, pos_column = 'pos')

plot_fp(probs, map, '1')
