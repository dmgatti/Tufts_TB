library(tidyverse)
library(qtl2)

base_dir   = '/media/dmgatti/hdb/projects/TB'
input_file = file.path(base_dir, 'data', 'tufts_do_tb_qtl2_input.Rdata') 
sim_out_dir = file.path(base_dir, 'simulations', 'qtl_sims_1')
fig_dir     = file.path(sim_out_dir, 'figures')

source(file.path(base_dir, 'scripts', 'phenotype_simulation_pipeline.R'))

# Read in the DOQTL formatted data.
load(input_file)

# CC SNP and gene database files (from Karl).
ccsnpdb = "/media/dmgatti/hda/data/MUGA/cc_variants.sqlite"
mgidb   = "/media/dmgatti/hda/data/MUGA/mouse_genes_mgi.sqlite"
snp_func  = create_variant_query_func(dbfile = ccsnpdb)
gene_func = create_gene_query_func(dbfile = mgidb)

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]

# Run a loop with 1000 simulations.
n_sim  = 1000 # Number of simulations
offset = 5    # Distance +/- for QTL location. (2 * offset is distance between simulated QTL)

sim_pheno = matrix(0, nrow = nrow(pheno), ncol = n_sim,
                   dimnames = list(rownames(pheno_rz), paste0('sim', 1:n_sim)))
sim_loci  = NULL
sim_peaks = NULL
sims      = NULL

for(s in 1:n_sim) {
  
  print(paste('Simulation', s))
  t = proc.time()[3]
  
  chr     = NULL
  chrlen  = NULL
  chr_rng = NULL
  pos     = NULL
  
  # Keep trying different loci until we get a result. Sometimes, we can't find the
  # requested allele split at a given locus.
  sims = NULL
  while(is.null(sims)) {
    chr     = sample(names(probs), 1)
    chrlen  = max(map[[chr]])
    chr_rng = which(map[[chr]] > offset + min(map[[chr]]) & map[[chr]] < chrlen - offset - 1)   # This protects us from going off the end of the chromosome.
    pos     = sample(map[[chr]][chr_rng], 1)
    # Simulate 1 peak.
    effects = data.frame(chr = rep(chr, 1), pos = c(max(map[[chr]][1], pos)), 
                         split = sample(1:4, 1))  # This also protects  us from going off the beginning of chr.
    # Simulate 2 peaks.
#    effects = data.frame(chr = rep(chr, 2), pos = c(max(map[[chr]][1], pos - offset), pos + offset), 
#                         split = sample(1:4, 2))  # This also protects  us from going off the beginning of chr.
    sims    = simulate_phenotype(probs, Kcor, map, effects, hsq = 0.3, query_func = snp_func)
  } # while(!is.null(sims))
  
  sim_pheno[,s] = sims$pheno[,1]
  sim_geno = data.frame(sim = rep(s, nrow(sims$sim_geno)), hsq =  rep(attr(sims$pheno, 'hsq'), nrow(sims$sim_geno)), sims$sim_geno)
  if(is.null(sim_loci)) {
    sim_loci = sim_geno
  } else {
    sim_loci = rbind(sim_loci, sim_geno)
  } # else
  
  lod = scan1(genoprobs = probs[,chr], 
              pheno     = sims$pheno,
              kinship   = K[[chr]],
              addcovar  = addcovar,
              cores     = 4)
  
  peaks = find_peaks(lod, map, threshold = 3, prob = 0.95)
  
  png(file.path(sim_out_dir, 'lod_plots', paste0('sim', s, '_lod.png')),
      width = 1000, height = 800, res = 128)
  plot_scan1(lod, map, main = paste('Sim', s, paste(sims$sim_geno$pos, collapse = ' & ')),
             col = 'black')
  gr_params = par()
  rect(peaks$ci_lo, gr_params$usr[3], peaks$ci_hi, gr_params$usr[4], col = rgb(0, 0, 1, 0.1))
  abline(v = sims$sim_geno$pos, col = 'red', lty = 3, lwd = 2)
  dev.off()
  
  # Save LOD data.
  saveRDS(lod, file = file.path(sim_out_dir, 'lod_plots', paste0('sim', s, '_lod.rds')))
  
  # Save current peaks to overall peaks.
  peaks = data.frame(sim = s, peaks)
  if(is.null(sim_peaks)) {
    sim_peaks = peaks
  } else {
    sim_peaks = rbind(sim_peaks, peaks)
  } # else
  
  pr_tmp  = probs[,chr]
  mkr_rng = 0
  if(nrow(sims$sim_geno) == 1) {
    mkr_rng = which(map[[chr]] >= sims$sim_geno$pos[1] - 5 & map[[chr]] <= sims$sim_geno$pos[1] + 5)
  } else {
    mkr_rng = which(map[[chr]] >= sims$sim_geno$pos[1] - 5 & map[[chr]] <= sims$sim_geno$pos[2] + 5)
  } # else
  pr_tmp[[chr]] = pr_tmp[[chr]][,,names(mkr_rng)]
  
  blup = scan1blup(genoprobs = pr_tmp, 
                   pheno     = sims$pheno,
                   kinship   = K[[chr]],
                   addcovar  = addcovar,
                   cores     = 8)
  
  png(file.path(sim_out_dir, 'blup_plots', paste0('sim', s, '_blup.png')),
      width = 1000, height = 800, res = 128)
  plot_coefCC(blup, map, main = paste('Sim', s, paste(sims$sim_geno$pos, collapse = ' & ')),
              scan1_output = lod)
  gr_params = par()
  rect(peaks$ci_lo, -1000, peaks$ci_hi, gr_params$usr[4], col = rgb(0, 0, 1, 0.1))
  # Had to use lines() instead of abline() because of a weird plotting artifact that only
  # extended the line through part of the top panel.
  for(j in 1:nrow(sims$sim_geno)) {
    lines(rep(sims$sim_geno$pos[j], 2), c(-1000, 1000), col = 'red', lty = 3, lwd = 2)
  } # for(j)
  dev.off()
  
  # Save BLUP data.
  saveRDS(blup, file = file.path(sim_out_dir, 'blup_plots', paste0('sim', s, '_blup.rds')))
  
  assoc = scan1snps(genoprobs = probs[,chr],
                    map       = map,
                    pheno     = sims$pheno,
                    kinship   = K[[chr]],
                    addcovar  = addcovar,
                    chr       = chr,
                    start     = sims$sim_geno$pos[1] - 1,
                    end       = ifelse(nrow(sims$sim_geno) == 1, sims$sim_geno$pos[1] + 1, sims$sim_geno$pos[2] + 1),
                    query_func = snp_func,
                    keep_all_snps = TRUE,
                    cores     = 4)
  
  png(file.path(sim_out_dir, 'assoc_plots', paste0('sim', s, '_assoc.png')),
      width = 1000, height = 800, res = 128)
  plot_snpasso(assoc$lod, assoc$snpinfo, show_all_snps = TRUE, col = 'black', 
               main = paste('Sim', s, paste(sims$sim_geno$pos, collapse = ' & ')),
               minlod = 1)
  abline(v = sims$sim_geno$pos, col = 'red', lty = 3, lwd = 2)
  gr_params = par()
  rect(peaks$ci_lo, gr_params$usr[3], peaks$ci_hi, gr_params$usr[4], col = rgb(0, 0, 1, 0.1))
  dev.off()
  
  # Save association mapping data.
  saveRDS(assoc, file = file.path(sim_out_dir, 'assoc_plots', paste0('sim', s, '_assoc.rds')))
  
  # Get the top 0.5% of SNPs.
  max_lod = maxlod(assoc$lod)
  thr     = quantile(assoc$lod, probs = 0.995)
  # Careful here: the gene column contains comma's, so use TAB delimiter.
  write.table(top_snps(assoc$lod, assoc$snpinfo, drop = max_lod - thr), 
              file = file.path(sim_out_dir, 'assoc_plots', paste0('sim', s, '_top_snps.txt')),
              quote = FALSE, sep = '\t')
  
  write.csv(sim_pheno, file.path(sim_out_dir, 'simulated_phenotypes.csv'), quote = F)
  write.csv(sim_loci,  file.path(sim_out_dir, 'simulated_loci.csv'), quote = F)
  write.csv(sim_peaks, file.path(sim_out_dir, 'found_peaks.csv'), quote = F)
  
  print(paste('   Time:', proc.time()[3] - t))
  
} # for(s)

##########################################################
# Read in the simulated and found peaks and compare them.
sim_peaks   = read.csv(file.path(sim_out_dir, 'simulated_loci.csv'))
found_peaks = read.csv(file.path(sim_out_dir, 'found_peaks.csv'))

comp = compare_peaks(sim_peaks, found_peaks)

png(file.path(fig_dir, 'hsq_hist.png'), res = 128)
print(ggplot(comp) +
        geom_histogram(aes(hsq), color = 'black') +
        labs(title = 'Distribution of h^2'))
dev.off()

png(file.path(fig_dir, 'lod_hist.png'), res = 128)
print(ggplot(comp) +
        geom_histogram(aes(lod), color = 'black') +
        labs(title = 'Distribution of LODs', x = 'LOD') )
dev.off()

png(file.path(fig_dir, 'num_found_vs_hsq.png'), width = 1000, height = 800, res = 164)
print(comp %>% 
        filter(!is.na(num_found)) %>% 
        ggplot(aes(num_found, hsq)) +
        geom_beeswarm(alpha = 0.3) +
        geom_smooth(method = 'lm')) +
        labs(title = 'Number of peaks found vs. h^2', x = 'Number of Peaks Found',
             y = 'h^2')
dev.off()


png(file.path(fig_dir, 'num_found_vs_lod.png'), width = 1000, height = 800, res = 164)
print(comp %>% 
        filter(!is.na(num_found)) %>% 
        ggplot(aes(num_found, lod)) +
        geom_beeswarm(alpha = 0.3) +
        geom_smooth(method = 'lm')) +
       labs(title = 'Number of peaks found vs. LOD', x = 'Number of Peaks Found',
            y = 'LOD')
dev.off()

# How frequently does the CI cover at least one of the peaks.
comp %>% 
  summarise(mean(num_found > 0))

# About 90% of the time. Lower than 95%.

# Add peak CI width column.
comp = comp %>% 
         mutate(width = ci_hi - ci_lo) %>% 
         arrange(num_found)

# LOD vs. peak width.
png(file.path(fig_dir, 'lod_vs_peak_width.png'), width = 1000, height = 800, res = 164)
comp %>% 
  mutate(num_found = as.character(num_found)) %>% 
  ggplot() +
    geom_point(aes(lod, width + 0.01, color = num_found), alpha = 0.5) +
    scale_color_brewer(palette = 'Dark2') +
    scale_y_log10() +
    labs(title = 'LOD vs. C.I. width', x = 'LOD', y = '95% Baysian C.I. Width (Mb)')
dev.off()


######################################################################
# Read in simulated peak locations and compare with top SNPs from 
# association mapping.
assoc_dir = file.path(sim_out_dir, 'assoc_plots')
top_snps_files = list.files(path = assoc_dir, pattern = '_top_snps.txt$', full.names = TRUE)

# Get simulation number for each file.
sim_num = strsplit(top_snps_files, split = '/')
sim_num = sapply(sim_num, function(z) { z[length(z)] })
sim_num = parse_number(sim_num)

sim_loci = read.csv(file.path(sim_out_dir, 'simulated_loci.csv'))[,-1]

found_snps = NULL

for(i in 1:length(top_snps_files)) {
  
  # Read in top_snps file and fix column names.
  top_snps = read.delim(top_snps_files[i], row.names = NULL, sep = '\t')
  
  found_snps = bind_rows(found_snps, compare_top_snps(sim_loci, sim_num[i], top_snps))
  
} # for(i)

# How frequently do we find one of the peaks?
found_snps %>% 
  summarize(pct_found = mean(!is.na(lod)))

# ~45.5 % of the time.

# How frequently do we find one of the peaks with the highest ranked LOD?
png(file.path(fig_dir, 'lod_rank_top_peak_hist.png'), width = 1000, height = 800, res = 164)
print(found_snps %>% 
        mutate(lod_rank = if_else(is.na(lod_rank), 0L, lod_rank)) %>% 
        count(lod_rank) %>% 
        mutate(n = n / 2) %>% 
  ggplot() +
    geom_col(aes(lod_rank, n)) +
    labs(title = 'Rank of simulated SNP among top LOD scores',
         x = 'Rank of simulated SNP', y = 'Count'))
dev.off()


found_snps %>% 
  dplyr::select(lod_rank) %>%
  mutate(lod_rank = if_else(is.na(lod_rank), 0L, lod_rank)) %>% 
  mutate(prop = cume_dist(lod_rank)) %>% 
  ggplot(aes(lod_rank, prop)) +
    geom_line() +
    scale_x_continuous(breaks = 0:8 * 5) +
    labs(title = 'Proportion of simulations causal SNP is in top n LOD',
         x     = 'Rank of LOD score',
         y     = 'Proportion of simulations')

found_snps %>% 
  count(lod_rank) %>% 
  mutate(lr = lod_rank == 1) %>% 
  group_by(lr) %>% 
  summarise(num_found = sum(n)) %>% 
  pivot_wider(names_from = lr, values_from = num_found) %>% 
  summarize(pct_found = `TRUE` / (`TRUE` + `FALSE` + `NA`))

# Only 7% of the time...


############################################################################
# Go through the previous simulations and regress out the highest SNP.
# Then see if we identify the simulated SNP at the other locus.
assoc_in_dir   = file.path(sim_out_dir, 'assoc_plots')
assoc_regr_dir = file.path(sim_out_dir, 'assoc_regr_plots')

sim_loci      = read.csv(file.path(sim_out_dir, 'simulated_loci.csv'))[,-1]
sim_pheno     = read.csv(file.path(sim_out_dir, 'simulated_phenotypes.csv'))
rownames(sim_pheno) = sim_pheno[,1]
sim_pheno    = sim_pheno[,-1]
top_snp_files = dir(assoc_in_dir, pattern = '_top_snps.txt')
sim_num = parse_number(top_snp_files)

for(i in seq_along(top_snp_files)) {

  print(i)
  top_snps = read.delim(file.path(assoc_in_dir, top_snp_files[i])) %>% 
               slice_max(lod, n = 1) %>% 
               slice(1)
  
  chr = top_snps$chr[1]
  pos = top_snps$pos[1]
  sim = filter(sim_loci, sim == sim_num[i])
  sim_id = paste0('sim', sim_num[i])
  
  rescan1(probs = probs, pheno = sim_pheno[,sim_id,drop = FALSE], K = K, addcovar = addcovar, 
          map = map, chr = chr, pos = pos, sim_loci = sim, query_func = snp_func)
    
} # for(f)

# Using 'found_snps' from the block above this one, see if we found any more SNPs 
# once we regressed out the top SNP.
assoc_dir = file.path(sim_out_dir, 'assoc_regr_plots')
top_snps_files = list.files(path = assoc_dir, pattern = '_top_snps.txt$', full.names = TRUE)

# Get simulation number for each file.
sim_num = strsplit(top_snps_files, split = '/')
sim_num = sapply(sim_num, function(z) { z[length(z)] })
sim_num = parse_number(sim_num)

found_snps_2 = NULL

for(i in 1:length(top_snps_files)) {
  
  # Read in top_snps file and fix column names.
  top_snps = read.delim(top_snps_files[i], row.names = NULL, sep = '\t')
   
  # The first SNPs found will be in lod.x and lod_rank.x. The second ones found
  # after regressing out the max SNP will be in lod.y and lod_rank.y.
  found_snps_2 = bind_rows(found_snps_2, compare_top_snps(found_snps, sim_num[i], top_snps))
  
} # for(i)

found_snps_2 = found_snps_2 %>% 
                 mutate(found_order = if_else(is.na(lod.x), NA_integer_, 1L),
                        found_order = if_else(is.na(lod.y), found_order, 2L),
                        lod.x       = if_else(is.na(lod.x) & !is.na(lod.y), lod.y, lod.x),
                        lod_rank.x  = if_else(is.na(lod_rank.x) & !is.na(lod_rank.y), lod_rank.y, lod_rank.x)) %>%
                 dplyr::select(-lod.y, -lod_rank.y) %>% 
                 dplyr::rename(lod      = lod.x,
                               lod_rank = lod_rank.x)

# How frequently do we find one of the peaks after regression?
found_snps_2 %>% 
  summarize(pct_found = mean(!is.na(lod)))

# Better: 72.8 % of the time.

# How frequently do we find one of the peaks with the highest ranked LOD?
png(file.path(fig_dir, 'lod_rank_top_peak_regr_hist.png'), width = 1000, height = 800, res = 164)
print(ggplot(found_snps_2) +
        geom_histogram(aes(lod_rank), binwidth = 1, fill = 'grey70', color = 'grey30') +
        labs(title = 'Histogram of LOD Rank containing simulated SNP (after regressing out top SNP)',
             x = 'rank(LOD)', y = 'Count'))
dev.off()

found_snps_2 %>% 
  count(lod_rank) %>% 
  mutate(lr = lod_rank == 1) %>% 
  group_by(lr) %>% 
  summarise(num_found = sum(n)) %>% 
  pivot_wider(names_from = lr, values_from = num_found) %>% 
  summarize(pct_found = `TRUE` / (`TRUE` + `FALSE` + `NA`))

# Simulated SNPs have the highest LOD score only 9.65% of the time.

write.csv(found_snps_2, file = file.path(sim_out_dir, 'assoc_found_peak_summary.csv'))

# How often do we find 0, 1 or 2 peaks?
found_snps_2 %>% 
  group_by(sim) %>% 
  summarize(num_found = sum(!is.na(lod))) %>% 
  count(num_found)

# We find 2 peaks 53.7%, 1 peak 38.2% and 0 peaks 8.1%.

################################################################################
# Make P x G plots for each simulated SNP.
fig_out_dir = file.path(sim_out_dir, 'pxg_plots')
sim_loci  = read.csv(file.path(sim_out_dir, 'simulated_loci.csv'))[,-1]
sim_pheno = read.csv(file.path(sim_out_dir, 'simulated_phenotypes.csv'))
rownames(sim_pheno) = sim_pheno[,1]
sim_pheno = sim_pheno[,-1]

sim_nums = sim_loci %>% 
             distinct(sim) %>% 
             arrange(sim) %>% 
             pull(sim)

# Also get MAF for each SNP.
sim_maf = NULL

for(s in sim_nums) {
  
  current_pos = subset(sim_loci, sim == s)
  
  png(file.path(fig_out_dir, str_c('sim_', s, '_PxG.png')), width = 800, height = 800, res = 164)
  print(plot_pxg_snps(genoprobs = probs, pheno = sim_pheno[,s,drop = FALSE], 
                    addcovar = addcovar, map = map, 
                    pos = current_pos, query_func = snp_func))
  dev.off()
  
  sim_maf = bind_rows(sim_maf, get_snp_maf(genoprobs = probs, map = map, 
                                           pos = current_pos, query_func = snp_func))
  
} # for(s)

# NOTE: Overwriting 'simulated_loci.csv' here.
write.csv(sim_maf, file.path(sim_out_dir, 'simulated_loci.csv'))

##########
# Read in the association mapping found peaks and the new simulated loci with MAFs and 
# compare whether higher MAF is associated with finding the peak.
sim_loci = read.csv(file.path(sim_out_dir, 'simulated_loci.csv'))
found_snps_2 = read.csv(file = file.path(sim_out_dir, 'assoc_found_peak_summary.csv'))

sim_loci = full_join(sim_loci, dplyr::select(found_snps_2, -X, -hsq, -(chr:WSB_EiJ)), by = c('sim', 'snp_id')) %>% 
             mutate(peak_found = !is.na(lod))

png(file.path(fig_dir, 'snp_found_vs_maf.png'), width = 800, height = 800, res = 164)
ggplot(sim_loci, aes(peak_found, maf)) +
  geom_beeswarm(col = 'grey80') +
  geom_boxplot(fill = NA) +
  labs(title = 'Simulated SNP Found vs. Minor Allele Frequency',
       x = 'Simulated SNP Found', y = 'Minor Allele Frequency')
dev.off()

##########
# Sim 390 found both peaks with the simulated SNP as the top ranked LOD.
# Sim 11 found neither peak.
# Read them in and compare them.
assoc390 = readRDS(file.path(sim_out_dir, 'assoc_plots', 'sim390_assoc.rds'))
assoc11  = readRDS(file.path(sim_out_dir, 'assoc_plots', 'sim11_assoc.rds'))

# Expand to all SNPs.
tmp = assoc11$snpinfo[match(rownames(assoc11$lod), assoc11$snpinfo$snp_id),]
tmp = data.frame(snp_id = tmp$snp_id, pos = tmp$pos, lod = assoc11$lod[,1])

#snpmap = qtl2:::snpinfo_to_map(assoc11$snpinfo)
#tmp    = qtl2:::expand_snp_results(assoc11$lod, snpmap, assoc11$snpinfo)

tmp = subset(tmp, pos >= 100 & pos <= 101)

tmp_cor = get_snp_cor(genoprobs, map, chr = '5', start = min(tmp$pos), end = max(tmp$pos), query_func)
tmp_cor = tmp_cor[tmp$snp_id, tmp$snp_id]






