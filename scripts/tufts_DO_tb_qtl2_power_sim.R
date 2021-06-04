################################################################################
# Tufts TB QTL Power simulations.
# DMG
# Mar 12, 2021
##############################
library(survival)
library(tidyverse)
library(qtl2)

base_dir   = '/media/dmgatti/hdb/projects/TB'
input_file = file.path(base_dir, 'data', 'tufts_do_tb_qtl2_input.Rdata') 
qtl_dir    = file.path(base_dir, 'results', 'qtl2', 'gen_factor2')
result_dir = file.path(base_dir, 'results', 'power_sim')

# Read in the DOQTL formatted data.
load(input_file)

peaks = read_csv(file.path(qtl_dir, 'tb_qtl_peaks.csv'))

##########
# Heritability.

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]

herit = data.frame(pheno = colnames(pheno_rz),
                   herit = 0)

K_all = calc_kinship(probs, type = 'overall', cores = 1)

for(i in 1:ncol(pheno_rz)) {
  
  herit$herit[i] = qtl2::est_herit(pheno = pheno_rz[,i,drop = FALSE], 
                                   kinship = K_all, addcovar = addcovar)
  
} # for(i)

write_csv(herit, file.path(result_dir, 'tufts_do_tb_herit.csv'))

##########
# Load in perms for thresholds.
perms = readRDS(file.path(qtl_dir, 'perms.rds'))
thr = quantile(perms, probs = 0.95)

# Go through each phenotype's QTL results and get the highest peak.

rdata_files = dir(path = qtl_dir, pattern = 'Rdata$')
# Don't keep the principal components.
rdata_files = rdata_files[!startsWith(rdata_files, 'PC')]
rdata_files = rdata_files[-grep('60day|assoc', rdata_files)]

peaks = data.frame(lodindex = 0,
                   lodcolumn = '',
                   chr = '0',
                   pos = 0,
                   lod = 0)

for(f in rdata_files) {

  load(file.path(qtl_dir, f))
  
  p = find_peaks(lod, map, threshold = 5)
  p = subset(p, lod == max(p$lod))
  
  peaks = rbind(peaks, p)

  rm(lod, coefs, assocs)
  gc()
     
} # for(f)

peaks = peaks[-1,]
colnames(peaks)[2] = 'phenotype'

# For each peak, run 100 simulations of 200, 400 and 600 samples each.
# Only map the location of the maximum QTL peak & save the LODs.

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]

sample_sizes = c(100, 200, 300, 400, 500, 600)
nperm = 100

for(i in 1:nrow(peaks)) {
  
  pheno_name = peaks$phenotype[i]
  chr        = peaks$chr[i]
  pos        = peaks$pos[i]
  
  samples2use = which(!is.na(pheno_rz[,pheno_name]))
  
  print(str_c(pheno_name, " : ", length(samples2use)))
  
  # Only work with phenotypes that have at least 600 samples.
  if(length(samples2use) < 600) {
    next
  } # if(samples2use < 600)
  
  # Get the subset of complete data for this phenotype.
  tmp_pheno = pheno_rz[samples2use,,drop = FALSE]
  # Pull genoprobs for current marker.
  pr = pull_genoprobpos(probs, map, chr, pos)
  
  # Just get K for the current chr.
  tmpK = K[[chr]][samples2use, samples2use]
  tmp_covar = addcovar[samples2use,,drop = FALSE]
  tmp_covar = tmp_covar[,colSums(tmp_covar) > 0, drop = FALSE]
  
  sim_result = data.frame(chr = chr, pos = pos, 
                          samples = rep(sample_sizes, each = nperm),
                          lod = rep(0, length(sample_sizes) * nperm))
  sim_index = 1
  
  for(ss in sample_sizes) {
    
    print(paste(ss, 'samples'))
    
    for(p in 1:nperm) {
  
      samples = sample(rownames(tmp_pheno), size = ss)
      
      print(p)
      lod = qtl2::fit1(genoprobs  = pr[samples,], 
                        pheno     = pheno_rz[samples, pheno_name, drop = FALSE], 
                        kinship   = tmpK[samples, samples], 
                        addcovar  = tmp_covar[samples,,drop = F])
      sim_result$lod[sim_index] = lod$lod
      
      write.csv(sim_result, file = file.path(result_dir, str_c(pheno_name, '_power_sim.csv')),
                row.names = FALSE, quote = FALSE)
      
      sim_index = sim_index + 1
      
    } # for(p)
  } # for(ss)
} # for(i)

# Make plots of LOD vs sample size with 0.05 threshold.
files = dir(result_dir, pattern = '_sim.csv$')
for(f in files) {

  pheno_name = str_replace(f, '_power_sim\\.csv$', '') %>% 
               str_replace('_', ' ') %>% 
               str_to_title()
  
  res = read.csv(file.path(result_dir, f))
    
  png(file.path(result_dir, sub('csv', 'png', f)), width = 800, height = 600, res = 124)
  print(ggplot(res, aes(samples, lod)) +
          geom_boxplot(aes(group = samples)) +
          geom_hline(yintercept = thr, color = 'red', size = 1) +
          scale_x_continuous(breaks = 1:6 * 100) +
          labs(title = pheno_name, x = 'Num. Samples', y = 'LOD'))
  dev.off()

} # for(f)





