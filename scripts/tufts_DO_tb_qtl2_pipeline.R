################################################################################
# Tufts QTL mapping using new pipeline scripts.
# May 25, 2020
# Daniel Gatti
# dmgatti@coa.edu
################################################################################
options(stringsAsFactors = FALSE)

# Libraries
library(pcaMethods)
library(readxl)
library(tidyverse)
library(qtl2)

# Directories & Files
base_dir   = '/media/dmgatti/hdb/projects/TB'
input_file = file.path(base_dir, 'data', 'tufts_do_tb_qtl2_input.Rdata')
pheno_file = file.path(base_dir, 'data', 'phenotypes', '2020-02-17 JDO TB lung phenotypes.xlsx')
snp_file   = '/media/dmgatti/hda/data/MUGA/cc_variants.sqlite'
result_dir = file.path(base_dir, 'results', 'qtl2', 'gen_factor')

source(file.path(base_dir, 'scripts', 'gwas_pipeline.R'))

load(input_file)

# Create generation as factor covariates. Dose is confounded with this and
# doesn't help the results.
covar$gen = factor(covar$gen)

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]

# GWAS
run_gwas(genoprobs = probs, 
         map       = map, 
         pheno     = pheno, 
         kinship   = K, 
         addcovar  = addcovar, 
         snp_file  = snp_file, 
         cores     = 4, 
         out_dir   = result_dir,
         verbose   = TRUE)

##############
# Permutations.

pheno_name = "cxcl5"

samples2use = which(!is.na(pheno[,pheno_name]))

print(paste0(pheno_name, " : ", length(samples2use)))

perms =
if(!file.exists(file.path(result_dir, "gwas_perms.rds"))) {
   tmpK = K
   for(j in 1:length(K)) {
     tmpK[[j]] = K[[j]][samples2use, samples2use]
   }
   tmp_covar = addcovar[samples2use,,drop = FALSE]
   tmp_covar = tmp_covar[,colSums(tmp_covar) > 0,drop = FALSE]

   perms = qtl2::scan1perm(genoprobs = probs[samples2use,], 
                           pheno     = pheno[samples2use,pheno_name, drop = FALSE], 
                           kinship   = tmpK, 
                           addcovar  = tmp_covar, 
                           cores     = 10,
                           n_perm    = 1000)
   saveRDS(perms, file = file.path(result_dir, "gwas_perms.rds"))
} else {
  perms = readRDS(file = file.path(result_dir, "gwas_perms.rds"))
} # else

# GWAS plots.
h5files     = dir(result_dir, pattern = 'h5$')
pheno_names = str_replace(h5files, '_gwas.h5', '')

gwas_files = data.frame(h5files, pheno_names)

for(i in 1:nrow(gwas_files)) {

  print(gwas_files$pheno_names[i])
  
  png(file.path(result_dir, str_replace(gwas_files$h5files[i], 'h5$', 'png')),
      width = 2000, height = 1200, res = 128)
  print(plot_gwas(h5filename = file.path(result_dir, gwas_files$h5files[i]), 
            thr        = quantile(perms, probs = 0.95), 
            title      = gwas_files$pheno_names[i]))
  dev.off()

} # for(i)

# Harvest peaks.
for(i in 1:nrow(gwas_files)) {
  
  print(gwas_files$pheno_names[i])
  
  top_snps = harvest_gwas(file.path(result_dir, gwas_files$h5files[i]), thr = quantile(perms, probs = 0.95), ensembl_ver = 93)
  write_csv(top_snps, path = file.path(result_dir, str_replace(gwas_files$h5files[i], '\\.h5', '_top_snps.csv')))
  
} # for(i)

# Gather them all up into one file.
files = dir(result_dir, pattern = 'gwas_top_snps.csv$')
files = files[files != 'do_tb_gwas_top_snps.csv']

results = NULL

for(f in files) {
  
  pheno_name = str_replace(f, '_gwas_top_snps.csv$', '')
  res = read_csv(file.path(result_dir, f))
  
  if(nrow(res) > 0) {
    
    res = data.frame(phenotype = pheno_name, res) %>%
            select(-index, -interval)
    results = bind_rows(results, res)
    
  } # if(nrow(res) > 0)
  
} # for(f)

results$pos = as.numeric(results$pos)
results$lod = as.numeric(results$lod)

write_csv(results, path = file.path(result_dir, 'do_tb_gwas_top_snps.csv'))

results_summ = summarize_gwas(results, thr = quantile(perms, probs = 0.95))

write_csv(results_summ, path = file.path(result_dir, 'do_tb_gwas_summary.csv'))

#####################################
# GWAS of PCA.
# We don't want too much missing data here. 
# Remove colums with > 40% missing data.
pheno4pca = pheno[,colMeans(is.na(pheno)) < 0.4]
pc_pheno  = pca(pheno4pca, method = 'ppca', nPcs = ncol(pheno4pca))

png(file.path(result_dir, 'pca_pct_var.png'))
pctvar = sDev(pc_pheno)
pctvar = pctvar^1 / sum(pctvar^2)
barplot(pctvar, names = factor(names(pctvar)), las = 1)
dev.off()

png(file.path(result_dir, 'pheno_cor.png'))
corrplot(cor(pheno4pca, use = 'pair'), order = 'hclust')
dev.off()


run_gwas(genoprobs = probs, 
         map       = map, 
         pheno     = scores(pc_pheno),
         kinship   = K, 
         addcovar  = addcovar, 
         snp_file  = '/media/dmgatti/hda/data/MUGA/cc_variants.sqlite', 
         cores     = 4,
         out_dir   = result_dir,
         verbose   = TRUE)

# PCA GWAS plots.
h5files     = dir(result_dir, pattern = '^PC(.)+h5$')
pheno_names = str_replace(h5files, '_gwas.h5', '')

gwas_files = data.frame(h5files, pheno_names)

for(i in 1:nrow(gwas_files)) {
  
  print(gwas_files$pheno_names[i])
  
  png(file.path(result_dir, str_replace(gwas_files$h5files[i], 'h5$', 'png')),
      width = 2000, height = 1200, res = 128)
  print(plot_gwas(h5filename = file.path(result_dir, gwas_files$h5files[i]), 
                  thr        = quantile(perms, probs = 0.95), 
                  title      = gwas_files$pheno_names[i]))
  dev.off()
  
} # for(i)

# Harvest PCA GWAS peaks.
for(i in 1:nrow(gwas_files)) {
  
  print(gwas_files$pheno_names[i])
  
  top_snps = harvest_gwas(file.path(result_dir, gwas_files$h5files[i]), thr = quantile(perms, probs = 0.95), ensembl_ver = 93)
  write_csv(top_snps, path = file.path(result_dir, str_replace(gwas_files$h5files[i], '\\.h5', '_top_snps.csv')))
  
} # for(i)

# Gather them all up into one file.
files = dir(result_dir, pattern = 'PC[0-9]_gwas_top_snps.csv$')
files = files[files != 'do_tb_pca_gwas_top_snps.csv']

results = NULL

for(f in files) {
  
  pheno_name = str_replace(f, '_gwas_top_snps.csv$', '')
  res = read_csv(file.path(result_dir, f))
  
  if(nrow(res) > 0) {
    
    res = data.frame(phenotype = pheno_name, res) %>%
      select(-index, -interval)
    results = bind_rows(results, res)
    
  } # if(nrow(res) > 0)
  
} # for(f)

results$pos = as.numeric(results$pos)
results$lod = as.numeric(results$lod)

write_csv(results, path = file.path(result_dir, 'do_tb_pca_gwas_top_snps.csv'))

results_summ = summarize_gwas(results, thr = quantile(perms, probs = 0.95))

write_csv(results_summ, path = file.path(result_dir, 'do_tb_pca_gwas_summary.csv'))

################################################################################
# Survival GWAS.
# Load in phenotypes and genoprobs.
pheno_orig = read_xlsx(pheno_file)
pheno_orig$`Mouse #`[1071] = 1091  # Fixing sample ID.
pheno_surv = pheno_orig %>%
             filter(`Mouse strain` == 'J:DO') %>%
             select(mouse = `Mouse #`, gen = `DOB or generation`, 
                    survival = `Day of Euthanasia`, event = `Euthanized due to pulmonary TB`) %>%
             mutate(mouse    = as.character(mouse),
                    survival = as.numeric(survival),
                    event    = (event == 'y') * 1,
                    gen      = if_else(gen == 'na', NA_character_, gen)) %>%
             filter(!is.na(survival)) %>%
             as.data.frame(pheno_surv)
rownames(pheno_surv) = pheno_surv$mouse

addcovar = model.matrix(~gen, data = pheno_surv)[, -1, drop = FALSE]

probs = readRDS('/media/dmgatti/hdb/projects/TB/haplo_reconstr/tufts_do_alleleprobs.rds')
probs = probs[pheno_surv$mouse,]

# Synch samples.
samples = intersect(rownames(pheno_surv), rownames(probs[[1]]))
pheno_surv = pheno_surv[samples,]
addcovar   = addcovar[samples,]
probs      = probs[samples,]
  
# Run the survival model (takes a while).
# genoprobs = probs; map = map; pheno = pheno_surv; kinship = K; addcovar = addcovar; snp_file = snp_file; cores = 1; out_dir = result_dir; verbose = FALSE
gwas_surv(genoprobs = probs, map = map, pheno = pheno_surv, kinship = K, addcovar = addcovar, 
          snp_file = snp_file, cores = 1, out_dir = result_dir, verbose = FALSE)

# Survival GWAS plot.
h5filename = file.path(result_dir, 'survival_gwas.h5')
png(file.path(result_dir, str_replace(h5filename, 'h5', 'png')), width = 2000, height = 1200, res = 128)
print(plot_gwas(h5filename = h5filename, 
                thr = 6, 
                title = 'Tufts DO TB Survival GWAS', 
                ensembl_ver = 93))
dev.off()

# Zoom in on chr 1 and plot with genes.
# Get the chromosomes from the file. (Called 'groups' in hdf5 lingo)
grps    = h5ls(file = h5filename, recursive = FALSE)

lod_all  = h5read(h5filename, paste('chr1', 'lod',     sep = '/'))
info_all = h5read(h5filename, paste('chr1', 'snpinfo', sep = '/'))

genes = get_ensembl(93)
chr   = 1
st    = 135
end   = 137
genes = subset(genes, chr == 1 & stop >= st & start <= end)

png(file.path(result_dir, 'survival_gwas_chr1.png'), width = 2000, height = 1800, res = 200)
qtl2::plot_snpasso(scan1output = lod_all, snpinfo = info_all, genes = genes, xlim = c(st, end), 
                   show_all_snps = TRUE, top_panel_prop = 0.5, color = 'black')
dev.off()

# Get top SNPs.
snp_map  = qtl2:::snpinfo_to_map(info_all)
snp_exp  = qtl2:::expand_snp_results(lod_all, snp_map, info_all)
info_all = data.frame(info_all, lod = snp_exp$lod)
  
top_snps_chr1 = subset(info_all, lod > 6.5)
top_snps_chr1 = top_snps_chr1 %>%
                  left_join(genes, by = c('ensembl_gene' = 'ensembl'))
count(top_snps_chr1, Name, consequence)

