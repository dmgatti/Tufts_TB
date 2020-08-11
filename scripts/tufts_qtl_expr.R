# Look at each QTL interval and see if any of the differentially expressed genes
# from the microarray core are in the intervals.
options(stringsAsFactors = FALSE)
library(sva)
library(pcaMethods)
library(biomaRt)
library(readxl)
library(tidyverse)
library(qtl2)

mart = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')

base_dir    = '/media/dmgatti/hdb/projects/TB'
expr_dir    = file.path(base_dir, 'data', 'expression', 'excel')
results_dir = file.path(base_dir, 'results', 'qtl2', 'eqtl')
expr_files  = c('2020_Beamer_GB07.2_analysis_sans_snps.xlsx', '2020_Beamer_GB07.3_analysis_sans_snps.xlsx')
expr_sheets = c('2020_Beamer_GB07.2', '2020_Beamer_GB07.3')
expr_ranges = c('A2:CI26305', 'A2:CW26305')
batch_file  = file.path(base_dir, 'data', 'expression', 'microarray_batches.csv')
qtl_file    = file.path(base_dir, 'results/qtl2/gen_factor/tb_qtl_peaks.csv')
pheno_file  = file.path(base_dir, 'data', 'phenotypes', '2020-02-17 JDO TB lung phenotypes.xlsx')

source(file.path(base_dir, 'scripts', 'cis_eqtl_pipeline.R'))

# Read in QTL data, but only keep phenotypes and probs.
load(file = file.path(base_dir, "data/tufts_do_tb_qtl2_input.Rdata"))
pheno = read_xlsx(pheno_file) %>%
          dplyr::select(sample = `Mouse #`, 
                        dose   = `Mtb initial dose`,
                        class  = `Susceptibility Class`,
                        gen    = `DOB or generation`) %>%
          mutate(sample = as.character(sample))
pheno = as.data.frame(pheno)
pheno$sample[duplicated(pheno$sample)] = '1091'
rownames(pheno) = pheno$sample

# Read in the QTL peaks
peaks = readr::read_csv(qtl_file)

# Read in the expression data.
expr1 = read_xlsx(file.path(expr_dir, expr_files[1]), sheet = expr_sheets[1],
                 range = expr_ranges[1]) %>%
        select(`Brainarray probeset ID`:`KEGG Pathway(s)`, starts_with(as.character(0:9)))
expr2 = read_xlsx(file.path(expr_dir, expr_files[2]), sheet = expr_sheets[2],
                  range = expr_ranges[2]) %>%
        select(`Brainarray probeset ID`:`KEGG Pathway(s)`, starts_with(as.character(0:9)))

# It looks like Adam collapsed the replicates already. This is a reasonable decision, but
# it means that I can't do batch correction. I'm removing duplicates here.
batch = read.csv(batch_file) %>%
          mutate_all(as.character) %>%
          filter(!duplicated(sample))

expr = full_join(expr1, expr2)
rm(expr1, expr2)

# Reshape data & add in batch.
expr = expr %>%
         dplyr::rename(probeset = `Brainarray probeset ID`,
                       entrez   = `Mouse Entrez Gene ID`,
                       human_entrez = `Human Entrez Gene ID(s)`,
                       symbol   =  Symbol,
                       desc     = `Description`,
                       go       = `GO Term(s)`,
                       kegg     = `KEGG Pathway(s)`)

annot = expr %>%
          dplyr::select(probeset:kegg)
expr = expr %>%
         dplyr::select(entrez, starts_with(as.character(0:9))) %>%
         gather(sample, expr, -entrez) %>% 
         left_join(batch, by = 'sample') %>%
         left_join(select(pheno, sample, dose), by = 'sample')

# PCA of expression.
covar = expr %>%
          select(sample, batch, expt, dose) %>%
          mutate(dose = as.character(dose)) %>%
          distinct() %>%
          arrange(sample) %>%
          as.data.frame()
expr_mat = expr %>%
             select(-batch, -expt, -dose) %>%
             spread(sample, expr)

# Write out expression matrix for QTL mapping.
expr_mat = as.data.frame(expr_mat)
rn       = expr_mat$entrez
expr_mat = as.matrix(expr_mat[,-1])
rownames(expr_mat) = rn
saveRDS(expr_mat, file = file.path(base_dir, 'data', 'expression', 'tb_expr.rds'))
saveRDS(annot, file = file.path(base_dir, 'data', 'expression', 'tb_expr_annot.rds'))
saveRDS(covar, file = file.path(base_dir, 'data', 'expression', 'tb_expr_covar.rds'))

all(covar$sample == colnames(expr_mat))

pc_expr = pcaMethods::pca(expr_mat, nPcs = 50)
plot(loadings(pc_expr), pch = 16, col = covar$batch)
plot(loadings(pc_expr)[,2:3], pch = 16, col = covar$batch)
plot(loadings(pc_expr)[,2:3], pch = 16, col = as.numeric(factor(covar$dose)))

# Argh! Dose and batch are completely counfounded! There's no way to disentagle it.
table(covar$batch, covar$dose)

# I may use ComBat with the 0 dose mice and see if I can adjust for array batch.
# This assumes little difference between dose = 16 & 28.
mod = model.matrix(~dose, data = covar)
expr_cb = ComBat(dat = expr_mat, batch = covar$batch, mod = NULL, par.prior = TRUE)

pc_expr = pcaMethods::pca(expr_cb, nPcs = 50)
plot(loadings(pc_expr), pch = 16, col = covar$batch)
plot(loadings(pc_expr)[,2:3], pch = 16, col = covar$batch)
plot(loadings(pc_expr)[,2:3], pch = 16, col = as.numeric(factor(covar$dose)))

saveRDS(expr_cb, file = file.path(base_dir, 'data', 'expression', 'tb_expr_combat.rds'))

# Get cis-eQTL for non-normalized and ComBat normalized data. I'm going to assume that
# the one with with higher cis-eQTL values is the better normalization.

samples = intersect(colnames(expr_mat), rownames(probs[[1]]))
expr_mat = expr_mat[,samples]
expr_cb  = expr_cb[,samples]
addcovar = model.matrix(~gen, data = pheno)[,-1, drop = FALSE]
addcovar = addcovar[samples,]
addcovar = addcovar[,colSums(addcovar) > 0]

probs = probs[samples,]
for(i in seq_along(K)) {
  K[[i]] = K[[i]][samples, samples]
}

res_non = cis_eqtl(genoprobs = probs, map = map, expr = t(expr_mat), addcovar = addcovar, K = K, cores = 5)
save(res_non, file = file.path(results_dir, 'tufts_do_tb_cis_eqtl.Rdata'))

res_cb  = cis_eqtl(genoprobs = probs, map = map, expr = t(expr_cb),  addcovar = addcovar, K = K, cores = 5)
save(res_cb, file = file.path(results_dir, 'tufts_do_tb_cis_eqtl_combat.Rdata'))

# Plot the difference between the two.
plot(res_non$qtl_lod - res_cb$qtl_lod)
abline(h = mean(res_non$qtl_lod - res_cb$qtl_lod), col = 2)
# It doesn't look like ComBat normalization is making a difference. 
# Use the non batch normalized data.
rm(expr_cb)


# Get all of the genes, with entrez and ensembl IDs, in a given interval.
get_qtl_genes = function(chr, start, end) {
  
  filters = c('chromosome_name', 'start', 'end')
  values = list(chr, start, end)
  attributes = c('chromosome_name', 'start_position',  'end_position', 
                 'entrezgene_id',   'ensembl_gene_id', 'mgi_symbol')
  result = getBM(attributes = attributes, filters = filters, values = values,
                 mart = mart)
  return(result)
  
} # get_qtl_genes()

# Get the maximum marginal genotype.
get_max_geno = function(pr) {
  
  r  = apply(pr, 1, rank)
  wh = t(apply(apply(r, 2, '>=', 7), 2, which))
  mat = matrix(colnames(pr)[wh], ncol = 2, dimnames = list(rownames(pr), NULL))
  return(apply(mat, 1, str_c, collapse = ''))
  
} # get_max_geno()

# For each phenotype QTL, look for genes with cis-eQTL within the QTL interval.
for(i in 1:nrow(peaks)) {
  
  
  
} # for(i)

