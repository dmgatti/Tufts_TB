# Fit the survival model censoring at different dates.
# I'm trying to see if the Chr 1 signal changes if we change the censoring time.
# If the signal gets stronger when we censor closer to 50 days, then
# the peak may be due to the early "death step".

library(survival)
library(tidyverse)
library(qtl2)

base_dir   = '/media/dmgatti/hdb/projects/TB'
input_file = file.path(base_dir, 'data', 'tufts_do_tb_qtl2_input.Rdata') 
result_dir = file.path(base_dir, 'results', 'qtl2', 'survival_test')

# Read in the DOQTL formatted data.
load(input_file)

# CC SNP and gene database files (from Karl).
ccsnpdb = "/media/dmgatti/hda/data/MUGA/cc_variants.sqlite"
mgidb   = "/media/dmgatti/hda/data/MUGA/mouse_genes_mgi.sqlite"
snp_func  = create_variant_query_func(dbfile = ccsnpdb)
gene_func = create_gene_query_func(dbfile = mgidb)
ensembl = 93

peaks = NULL

##########
# QTL Mapping with coef calculations.

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]

################################################################################
# Survival mapping
# Set up survival model.
pheno_surv = data.frame(mouse = as.character(covar$mouse),
                        surv = Surv(time = covar$euth_day, event = (covar$euth == 'y') * 1))
rownames(pheno_surv) = pheno_surv$mouse
pheno_surv = pheno_surv[complete.cases(pheno_surv),]

# Subset to keep only samples with data.
samples = intersect(rownames(pheno_surv), rownames(addcovar))
pheno_surv = pheno_surv[samples,]
probs = probs[samples,]
for(i in seq_along(K)) {
  K[[i]] = K[[i]][samples,samples]
}

# Map survival censoring at 

