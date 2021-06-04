################################################################################
# Tufts TB QTL Mapping. Image data from Metin's group.
# DMG
# May 20, 2020
##############################
library(biomaRt)
library(survival)
library(tidyverse)
library(qtl2)

base_dir   = '/media/dmgatti/hdb/projects/TB'
input_file = file.path(base_dir, 'data', 'tufts_do_tb_qtl2_input.Rdata') 
result_dir = file.path(base_dir, 'results', 'qtl2', 'gen_factor_all_coefs')

# Set up biomaRt.
ensembl = useEnsembl(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
#bmfilt = listFilters(ensembl)
#bmattr = listAttributes(ensembl)

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
# Heritability.

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]

herit = data.frame(pheno = colnames(pheno_rz),
                   herit = 0)

K_all = calc_kinship(probs, type = 'overall', cores = 4)

for(i in 1:ncol(pheno_rz)) {
  
  herit$herit[i] = qtl2::est_herit(pheno = pheno_rz[,i,drop = FALSE], 
                                   kinship = K_all, addcovar = addcovar)
  
} # for(i)

write_csv(herit, file.path(result_dir, 'tufts_do_tb_herit.csv'))

##########
# QTL Mapping with coef calculations.

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]

for(i in 1:ncol(pheno_rz)) {
  
  pheno_name = colnames(pheno_rz)[i]
  
  samples2use = which(!is.na(pheno_rz[,pheno_name]))
  
  print(str_c(pheno_name, " : ", length(samples2use)))
  
  tmpK = K
  for(j in 1:length(K)) {
    tmpK[[j]] = K[[j]][samples2use, samples2use]
  }
  tmp_covar = addcovar[samples2use,,drop = FALSE]
  tmp_covar = tmp_covar[,colSums(tmp_covar) > 0, drop = FALSE]
  
  print("    QTL")  
  lod = qtl2::scan1(genoprobs = probs[samples2use,], 
                    pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE], 
                    kinship   = tmpK, 
                    addcovar  = tmp_covar, 
                    cores     = 10)
  
  png(file.path(result_dir, str_c(pheno_name, "_qtl.png")), width = 1000, height = 800, res = 128)
  qtl2::plot_scan1(lod, map = map, main = pheno_name)
  dev.off()
  
  # Estimate BLUPS on all chromosomes.
  coefs = vector('list', length(map))
  names(coefs) = names(map)

  # For each peak, plot coefs and association mapping.
  for(curr_chr in names(map)) {
      
    # BLUP mapping on one chromosome.
    print(paste("    coef", curr_chr))
    coefs[[curr_chr]] = qtl2::scan1blup(genoprobs = probs[samples2use, curr_chr], 
                                        pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE], 
                                        kinship   = K[[curr_chr]], 
                                        addcovar  = tmp_covar,
                                        se = TRUE,
                                        cores = 10)
    stopifnot(!is.nan(coefs[[j]]))
      
    png(file.path(result_dir, str_c(pheno_name, "_coef_chr", curr_chr, ".png")), width = 1000, height = 800, res = 128)
    plot_coefCC(coefs[[curr_chr]], map, scan1_output = lod, top_panel_prop = 0.6, main = pheno_name)
    dev.off()
      
  } # for(curr_chr)
  
  save(lod, coefs, file = file.path(result_dir, str_c(pheno_name, "_qtl2.Rdata")))
  
} # for(i)

