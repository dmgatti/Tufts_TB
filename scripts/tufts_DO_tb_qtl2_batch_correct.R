##############################
# Tufts TB QTL Mapping.
# DMG
# Mar. 20, 2019
##############################
library(tidyverse)
library(qtl2)
library(qtl2convert)

# CC SNP and gene database files (from Karl).
ccsnpdb = "/media/dmgatti/hda/data/MUGA/cc_variants.sqlite"
mgidb   = "/media/dmgatti/hda/data/MUGA/mouse_genes_mgi.sqlite"
snp_func  = create_variant_query_func(dbfile = ccsnpdb)
gene_func = create_gene_query_func(dbfile = mgidb)
ensembl = 93

base_dir = "/media/dmgatti/hdb/projects/TB/"

# Read in the DOQTL formatted data.
load(stringr::str_c(base_dir, "data/phenotypes/GB_Tufts_mapping_input.Rdata"))
markers = snps[rownames(snps) %in% dimnames(probs)[[3]],]

rm(pheno, pheno.rz, covar)

pheno    = read.csv(str_c(base_dir, "data/phenotypes/tb_pheno_log_batch_correct.csv"))
#pheno_rz = read.csv(str_c(base_dir, "data/phenotypes/tb_pheno_rz.csv"))

# Make dose a factor.
pheno$dose = factor(pheno$dose)

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

pheno_rz = pheno
for(i in 5:ncol(pheno_rz)) {
  pheno_rz[[i]] = rankZ(pheno_rz[[i]])
} # for(i)

rownames(pheno)    = pheno$mouse
rownames(pheno_rz) = pheno_rz$mouse

# Synch samples.
samples2keep = intersect(rownames(pheno), rownames(probs))
pheno    = pheno[samples2keep,]
pheno_rz = pheno_rz[samples2keep,]
probs    = probs[samples2keep,,]
pheno$dose    = factor(pheno$dose,    levels = c(16, 28, 97, 127))
pheno_rz$dose = factor(pheno_rz$dose, levels = c(16, 28, 97, 127))

# Convert to qtl2 format and resave.
probs = qtl2convert::probs_doqtl_to_qtl2(probs = probs, map = markers, pos_column = "pos")
map   = qtl2convert::map_df_to_list(map = markers, pos_column = "pos")
K = qtl2::calc_kinship(probs = probs, type = "loco", cores = 8, quiet = FALSE)
covar = model.matrix(~expt+dose, data = pheno)[,-1,drop = FALSE]

save(pheno, pheno_rz, covar, probs, map, K, markers, file = str_c(base_dir, "data/phenotypes/GB_Tufts_qtl2_input_batch_correct.Rdata"))

#############################
# Start here once the above has been run.
library(tidyverse)
library(qtl2)
library(qtl2convert)

base_dir = "/media/dmgatti/hdb/projects/TB/"

load(file = str_c(base_dir, "data/phenotypes/GB_Tufts_qtl2_input_batch_correct.Rdata"))

##########
# QTL Mapping with coef calculations.

result_dir = str_c(base_dir, "results/qtl2/dose_as_factor_batch_correct/")
peaks = NULL

for(i in 6:ncol(pheno_rz)) {
  
  pheno_name = colnames(pheno_rz)[i]
  
  samples2use = which(!is.na(pheno_rz[,pheno_name]))
  
  print(str_c(pheno_name, " : ", length(samples2use)))
  
  tmpK = K
  for(j in 1:length(K)) {
    tmpK[[j]] = K[[j]][samples2use, samples2use]
  }
  tmp_covar = covar[samples2use,,drop = FALSE]
  tmp_covar = tmp_covar[,colSums(tmp_covar) > 0, drop = FALSE]
  
  lod = qtl2::scan1(genoprobs = probs[samples2use,], 
                    pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE], 
                    kinship   = tmpK, 
                    addcovar  = tmp_covar, 
                    cores = 10)
  
  png(str_c(result_dir, pheno_name, "_qtl.png"), width = 1000, height = 800, res = 128)
  qtl2::plot_scan1(lod, map = map, main = pheno_name)
  dev.off()
  
  current_peaks = find_peaks(lod, map, threshold = 6, prob = 0.95)
  peaks = rbind(peaks, current_peaks)

  for(j in 1:nrow(current_peaks)) {
    
    curr_chr = current_peaks$chr[j]
    print(str_c("    chr ", curr_chr))
    start = current_peaks$ci_lo[j]
    end   = current_peaks$ci_hi[j]
    
    coefs = qtl2::scan1blup(genoprobs = probs[samples2use, curr_chr], 
                                 pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE], 
                                 kinship   = K[[curr_chr]], 
                                 addcovar  = tmp_covar,
                                 se = TRUE,
                                 cores = 10)
    stopifnot(!is.nan(coefs))

    png(str_c(result_dir, pheno_name, "_coef_chr", curr_chr, ".png"), width = 1000, height = 800, res = 128)
    plot_coefCC(coefs, map, scan1_output = lod, top_panel_prop = 0.6)
    dev.off()
    
    assoc = scan1snps(genoprobs = probs[,curr_chr], pheno = pheno_rz[,pheno_name, drop = FALSE], kinship = tmpK[[curr_chr]], addcovar = tmp_covar,
                      map = map, chr = curr_chr, start = start - 1, end = end + 1, query_func = snp_func, 
                      keep_all_snps = TRUE, cores = 10)
    
    genes = gene_func(chr = curr_chr, start = start - 1, end = end + 1)
    png(str_c(result_dir, pheno_name, "_assoc_chr", curr_chr, ".png"), width = 1000, height = 800, res = 128)
    plot_snpasso(scan1output = assoc$lod, snpinfo = assoc$snpinfo, genes = genes, drop_hilit = 1, top_panel_prop = 0.5, main = pheno_name)
    dev.off()
  
  } # for(j)
  
  save(lod, coefs, file = str_c(result_dir, pheno_name, "_qtl2.Rdata"))
  
} # for(i)

write_csv(peaks, path = str_c(result_dir, "tb_qtl_peaks.csv"))

##############
# Permutations.

pheno_name = "cxcl5"

samples2use = which(!is.na(pheno_rz[,pheno_name]))

print(str_c(pheno_name, " : ", length(samples2use)))

tmpK = K
for(j in 1:length(K)) {
  tmpK[[j]] = K[[j]][samples2use, samples2use]
}
tmp_covar = covar[samples2use,,drop = FALSE]
tmp_covar = tmp_covar[,colSums(tmp_covar) > 0,drop = FALSE]

perms = qtl2::scan1perm(genoprobs = probs[samples2use,], 
                        pheno = pheno_rz[samples2use,pheno_name, drop = FALSE], 
                        kinship = tmpK, 
                        addcovar = tmp_covar, 
                        cores = 10,
                        n_perm = 1000)
saveRDS(perms, file = str_c(result_dir, "perms.rds"))

