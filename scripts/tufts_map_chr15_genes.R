# Map the traits that have QTL on Chr 5.
# cxcl2, tnf, s100a8, cxcl1, cxcl5
library(tidyverse)
library(qtl2)
library(qtl2convert)

# CC SNP and gene database files (from Karl).
ccsnpdb = "/home/dmgatti/Documents/data/qtl2/cc_variants.sqlite"
mgidb   = "/home/dmgatti/Documents/data/qtl2/mouse_genes_mgi.sqlite"
snp_func  = create_variant_query_func(dbfile = ccsnpdb)
gene_func = create_gene_query_func(dbfile = mgidb)
ensembl = 93

base_dir = "/home/dmgatti/Documents/data/TB/"

load(file = str_c(base_dir, "data/phenotypes/GB_Tufts_qtl2_input.Rdata"))

# Make a single phenotype that is the mean of all 5.
pheno_names = c('cxcl2', 'tnf', 's100a8', 'cxcl1', 'cxcl5')
pheno_chr5 = rowMeans(pheno_rz[,pheno_names], na.rm = TRUE)
pheno_chr5 = matrix(pheno_chr5, ncol = 1, dimnames = list(rownames(pheno_rz), 'pheno_means'))

lod = qtl2::scan1(genoprobs = probs, 
                  pheno     = pheno_chr5[, 'pheno_means', drop = FALSE], 
                  kinship   = K, 
                  addcovar  = covar, 
                  cores = 10)
plot_scan1(lod, map)

# Add in the haplotypes on Chr 17.
peaks = find_peaks(lod, map, 7)

map17 = map[[peaks$chr[2]]]
mkr = names(map17)[which.min(abs(map17 - peaks$pos[2]))]
rm(map17)

gt17 = probs[[17]][,,mkr]

lod17 = qtl2::scan1(genoprobs = probs, 
                    pheno     = pheno_chr5[, 'pheno_means', drop = FALSE], 
                    kinship   = K, 
                    addcovar  = cbind(covar, gt17), 
                    cores = 10)
plot_scan1(lod17, map)
find_peaks(lod17, map, 7)

plot_scan1(lod, map)
plot_scan1(lod17, map, add = TRUE, col = 3)

# The end result is that the genotype of mice on chr 17 does NOT affect the 
# peak on Chr 15.

###########################################################################
# Look at Mtb_burden
lod = qtl2::scan1(genoprobs = probs, 
                  pheno     = pheno_rz[, 'mtb_burden', drop = FALSE], 
                  kinship   = K, 
                  addcovar  = covar, 
                  cores = 10)
plot_scan1(lod, map)

# Add in the haplotypes on Chr 17.
peaks = find_peaks(lod, map, 4.8)

map17 = map[[peaks$chr[7]]]
mkr = names(map17)[which.min(abs(map17 - peaks$pos[7]))]
rm(map17)

gt17 = probs[[17]][,,mkr]

lod17 = qtl2::scan1(genoprobs = probs, 
                    pheno     = pheno_rz[, 'mtb_burden', drop = FALSE], 
                    kinship   = K, 
                    addcovar  = cbind(covar, gt17), 
                    cores = 10)
plot_scan1(lod17, map)
find_peaks(lod17, map, 7)

plot_scan1(lod, map, ylim = c(0, 8))
plot_scan1(lod17, map, add = TRUE, col = 3)

# Slight LOD increase from 7.58 to 7.89 on Chr 5. Too small to really get excited about.

# The end result is that the genotype of mice on chr 17 does NOT affect the 
# peak on Chr 5.





