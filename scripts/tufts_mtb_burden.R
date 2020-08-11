# Mtb burden on Chr 5.
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

# Go wide on Chr 5.
curr_chr = 5
start = 15
end   = 25
pheno_name = 'mtb_burden'
assoc = scan1snps(genoprobs = probs[,curr_chr], pheno = pheno_rz[,pheno_name, drop = FALSE], 
                  kinship = K[[curr_chr]], addcovar = covar, map = map, chr = curr_chr, 
                  start = max(1, start - 1), end = end + 1, query_func = snp_func, 
                  keep_all_snps = TRUE, cores = 10)

# Plot peaks and genes.
genes = gene_func(chr = curr_chr, start = start - 1, end = end + 1)
plot_snpasso(scan1output = assoc$lod, snpinfo = assoc$snpinfo, genes = genes, drop_hilit = 1, 
             top_panel_prop = 0.35, main = pheno_name)

top = top_snps(assoc$lod, assoc$snpinfo, drop = 0.5)

# Get SNPs in QTL interval and look for 129 private SNPs.
snps = snp_func(curr_chr, 18, 22)
snps = subset(snps, rowSums(snps[,10] != snps[,8:15]) > 6)

