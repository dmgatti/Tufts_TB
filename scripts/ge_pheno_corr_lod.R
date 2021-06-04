################################################################################
# Data for Anna & Matt's functional QTL analysis.
# Correlate gene expression with physiological phenotypes and otuput adjusted
# p-values. 
# For each gene in genome, get GWAS LOD score for each phenotype.
# April 18, 2021
# Daniel Gatti
# dmgatti@coa.edu
################################################################################
options(stringsAsFactors = FALSE)
options(connectionObserver = NULL)

library(readxl)
library(tidyverse)
library(Mus.musculus)
library(biomaRt)
library(GenomicRanges)
library(qtl2)

base_dir  = '/media/dmgatti/hdb/projects/TB'
pheno_dir = file.path(base_dir, 'data', 'phenotypes')
expr_dir  = file.path(base_dir, 'data', 'expression')
result_dir = file.path(base_dir, 'results', 'expr')
input_file = file.path(base_dir, 'data', 'tufts_do_tb_qtl2_input.Rdata')
cc_snp_file = '/media/dmgatti/hda/data/MUGA/cc_variants.sqlite'
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

snp_func = qtl2::create_variant_query_func(cc_snp_file)

# Set up Biomart. Using Ensembl 102, GRCm38
mart = useEnsembl(biomart = 'genes', dataset = 'mmusculus_gene_ensembl', version = 102)
bmattr = listAttributes(mart)
bmfilt = listFilters(mart)

# Read in phenotypes.
pheno = read.csv(file.path(pheno_dir, 'tufts_do_tb_pheno_cleaned.csv'))

# Read in gene expression. Note that they used Entrez IDs.
expr1 = readxl::read_xlsx(file.path(expr_dir, 'excel', '2020_Beamer_GB07.2_analysis_sans_snps.xlsx'),
                          range = 'A2:CI26305') %>% 
        dplyr::select(entrez = `Mouse Entrez Gene ID`, symbol = Symbol, name = Description, starts_with('1'), starts_with('2'))
expr2 = readxl::read_xlsx(file.path(expr_dir, 'excel', '2020_Beamer_GB07.3_analysis_sans_snps.xlsx'),
                          range = 'A2:CW26305') %>% 
        dplyr::select(entrez = `Mouse Entrez Gene ID`, symbol = Symbol, name = Description, starts_with('2'), starts_with('3'), starts_with('4'))

expr = full_join(expr1, expr2) %>% 
          mutate(entrez = as.integer(entrez))
rm(expr1, expr2)

# Subset and align the samples in each data sets.
expr_annot = expr[,1:6]
expr_mat = as.matrix(t(expr[,-(1:6)]))
colnames(expr_mat) = expr_annot$entrez
samples = rownames(expr_mat)

pheno_mat = dplyr::select(pheno, pct_wt_loss:necr_ratio) %>% 
              as.matrix()

# Log transform the phenotypes.
pheno_mat = apply(pheno_mat, 2, log1p)
rownames(pheno_mat) = pheno$mouse

samples = intersect(samples, rownames(pheno_mat))

expr_mat  = expr_mat[samples,]
pheno_mat = pheno_mat[samples,]

stopifnot(all(nrow(expr_mat) == nrow(pheno_mat)))
stopifnot(all(rownames(expr_mat) == rownames(pheno_mat)))

covar = model.matrix(~gen, data = pheno)[,-1,drop = FALSE]
covar = covar[samples,]
covar = covar[,colSums(covar) > 0]

# Regress out batch effects.
for(i in 1:ncol(pheno_mat)) {
  mod = lm(pheno_mat[,i] ~ covar, na.action = na.exclude)
  pheno_mat[,i] = residuals(mod)
} # for(i)

for(i in 1:ncol(expr_mat)) {
  mod = lm(expr_mat[,i] ~ covar, na.action = na.exclude)
  expr_mat[,i] = residuals(mod)
} # for(i)

# Calculate gene/phenotype correlations.
n = nrow(expr_mat)
gp_cor = cor(pheno_mat, expr_mat, use = 'pairwise')
gp_t   = gp_cor * sqrt(n - 2) / sqrt(1.0 - gp_cor^2)
gp_p   = 2 * pt(abs(gp_t), df = n - 2, lower.tail = FALSE)
gp_p_adj = matrix(p.adjust(gp_p, method = 'holm'), nrow = nrow(gp_p), 
                  ncol = ncol(gp_p), dimnames = dimnames(gp_p))
gp_p_adj = t(gp_p_adj)
gp_p_adj = data.frame(entrez = rownames(gp_p_adj), gp_p_adj)
rownames(gp_p_adj) = NULL

write.csv(gp_p_adj, file.path(result_dir, 'tufts_do_tb_gene_phenotype_cor_pvalue_adj.csv'))

rm(bmattr, bmfilt, gp_cor, gp_p, gp_t, gp_p_adj)

# Need to remap each phenotype, one chr at a time to get the LOD for ALL SNPs.
# Then intersect with ALL gene locations.
load(input_file)

# Get genes from bioconductor.
tx = transcriptsBy(txdb, by = 'gene')
tx = lapply(tx, function(z) { GRanges(seqnames = seqnames(z)[1], 
                                      ranges = IRanges(start = min(start(z)), end = max(end(z)))) })
# unlist(tx) doesn't work here for some reason. 
seqs  = sapply(tx, seqnames)
seqs  = sapply(seqs, as.character)
seqs  = sub('^chr', '', seqs)
start = sapply(tx, start)
end   = sapply(tx, end)

gene_gr = GRanges(seqnames = seqs, 
                  ranges = IRanges(start, end),
                  entrez = names(tx))

addcovar = model.matrix(~gen, data = covar)[,-1]

for(chr in names(probs)) {
  
    print(paste('CHR:', chr))
  
    # Map all phenotypes on this chromosome.
    gwas = scan1snps(genoprobs = probs, map = map, pheno = pheno_rz, kinship = K, 
                     addcovar = addcovar, chr = chr, start = 0, end = 200, query_func = snp_func,
                     keep_all_snps = TRUE)
    
    # Expand to get all SNPs.
    snp_map = qtl2:::snpinfo_to_map(gwas$snpinfo)
    gwas    = qtl2:::expand_snp_results(gwas$lod, snp_map, gwas$snpinfo)
    
    lod_gr = GRanges(seqnames = chr, 
                     ranges   = IRanges(start = gwas$map[[1]] * 1e6, width = 1),
                     mcols    = gwas$lod)

    ol  = findOverlaps(gene_gr, lod_gr)
    tmp = cbind(gene_gr$entrez[queryHits(ol)], mcols(lod_gr)[subjectHits(ol),])
    colnames(tmp)[1] = 'entrez'
    colnames(tmp) = sub('^mcols\\.', '', colnames(tmp))
    tmp = split(tmp, tmp$entrez)
    tmp = lapply(tmp, apply, 2, max)
    tmp = matrix(unlist(tmp), nrow = length(tmp[[1]]), ncol = length(tmp),
                 dimnames = list(names(tmp[[1]]), names(tmp)))
    tmp = t(tmp[-1,])
    
    max_lod = tibble(entrez = rownames(tmp), as.data.frame(tmp))
    write_csv(max_lod, file.path(result_dir, 'tufts_do_tb_gene_lod.csv'), append = chr != '1')
    
    rm(gwas, tmp, snp_map, lod_gr, ol, max_lod)
    gc()
  
} # for(chr)
  

