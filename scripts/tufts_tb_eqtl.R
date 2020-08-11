# Start here once the above has been run.
options(stringsAsFactors = FALSE)
library(GenomicRanges)
library(biomaRt)
library(tidyverse)
library(qtl2)
library(qtl2convert)

# CC SNP and gene database files (from Karl).
ccsnpdb = "/home/dmgatti/Documents/data/qtl2/cc_variants.sqlite"
mgidb   = "/home/dmgatti/Documents/data/qtl2/mouse_genes_mgi.sqlite"
snp_func  = create_variant_query_func(dbfile = ccsnpdb)
gene_func = create_gene_query_func(dbfile = mgidb)
ensembl = 93

mart = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')

base_dir = "/media/dmgatti/hdb/projects/TB/"

load(file = str_c(base_dir, "data/phenotypes/GB_Tufts_qtl2_input.Rdata"))
peaks = readr::read_csv(str_c(base_dir, 'results/qtl2/dose_as_factor/tb_qtl_peaks_coefs.csv'))

# Load in expression data.
expr = read.csv(str_c(base_dir, 'data/expression/tb_expr.csv'))
rownames(expr) = expr[,1]
expr = as.matrix(expr[,-1])

# Rank Z transform expression.
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}
expr_rz = apply(expr, 2, rankZ)

# Synch up sample IDs
samples = intersect(rownames(expr_rz), rownames(pheno))
expr_rz = expr_rz[samples,]
covar = covar[samples,]
for(i in 1:length(probs)) {
  probs[[i]] = probs[[i]][samples,,]
  K[[i]]     = K[[i]][samples,samples]
}

# Filter to keep genes with the highest variance. Note that this data hasn't been
# de-SNPed since we can't yet get probe level data from affymetrix.
expr_mean = colMeans(expr)
expr_sd   = apply(expr, 2, sd)

plot(expr_mean, expr_sd)
hist(expr_sd, breaks = 100)

expr_rz = expr_rz[,expr_mean > 5]

##########
# QTL Mapping with just LOD.

result_dir = str_c(base_dir, "results/qtl2/eqtl/")

for(i in 1:11) {
  
  rng = ((i-1) * 1000 + 1):(i * 1000)
  if(max(rng) > ncol(expr_rz)) {
    rng = ((i-1) * 1000 + 1):ncol(expr_rz)
  }
  print(range(rng))
  
  eqtl = scan1(genoprobs = probs, pheno = expr_rz[,rng], kinship = K, 
               addcovar = covar, cores = 4)
  
  #saveRDS(eqtl, file = str_c(result_dir, 'tb_expr_qtl2_lod', i,'.rds'))
  rm(eqtl)
  gc()

} # for(i)

##########
# Summarize eQTL results.

# Read the data back in and keep only the maximum QTL location.
rds_files = dir(result_dir, pattern = 'tb_expr_qtl2_lod', full.names = TRUE)
qtl = data.frame(entrez = NA,
                 marker = NA,
                 qtl_chr = NA,
                 qtl_pos = NA,
                 lod    = NA)
for(f in rds_files) {

  print(f)
  q = readRDS(f)
  max_lod = apply(q, 2, max)
  max_mkr = apply(q, 2, which.max)
  
  qtl = rbind(qtl, data.frame(entrez = sub('^X', '', colnames(q)),
                              marker = rownames(q)[max_mkr],
                              qtl_chr = markers$chr[max_mkr],
                              qtl_pos = markers$pos[max_mkr],
                              lod    = max_lod))
  
} # for(f)

qtl = qtl[-1,]

# Add in gene locations.
attributes = c('entrezgene_id', 'ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position')
filters = 'entrezgene_id'
values = qtl$entrez
genes = getBM(attributes = attributes, filters = filters, values = values, mart = mart)
genes$entrezgene_id = as.character(genes$entrezgene_id)

qtl = full_join(qtl, genes, by = c('entrez' = 'entrezgene_id')) %>%
      dplyr::select(entrez, ensembl = ensembl_gene_id, symbol = external_gene_name, 
                    gene_chr = chromosome_name, gene_start = start_position,
                    marker:lod) %>%
      mutate(gene_start = gene_start * 1e-6) %>%
      filter(gene_chr == qtl_chr)

write_csv(qtl, path = str_c(result_dir, 'tb_eqtl_summary.csv'))

# Intersect phenotype QTL peaks with eQTL peaks with LOD > 8.
qtl = subset(qtl, qtl$lod >= 8.0)

eqtl_gr = GRanges(seqnames = qtl$qtl_chr, ranges = IRanges(start = qtl$qtl_pos * 1e6, width = 1), mcols = qtl[,c(1:3, 6:9)])
qtl_gr  = GRanges(seqnames = peaks$chr,   ranges = IRanges(start = (peaks$pos - 5) * 1e6, end = (peaks$pos + 5) * 1e6),
                  mcols = peaks[,c(1, 6, 4, 5)])

ol = findOverlaps(qtl_gr, eqtl_gr)

result = data.frame(qtl_gr[queryHits(ol),], eqtl_gr[subjectHits(ol),])
result = result %>%
           dplyr::select(phenotype = mcols.phenotype,
                  qtl_chr   = seqnames,
                  qtl_start = start,
                  qtl_end   = end,
                  qtl_lod   = mcols.lod,
                  qtl_mkr   = mcols.mkr,
                  entrez    = mcols.entrez,
                  ensmusg   = mcols.ensembl,
                  symbol    = mcols.symbol,
                  gene_chr  = seqnames.1,
                  gene_start = start.1,
                  eqtl_lod  = mcols.lod.1,
                  eqtl_mkr  = mcols.marker)
result$qtl_chr = as.character(result$qtl_chr)
result$gene_chr = as.character(result$gene_chr)

write_csv(result, path = str_c(result_dir, 'tb_eqtl_pheno_qtl_intersection.csv'))

##########
# Calculate BLUPs for the eQTL that intersect with the phenotype QTL.
expr_rz = expr_rz[,colnames(expr_rz) %in% paste0('X', result$entrez)]

for(j in 1:ncol(expr_rz)) {
  
  entrez = sub('^X', '', colnames(expr_rz)[j])
  curr_row = which(result$entrez == entrez)[1]
  ensembl = result$ensmusg[curr_row]
  symbol  = result$symbol[curr_row]
  curr_chr = result$chr[curr_row]
  
  print(str_c(entrez, ' chr ', curr_chr))

  coefs = qtl2::scan1blup(genoprobs = probs[, curr_chr],
                          pheno     = expr_rz[, j, drop = FALSE],
                          kinship   = K[[curr_chr]], 
                          addcovar  = covar,
                          se = TRUE,
                          cores = 8)
  stopifnot(!is.nan(coefs))
  
  png(str_c(result_dir, 'blup/', entrez, "_coef_chr", curr_chr, ".png"), width = 1000, height = 800, res = 128)
  plot_coefCC(coefs, map, scan1_output = lod, top_panel_prop = 0.6, main = str_c(entrez, ensembl, symbol, sep = ' '))
  dev.off()
} # for(j)





