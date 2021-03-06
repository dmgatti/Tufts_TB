################################################################################
# cis-eQTL pipeline.
# Given a gene expression and genoprobs for a set of samples, map the cis-eqtl
# for each gene.
# Daniel Gatti
# June 30, 2020
# dmgatti@coa.edu
################################################################################
options(stringsAsFactors = FALSE)

library(biomaRt)
library(qtl2)

# Arguments:
# genoprobs: list of qtl2-style genoprobs objects.
# map:       list of qtl2-style marker locations (in Mb).
# expr:     numeric matrix containing gene expression. Must have Ensembl IDs as colnames, samples in rownames.
# addcovar: numeric matrix containing additive covariates.
# K:     list of qtl2 style kinship matrices.
# cores: integer that is the number of cores to use.
# Returns: data.frame with cis-eQTL positions and LOD scores.
cis_eqtl = function(genoprobs, map, expr, addcovar, K, cores = 1) {
  
  # Get gene locations from Biomart.
  mart  = biomaRt::useEnsembl(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
  genes = get_gene_info(entrez = colnames(expr), mart = mart)
  
  # For now, remove genes with multiple locations.
  genes = subset(genes, !grepl(',', genes$start_position))
  
  # For cis-eQTL, remove genes on Chr M or Y.
  genes = subset(genes, !chromosome_name %in% c('M', 'Y'))
  
  # Intersect genes and Ensembl data.
  common_genes = intersect(colnames(expr), genes$entrezgene_id)
  common_genes = as.character(common_genes)
  
  print(paste('Keeping', length(common_genes), 'out of', ncol(expr)))
  
  expr  = expr[,common_genes]
  genes = genes[common_genes,]
  stopifnot(colnames(expr) == rownames(genes))
  
  genes = genes[order(genes$chromosome_name),]
  expr  = expr[,rownames(genes)]
  
  # Convert from bp to Mb.
  genes$start_position = as.numeric(genes$start_position) * 1e-6
  genes$end_position   = as.numeric(genes$end_position)   * 1e-6
  
  # Combine genes and expr into one data structure.
  ge = merge(genes, t(expr), by = 'row.names')
  ge = split(ge, ge$chromosome_name)
  
  results = NULL
  
  for(chr in names(ge)) {
    
    print(paste('Chr', chr))
    
    closest_mkr = abs(outer(map[[chr]], ge[[chr]]$start_position, '-'))
    closest_mkr = apply(closest_mkr, 2, which.min)
    closest_mkr = names(map[[chr]])[closest_mkr]
    
    # Make expression matrix.
    e = t(as.matrix(ge[[chr]][,-(1:9)]))
    colnames(e) = ge[[chr]]$entrezgene_id
    
    # Scan genes on current chromosome.
    t1 = proc.time()[3]
    qtl = scan1(genoprobs = genoprobs[,chr], pheno = e, kinship = K[[chr]], addcovar = addcovar, cores = cores)
    print(proc.time()[3] - t1)
    
    # Make results.
    results = rbind(results, data.frame(ge[[chr]][,2:9], 
                                        qtl_chr = chr,
                                        qtl_pos = map[[chr]][closest_mkr],
                                        qtl_lod = qtl[matrix(c(closest_mkr, colnames(qtl)), ncol = 2)]))
    
  } # for(chr)
  
  return(results)
  
} # cis_eqtl()


# Arguments:
# genoprobs: list of qtl2-style genoprobs objects.
# map:       list of qtl2-style marker locations (in Mb).
# expr:     numeric matrix containing gene expression. Must have Ensembl IDs as colnames, samples in rownames.
# addcovar: numeric matrix containing additive covariates.
# K:        list of qtl2 style kinship matrices.
# lod:     data.frame containing maximum cis-QTL LOD scores for each gene as produced by cis_eqtl() above.
# cores:   integer that is the number of cores to use.
# Returns: data.frame with cis-eQTL BLUP positions and LOD scores.
cis_blup = function(genoprobs, map, expr, addcovar, K, lod, cores = 1) {
  
  # Get gene locations from Biomart.
  mart  = biomaRt::useEnsembl(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
  genes = get_gene_info(entrez = colnames(expr), mart = mart)
  
  # For now, remove genes with multiple locations.
  genes = subset(genes, !grepl(',', genes$start_position))
  
  # For cis-eQTL, remove genes on Chr M or Y.
  genes = subset(genes, !chromosome_name %in% c('M', 'Y'))
  
  # Intersect genes and Ensembl data.
  common_genes = intersect(colnames(expr), genes$entrezgene_id)
  common_genes = as.character(common_genes)
  
  print(paste('Keeping', length(common_genes), 'out of', ncol(expr)))
  
  expr  = expr[,common_genes]
  genes = genes[common_genes,]
  stopifnot(colnames(expr) == rownames(genes))
  
  genes = genes[order(genes$chromosome_name),]
  expr  = expr[,rownames(genes)]
  
  # Convert from bp to Mb.
  genes$start_position = as.numeric(genes$start_position) * 1e-6
  genes$end_position   = as.numeric(genes$end_position)   * 1e-6
  
  # Combine genes and expr into one data structure.
  ge = merge(genes, t(expr), by = 'row.names')
  rownames(ge) = ge$Row.names
  ge = ge[,-1]
  
  results = data.frame(lod, matrix(0, nrow(lod), 8, dimnames = list(NULL, LETTERS[1:8])))
  
  ge = ge[lod$entrezgene_id,]
  stopifnot(all(lod$entrezgene_id == ge$entrezgene_id))
  
  for(i in 1:nrow(lod)) {
    
    print(paste('gene', i))
    
    # Get local genoprobs.
    chr    = lod$qtl_chr[i]
    pr_max = pull_genoprobpos(probs, map, chr = chr, pos = lod$qtl_pos[i])
    
    # Make expression matrix.
    e = t(as.matrix(ge[i,-(1:9),drop = FALSE]))
    colnames(e) = ge$entrezgene_id[i]

    # Get BLUP at closest marker.
    t1 = proc.time()[3]
    blup = fit1(genoprobs = pr_max, pheno = e, kinship = K[[chr]], addcovar = addcovar, cores = cores)
    print(proc.time()[3] - t1)

    # Make results.
    results[i, names(blup$coef[1:8])] = blup$coef[1:8]
    
  } # for(i)
  
  return(results)
  
} # cis_blup()


# Arguments:
# genoprobs: list of qtl2-style genoprobs objects.
# map:       list of qtl2-style marker locations (in Mb).
# expr:      numeric matrix containing gene expression. Must have Ensembl IDs as colnames, samples in rownames.
# addcovar:  numeric matrix containing additive covariates.
# K:         list of qtl2 style kinship matrices.
# snp_file:  full path to the Sanger VCF.
# cores:     integer that is the number of cores to use.
# Returns:   data.frame with cis-eQTL SNP positions and LOD scores.
cis_assoc = function() {
  
} # cis_assoc()



# Get info about genes using Entrez IDs from Biomart.
get_gene_info = function(entrez, mart = NA) {
  
  if(is.na(mart)) {
    mart = biomaRt::useEnsembl(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
  }
  filters = 'entrezgene_id'
  attributes = c('entrezgene_id', 'ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position', 'end_position', 'strand', 'gene_biotype')
  # , 'go_id', 'name_1006', 'namespace_1003'
  bm = getBM(attr = attributes, filt = filters, values = colnames(expr), mart = mart)
  
  # Keep only canonical chromosomes.
  bm = subset(bm, chromosome_name %in% c(1:19, 'X', 'Y', 'M'))
  
  # Handle one-to-many mappings.
  dupl = which(duplicated(bm$entrezgene_id))
  if(length(dupl) > 0) {
    dupl = bm$entrezgene_id[dupl]
    
    tmp = subset(bm, entrezgene_id %in% dupl)
    bm  = subset(bm, !entrezgene_id %in% dupl)
    spl = split(tmp, tmp$entrezgene_id)
    tmp = data.frame(t(sapply(spl, apply, 2, paste, collapse = ',')))
    tmp$entrezgene_id = sub(pat = ',[0-9]+$', repl = '', tmp$entrezgene_id)
    bm = rbind(bm, tmp)
  } # if(length(dupl) > 0)
  
  rownames(bm) = bm$entrezgene_id
  return(bm)
  
} # get_gene_info()


# Test code
load('/media/dmgatti/hdb/projects/TB/data/tufts_do_tb_qtl2_input.Rdata')

expr = read.csv(file = '/media/dmgatti/hdb/projects/TB/data/expression/tb_expr.csv')
rownames(expr) = expr[,1]
expr = as.matrix(expr[,-1])
colnames(expr) = sub('^X', '', colnames(expr))
expr = t(expr)

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]
samples  = intersect(rownames(expr), rownames(addcovar))
samples  = intersect(samples, rownames(probs[[1]]))
expr     = expr[samples,]
addcovar = addcovar[samples,]
for(i in seq_along(probs)) {
  K[[i]] = K[[i]][samples, samples]
}
probs = probs[samples,]

results = cis_eqtl(genoprobs = probs, map = map, expr = expr, addcovar = addcovar, K = K, cores = 4)
write.csv(results, file = '/media/dmgatti/hdb/projects/TB/results/qtl2/eqtl/tufts_cis_eqtl.csv', quote = F, row.names = FALSE)

blup = cis_blup(genoprobs = probs, map = map, expr = expr, addcovar = addcovar, K = K, lod = results, cores = 4)
write.csv(blup, file = '/media/dmgatti/hdb/projects/TB/results/qtl2/eqtl/tufts_cis_eqtl_blup.csv', quote = F, row.names = FALSE)

# Intersect eQTL with phenotype QTL.
peaks = readr::read_csv('/media/dmgatti/hdb/projects/TB/results/qtl2/gen_factor2/tb_qtl_peaks.csv')
blup  = readr::read_csv('/media/dmgatti/hdb/projects/TB/results/qtl2/eqtl/tufts_cis_eqtl_blup.csv',
                        guess_max = 20000)
blup = dplyr::filter(blup, qtl_lod >= 7.5)

res = NULL

for(i in 1:nrow(peaks)) {
  
  wh = which(blup$qtl_chr == peaks$chr[i] &  blup$qtl_pos >= peaks$ci_lo[i] & blup$qtl_pos <= peaks$ci_hi[i])
  if(length(wh) > 0) {
    print(paste(peaks$lodcolumn[i], peaks$chr[i]))
    print(paste(' ', blup$external_gene_name[wh]))
  } # if(length(wh) > 0)

} # for(i)

